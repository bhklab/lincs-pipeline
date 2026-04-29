# ruff: noqa: ANN001, ANN201, EM101, EM102, F821, T201, TRY003

from pathlib import Path

import numpy as np
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse

SIGNATURE_COLUMNS = [
	'compound_id',
	'feature_index',
	'feature_name',
	'feature_type',
	'feature_value',
	'coefficient',
]
METADATA_COLUMNS = [
	'compound_id',
	'cid',
	'cmap_name',
	'annotationdb_name',
	'annotationdb_smiles',
	'inchikey',
	'lincs_moa',
	'lincs_targets',
	'lincs_aliases',
	'annotationdb_aliases',
	'n_signatures',
	'n_genes',
	'n_features',
]


def parse_bool(value):
	if isinstance(value, bool):
		return value
	if value is None:
		return False
	return str(value).strip().lower() in {'1', 'true', 'yes', 'on'}


def build_design_matrix(compound_input, batch_categories, cell_categories):
	pert_time = compound_input['pert_time'].to_numpy(dtype=np.float64).reshape(-1, 1)
	pert_dose = compound_input['pert_dose'].to_numpy(dtype=np.float64).reshape(-1, 1)

	batch_values = compound_input['bead_batch'].astype(str).to_numpy()
	cell_values = compound_input['cell_iname'].astype(str).to_numpy()

	batch_features = (batch_values[:, None] == batch_categories[None, :]).astype(
		np.float64
	)
	cell_features = (cell_values[:, None] == cell_categories[None, :]).astype(
		np.float64
	)

	return np.hstack([pert_time, pert_dose, batch_features, cell_features])


def fit_coefficients(design_with_intercept, response_matrix):
	missing_mask = np.isnan(response_matrix)
	if not missing_mask.any():
		coefficient_matrix, *_ = np.linalg.lstsq(
			design_with_intercept,
			response_matrix,
			rcond=None,
		)
		return coefficient_matrix[1:, :], 0, 0, 0

	all_missing_genes = missing_mask.all(axis=0)
	if all_missing_genes.any():
		raise ValueError(
			'Some selected genes have no expression values for this compound'
		)

	coefficient_matrix = np.empty(
		(design_with_intercept.shape[1] - 1, response_matrix.shape[1]),
		dtype=np.float64,
	)
	for gene_index in range(response_matrix.shape[1]):
		response = response_matrix[:, gene_index]
		observed = ~np.isnan(response)
		gene_coefficients, *_ = np.linalg.lstsq(
			design_with_intercept[observed, :],
			response[observed],
			rcond=None,
		)
		coefficient_matrix[:, gene_index] = gene_coefficients[1:]

	return (
		coefficient_matrix,
		int(missing_mask.sum()),
		int(missing_mask.any(axis=1).sum()),
		int(missing_mask.any(axis=0).sum()),
	)


compound_id = snakemake.wildcards.compound_id
compound_input_path = Path(snakemake.params.compound_input_path)
manifest_df = pd.read_csv(snakemake.input.manifest_tsv, sep='\t')
manifest_row = manifest_df.loc[manifest_df['compound_id'] == compound_id]
if manifest_row.empty:
	raise ValueError(f'Compound {compound_id} is missing from the manifest')
manifest_row = manifest_row.iloc[0]

feature_metadata = pd.read_csv(
	snakemake.input.feature_metadata_tsv, sep='\t'
).sort_values('feature_index')
selected_gene_metadata = pd.read_csv(snakemake.input.gene_metadata_tsv, sep='\t')
compound_input = pd.read_csv(compound_input_path, sep='\t').sort_values('sig_id')

batch_categories = (
	feature_metadata.loc[
		feature_metadata['feature_type'] == 'bead_batch', 'feature_value'
	]
	.astype(str)
	.to_numpy()
)
cell_categories = (
	feature_metadata.loc[
		feature_metadata['feature_type'] == 'cell_iname', 'feature_value'
	]
	.astype(str)
	.to_numpy()
)

design_matrix = build_design_matrix(
	compound_input=compound_input,
	batch_categories=batch_categories,
	cell_categories=cell_categories,
)
sig_ids = compound_input['sig_id'].tolist()
gene_ids = [str(gene_id) for gene_id in selected_gene_metadata['gene_id'].tolist()]

gctx = parse(
	str(snakemake.input.expression_matrix),
	cid=sig_ids,
	rid=gene_ids,
)
expression_df = gctx.data_df
expression_df.index = expression_df.index.astype(str)
expression_df = expression_df.reindex(index=gene_ids, columns=sig_ids)

response_matrix = expression_df.to_numpy(dtype=np.float64).T
design_with_intercept = np.column_stack(
	[np.ones(design_matrix.shape[0], dtype=np.float64), design_matrix]
)
(
	coefficient_matrix,
	missing_expression_values,
	signatures_with_missing_expression,
	genes_with_missing_expression,
) = fit_coefficients(
	design_with_intercept=design_with_intercept,
	response_matrix=response_matrix,
)
if parse_bool(snakemake.params.only_nonnegative_signatures):
	coefficient_matrix = np.abs(coefficient_matrix)

mean_coefficients = coefficient_matrix.mean(axis=1)
if mean_coefficients.shape[0] != feature_metadata.shape[0]:
	raise ValueError(
		'Computed coefficient length does not match the feature metadata table'
	)

signature_output = feature_metadata.copy()
signature_output.insert(0, 'compound_id', compound_id)
signature_output['coefficient'] = mean_coefficients
signature_output = signature_output.loc[:, SIGNATURE_COLUMNS]

compound_metadata = pd.DataFrame(
	[
		{
			'compound_id': compound_id,
			'cid': manifest_row['cid'],
			'cmap_name': manifest_row['cmap_name'],
			'annotationdb_name': manifest_row['annotationdb_name'],
			'annotationdb_smiles': manifest_row['annotationdb_smiles'],
			'inchikey': manifest_row['inchikey'],
			'lincs_moa': manifest_row['lincs_moa'],
			'lincs_targets': manifest_row['lincs_targets'],
			'lincs_aliases': manifest_row['lincs_aliases'],
			'annotationdb_aliases': manifest_row['annotationdb_aliases'],
			'n_signatures': int(compound_input.shape[0]),
			'n_genes': int(len(gene_ids)),
			'n_features': int(feature_metadata.shape[0]),
		}
	],
	columns=METADATA_COLUMNS,
)

signature_output_path = Path(snakemake.output.signature_tsv)
signature_output_path.parent.mkdir(parents=True, exist_ok=True)
signature_output.to_csv(signature_output_path, sep='\t', index=False)

metadata_output_path = Path(snakemake.output.metadata_tsv)
metadata_output_path.parent.mkdir(parents=True, exist_ok=True)
compound_metadata.to_csv(metadata_output_path, sep='\t', index=False)

print(
	'[compute_compound_signature] '
	f'compound_id={compound_id} '
	f'n_signatures={compound_input.shape[0]} '
	f'n_genes={len(gene_ids)} '
	f'n_features={feature_metadata.shape[0]} '
	f'missing_expression_values={missing_expression_values} '
	f'signatures_with_missing_expression={signatures_with_missing_expression} '
	f'genes_with_missing_expression={genes_with_missing_expression}',
	flush=True,
)
