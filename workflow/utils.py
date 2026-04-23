from typing import List, Tuple

import numpy as np
import pandas as pd
from damply import dirs
from sklearn.linear_model import LinearRegression


def featurize_signatures(
	expression_vector: pd.Series,
	signature_info: pd.DataFrame,
	cell_encoder,
	batch_encoder,
):

	times = []
	batches = []
	doses = []
	cells = []
	for sig_id in expression_vector.index:
		sig_paramters = signature_info[signature_info['sig_id'] == sig_id]

		sig_time = float(sig_paramters['pert_time'].item())
		sig_dose = float(sig_paramters['pert_dose'].item())
		sig_batch = sig_paramters['bead_batch'].item()
		sig_cell = sig_paramters['cell_iname'].item()
		times.append(sig_time)
		doses.append(sig_dose)

		enc_batch = batch_encoder.transform(
			np.array(sig_batch).reshape(-1, 1)
		).toarray()
		enc_cell = cell_encoder.transform(np.array(sig_cell).reshape(-1, 1)).toarray()
		batches.append(enc_batch)
		cells.append(enc_cell)

	doses = np.array(doses).reshape(-1, 1)
	times = np.array(times).reshape(-1, 1)

	X_scalars = np.concatenate((times, doses), axis=1)
	X_batch = np.concatenate(batches, axis=0)
	X_cell = np.concatenate(cells, axis=0)
	X = np.concatenate((X_scalars, X_batch, X_cell), axis=1)

	return X


def compute_gene_specific_signature(
	expression_vector: pd.Series,
	signature_info: pd.DataFrame,
	cell_encoder,
	batch_encoder,
	only_nonneg_signatures: bool = True,
):
	"""
	Compute the gene specific signature from the data.

	Expression Vector: their proprietery. Should be a 1xN vector. Each entry has signature value. Si

	"""
	y = [expression_vector[x] for x in expression_vector.index]
	X = featurize_signatures(
		expression_vector=expression_vector,
		signature_info=signature_info,
		cell_encoder=cell_encoder,
		batch_encoder=batch_encoder,
	)

	model = LinearRegression()
	model.fit(X, y)

	if only_nonneg_signatures:
		coefficient_signature = np.abs(model.coef_)
	else:
		coefficient_signature = model.coef_

	return coefficient_signature


def load_metadata(
	gene_file: str = 'geneinfo_beta.txt',
	cell_file: str = 'cellinfo_beta.txt',
	cpd_file: str = 'compoundinfo_beta.txt',
	sig_file: str = 'siginfo_beta.txt',
	only_named_cells: bool = True,
	only_cancer_cells: bool = True,
	valid_gene_types: List[str] = ['landmark'],
	valid_dose_units: List[str] = ['uM'],
) -> Tuple[pd.Dataframe, pd.DataFrame]:

	cell_line_metadata = pd.read_csv(dirs.RAWDATA / cell_file, sep='\t')

	if only_named_cells:
		# filter only those cell lines with either a CCLE Name or Cellosaurus ID
		cell_line_metadata = cell_line_metadata.dropna(
			how='any', subset=['cellosaurus_id', 'ccle_name']
		)

	if only_cancer_cells:
		# filter for only tumor cells
		cell_line_metadata = cell_line_metadata[
			cell_line_metadata['cell_type'] == 'tumor'
		]

	keep_cells = list(cell_line_metadata['cell_iname'])
	gene_metadata = pd.read_csv(dirs.RAWDATA / gene_file, sep='\t')

	gene_metadata = gene_metadata[gene_metadata['feature_space'].isin(valid_gene_types)]

	signature_data = pd.read_csv(dirs.RAWDATA / sig_file, sep='\t')

	signature_data = signature_data[
		signature_data['pert_dose_unit'].isin(valid_dose_units)
	]

	cpd_metadata = pd.read_csv(dirs.RAWDATA / cpd_file, sep='\t')

	# select only those compounds which have at least one computed signature AND metadata
	common_cpds = set(signature_data['cmap_name']).intersection(
		set(cpd_metadata['cmap_name'])
	)

	cpd_metadata = cpd_metadata[cpd_metadata['cmap_name'].isin(common_cpds)]

	signature_data = signature_data[signature_data['cell_iname'].isin(keep_cells)]
	signature_data = signature_data[signature_data['cmap_name'].isin(common_cpds)]

	# gene_order = [x.item() for x in sorted(list(pd.unique(gene_metadata['gene_id'])))]
	signature_data = signature_data.dropna(
		subset=['pert_time', 'pert_dose', 'bead_batch', 'cell_iname']
	)
	# cpd_metadata = cpd_metadata.dropna(subset=['inchi_key'])
	return cpd_metadata, signature_data, gene_metadata
