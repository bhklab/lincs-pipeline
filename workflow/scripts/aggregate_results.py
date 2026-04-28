# ruff: noqa: EM101, EM102, F821, T201, TRY003

from pathlib import Path

import numpy as np
import pandas as pd

feature_metadata = pd.read_csv(
	snakemake.input.feature_metadata_tsv, sep='\t'
).sort_values('feature_index')

metadata_frames = []
for metadata_path in snakemake.input.metadata_tsvs:
	metadata_frame = pd.read_csv(metadata_path, sep='\t')
	if metadata_frame.shape[0] != 1:
		raise ValueError(f'Expected exactly one metadata row in {metadata_path}')
	metadata_frames.append(metadata_frame)

if not metadata_frames:
	raise ValueError('No compound metadata files were produced')

compound_metadata = pd.concat(metadata_frames, ignore_index=True).sort_values(
	'cmap_name'
)
if (
	compound_metadata['cmap_name'].isna().any()
	or (compound_metadata['cmap_name'].astype(str).str.strip() == '').any()
):
	raise ValueError('LINCS cmap_name values must be non-missing')
if compound_metadata['cmap_name'].duplicated().any():
	raise ValueError('LINCS cmap_name values must be unique')

signature_columns = {}
for signature_path in snakemake.input.signature_tsvs:
	signature_frame = pd.read_csv(signature_path, sep='\t').sort_values('feature_index')
	if signature_frame.empty:
		raise ValueError(f'Signature file is empty: {signature_path}')

	compound_id = signature_frame['compound_id'].iloc[0]
	observed_feature_index = signature_frame['feature_index'].tolist()
	expected_feature_index = feature_metadata['feature_index'].tolist()
	if observed_feature_index != expected_feature_index:
		raise ValueError(
			f'Feature order mismatch while aggregating signatures for {compound_id}'
		)

	signature_columns[compound_id] = signature_frame['coefficient'].to_numpy(
		dtype=np.float64
	)

ordered_compound_ids = compound_metadata['compound_id'].tolist()
ordered_cmap_names = compound_metadata['cmap_name'].tolist()
missing_signature_columns = [
	compound_id
	for compound_id in ordered_compound_ids
	if compound_id not in signature_columns
]
if missing_signature_columns:
	raise ValueError(
		'Missing signature files for compounds: '
		+ ', '.join(missing_signature_columns[:10])
	)

signature_matrix = np.column_stack(
	[signature_columns[compound_id] for compound_id in ordered_compound_ids]
)
signatures_df = pd.DataFrame(
	signature_matrix,
	index=feature_metadata['feature_name'].tolist(),
	columns=ordered_cmap_names,
)
signatures_df.index.name = 'feature_name'

signatures_output_path = Path(snakemake.output.signatures_csv)
signatures_output_path.parent.mkdir(parents=True, exist_ok=True)
signatures_df.to_csv(signatures_output_path)

coldata_output_path = Path(snakemake.output.coldata_csv)
coldata_output_path.parent.mkdir(parents=True, exist_ok=True)
compound_metadata.to_csv(coldata_output_path, index=False)

compound_metadata_output_path = Path(snakemake.output.compound_metadata_tsv)
compound_metadata_output_path.parent.mkdir(parents=True, exist_ok=True)
compound_metadata.to_csv(compound_metadata_output_path, sep='\t', index=False)

print(
	'[aggregate_results] '
	f'compound_rows={compound_metadata.shape[0]} '
	f'signature_shape={signatures_df.shape} '
	f'signatures={signatures_output_path}',
	flush=True,
)
