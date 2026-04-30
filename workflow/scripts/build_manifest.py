# ruff: noqa: ANN001, ANN201, EM101, F821, T201, TRY003

import csv
import hashlib
import re
from pathlib import Path

import pandas as pd

SIGNATURE_INPUT_COLUMNS = [
	'sig_id',
	'pert_time',
	'pert_dose',
	'bead_batch',
	'cell_iname',
]
MANIFEST_COLUMNS = [
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
]
FEATURE_METADATA_COLUMNS = [
	'feature_index',
	'feature_name',
	'feature_type',
	'feature_value',
]
SOURCE_CELL_METADATA_COLUMNS = [
	'cell_iname',
	'cellosaurus_id',
	'ccle_name',
	'cell_type',
	'primary_disease',
	'subtype',
	'cell_lineage',
	'donor_age',
	'donor_sex',
	'donor_ethnicity',
	'donor_tumor_phase',
	'growth_pattern',
	'provider_name',
	'provider_catalog_id',
	'cell_alias',
	'doubling_time',
]
CELL_LINE_METADATA_COLUMNS = [
	'LINCS.Cell.Name',
	'Cellosaurus.ID',
	'In.AnnotationDB.CellLine',
	'Cell.Line.Name',
	'LINCS.CCLE.Name',
	'LINCS.Cell.Type',
	'LINCS.Primary.Disease',
	'LINCS.Subtype',
	'LINCS.Cell.Lineage',
	'LINCS.Donor.Age',
	'LINCS.Donor.Sex',
	'LINCS.Donor.Ethnicity',
	'LINCS.Donor.Tumor.Phase',
	'LINCS.Growth.Pattern',
	'LINCS.Provider.Name',
	'LINCS.Provider.Catalog.ID',
	'LINCS.Cell.Alias',
	'LINCS.Doubling.Time',
]


def parse_bool(value):
	if isinstance(value, bool):
		return value
	if value is None:
		return False
	return str(value).strip().lower() in {'1', 'true', 'yes', 'on'}


def parse_optional_int(value):
	if value is None:
		return None
	if isinstance(value, str) and value.strip().lower() in {'', 'none', 'null'}:
		return None
	return int(value)


def parse_selection_list(values):
	if values is None:
		return []
	if isinstance(values, str):
		values = [item.strip() for item in values.split(',')]

	normalized = []
	seen = set()
	for value in values:
		item = str(value).strip()
		if not item or item in seen:
			continue
		normalized.append(item)
		seen.add(item)
	return normalized


def clean_optional_string(value):
	if pd.isna(value):
		return None
	text = str(value).strip()
	if not text:
		return None
	if text.upper() in {'NA', 'N/A', 'NULL', 'NONE', 'NAN'}:
		return None
	return text


def normalize_pipe_list(values):
	cleaned = []
	seen = set()
	for value in values:
		text = clean_optional_string(value)
		if text is None or text in seen:
			continue
		cleaned.append(text)
		seen.add(text)
	return cleaned


def stable_compound_id(cmap_name, cid=None, inchikey=None):
	slug = re.sub(r'[^A-Za-z0-9]+', '-', str(cmap_name).strip().lower()).strip('-')
	slug = slug or 'compound'
	identity = cid if cid not in (None, '') else inchikey or cmap_name
	digest = hashlib.sha1(f'{cmap_name}\t{identity}'.encode('utf-8')).hexdigest()[:10]
	if cid not in (None, '') and not pd.isna(cid):
		return f'cid_{int(cid)}_{slug[:48]}_{digest}'
	return f'lincs_{slug[:48]}_{digest}'


def normalize_cid(value):
	if pd.isna(value):
		return ''
	return int(float(value))


def normalize_accession_key(value):
	text = clean_optional_string(value)
	if text is None:
		return ''
	return text.upper()


def normalize_cell_name_key(value):
	text = clean_optional_string(value)
	if text is None:
		return ''
	return re.sub(r'[^A-Z0-9]+', '', text.upper())


def clean_string_columns(frame):
	for column in frame.columns:
		if pd.api.types.is_object_dtype(frame[column]):
			frame[column] = frame[column].map(clean_optional_string)
	return frame


def write_feature_metadata(path, batch_categories, cell_categories):
	rows = [
		{
			'feature_index': 0,
			'feature_name': 'pert_time',
			'feature_type': 'scalar',
			'feature_value': 'pert_time',
		},
		{
			'feature_index': 1,
			'feature_name': 'pert_dose',
			'feature_type': 'scalar',
			'feature_value': 'pert_dose',
		},
	]

	feature_index = len(rows)
	for batch in batch_categories:
		rows.append(
			{
				'feature_index': feature_index,
				'feature_name': f'bead_batch:{batch}',
				'feature_type': 'bead_batch',
				'feature_value': batch,
			}
		)
		feature_index += 1

	for cell in cell_categories:
		rows.append(
			{
				'feature_index': feature_index,
				'feature_name': f'cell_iname:{cell}',
				'feature_type': 'cell_iname',
				'feature_value': cell,
			}
		)
		feature_index += 1

	feature_frame = pd.DataFrame(rows, columns=FEATURE_METADATA_COLUMNS)
	feature_frame.to_csv(path, sep='\t', index=False)
	return feature_frame


def build_cell_line_metadata(cell_metadata, annotationdb_cell_lines, cell_categories):
	cell_source = cell_metadata.drop_duplicates(subset=['cell_iname']).copy()
	cell_source = (
		cell_source.set_index('cell_iname').reindex(cell_categories).reset_index()
	)

	annotationdb_cell_lines = annotationdb_cell_lines.rename(
		columns={
			'name': 'annotationdb_cell_line_name',
			'accession': 'annotationdb_cellosaurus_id',
		}
	)
	annotationdb_cell_lines = clean_string_columns(annotationdb_cell_lines)
	annotationdb_cell_lines['accession_key'] = annotationdb_cell_lines[
		'annotationdb_cellosaurus_id'
	].map(normalize_accession_key)
	annotationdb_cell_lines['name_key'] = annotationdb_cell_lines[
		'annotationdb_cell_line_name'
	].map(normalize_cell_name_key)

	by_accession = (
		annotationdb_cell_lines[annotationdb_cell_lines['accession_key'] != '']
		.drop_duplicates(subset=['accession_key'])
		.set_index('accession_key')
	)
	by_name = (
		annotationdb_cell_lines[annotationdb_cell_lines['name_key'] != '']
		.drop_duplicates(subset=['name_key'])
		.set_index('name_key')
	)

	cell_source['source_accession_key'] = cell_source['cellosaurus_id'].map(
		normalize_accession_key
	)
	cell_source['source_name_key'] = cell_source['cell_iname'].map(
		normalize_cell_name_key
	)
	cell_source['matched_name'] = None
	cell_source['matched_accession'] = None

	accession_matches = cell_source['source_accession_key'].isin(by_accession.index)
	if accession_matches.any():
		matched_rows = by_accession.loc[
			cell_source.loc[accession_matches, 'source_accession_key']
		]
		cell_source.loc[accession_matches, 'matched_name'] = matched_rows[
			'annotationdb_cell_line_name'
		].to_numpy()
		cell_source.loc[accession_matches, 'matched_accession'] = matched_rows[
			'annotationdb_cellosaurus_id'
		].to_numpy()

	name_matches = cell_source['matched_accession'].isna() & cell_source[
		'source_name_key'
	].isin(by_name.index)
	if name_matches.any():
		matched_rows = by_name.loc[cell_source.loc[name_matches, 'source_name_key']]
		cell_source.loc[name_matches, 'matched_name'] = matched_rows[
			'annotationdb_cell_line_name'
		].to_numpy()
		cell_source.loc[name_matches, 'matched_accession'] = matched_rows[
			'annotationdb_cellosaurus_id'
		].to_numpy()

	out = pd.DataFrame(
		{
			'LINCS.Cell.Name': cell_source['cell_iname'],
			'Cellosaurus.ID': cell_source['cellosaurus_id'].combine_first(
				cell_source['matched_accession']
			),
			'In.AnnotationDB.CellLine': cell_source['matched_accession'].notna(),
			'Cell.Line.Name': cell_source['matched_name'],
			'LINCS.CCLE.Name': cell_source['ccle_name'],
			'LINCS.Cell.Type': cell_source['cell_type'],
			'LINCS.Primary.Disease': cell_source['primary_disease'],
			'LINCS.Subtype': cell_source['subtype'],
			'LINCS.Cell.Lineage': cell_source['cell_lineage'],
			'LINCS.Donor.Age': cell_source['donor_age'],
			'LINCS.Donor.Sex': cell_source['donor_sex'],
			'LINCS.Donor.Ethnicity': cell_source['donor_ethnicity'],
			'LINCS.Donor.Tumor.Phase': cell_source['donor_tumor_phase'],
			'LINCS.Growth.Pattern': cell_source['growth_pattern'],
			'LINCS.Provider.Name': cell_source['provider_name'],
			'LINCS.Provider.Catalog.ID': cell_source['provider_catalog_id'],
			'LINCS.Cell.Alias': cell_source['cell_alias'],
			'LINCS.Doubling.Time': cell_source['doubling_time'],
		},
		columns=CELL_LINE_METADATA_COLUMNS,
	)
	return out.sort_values('LINCS.Cell.Name')


cfg = snakemake.config
selection_cfg = cfg['selection']

only_named_cells = parse_bool(selection_cfg.get('only_named_cells', True))
only_cancer_cells = parse_bool(selection_cfg.get('only_cancer_cells', True))
valid_gene_types = parse_selection_list(selection_cfg.get('valid_gene_types'))
valid_dose_units = parse_selection_list(selection_cfg.get('valid_dose_units'))
valid_pert_types = parse_selection_list(selection_cfg.get('valid_pert_types'))
include_compounds = set(parse_selection_list(selection_cfg.get('include_compounds')))
exclude_compounds = set(parse_selection_list(selection_cfg.get('exclude_compounds')))
compound_limit = parse_optional_int(selection_cfg.get('compound_limit'))

cell_metadata = pd.read_csv(
	snakemake.input.cell_metadata,
	sep='\t',
	usecols=SOURCE_CELL_METADATA_COLUMNS,
)
cell_metadata = clean_string_columns(cell_metadata)
if only_named_cells:
	cell_metadata = cell_metadata.dropna(subset=['cellosaurus_id', 'ccle_name'])
if only_cancer_cells:
	cell_metadata = cell_metadata[cell_metadata['cell_type'] == 'tumor']
keep_cells = set(cell_metadata['cell_iname'].tolist())

gene_metadata = pd.read_csv(
	snakemake.input.gene_metadata,
	sep='\t',
	usecols=['gene_id', 'gene_symbol', 'ensembl_id', 'gene_title', 'feature_space'],
)
if valid_gene_types:
	gene_metadata = gene_metadata[gene_metadata['feature_space'].isin(valid_gene_types)]
gene_metadata = gene_metadata.drop_duplicates(subset=['gene_id']).sort_values('gene_id')

signature_data = pd.read_csv(
	snakemake.input.signature_metadata,
	sep='\t',
	usecols=[
		'sig_id',
		'cmap_name',
		'pert_time',
		'pert_dose',
		'pert_dose_unit',
		'cell_iname',
		'bead_batch',
		'pert_type',
	],
	low_memory=False,
)
if valid_dose_units:
	signature_data = signature_data[
		signature_data['pert_dose_unit'].isin(valid_dose_units)
	]
if valid_pert_types:
	signature_data = signature_data[signature_data['pert_type'].isin(valid_pert_types)]
signature_data = signature_data[signature_data['cell_iname'].isin(keep_cells)]
signature_data = signature_data.dropna(
	subset=['sig_id', 'cmap_name', 'pert_time', 'pert_dose', 'bead_batch', 'cell_iname']
)

compound_metadata = pd.read_csv(
	snakemake.input.compound_metadata,
	sep='\t',
	usecols=['cmap_name', 'pert_id', 'target', 'moa', 'inchi_key', 'compound_aliases'],
)
common_compounds = set(signature_data['cmap_name']).intersection(
	set(compound_metadata['cmap_name'])
)
signature_data = signature_data[signature_data['cmap_name'].isin(common_compounds)]
compound_metadata = compound_metadata[
	compound_metadata['cmap_name'].isin(common_compounds)
]

if include_compounds:
	signature_data = signature_data[signature_data['cmap_name'].isin(include_compounds)]
	compound_metadata = compound_metadata[
		compound_metadata['cmap_name'].isin(include_compounds)
	]
if exclude_compounds:
	signature_data = signature_data[
		~signature_data['cmap_name'].isin(exclude_compounds)
	]
	compound_metadata = compound_metadata[
		~compound_metadata['cmap_name'].isin(exclude_compounds)
	]

signature_data = signature_data.sort_values(['cmap_name', 'sig_id']).reset_index(
	drop=True
)
compound_metadata = compound_metadata.sort_values(['cmap_name', 'pert_id']).reset_index(
	drop=True
)

annotationdb = pd.read_csv(snakemake.input.annotationdb_cache)
if 'Unnamed: 0' in annotationdb.columns:
	annotationdb = annotationdb.drop(columns=['Unnamed: 0'])
annotationdb = annotationdb.rename(
	columns={
		'name': 'annotationdb_name',
		'cid': 'cid',
		'smiles': 'annotationdb_smiles',
		'inchikey': 'inchikey',
		'mapped_name': 'mapped_name',
	}
)
annotationdb = annotationdb.dropna(subset=['inchikey']).copy()

annotationdb_cell_lines = pd.read_csv(snakemake.input.annotationdb_cell_line_cache)
if 'Unnamed: 0' in annotationdb_cell_lines.columns:
	annotationdb_cell_lines = annotationdb_cell_lines.drop(columns=['Unnamed: 0'])

manifest_rows = []
summary_counts = {
	'missing_inchikey': 0,
	'unmapped_annotationdb': 0,
	'unmapped_cell_lines': 0,
}
compound_input_dir = Path(snakemake.params.compound_input_dir)
compound_input_dir.mkdir(parents=True, exist_ok=True)
for existing_file in compound_input_dir.glob('*.tsv'):
	existing_file.unlink()

candidate_compounds = sorted(signature_data['cmap_name'].drop_duplicates().tolist())
if compound_limit is not None:
	candidate_compounds = candidate_compounds[:compound_limit]

for cmap_name in candidate_compounds:
	compound_signature_rows = signature_data[signature_data['cmap_name'] == cmap_name]
	compound_rows = compound_metadata[compound_metadata['cmap_name'] == cmap_name]

	inchikeys = normalize_pipe_list(compound_rows['inchi_key'].tolist())
	if not inchikeys:
		summary_counts['missing_inchikey'] += 1
		continue

	annotationdb_rows = annotationdb[annotationdb['inchikey'].isin(inchikeys)].copy()
	if annotationdb_rows.empty:
		summary_counts['unmapped_annotationdb'] += 1
		selected_annotationdb_row = {
			'cid': '',
			'annotationdb_name': '',
			'annotationdb_smiles': '',
			'inchikey': inchikeys[0],
			'mapped_name': '',
		}
		cid = ''
	else:
		annotationdb_rows = annotationdb_rows.sort_values(
			by=['cid', 'annotationdb_name', 'mapped_name', 'inchikey'],
			na_position='last',
		)
		selected_annotationdb_row = annotationdb_rows.iloc[0]
		cid = normalize_cid(selected_annotationdb_row['cid'])

	compound_id = stable_compound_id(
		cmap_name=cmap_name,
		cid=cid,
		inchikey=clean_optional_string(selected_annotationdb_row['inchikey']),
	)

	compound_input_path = compound_input_dir / f'{compound_id}.tsv'
	compound_signature_rows.loc[:, SIGNATURE_INPUT_COLUMNS].to_csv(
		compound_input_path,
		sep='\t',
		index=False,
	)

	manifest_rows.append(
		{
			'compound_id': compound_id,
			'cid': cid,
			'cmap_name': cmap_name,
			'annotationdb_name': clean_optional_string(
				selected_annotationdb_row['annotationdb_name']
			),
			'annotationdb_smiles': clean_optional_string(
				selected_annotationdb_row['annotationdb_smiles']
			),
			'inchikey': clean_optional_string(selected_annotationdb_row['inchikey']),
			'lincs_moa': '|'.join(normalize_pipe_list(compound_rows['moa'].tolist())),
			'lincs_targets': '|'.join(
				normalize_pipe_list(compound_rows['target'].tolist())
			),
			'lincs_aliases': '|'.join(
				normalize_pipe_list(compound_rows['compound_aliases'].tolist())
			),
			'annotationdb_aliases': '|'.join(
				normalize_pipe_list(
					[]
					if annotationdb_rows.empty
					else annotationdb_rows['mapped_name'].tolist()
				)
			),
			'n_signatures': int(compound_signature_rows.shape[0]),
		}
	)

if not manifest_rows:
	raise ValueError('No compounds matched the configured selection filters.')

manifest_df = pd.DataFrame(manifest_rows, columns=MANIFEST_COLUMNS).sort_values(
	['compound_id']
)
manifest_path = Path(snakemake.output.manifest_tsv)
manifest_path.parent.mkdir(parents=True, exist_ok=True)
manifest_df.to_csv(manifest_path, sep='\t', index=False)

batch_categories = sorted(signature_data['bead_batch'].drop_duplicates().tolist())
cell_categories = sorted(signature_data['cell_iname'].drop_duplicates().tolist())
feature_metadata = write_feature_metadata(
	path=Path(snakemake.output.feature_metadata_tsv),
	batch_categories=batch_categories,
	cell_categories=cell_categories,
)

cell_line_metadata = build_cell_line_metadata(
	cell_metadata=cell_metadata,
	annotationdb_cell_lines=annotationdb_cell_lines,
	cell_categories=cell_categories,
)
summary_counts['unmapped_cell_lines'] = int(
	(~cell_line_metadata['In.AnnotationDB.CellLine']).sum()
)
cell_line_metadata_path = Path(snakemake.output.cell_line_metadata_tsv)
cell_line_metadata_path.parent.mkdir(parents=True, exist_ok=True)
cell_line_metadata.to_csv(cell_line_metadata_path, sep='\t', index=False)

gene_metadata_path = Path(snakemake.output.gene_metadata_tsv)
gene_metadata_path.parent.mkdir(parents=True, exist_ok=True)
gene_metadata.to_csv(gene_metadata_path, sep='\t', index=False)

summary_rows = [
	{'metric': 'manifest_rows', 'value': int(manifest_df.shape[0])},
	{'metric': 'filtered_signature_rows', 'value': int(signature_data.shape[0])},
	{'metric': 'selected_gene_rows', 'value': int(gene_metadata.shape[0])},
	{'metric': 'feature_rows', 'value': int(feature_metadata.shape[0])},
	{'metric': 'batch_categories', 'value': int(len(batch_categories))},
	{'metric': 'cell_categories', 'value': int(len(cell_categories))},
	{
		'metric': 'unmapped_cell_lines',
		'value': int(summary_counts['unmapped_cell_lines']),
	},
	{
		'metric': 'skipped_missing_inchikey',
		'value': int(summary_counts['missing_inchikey']),
	},
	{
		'metric': 'unmapped_annotationdb',
		'value': int(summary_counts['unmapped_annotationdb']),
	},
]
summary_path = Path(snakemake.output.summary_tsv)
summary_path.parent.mkdir(parents=True, exist_ok=True)
with summary_path.open('w', newline='', encoding='utf-8') as handle:
	writer = csv.DictWriter(handle, fieldnames=['metric', 'value'], delimiter='\t')
	writer.writeheader()
	writer.writerows(summary_rows)

print(
	'[build_manifest] '
	f'manifest_rows={manifest_df.shape[0]} '
	f'filtered_signature_rows={signature_data.shape[0]} '
	f'selected_genes={gene_metadata.shape[0]} '
	f'features={feature_metadata.shape[0]} '
	f'cell_lines={cell_line_metadata.shape[0]}',
	flush=True,
)
