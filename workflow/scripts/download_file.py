# ruff: noqa: ANN001, ANN201, EM102, F821, T201, TRY003

import shutil
from pathlib import Path
from urllib.request import Request, urlopen


def parse_bool(value):
	if isinstance(value, bool):
		return value
	if value is None:
		return False
	return str(value).strip().lower() in {'1', 'true', 'yes', 'on'}


def download_http(source_uri, destination_tmp):
	request = Request(source_uri, headers={'User-Agent': 'lincs-pipeline/1.0'})
	with urlopen(request) as response, destination_tmp.open('wb') as handle:
		shutil.copyfileobj(response, handle)


source_uri = snakemake.params.source_uri
destination = Path(snakemake.output[0])
destination.parent.mkdir(parents=True, exist_ok=True)
overwrite = parse_bool(snakemake.params.overwrite)

if destination.exists() and destination.stat().st_size > 0 and not overwrite:
	print(f'[download_file] skip existing {destination}', flush=True)
else:
	tmp_path = destination.with_suffix(destination.suffix + '.tmp')
	if tmp_path.exists():
		tmp_path.unlink()

	print(f'[download_file] {source_uri} -> {destination}', flush=True)
	if source_uri.startswith(('http://', 'https://')):
		download_http(source_uri=source_uri, destination_tmp=tmp_path)
	else:
		raise ValueError(f'Unsupported source URI: {source_uri}')

	tmp_path.replace(destination)

if not destination.exists() or destination.stat().st_size <= 0:
	raise ValueError(f'Downloaded file is missing or empty: {destination}')
