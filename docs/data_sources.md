# Data Sources

## LINCS 2020

The workflow downloads LINCS 2020 metadata and Level 5 treatment expression
matrices from the public CLUE S3 bucket configured in `config/pipeline.yaml`.
Downloaded files are written under `data/rawdata/`.

Configured source files:

- `cellinfo_beta.txt`
- `compoundinfo_beta.txt`
- `geneinfo_beta.txt`
- `siginfo_beta.txt`
- `level5_beta_trt_cp_n720216x12328.gctx`

The control matrix URL is retained in the config but is disabled by default.

## AnnotationDB

AnnotationDB enrichment uses the configurable API base URL at
`annotationdb_api`. The workflow appends `/compound/all` for compound metadata
and `/cell_line/all` for Cellosaurus cell-line identifiers. Downloaded caches
are stored under `data/rawdata/metadata/`, while generated manifests and
processed compound and cell-line metadata are stored under
`data/procdata/metadata/`.
