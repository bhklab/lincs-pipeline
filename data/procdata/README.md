# procdata

Processed LINCS intermediates live here. The pipeline assumes generated metadata under `procdata/metadata/`, per-compound intermediates under processing subdirectories, and aggregate signature tables before MAE construction.

Expected usage depends on filters and compound limits. Full runs can require tens of GB; small subset runs are much smaller. Files here are ignored by git.
