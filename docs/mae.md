# MAE Structure

The final object is written to `data/results/lincs_MultiAssayExperiment.rds`. It is compound-level: signatures are aggregated before MAE construction.

| Component | Structure |
| --- | --- |
| Assays | `signatures`, containing one compound-level expression signature matrix. |
| Assay columns | `LINCS.CMap.Name`, the source LINCS compound name. |
| Assay rows | Selected gene/signature features from the processed LINCS matrices. |
| sampleMap | One row per compound with `assay`, `primary`, and `colname`; both `primary` and `colname` are `LINCS.CMap.Name`. |
| colData | Compound-level metadata. Each row is keyed by `LINCS.CMap.Name`. |
| rowData | Selected feature metadata, with source `feature_*` fields renamed to `LINCS.Feature.*` and Cellosaurus fields populated for `cell_iname` rows. |

## metadata(mae)

| Object | Source | Purpose |
| --- | --- | --- |
| `Pipeline` | `config/pipeline.yaml` and Snakefile params | Named list containing `ID`, `Version`, and parsed run `Config`. |
| `Selected.Gene.Metadata` | `data/procdata/metadata/lincs_selected_gene_metadata.tsv` | Gene metadata for features retained in the assay. |
| `Cell.Line.Metadata` | `data/procdata/metadata/lincs_cell_line_metadata.tsv` plus AnnotationDB `/cell_line/all` | Selected LINCS cell-line metadata keyed by `LINCS.Cell.Name`. |
| `Drug.Metadata` | `data/procdata/metadata/lincs_compound_metadata.tsv` plus AnnotationDB cache | Harmonized compound metadata keyed by `LINCS.CMap.Name`. |
