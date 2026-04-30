# LINCS Pipeline

This pipeline curates LINCS2020 compound perturbation signatures into processed compound-level signature tables, a `MultiAssayExperiment`, and an MAE-derived tabular archive. Sources include CLUE-hosted LINCS expression and metadata files plus AnnotationDB compound and cell-line metadata.

Run from the repository root:

```bash
pixi run snakemake --cores <n>
```

Edit `config/pipeline.yaml` to change LINCS filters, compound limits, source URLs, paths, and the AnnotationDB API base URL.
