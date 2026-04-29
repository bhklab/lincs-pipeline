# Developer Notes

## Scope Decisions

- The default curation uses LINCS2020 treatment signatures and compound, gene, cell, and signature metadata from CLUE-hosted source files.
- The default signature selection keeps named cancer cell lines, landmark genes, micromolar doses, compound perturbations, and nonnegative signature values.
- Optional `include_compounds`, `exclude_compounds`, and `compound_limit` settings intentionally subset the source data when configured.
- The control matrix is configured as an optional download but is not part of the default MAE output.

## Identifier Decisions

- The compound-level MAE key is `LINCS.CMap.Name`, using the LINCS source compound name.
- `LINCS.CMap.Name` must be unique and non-missing before MAE construction. Future data that violates this should fail clearly rather than creating a public replacement key.
- Internal file-safe stems may still be used for intermediate filenames, but they are not public identifiers.
- Public columns with HDD-shared names are intended to be directly joinable to the base HDD; source-specific fields use source-specific prefixes.

## Metadata Decisions

- AnnotationDB enrichment is intentionally minimal: PubChem CID, AnnotationDB name, AnnotationDB SMILES, AnnotationDB aliases, and a match flag.
- The AnnotationDB `/compound/all` response is cached once at `data/rawdata/metadata/all_adb_compounds.csv` and downstream joins read the cache.
- Assay-valid compounds are retained when PubChem CID is missing if enough source metadata exists to build signatures.
- Derived display-name columns are omitted from public outputs; users can choose LINCS or AnnotationDB names downstream.

## Assay Decisions

- The `signatures` assay is derived by computing compound-level signatures from selected LINCS level-5 expression data.
- The assay is not a direct one-file copy from source metadata; it is a pipeline-derived summary intended to align one column per curated compound.
- Sparse missing expression values in the source GCTX are handled per gene by fitting the signature model on non-missing signatures for that gene. A compound still fails if any selected gene has no expression values for that compound.
