# colData Columns

Final `colData` is built in `workflow/scripts/build_mae.R` from processed LINCS compound metadata, computed compound signature summaries, and the cached AnnotationDB table. Columns are exported from the MAE to `data/results/lincs_tables/colData.tsv`.

| Column | Type | Description | Computed from / origin |
| --- | --- | --- | --- |
| `LINCS.CMap.Name` | character | Compound-level MAE primary key using the LINCS CMap compound name. | Source compound metadata `cmap_name`; values must be unique and non-missing before MAE construction. |
| `Pubchem.CID` | integer | PubChem compound identifier used to join back to the base HDD. | Processed compound metadata `cid`, after AnnotationDB/source metadata reconciliation. |
| `InChIKey` | character | Compound InChIKey used for metadata joining and no-CID retention. | Processed compound metadata `inchikey`. |
| `In.AnnotationDB` | logical | Whether a PubChem CID was mapped through AnnotationDB/source reconciliation. | Non-missing `Pubchem.CID` in processed compound metadata. |
| `AnnotationDB.Name` | character | AnnotationDB compound name. | Processed compound metadata `annotationdb_name`. |
| `AnnotationDB.SMILES` | character | AnnotationDB SMILES string. | Processed compound metadata `annotationdb_smiles`. |
| `LINCS.MOA` | character | LINCS mechanism-of-action annotation. | Source compound metadata `moa`, written as `lincs_moa` in processed metadata. |
| `LINCS.Targets` | character | LINCS target annotation. | Source compound metadata `target`, written as `lincs_targets` in processed metadata. |
| `LINCS.Aliases` | character | Source LINCS compound aliases. | Source compound metadata aliases, written as `lincs_aliases` in processed metadata. |
| `AnnotationDB.Aliases` | character | AnnotationDB aliases retained for review and matching context. | Processed AnnotationDB join field `annotationdb_aliases`. |
| `LINCS.Signature.Count` | integer | Number of LINCS signatures contributing to the compound-level assay column. | Processed compound metadata `n_signatures`, integer cast in MAE construction. |
| `LINCS.Gene.Count` | integer | Number of selected LINCS genes used for the compound signature. | Processed compound metadata `n_genes`, integer cast in MAE construction. |
| `LINCS.Feature.Count` | integer | Number of final assay features for the compound. | Processed compound metadata `n_features`, integer cast in MAE construction. |
