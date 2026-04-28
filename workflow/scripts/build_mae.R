suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(SummarizedExperiment)
  library(yaml)
})

set_vector_memory_limit <- function() {
  if (!exists("mem.maxVSize", mode = "function")) {
    return(invisible(NULL))
  }
  target_gb <- Sys.getenv("PIPELINE_R_MAX_VSIZE_GB", "64")
  target_mb <- suppressWarnings(as.numeric(target_gb) * 1024)
  if (tolower(target_gb) %in% c("inf", "infinite")) {
    target_mb <- Inf
  }
  if (is.na(target_mb)) {
    return(invisible(NULL))
  }
  current_mb <- mem.maxVSize()
  if (
    !is.finite(current_mb) || is.infinite(target_mb) || current_mb < target_mb
  ) {
    try(mem.maxVSize(target_mb), silent = TRUE)
  }
  invisible(NULL)
}

set_vector_memory_limit()

stopf <- function(message, ...) {
  stop(sprintf(message, ...), call. = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(
    path,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    na.strings = c("", "NA", "nan", "NaN"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

read_csv <- function(path, row.names = NULL) {
  utils::read.csv(
    path,
    row.names = row.names,
    check.names = FALSE,
    na.strings = c("", "NA", "nan", "NaN"),
    stringsAsFactors = FALSE
  )
}

normalize_blank_to_na <- function(values) {
  values <- as.character(values)
  blank <- !is.na(values) & !nzchar(trimws(values))
  values[blank] <- NA_character_
  values
}

first_non_missing <- function(...) {
  values <- list(...)
  out <- rep(NA_character_, length(values[[1L]]))
  for (value in values) {
    value <- normalize_blank_to_na(value)
    replace <- is.na(out) & !is.na(value)
    out[replace] <- value[replace]
  }
  out
}

as_integer_or_na <- function(values) {
  suppressWarnings(as.integer(as.character(values)))
}

build_public_drug_metadata <- function(compound_metadata) {
  pubchem_cid <- as_integer_or_na(compound_metadata$cid)
  lincs_key <- as.character(compound_metadata$cmap_name)
  if (any(is.na(lincs_key) | !nzchar(trimws(lincs_key)))) {
    stopf("LINCS cmap_name values must be non-missing before MAE construction")
  }
  if (anyDuplicated(lincs_key) > 0L) {
    stopf("LINCS cmap_name values must be unique before MAE construction")
  }
  public <- data.frame(
    LINCS.CMap.Name = lincs_key,
    Pubchem.CID = pubchem_cid,
    InChIKey = as.character(compound_metadata$inchikey),
    In.AnnotationDB = !is.na(pubchem_cid),
    AnnotationDB.Name = as.character(compound_metadata$annotationdb_name),
    AnnotationDB.SMILES = as.character(compound_metadata$annotationdb_smiles),
    LINCS.MOA = as.character(compound_metadata$lincs_moa),
    LINCS.Targets = as.character(compound_metadata$lincs_targets),
    LINCS.Aliases = as.character(compound_metadata$lincs_aliases),
    AnnotationDB.Aliases = as.character(compound_metadata$annotationdb_aliases),
    LINCS.Signature.Count = as_integer_or_na(compound_metadata$n_signatures),
    LINCS.Gene.Count = as_integer_or_na(compound_metadata$n_genes),
    LINCS.Feature.Count = as_integer_or_na(compound_metadata$n_features),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  row.names(public) <- public$LINCS.CMap.Name
  public
}

compound_metadata <- read_tsv(snakemake@input[["compound_metadata_tsv"]])
drug_metadata <- build_public_drug_metadata(compound_metadata)

signatures <- read_csv(snakemake@input[["signatures_csv"]], row.names = 1)
missing_signature_columns <- setdiff(
  drug_metadata$LINCS.CMap.Name,
  colnames(signatures)
)
if (length(missing_signature_columns) > 0L) {
  stopf(
    "Signature matrix is missing %d compound columns",
    length(missing_signature_columns)
  )
}
signatures <- signatures[, drug_metadata$LINCS.CMap.Name, drop = FALSE]
assay_matrix <- as.matrix(signatures)
storage.mode(assay_matrix) <- "numeric"

feature_metadata <- read_tsv(snakemake@input[["feature_metadata_tsv"]])
row_data <- data.frame(
  Feature.ID = rownames(assay_matrix),
  feature_metadata[
    match(rownames(assay_matrix), feature_metadata$feature_name),
  ],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
names(row_data) <- sub("^feature_", "LINCS.Feature.", names(row_data))
row.names(row_data) <- row_data$Feature.ID

signature_experiment <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(values = assay_matrix),
  rowData = S4Vectors::DataFrame(row_data, row.names = row.names(row_data)),
  colData = S4Vectors::DataFrame(row.names = colnames(assay_matrix))
)

sample_map <- data.frame(
  assay = "signatures",
  primary = colnames(assay_matrix),
  colname = colnames(assay_matrix),
  stringsAsFactors = FALSE
)

col_data <- drug_metadata
row.names(col_data) <- col_data$LINCS.CMap.Name

mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = MultiAssayExperiment::ExperimentList(
    signatures = signature_experiment
  ),
  colData = S4Vectors::DataFrame(col_data, row.names = row.names(col_data)),
  sampleMap = S4Vectors::DataFrame(sample_map)
)

metadata(mae) <- list(
  Pipeline = list(
    ID = snakemake@params[["dataset_id"]],
    Version = snakemake@params[["dataset_version"]],
    Config = yaml::read_yaml(
      snakemake@input[["configfile"]],
      eval.expr = FALSE
    )
  ),
  Selected.Gene.Metadata = read_tsv(snakemake@input[["gene_metadata_tsv"]]),
  Drug.Metadata = S4Vectors::DataFrame(
    drug_metadata,
    row.names = drug_metadata$LINCS.CMap.Name
  )
)

output_path <- snakemake@output[["mae_rds"]]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(mae, output_path)
