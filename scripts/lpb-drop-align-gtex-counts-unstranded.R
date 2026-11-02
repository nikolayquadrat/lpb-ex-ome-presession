#!/usr/bin/env Rscript
# Base R only for TSVs; reading DROP's count_ranges.Rds may require its
# Bioconductor packages. Run in the DROP R environment.
# Optional arguments: annotation.tsv local.tsv gtex.tsv count_ranges.Rds output_dir

project <- "/mnt/data/exome/RNASEQ/drop/ba9_gtex_nov_counts"
annotation_file <- file.path(project, "Data/sample_annotation_ba9_gtex.tsv")
local_file <- "/mnt/data/exome/RNASEQ/rnaseq-drop/10_rnaseqc/00_drop_gene_counts.unstranded.tsv"
external_file <- file.path(project, "Data/gene_reads_v11_brain_frontal_cortex_ba9.txt")
ranges_file <- file.path(project, "Output/processed_data/aberrant_expression/v47_rnaseqc/count_ranges.Rds")
output_dir <- file.path(project, "Data/combined_counts")

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
    if (length(args) != 5L) stop("Supply exactly five arguments: annotation local external ranges output_dir")
    annotation_file <- args[1]; local_file <- args[2]; external_file <- args[3]
    ranges_file <- args[4]; output_dir <- args[5]
}

read_tsv <- function(path, ...) {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
               quote = "", comment.char = "", ...)
}
anno <- read_tsv(annotation_file, colClasses = "character", na.strings = c("", "NA"))
required <- c("RNA_ID", "batch", "GENE_COUNTS_FILE", "GENE_ANNOTATION", "DROP_GROUP")
if (!all(required %in% names(anno))) stop("Sample annotation lacks required columns")
if (!nrow(anno) || anyNA(anno$RNA_ID) || anyDuplicated(anno$RNA_ID)) {
    stop("Sample annotation must contain unique, nonmissing RNA_IDs")
}
if (anyNA(anno$batch)) stop("Missing batch values: cannot distinguish GTEx and local samples")
if (anyNA(anno$GENE_ANNOTATION) || length(unique(anno$GENE_ANNOTATION)) != 1L) {
    stop("This script expects one shared GENE_ANNOTATION")
}
local_ids <- anno$RNA_ID[anno$batch != "GTEx"]
external_ids <- anno$RNA_ID[anno$batch == "GTEx"]
if (!length(local_ids) || !length(external_ids)) stop("Both local and GTEx samples are required")

# Supports ordinary TSVs and GCT 1.2, with geneID/Name in the FIRST column.
# Select samples before loading values; extra samples and Description are ignored.
read_counts <- function(path, ids, label) {
    first <- readLines(path, n = 1L, warn = FALSE)
    skip <- if (identical(first, "#1.2")) 2L else 0L
    header <- names(read_tsv(path, skip = skip, nrows = 0L))
    if (!length(header) || !header[1] %in% c("geneID", "gene_id", "Name", "Geneid")) {
        stop(label, ": first column must be geneID, gene_id, Name or Geneid")
    }
    if (anyDuplicated(header)) stop(label, ": duplicate column names")
    missing <- setdiff(ids, header[-1])
    if (length(missing)) stop(label, ": required samples missing:\n", paste(missing, collapse = "\n"))
    classes <- rep("NULL", length(header))
    classes[1] <- "character"
    classes[match(ids, header)] <- "numeric"
    d <- read_tsv(path, skip = skip, colClasses = classes, na.strings = c("NA", ""))
    genes <- d[[1]]
    if (!length(genes) || anyNA(genes) || any(genes == "") || anyDuplicated(genes)) {
        stop(label, ": missing or duplicate gene IDs")
    }
    x <- as.matrix(d[, ids, drop = FALSE])
    if (any(!is.finite(x)) || any(x < 0) || any(x != floor(x))) {
        stop(label, ": counts must be finite, nonnegative integers (not TPMs)")
    }
    rownames(x) <- genes
    message(label, ": ", nrow(x), " genes, ", ncol(x), " selected samples")
    x
}
local <- read_counts(local_file, local_ids, "Local matrix")
external <- read_counts(external_file, external_ids, "GTEx matrix")

# DROP requires exactly these gene IDs AND this order. Do not strip Ensembl
# versions, silently take an intersection, or fabricate missing genes as zeros.
suppressPackageStartupMessages(library(GenomicRanges))
ranges <- readRDS(ranges_file)
genes <- names(ranges)

if (is.null(genes) || !length(genes) || anyNA(genes) || any(genes == "") || anyDuplicated(genes)) {
    stop("count_ranges.Rds must have unique, nonmissing gene names")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
write_tsv <- function(x, name) {
    write.table(x, file.path(output_dir, name), sep = "\t", quote = FALSE,
                row.names = FALSE, na = "NA")
}
missing <- data.frame(geneID = genes,
                      missing_local = !genes %in% rownames(local),
                      missing_GTEx = !genes %in% rownames(external))
missing <- missing[missing$missing_local | missing$missing_GTEx, , drop = FALSE]
write_tsv(missing, "missing_annotation_genes.tsv")
if (nrow(missing)) stop(nrow(missing), " annotation genes missing; see missing_annotation_genes.tsv. No combined matrix written.")
extra_local <- setdiff(rownames(local), genes)
extra_external <- setdiff(rownames(external), genes)
write_tsv(rbind(data.frame(source = rep("local", length(extra_local)), geneID = extra_local),
                data.frame(source = rep("GTEx", length(extra_external)), geneID = extra_external)),
          "excluded_nonannotation_genes.tsv")
local <- local[genes, , drop = FALSE]
external <- external[genes, , drop = FALSE]

# SD is on RAW counts and describes within-cohort variation, not expression
# differences after library-size adjustment. Constant nonzero genes also flag.
stats <- function(x) {
    data.frame(n = ncol(x), mean = rowMeans(x),
               sd = if (ncol(x) > 1L) apply(x, 1L, sd) else rep(NA_real_, nrow(x)),
               n_nonzero = rowSums(x > 0))
}
ls <- stats(local); es <- stats(external)
qc <- data.frame(geneID = genes, local_n = ls$n, GTEx_n = es$n,
                 local_mean = ls$mean, GTEx_mean = es$mean,
                 local_sd = ls$sd, GTEx_sd = es$sd,
                 local_n_nonzero = ls$n_nonzero, GTEx_n_nonzero = es$n_nonzero)
qc$local_all_zero <- qc$local_n_nonzero == 0
qc$GTEx_all_zero <- qc$GTEx_n_nonzero == 0
qc$local_zero_sd <- !is.na(qc$local_sd) & qc$local_sd == 0
qc$GTEx_zero_sd <- !is.na(qc$GTEx_sd) & qc$GTEx_sd == 0
qc$local_zero_GTEx_variable <- qc$local_all_zero & !is.na(qc$GTEx_sd) & qc$GTEx_sd > 0
qc$GTEx_zero_local_variable <- qc$GTEx_all_zero & !is.na(qc$local_sd) & qc$local_sd > 0
write_tsv(qc, "gene_qc_all.tsv")
flag <- qc$local_zero_sd | qc$GTEx_zero_sd | qc$local_all_zero | qc$GTEx_all_zero
write_tsv(qc[flag, , drop = FALSE], "gene_qc_flagged.tsv")

combined <- cbind(local, external)[, anno$RNA_ID, drop = FALSE]
stopifnot(identical(rownames(combined), genes), identical(colnames(combined), anno$RNA_ID))
write_tsv(data.frame(geneID = genes, combined, check.names = FALSE), "combined_counts.tsv")
anno$GENE_COUNTS_FILE <- file.path(output_dir, "combined_counts.tsv")
write_tsv(anno, "sample_annotation_combined.tsv")
message("Written ", nrow(combined), " genes x ", ncol(combined), " samples to ", output_dir)
message(sum(flag), " genes flagged; none removed. Local all-zero / GTEx variable: ",
        sum(qc$local_zero_GTEx_variable), "; reverse: ", sum(qc$GTEx_zero_local_variable))
