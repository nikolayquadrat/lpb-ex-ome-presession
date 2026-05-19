#!/usr/bin/env Rscript

# Align a GTEx (or any external) counts file to the gene order DROP will use,
# WITHOUT requiring DROP to have run preprocessing first.
#
# DROP builds count_ranges.Rds from the GTF via:
#   txdb <- makeTxDbFromGFF(gtf)
#   exonsBy(txdb, by = "gene")
# which produces a GRangesList with genes in alphabetical order by gene_id.
#
# This script replicates that ordering directly from the GTF, then aligns
# the external counts file so it matches.
#
# Usage:
#   Rscript align_external_counts.R <gtf> <counts_in> <counts_out>
#
# Output:
#   - <counts_out>: counts file with rows in DROP's canonical order,
#     using "geneID" as the first column header (matching DROP's expectation).
#     Missing genes are zero-padded.

suppressPackageStartupMessages({
    library(data.table)
    library(txdbmaker)
    library(GenomicFeatures)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    cat("Usage: Rscript align_external_counts.R <gtf> <counts_in> <counts_out>\n")
    cat("\nExample:\n")
    cat("  Rscript align_external_counts.R \\\n")
    cat("    /path/to/gencode.v47.GRCh38.annotation.gtf \\\n")
    cat("    /path/to/gene_reads_v11_brain_frontal_cortex_ba9.txt \\\n")
    cat("    /path/to/gene_reads_v11_brain_frontal_cortex_ba9.aligned.txt\n")
    quit(status = 1)
}

GTF_PATH       <- args[1]
COUNTS_IN_PATH <- args[2]
COUNTS_OUT_PATH <- args[3]

# ---- Step 1: derive the canonical gene order from the GTF ----
cat("Building TxDb from GTF:\n  ", GTF_PATH, "\n", sep = "")
txdb <- suppressWarnings(makeTxDbFromGFF(GTF_PATH, format = "gtf"))

cat("Extracting exonsBy(gene)...\n")
exons_by_gene <- exonsBy(txdb, by = "gene")
canonical_order <- names(exons_by_gene)
cat("  ", length(canonical_order), " genes in canonical (DROP) order\n", sep = "")
cat("  first 3: ", paste(head(canonical_order, 3), collapse = ", "), "\n", sep = "")

# ---- Step 2: read the external counts file ----
cat("\nReading external counts:\n  ", COUNTS_IN_PATH, "\n", sep = "")
counts <- fread(COUNTS_IN_PATH)

# Determine the gene-ID column name (could be "gene_id" or "geneID")
id_col <- NULL
for (candidate in c("geneID", "gene_id", "Name", "Gene")) {
    if (candidate %in% colnames(counts)) {
        id_col <- candidate
        break
    }
}
if (is.null(id_col)) {
    stop("Could not find a gene-ID column. Expected one of: geneID, gene_id, Name, Gene.\n",
         "Columns present: ", paste(head(colnames(counts), 5), collapse = ", "), "...")
}
cat("  gene-ID column detected: '", id_col, "'\n", sep = "")

setnames(counts, id_col, "geneID")  # normalize to geneID (DROP's expected name)
sample_cols <- setdiff(colnames(counts), "geneID")
cat("  ", nrow(counts), " genes  x  ", length(sample_cols), " samples\n", sep = "")

# ---- Step 3: align ----
counts_gene_set <- counts$geneID
missing_from_external <- setdiff(canonical_order, counts_gene_set)
extra_in_external     <- setdiff(counts_gene_set, canonical_order)

cat("\nOverlap analysis:\n")
cat("  Genes in canonical & external: ", length(intersect(canonical_order, counts_gene_set)), "\n", sep = "")
cat("  Genes in canonical only:       ", length(missing_from_external),
    " (will be zero-padded)\n", sep = "")
cat("  Genes in external only:        ", length(extra_in_external),
    " (will be dropped)\n", sep = "")

if (length(extra_in_external) > 0) {
    cat("    first 3 dropped: ", paste(head(extra_in_external, 3), collapse = ", "), "\n", sep = "")
    cat("    (these are genes the external file has but the GTF does not -- they cannot be used)\n")
}

# Build the aligned counts data.table
cat("\nReordering and padding...\n")
setkey(counts, geneID)
result <- data.table(geneID = canonical_order)
for (s in sample_cols) {
    # Lookup; for genes not in counts, this returns NA
    v <- counts[result$geneID, ..s][[1]]
    v[is.na(v)] <- 0L
    result[[s]] <- v
}

# Verify
stopifnot(identical(result$geneID, canonical_order))
cat("  result dimensions: ", nrow(result), " x ", ncol(result), "\n", sep = "")
cat("  order matches canonical: TRUE\n")

# ---- Step 4: write ----
cat("\nWriting aligned counts to:\n  ", COUNTS_OUT_PATH, "\n", sep = "")
fwrite(result, COUNTS_OUT_PATH, sep = "\t", quote = FALSE)

cat("\nDone.\n")
