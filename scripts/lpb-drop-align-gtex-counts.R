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
# Supported input formats:
#   - Plain TSV / TXT with a gene-ID column named geneID / gene_id / Name / Gene
#     (with or without .gz compression)
#   - GTEx GCT (v1.2 and v1.3, .gct or .gct.gz). These ship a 2-line preamble
#     (#1.2/#1.3 + dimensions) and a "Description" column the script drops.
#
# Ensembl version handling: GTEx v11 uses GENCODE 39 (gene IDs like
# ENSG00000000003.14), while the GTF you align against (e.g., GENCODE 47)
# typically has different version suffixes. By default the script matches on
# the unversioned Ensembl ID (everything before the last "."), which gives
# very high overlap. Use --keep-versions to disable this and require exact
# versioned matches.
#
# Usage:
#   Rscript align_external_counts.R <gtf> <counts_in> <counts_out> [--keep-versions]
#
# Output:
#   - <counts_out>: counts file with rows in DROP's canonical order, using
#     "geneID" as the first column header (matching DROP's expectation).
#     Missing genes are zero-padded. The geneID column carries the GTF's
#     versioned IDs so the output is interchangeable with DROP's own output.

suppressPackageStartupMessages({
    library(data.table)
    library(txdbmaker)
    library(GenomicFeatures)
})

# ---- argument parsing ----
args <- commandArgs(trailingOnly = TRUE)
keep_versions <- "--keep-versions" %in% args
args <- args[args != "--keep-versions"]

if (length(args) != 3) {
    cat("Usage: Rscript align_external_counts.R <gtf> <counts_in> <counts_out> [--keep-versions]\n")
    cat("\nExamples:\n")
    cat("  # Plain TSV input:\n")
    cat("  Rscript align_external_counts.R \\\n")
    cat("    /path/to/gencode.v47.GRCh38.annotation.gtf \\\n")
    cat("    /path/to/gene_reads_v11_brain_frontal_cortex_ba9.txt \\\n")
    cat("    /path/to/gene_reads_v11_brain_frontal_cortex_ba9.aligned.tsv\n")
    cat("\n  # GTEx .gct.gz input (auto-detected):\n")
    cat("  Rscript align_external_counts.R \\\n")
    cat("    /path/to/gencode.v47.GRCh38.annotation.gtf \\\n")
    cat("    /path/to/gene_reads_adult_gtex_v11_brain_cortex.gct.gz \\\n")
    cat("    /path/to/brain_cortex_aligned.tsv\n")
    quit(status = 1)
}

GTF_PATH        <- args[1]
COUNTS_IN_PATH  <- args[2]
COUNTS_OUT_PATH <- args[3]

# ---- helper: strip Ensembl version suffix ----
# "ENSG00000000003.14"     -> "ENSG00000000003"
# "ENSG00000000003.14_PAR" -> "ENSG00000000003" (PAR_Y duplicates are stripped too)
# anything without a "."   -> returned unchanged
strip_version <- function(ids) {
    sub("\\..*$", "", ids)
}

# ---- helper: detect GCT format ----
# GCT v1.2 and v1.3 both start with a "#1.x" magic line. We check the first
# non-blank line of the file (handles .gz transparently via gzfile()).
is_gct_file <- function(path) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, "rt")
    } else {
        file(path, "rt")
    }
    on.exit(close(con))
    first_line <- ""
    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
        if (nzchar(trimws(line))) {
            first_line <- line
            break
        }
    }
    grepl("^#1\\.[0-9]", first_line)
}

# ---- helper: read a GCT file into a data.table ----
# GCT v1.2 layout:
#   line 1: #1.2
#   line 2: <n_rows> <TAB> <n_cols>
#   line 3: Name <TAB> Description <TAB> SAMPLE1 <TAB> SAMPLE2 ...
#   data rows below
#
# GCT v1.3 adds row/col metadata-counts to line 2 and may have extra
# metadata columns AFTER the ID column before the sample columns:
#   line 1: #1.3
#   line 2: <n_rows> <TAB> <n_cols> <TAB> <n_row_meta> <TAB> <n_col_meta>
#   line 3: id <TAB> Description <TAB> [extra row-meta cols] <TAB> SAMPLE1 ...
#
# We use line 2 to find n_cols (the number of sample columns), then take the
# LAST n_cols columns of the data as the samples. The first column is always
# the gene ID; everything between it and the sample block (Description, extra
# metadata) is discarded.
read_gct <- function(path) {
    open_con <- function() {
        if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt")
        else file(path, "rt")
    }

    # Read the two preamble lines from a fresh connection
    con <- open_con()
    version_line <- readLines(con, n = 1, warn = FALSE)
    dims_line    <- readLines(con, n = 1, warn = FALSE)
    close(con)
    cat("  GCT version line: ", version_line, "\n", sep = "")
    cat("  GCT dimensions  : ", dims_line, "\n", sep = "")

    dims <- as.integer(strsplit(trimws(dims_line), "\t")[[1]])
    if (length(dims) < 2 || any(is.na(dims[1:2]))) {
        stop("Could not parse GCT dimensions line: '", dims_line, "'")
    }
    n_rows <- dims[1]
    n_cols <- dims[2]

    # fread can skip the 2 preamble lines, but with .gz it needs R.utils
    # installed. To avoid that soft dependency, we instead read the file
    # through a gzfile() connection ourselves: read all lines, drop the
    # first 2, and feed the rest to fread via text=.
    con <- open_con()
    all_lines <- readLines(con, warn = FALSE)
    close(con)
    if (length(all_lines) < 3) {
        stop("GCT file has fewer than 3 lines after open: ", path)
    }
    dt <- fread(text = all_lines[-c(1, 2)], sep = "\t", header = TRUE)

    if (nrow(dt) != n_rows) {
        warning(sprintf(
            "GCT row-count mismatch: header says %d rows, parsed %d. Continuing with the parsed value.",
            n_rows, nrow(dt)))
    }

    # Last n_cols columns are samples; first column is the gene ID; anything
    # in between (Description, extra row metadata in v1.3) is discarded.
    n_total <- ncol(dt)
    if (n_total < n_cols + 1) {
        stop(sprintf("GCT has %d columns total but header says %d samples + at least 1 ID column expected",
                     n_total, n_cols))
    }
    id_colname      <- colnames(dt)[1]
    sample_colnames <- colnames(dt)[(n_total - n_cols + 1):n_total]
    discarded_meta  <- setdiff(colnames(dt), c(id_colname, sample_colnames))
    if (length(discarded_meta) > 0) {
        cat("  GCT metadata columns dropped: ",
            paste(discarded_meta, collapse = ", "), "\n", sep = "")
    }

    out <- dt[, c(id_colname, sample_colnames), with = FALSE]
    setnames(out, id_colname, "geneID")
    out
}

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
if (is_gct_file(COUNTS_IN_PATH)) {
    cat("  format: GCT (auto-detected)\n")
    counts <- read_gct(COUNTS_IN_PATH)
} else {
    cat("  format: plain TSV\n")
    counts <- fread(COUNTS_IN_PATH)

    id_col <- NULL
    for (candidate in c("geneID", "gene_id", "Name", "Gene", "id")) {
        if (candidate %in% colnames(counts)) {
            id_col <- candidate
            break
        }
    }
    if (is.null(id_col)) {
        stop("Could not find a gene-ID column. Expected one of: geneID, gene_id, Name, Gene, id.\n",
             "Columns present: ", paste(head(colnames(counts), 5), collapse = ", "), "...")
    }
    cat("  gene-ID column detected: '", id_col, "'\n", sep = "")
    setnames(counts, id_col, "geneID")
}

sample_cols <- setdiff(colnames(counts), "geneID")
cat("  ", nrow(counts), " genes  x  ", length(sample_cols), " samples\n", sep = "")

# ---- Step 3: align (optionally on unversioned IDs) ----
if (keep_versions) {
    cat("\n--keep-versions: matching on full Ensembl IDs (including version)\n")
    counts_lookup_ids <- counts$geneID
    target_lookup_ids <- canonical_order
} else {
    cat("\nMatching on unversioned Ensembl IDs (strip .NN suffix)\n")
    cat("  Use --keep-versions to require exact versioned matches.\n")
    counts_lookup_ids <- strip_version(counts$geneID)
    target_lookup_ids <- strip_version(canonical_order)

    pct_stripped_counts <- mean(counts_lookup_ids != counts$geneID) * 100
    pct_stripped_target <- mean(target_lookup_ids != canonical_order) * 100
    cat(sprintf("  counts file: %.1f%% of IDs had a version suffix\n",
                pct_stripped_counts))
    cat(sprintf("  GTF        : %.1f%% of IDs had a version suffix\n",
                pct_stripped_target))

    n_dup_counts <- sum(duplicated(counts_lookup_ids))
    n_dup_target <- sum(duplicated(target_lookup_ids))
    if (n_dup_counts > 0) {
        cat(sprintf("  WARNING: %d duplicate unversioned IDs in counts file (typically PAR_Y). Keeping first occurrence.\n",
                    n_dup_counts))
    }
    if (n_dup_target > 0) {
        cat(sprintf("  WARNING: %d duplicate unversioned IDs in GTF (typically PAR_Y). Affected genes will be zero-padded.\n",
                    n_dup_target))
    }
}

# Build a lookup keyed on the (possibly unversioned) ID, keeping first on dups.
counts_dt <- copy(counts)
counts_dt[, .lookup := counts_lookup_ids]
counts_dt <- counts_dt[!duplicated(.lookup)]
setkey(counts_dt, .lookup)

n_overlap <- sum(target_lookup_ids %in% counts_dt$.lookup)
n_missing <- length(target_lookup_ids) - n_overlap
n_extra   <- nrow(counts_dt) - n_overlap

cat("\nOverlap analysis:\n")
cat("  Genes in canonical & external: ", n_overlap, "\n", sep = "")
cat("  Genes in canonical only:       ", n_missing,
    " (will be zero-padded)\n", sep = "")
cat("  Genes in external only:        ", n_extra,
    " (will be dropped)\n", sep = "")

if (n_overlap < 0.5 * length(target_lookup_ids)) {
    cat("\n  WARNING: overlap is below 50%. ")
    if (keep_versions) {
        cat("This is often due to Ensembl version mismatch between the GTF and the external counts. ")
        cat("Try re-running without --keep-versions.\n")
    } else {
        cat("Even with versions stripped, overlap is low. ")
        cat("Check that the external counts file uses Ensembl IDs (ENSG...) rather than gene symbols.\n")
    }
}

# Reorder + zero-pad. The output keeps GTF-style versioned IDs in geneID.
cat("\nReordering and padding...\n")
result <- data.table(geneID = canonical_order)
match_idx <- counts_dt[J(target_lookup_ids), on = ".lookup", which = TRUE]

for (s in sample_cols) {
    v <- counts_dt[[s]][match_idx]
    v[is.na(v)] <- 0L
    result[[s]] <- v
}

stopifnot(identical(result$geneID, canonical_order))
cat("  result dimensions: ", nrow(result), " x ", ncol(result), "\n", sep = "")
cat("  order matches canonical: TRUE\n")

# ---- Step 4: write ----
cat("\nWriting aligned counts to:\n  ", COUNTS_OUT_PATH, "\n", sep = "")
fwrite(result, COUNTS_OUT_PATH, sep = "\t", quote = FALSE)

cat("\nDone.\n")
