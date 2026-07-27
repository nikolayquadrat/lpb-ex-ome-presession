#!/usr/bin/env Rscript
# r09e_cibersortx_prep.R -- write the two files a manual CIBERSORTx run needs,
# for ONE reference, into 09_deconv/_cibersortx/{ref}/:
#   refsample.txt : genes x cells, tab-delimited, COLUMN HEADERS = each cell's
#                   cell-type label (CIBERSORTx derives its signature from this),
#                   first column = gene symbol. Cells capped per class.
#   mixture.txt   : genes x samples, first column = gene symbol (the r09a bulk
#                   CPM matrix, just reformatted). Shared across references but
#                   copied per ref so each folder is a self-contained run dir.
#
# The pipeline only PREPARES these; you run CIBERSORTx yourself, e.g.:
#   docker run --rm -v $PWD:/src/data -v $PWD/out:/src/outdir cibersortx/fractions \
#     --username <you> --token <token> --single_cell TRUE \
#     --refsample refsample.txt --mixture mixture.txt --fraction 0 \
#     --rmbatchSmode TRUE --QN FALSE
# (S-mode batch correction is recommended: it targets the snRNA-ref -> bulk-mix
#  platform gap that biases marker-based methods.)
#
# argv: 1 matrix.mtx.gz  2 features.tsv.gz  3 cells.tsv.gz
#       4 mixture_symbol_cpm.tsv.gz (from r09a)
#       5 out_refsample.txt  6 out_mixture.txt  7 max_cells_per_class
suppressWarnings(suppressMessages({ library(Matrix); library(data.table) }))
args <- commandArgs(trailingOnly = TRUE)
mtx_f <- args[1]; feat_f <- args[2]; cells_f <- args[3]; mix_f <- args[4]
out_ref <- args[5]; out_mix <- args[6]; cap <- as.integer(args[7])

# ---- refsample from canonical reference ----
counts <- readMM(gzfile(mtx_f))                 # genes x cells
genes  <- readLines(gzfile(feat_f))
cells  <- read.delim(gzfile(cells_f), header = FALSE,
                     col.names = c("cell_id","subject_id","source_label","class"),
                     stringsAsFactors = FALSE)
rownames(counts) <- genes

# collapse duplicate gene symbols sparsely (never densify the full matrix)
if (any(duplicated(genes))) {
  ug  <- unique(genes)
  ind <- sparseMatrix(i = match(genes, ug), j = seq_along(genes), x = 1,
                      dims = c(length(ug), length(genes)))
  counts <- ind %*% counts; rownames(counts) <- ug; genes <- ug
}

# cap cells per class (reproducible)
set.seed(42)
keep <- unlist(lapply(split(seq_len(nrow(cells)), cells$class), function(idx)
  if (length(idx) > cap) sample(idx, cap) else idx), use.names = FALSE)
keep <- sort(keep)
counts <- counts[, keep, drop = FALSE]
labels <- cells$class[keep]

dt <- as.data.table(as.matrix(counts))          # genes x (capped) cells -- small
setnames(dt, labels)                            # headers = cell-type labels
dt <- cbind(GeneSymbol = genes, dt)
fwrite(dt, out_ref, sep = "\t")
message(sprintf("[r09e] refsample: %d genes x %d cells (%d classes) -> %s",
                length(genes), ncol(counts), length(unique(labels)), basename(out_ref)))

# ---- mixture: reformat the r09a symbol-space CPM matrix ----
mix <- fread(cmd = paste("gzip -dc", shQuote(mix_f)))
setnames(mix, 1, "GeneSymbol")
fwrite(mix, out_mix, sep = "\t")
message(sprintf("[r09e] mixture: %d genes x %d samples -> %s",
                nrow(mix), ncol(mix) - 1, basename(out_mix)))
