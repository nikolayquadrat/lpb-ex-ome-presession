#!/usr/bin/env Rscript
# r09b_reference.R -- turn a canonical single-cell reference into an hspe-ready
# reference object: pseudobulk profiles per (subject x class) restricted to the
# genes shared with the bulk mixture, plus hspe marker genes selected once per
# reference (so every sample of this reference uses identical markers).
#
# Inputs  (argv):
#   1 matrix.mtx.gz     genes x cells, integer counts (MatrixMarket)
#   2 features.tsv.gz   gene symbol per row (matches matrix rows)
#   3 cells.tsv.gz      cell_id \t subject_id \t source_label \t class
#   4 mixture_genes.txt one gene symbol per line (the bulk's symbol space)
#   5 n_markers         markers per class (e.g. 25)
#   6 out_rds           reference object for r09c
#   7 out_qc            marker-QC table (TSV)
suppressWarnings(suppressMessages({
  library(Matrix)
  library(hspe)
}))

args <- commandArgs(trailingOnly = TRUE)
mtx_f <- args[1]; feat_f <- args[2]; cells_f <- args[3]
mixgenes_f <- args[4]; n_markers <- as.integer(args[5])
out_rds <- args[6]; out_qc <- args[7]

# ---- load canonical reference ----
counts <- readMM(gzfile(mtx_f))                       # genes x cells
genes  <- readLines(gzfile(feat_f))
cells  <- read.delim(gzfile(cells_f), header = FALSE,
                     col.names = c("cell_id","subject_id","source_label","class"),
                     stringsAsFactors = FALSE)
stopifnot(nrow(counts) == length(genes), ncol(counts) == nrow(cells))
rownames(counts) <- genes

# collapse duplicate gene symbols (sum rows) so the gene axis is unique.
# Done sparsely: a (unique-gene x gene) 0/1 indicator matrix times counts sums
# rows sharing a symbol without ever densifying (as.matrix on a ~65k-cell
# reference would allocate ~19 GB and segfault).
if (any(duplicated(genes))) {
  ug  <- unique(genes)
  ind <- sparseMatrix(i = match(genes, ug), j = seq_along(genes), x = 1,
                      dims = c(length(ug), length(genes)))
  counts <- ind %*% counts                             # unique-genes x cells
  rownames(counts) <- ug
  genes  <- ug
}

# ---- restrict to genes shared with the bulk mixture ----
mix_genes <- readLines(mixgenes_f)
shared <- intersect(genes, mix_genes)
if (length(shared) < 200)
  stop(sprintf("only %d genes shared between reference and mixture -- gene-id ",
               length(shared)),
       "spaces likely mismatched (symbols vs Ensembl?).")
counts <- counts[shared, , drop = FALSE]

# ---- pseudobulk per (subject x class): sum counts, then CPM (linear) ----
# Sparse aggregation: (cells x group) indicator, so counts %*% ind sums each
# group's columns without densifying the full cells matrix.
grp     <- paste(cells$subject_id, cells$class, sep = "||")
grp_lvl <- unique(grp)
ind_c   <- sparseMatrix(i = seq_along(grp), j = match(grp, grp_lvl), x = 1,
                       dims = c(length(grp), length(grp_lvl)))
pb <- as.matrix(counts %*% ind_c)                     # genes x pseudobulk (small)
colnames(pb) <- grp_lvl
pb_class   <- sub(".*\\|\\|", "", colnames(pb))
# CPM (linear scale; hspe log-transforms internally)
libsize <- colSums(pb)
libsize[libsize == 0] <- 1
pb_cpm <- sweep(pb, 2, libsize, "/") * 1e6

# hspe wants samples x genes; pure_samples = per-class column indices
ref_mat <- t(pb_cpm)                                  # pseudobulk x genes
classes <- sort(unique(pb_class))
pure_samples <- lapply(classes, function(k) which(pb_class == k))
names(pure_samples) <- classes

# drop classes with no pseudobulk sample (shouldn't happen) and warn on thin ones
thin <- names(pure_samples)[vapply(pure_samples, length, 1L) < 2]
if (length(thin))
  message("[warn] classes with <2 pseudobulk samples (1 donor only): ",
          paste(thin, collapse = ", "))

# ---- marker QC on LOG scale (same transform hspe will use): confirm each
#      class's top ratio-markers are actually highest in that class. This is a
#      sanity check only -- hspe selects its own markers at run time via
#      n_markers, on the shared Y/references gene axis. ----
ref_log <- log2(ref_mat + 1)
mk <- find_markers(Y = ref_log, references = ref_log,
                   pure_samples = pure_samples, marker_method = "ratio")
qc <- do.call(rbind, lapply(classes, function(k) {
  mi <- head(mk$L[[k]], n_markers)
  if (length(mi) == 0) return(NULL)
  sub <- ref_log[, mi, drop = FALSE]
  class_mean <- colMeans(sub[pure_samples[[k]], , drop = FALSE])
  other_mean <- colMeans(sub[-pure_samples[[k]], , drop = FALSE])
  data.frame(class = k, n_markers = length(mi),
             frac_enriched = mean(class_mean > other_mean))
}))
write.table(qc, out_qc, sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(list(reference = ref_mat, pure_samples = pure_samples,
             n_markers = n_markers, classes = classes, genes = shared),
        out_rds)
message(sprintf("[r09b] reference ready: %d genes, %d classes, %d pseudobulk samples",
                length(shared), length(classes), nrow(ref_mat)))
