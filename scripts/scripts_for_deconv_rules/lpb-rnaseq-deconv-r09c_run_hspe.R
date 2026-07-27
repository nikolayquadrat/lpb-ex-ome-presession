#!/usr/bin/env Rscript
# r09c_run_hspe.R -- deconvolve ONE bulk sample against ONE prepared reference.
# Uses the reference's precomputed markers so every sample of a reference is
# estimated on identical marker genes. Writes tidy per-sample proportions.
#
# Inputs (argv):
#   1 reference_rds   from r09b (reference, pure_samples, markers, genes)
#   2 mixture_tsv_gz  genes(symbol) x samples CPM (from r09a)
#   3 sample_name     column to deconvolve
#   4 out_props       proportions.tsv (class \t proportion)
#   5 out_rds         full hspe result (for provenance/debug)
suppressWarnings(suppressMessages({
  library(data.table)
  library(hspe)
}))

args <- commandArgs(trailingOnly = TRUE)
ref <- readRDS(args[1])
# Read the gzipped mixture via a shell pipe -- fread's direct .gz support needs
# the R.utils package, which isn't in the container; gzip -dc is always present.
mix <- fread(cmd = paste("gzip -dc", shQuote(args[2])))
sample_name <- args[3]
out_props <- args[4]; out_rds <- args[5]

if (!sample_name %in% colnames(mix))
  stop(sprintf("sample '%s' not in mixture columns", sample_name))

# align mixture to the reference gene axis (already the shared space from r09b)
rownames_mix <- mix[[1]]
y <- mix[[sample_name]]
names(y) <- rownames_mix
y <- y[ref$genes]                          # order/subset to reference genes
y[is.na(y)] <- 0
Y <- matrix(y, nrow = 1, dimnames = list(sample_name, ref$genes))

# hspe (dtangle family) expects LOG-scale expression, and the argument is
# `references` (plural). Both Y and references are log2(CPM + 1) on the SAME
# gene axis (ref$genes). Markers are chosen by hspe itself on this shared axis
# via n_markers -- more robust than passing precomputed indices, which must be
# column indices of Y and are easy to misalign.
Y_log   <- log2(Y + 1)
ref_log <- log2(ref$reference + 1)         # ref$reference is pseudobulk CPM

res <- hspe(Y = Y_log, references = ref_log,
            pure_samples = ref$pure_samples,
            n_markers = ref$n_markers)

est <- res$estimates                       # 1 x n_class
props <- data.frame(sample = sample_name,
                    class = colnames(est),
                    proportion = as.numeric(est[1, ]),
                    row.names = NULL)
props <- props[order(-props$proportion), ]
fwrite(props, out_props, sep = "\t")
saveRDS(res, out_rds)
message(sprintf("[r09c] %s x %d classes -> %s",
                sample_name, ncol(est), basename(out_props)))
