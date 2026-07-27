#!/usr/bin/env Rscript
# r09a_mixture.R -- assemble the bulk mixture in the reference's SYMBOL space.
# Reads every sample's STAR ReadsPerGene.out.tab, picks the strand column, maps
# versioned GENCODE gene_id -> gene symbol via the GTF, sums duplicate symbols,
# and writes a genes(symbol) x samples CPM matrix plus the gene list that r09b
# intersects against.
#
# Inputs (argv):
#   1 gtf              gencode.v47 GTF (for gene_id -> gene_name)
#   2 strand           0|1|2  STAR ReadsPerGene column: 1=unstranded(col2),
#                      2=stranded-fwd(col3), 3=reverse(col4). Pass the DATA index
#                      (we add 1 internally). GTEx/TruSeq = reverse -> use 2.
#   3 out_matrix       mixture_symbol_cpm.tsv.gz (genes x samples)
#   4 out_genes        mixture_genes.txt (one symbol per line)
#   5..  pairs of: <sample_name> <ReadsPerGene.out.tab path>
suppressWarnings(suppressMessages({ library(data.table) }))

args <- commandArgs(trailingOnly = TRUE)
gtf_f <- args[1]; strand <- as.integer(args[2])
out_mat <- args[3]; out_genes <- args[4]
rest <- args[-(1:4)]
stopifnot(length(rest) %% 2 == 0)
samp_names <- rest[c(TRUE, FALSE)]
samp_paths <- rest[c(FALSE, TRUE)]
col_data <- strand + 1L                 # DATA column index within the 4-col tab

# ---- gene_id(versioned) -> symbol from GTF ----
g <- fread(cmd = sprintf("grep -P '\\tgene\\t' %s", shQuote(gtf_f)),
           header = FALSE, sep = "\t")
attr9 <- g$V9
gid  <- sub('.*gene_id "([^"]+)".*',   "\\1", attr9)   # ENSG...N (versioned)
gsym <- sub('.*gene_name "([^"]+)".*', "\\1", attr9)
id2sym <- setNames(gsym, gid)

# ---- read each ReadsPerGene.out.tab ----
read_star <- function(path) {
  d <- fread(path, header = FALSE)
  d <- d[!V1 %in% c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous")]
  # STAR tab columns: V1 gene_id, V2 unstranded, V3 fwd, V4 reverse
  cnt <- d[[col_data + 1L]]              # +1: V1 is gene_id, data starts at V2
  setNames(cnt, d$V1)
}
mats <- lapply(samp_paths, read_star)
all_ids <- Reduce(union, lapply(mats, names))
M <- sapply(mats, function(x) x[all_ids]); M[is.na(M)] <- 0
rownames(M) <- all_ids; colnames(M) <- samp_names

# ---- map to symbol, sum duplicate symbols, drop unmapped ----
sym <- id2sym[rownames(M)]
keep <- !is.na(sym)
M <- M[keep, , drop = FALSE]; sym <- sym[keep]
M <- rowsum(M, group = sym)              # symbol x samples (summed)

# ---- CPM (linear) ----
libsize <- colSums(M); libsize[libsize == 0] <- 1
cpm <- sweep(M, 2, libsize, "/") * 1e6

out <- data.frame(gene = rownames(cpm), cpm, check.names = FALSE)
# Write plain TSV then gzip in the shell -- avoids depending on whether this
# container's data.table was compiled with zlib (fwrite compress="gzip" fails
# otherwise). gzip is present in every image.
out_plain <- sub("\\.gz$", "", out_mat)
fwrite(out, out_plain, sep = "\t")
system2("gzip", c("-f", shQuote(out_plain)))          # -> out_mat (out_plain.gz)
writeLines(rownames(cpm), out_genes)
message(sprintf("[r09a] mixture: %d symbols x %d samples (strand col %d)",
                nrow(cpm), ncol(cpm), col_data))