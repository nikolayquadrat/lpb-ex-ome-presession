rm(list = ls())

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(purrr)
library(readxl)
library(writexl)
library(msigdbr) # v7.5.1
library(edgeR) # 4.2.2
library(limma) # 3.60.6
library(GSVA)  # 1.52.3

# Arguments ============
git_folder <- "C:/Users/Nikolay/Dropbox/Git/lpb-ex-ome-presession" # <-- SET THE GIT PATH
exome_name <- "e1-19"

# ENSG --> HGNC dictionary ==========
constraint <- fread(sprintf("%s/data/exome-pipe/data/gnomad_v4.1_constraint_metrics.tsv", git_folder))
constraint <- constraint[mane_select == TRUE]  # one row per gene
ensg_to_hgnc <- setNames(constraint$gene, sub("\\..*$", "", constraint$gene_id))

# Candidate mutations ======
tier_mutations <- read.delim(sprintf("%s/data/exome-pipe/fastq/12_tiered/%s-combined_master.tsv", git_folder, exome_name)) %>% 
    filter(tier %in% c("A","B","C"))
tier_genes <- names(table(tier_mutations$SYMBOL)) # here used only got visualisation

# GSEA signatures ===========
hallmark <- msigdbr(species="Homo sapiens", category="H")
reactome <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
biocarta <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:BIOCARTA")
kegg     <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:KEGG")
go_bp    <- msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP")
go_cc    <- msigdbr(species="Homo sapiens", category="C5", subcategory="GO:CC")
go_mf    <- msigdbr(species="Homo sapiens", category="C5", subcategory="GO:MF")
wiki     <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:WIKIPATHWAYS")
hpo      <- msigdbr(species="Homo sapiens", category="C5", subcategory="HPO")

pathway_list <- c(
    lapply(split(hallmark$gene_symbol, hallmark$gs_name), unique),
    lapply(split(reactome$gene_symbol, reactome$gs_name), unique),
    lapply(split(kegg$gene_symbol, kegg$gs_name), unique),
    lapply(split(biocarta$gene_symbol, biocarta$gs_name), unique),
    lapply(split(wiki$gene_symbol, wiki$gs_name), unique),
    lapply(split(go_bp$gene_symbol, go_bp$gs_name), unique),
    lapply(split(go_cc$gene_symbol, go_cc$gs_name), unique),
    lapply(split(go_mf$gene_symbol, go_mf$gs_name), unique),
    lapply(split(hpo$gene_symbol, hpo$gs_name), unique)
)

# Speificity from Siletti et al (2023) data =========
if (0) {
    cib_siletti_specificity <- read.delim(sprintf("%s/data/rnaseq-pipe/00_additional_files/deconv/reference_canonical/siletti_cortex/specificity.tsv", git_folder))
    spec_vec <- setNames(cib_siletti_specificity$delta_log2,
                         cib_siletti_specificity$gene)
    
    pathway_neuronal_score <- sapply(rownames(scores), function(y) {
        mean(spec_vec[pathway_list[[y]]], na.rm = TRUE)   # named lookup, vectorized
    })
    write.table(
        data.frame("pathway"=names(pathway_neuronal_score),
                   "neuronality"=pathway_neuronal_score),
        sprintf("%s/data/rnaseq-pipe/00_additional_files/deconv/reference_canonical/siletti_cortex/neuronality.tsv", git_folder),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
}
pathway_neuronal_score <- read.delim(sprintf("%s/data/rnaseq-pipe/00_additional_files/deconv/reference_canonical/siletti_cortex/neuronality.tsv", git_folder))

# RNA data =========
sample_data <- readxl::read_xlsx(sprintf("%s/data/supplementary_files/LPB_Supplementary_Tables.xlsx", git_folder), sheet = "S2. RNAseq samples", skip = 1)
counts_raw <- readxl::read_xlsx(sprintf("%s/data/rnaseq-pipe/04_qc/00_expression_summary.xlsx", git_folder), sheet = "counts") 

sym    <- as.character(counts_raw$gene_name)
counts <- as.matrix(counts_raw[, c(-1, -2)])
storage.mode(counts) <- "integer"          # or "double"; ensures numeric matrix
rownames(counts) <- NULL
counts <- counts[, colnames(counts) %in% sample_data$`sample name`]

# drop rows with no real HGNC symbol (unmapped -> kept as ENSG / empty)-
unmapped <- is.na(sym) | sym == "" | grepl("^ENSG[0-9]+", sym)
counts <- counts[!unmapped, , drop = FALSE]
sym    <- sym[!unmapped]
cat(sprintf("dropped %d unmapped rows; %d remain\n", sum(unmapped), length(sym)))

# collapse duplicate symbols by summing counts (one row per symbol) ----
counts <- rowsum(counts, group = sym)
cat(sprintf("after collapsing duplicates: %d unique genes\n", nrow(counts)))

# meta-data ========
coldata <- data.frame(
    row.names = colnames(counts),
    "condition"     = sapply(colnames(counts), function(x) {sample_data$schizophrenia[sample_data$`sample name` == x]}),
    "brain_id"      = sapply(colnames(counts), function(x) {sample_data$`brain ID`[sample_data$`sample name` == x]}),
    "brain_region"  = sapply(colnames(counts), function(x) {sample_data$`Brodmann area`[sample_data$`sample name` == x]}),
    "batch"         = sapply(colnames(counts), function(x) {sample_data$batch[sample_data$`sample name` == x]}),
    "sex"           = sapply(colnames(counts), function(x) {sample_data$sex[sample_data$`sample name` == x]}),
    "age"           = sapply(colnames(counts), function(x) {sample_data$`age at death`[sample_data$`sample name` == x]})
)
all(rownames(coldata) == colnames(counts))

# Helper functions ===========
source(sprintf("%s/scripts/lpb-post-drop-donor-specificity-table-gsva.R", git_folder))
source(sprintf("%s/scripts/lpb-post-drop-classify-pathways.R", git_folder))

# Main loop =========
dir.create(sprintf("%s/data/post-drop/trd-pass-gsva", git_folder), showWarnings = FALSE)
for (region in c("BA9", "BA22p", "BA4")) {
    # region <- "BA9"
    cat(region, "\n")
    
    coldata_reg <- coldata[coldata$brain_region == region, ]
    counts_reg  <- counts[, rownames(coldata_reg), drop = FALSE]
    stopifnot(identical(colnames(counts_reg), rownames(coldata_reg)))
    
    keep <- rowSums(counts_reg > 10) >= (ncol(counts_reg) / 2)
    counts_reg <- counts_reg[keep, , drop = FALSE]
    
    dge <- calcNormFactors(DGEList(counts_reg))
    logcpm <- cpm(dge, log = TRUE, prior.count = 1)
    
    # collapse within donor-region-batch
    grp1 <- paste(coldata_reg$brain_id, region, coldata_reg$batch, sep = "_")
    idx1 <- split(seq_len(ncol(logcpm)), grp1)
    logcpm_cb <- sapply(idx1, function(c) rowMeans(logcpm[, c, drop = FALSE]))
    coldata_cb <- data.frame(
        sample = names(idx1),
        brain_id = sapply(idx1, \(c) coldata_reg$brain_id[c][1]),
        batch = sapply(idx1, \(c) coldata_reg$batch[c][1]),
        condition = sapply(idx1, \(c) coldata_reg$condition[c][1]),
        sex = sapply(idx1, \(c) coldata_reg$sex[c][1]),
        age = sapply(idx1, \(c) coldata_reg$age[c][1]), row.names = NULL)
    stopifnot(identical(colnames(logcpm_cb), coldata_cb$sample))
    
    # batch-correct (skip if single batch), preserving condition/sex/age
    if (length(unique(coldata_cb$batch)) > 1) {
        design_cb <- model.matrix(~ condition + sex + age, data = coldata_cb)
        logcpm_cb_corr <- limma::removeBatchEffect(logcpm_cb, batch = coldata_cb$batch, design = design_cb)
    } else {
        logcpm_cb_corr <- logcpm_cb          # BA4: nothing to remove
    }
    
    # collapse over batches (donor-region)
    grp2 <- paste(coldata_cb$brain_id, region, sep = "_")
    idx2 <- split(seq_len(ncol(logcpm_cb_corr)), grp2)
    logcpm_collapsed <- sapply(idx2, function(c) rowMeans(logcpm_cb_corr[, c, drop = FALSE]))
    coldata_collapsed <- data.frame(
        sample = names(idx2),
        brain_id = sapply(idx2, \(c) coldata_cb$brain_id[c][1]),
        condition = sapply(idx2, \(c) coldata_cb$condition[c][1]),
        sex = sapply(idx2, \(c) coldata_cb$sex[c][1]),
        age = sapply(idx2, \(c) coldata_cb$age[c][1]), row.names = NULL)
    stopifnot(identical(colnames(logcpm_collapsed), coldata_collapsed$sample))
    
    # GSVA on the corrected log-matrix, Gaussian KCDF, SAME pathways/sizes as fgsea
    par <- gsvaParam(exprData = logcpm_collapsed,
                     geneSets = pathway_list,
                     minSize = 20, maxSize = 200,     # matches fgsea on OUTRIDER results
                     kcdf = "Gaussian")               # correct for log-continuous data
    scores <- gsva(par)                              # pathways x samples
    # SZ07's rank among all donors, per pathway (1 = lowest, 15 = highest)
    rank_sz07 <- apply(scores, 1, function(x) rank(x)[sprintf("SZ07_%s", region)])
    # pathways where SZ07 is the most extreme (rank 15 = top of all 15 donors)
    highest_rank_scores <- names(rank_sz07)[rank_sz07 == 15]
    
    # gsva_scores: pathways x donors (your 6847 x 15 matrix)
    others <- setdiff(colnames(scores), sprintf("SZ07_%s", region))

    gsva_res <- data.frame(
        pathway = rownames(scores),
        sz07_score = scores[, sprintf("SZ07_%s", region)],  # insted of NES                                
        rank_sz07  = apply(scores, 1, function(x) rank(x)[sprintf("SZ07_%s", region)]), # instead of padj
        robust_z   = (scores[,sprintf("SZ07_%s", region)] - apply(scores[,others],1,median)) /
            apply(scores[,others],1,mad),  
        donors_gt_sz07 = apply(scores, 1, function(x) sum(x[others] > x[sprintf("SZ07_%s", region)])),
        pathway_neuronal_score = sapply(rownames(scores), function(x) {
            pathway_neuronal_score$neuronality[pathway_neuronal_score$pathway == x][1]
        }),
        stringsAsFactors = FALSE
    )
    sz  <- scores[gsva_res$pathway, sprintf("SZ07_%s", region)] # SZ07 score per pathway
    M   <- scores[gsva_res$pathway, others, drop = FALSE]
    sgn <- sign(sz); sgn[sgn == 0] <- 1                     # SZ07 direction (+1/-1)
    gsva_res$runnerup_gap <- sgn * sz - apply(sgn * M, 1, max, na.rm = TRUE)

    pg <- setNames(vapply(row.names(scores), classify, character(1)), row.names(scores))
    
    writexl::write_xlsx(gsva_res %>%
                            mutate(classifier = pg[pathway]) %>%
                            select(pathway, classifier, pathway_neuronal_score, everything()) %>%
                            arrange(desc(runnerup_gap)) ,
                        sprintf("%s/data/post-drop/trd-pass-gsva/gsva_results_%s.xlsx", git_folder, region))
    
    
    
    plot_donor_specificity <- plot_gsva_specificity_table(
        scores       = scores,                 # pathways x donors
        sz07_id      = sprintf("SZ07_%s", region),
        donor_groups = setNames(ifelse(grepl("SZ", coldata_collapsed$sample[coldata_collapsed$sample != sprintf("SZ07_%s", region)]), "SZ", "HC"),
                                coldata_collapsed$sample[coldata_collapsed$sample != sprintf("SZ07_%s", region)]),
        pathway_groups = pg,  # same named vector you already built
        group_order    = c("Synaptic","RNA / Ribosome Biogenesis","Calcium / Ion Transport",
                           "TLR & Cytokine Signalling","Humoral/Complement Immunity",
                           "Neurodevelopment & Axonal Repair","Other"),
        top_n_per_group = c(6,5,5,7,5,3,3),
        rank_by      = "gap", # "gap"/"z"/"score" — NOT the fgsea "padj"
        spec_all     = setNames(gsva_res$pathway_neuronal_score, gsva_res$pathway),  # composition column
        label_fn     = strip_msigdb_prefix,
        label_truncate = 65,  # compress long labels to <=80 chars
        label_wrap     = NULL, 
        render       = FALSE)
    
    ggplot2::ggsave(sprintf("%s/data/post-drop/trd-pass-gsva/plot_donor_specificity_table_%s.png", git_folder, region),
                    plot = plot_donor_specificity$grob,
                    width = 8, height = 6.5, dpi = 300)
}
    