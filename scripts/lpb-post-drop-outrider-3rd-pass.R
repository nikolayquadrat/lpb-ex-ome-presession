rm(list = ls())

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(purrr)
library(readxl)
library(writexl)
library(OUTRIDER)   # v1.22.0
library(fgsea)      # v1.30.0
library(msigdbr)    # v7.5.1
library(enrichplot) # v1.24.4

# Arguments ============
experiments <- list(
    "ba9_gtex" = list(
        "outrider_tag" = "ba9_gtex",
        "samples" = list(
            "SZ07"  = c("KH16_SZ07_BA9_S68", "SZ07_BA9_14"),
            "HC91"  = c("KH3_HC91_BA9_S57",  "HC91_BA9_29"),
            "SZ06"  = c("KH15_SZ06_BA9_S67", "SZ06_BA9_5"),
            "SZ08"  = c("SZ08_BA9_20", "KH17_SZ08_BA9_S69"),
            "HC1M"  = "KH7_HC1M_BA9_S60",
            "HC2M"  = "KH8_HC2M_BA9_S61",
            "HC31"  = "KH2_HC31_BA9_S56",
            "HC318" = "KH5_HC318_BA9_S58",
            "HC3M"  = "KH9_HC3M_BA9_S62",
            "HC79"  = "KH6_HC79_BA9_S59",
            "SZ01"  = "KH11_SZ01_BA9_S63",
            "SZ04"  = "KH13_SZ04_BA9_S65",
            "SZ05"  = "SZ05_BA9_2",
            "SZ10"  = "KH18_SZ10_BA9_S70",
            "SZ11"  = "KH19_SZ11_BA9_S71"
        )
    )
)
git_folder <- "C:/Users/Nikolay/Dropbox/Git/lpb-ex-ome-presession"
exome_name <- "e1-19"

# Functions ========
# **** custom plotting functions -------
source(sprintf("%s/scripts/lpb-post-drop-fgsea-emap.R", git_folder))
source(sprintf("%s/scripts/lpb-post-drop-outrider-plot-gsea-highlighted.R", git_folder))
source(sprintf("%s/scripts/lpb-post-drop-donor-specificity-table.R", git_folder))

# ENSG --> HGNC dictionary ==========
constraint <- fread(sprintf("%s/data/exome-pipe/data/gnomad_v4.1_constraint_metrics.tsv", git_folder))
constraint <- constraint[mane_select == TRUE]  # one row per gene
ensg_to_hgnc <- setNames(constraint$gene, sub("\\..*$", "", constraint$gene_id))

# Candidate mutations ===========
significant_pathways_w_tiered_genes <- readxl::read_xlsx(sprintf("%s/data/post-drop/fst-pass/ba9_gtex_SZ07_fgsea_results.xlsx", git_folder), sheet = "tiered_genes_in_pathways")
tier_genes <- unique(significant_pathways_w_tiered_genes$tiered_gene)

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

# Main loop =========
dir.create(sprintf("%s/data/post-drop/trd-pass", git_folder), showWarnings = FALSE)
for(exp in names(experiments)) {
    # exp <- "ba9_gtex"
    cat("Experiment", exp, "\n")
    
    # **** Data =========
    outrider_tag <- experiments[[exp]][["outrider_tag"]]
    ods <- readRDS(sprintf("%s/data/drop-pipe/%s/Output/processed_results/aberrant_expression/v47/outrider/outrider/ods.Rds",
                           git_folder, outrider_tag))
    z <- zScore(ods) # z-scores; Genes x Samples matrix
    
    list_of_samples <- experiments[[exp]][["samples"]]
    list_of_ranks <- list()
    
    for (donor in names(list_of_samples)) {
        # donor <- "SZ07"
        donor_z <- zScore(ods)[, list_of_samples[[donor]], drop = FALSE]
        donor_combined <- apply(donor_z, 1, median)
        gene_ranks <- data.frame(
            ensembl_id = rownames(donor_z),
            hgnc       = ensg_to_hgnc[sub("\\..*$", "", rownames(donor_z))],
            combined_z = donor_combined
        ) %>% 
            drop_na() %>% 
            group_by(hgnc) %>%
            slice_max(abs(combined_z), n = 1, with_ties = FALSE) %>%
            ungroup()
        
        gene_ranks <- gene_ranks[order(-gene_ranks$combined_z, na.last=TRUE), ]
        gene_ranks_vector <- gene_ranks$combined_z
        names(gene_ranks_vector) <- gene_ranks$hgnc
        list_of_ranks[[donor]] <- gene_ranks_vector
    }
    
    # **** GSEA test with the OURIDER results ========
    if (!file.exists(sprintf("%s/data/post-drop/trd-pass/%s_list_of_fgsea_results.rds", git_folder, exp))) {
        list_of_fgsea_results <- list()
        for (donor in names(list_of_samples)) {
            # Run fgsea
            cat("fgsea:", donor, "\n")
            fgsea_res <- fgsea(
                pathways = pathway_list,
                stats = list_of_ranks[[donor]],
                minSize = 20,  # exclude too narrow categories
                maxSize = 200, # exclude too broad categories
                nPermSimple = 10000
            )
            list_of_fgsea_results[[donor]] <- fgsea_res[order(padj), ]
        }
        
        # save the fgsea run results
        saveRDS(list_of_fgsea_results, file = sprintf("%s/data/post-drop/trd-pass/%s_list_of_fgsea_results.rds",
                                                      git_folder, exp))
    }
    
    # **** Aggregate and collapse the results ============
    list_of_fgsea_results <- readRDS(sprintf("%s/data/post-drop/trd-pass/%s_list_of_fgsea_results.rds",
                                             git_folder, exp))
    
    for (donor in names(list_of_fgsea_results)) {
        # donor <- "SZ07"
        cat("collapse:", donor, "\n")
        fgsea_res <- list_of_fgsea_results[[donor]]
        fgsea_res_sig_uncorrected <- fgsea_res[order(pval)] # collapsePathways walks the pathways in the order given and greedily keeps the first (most significant) representative of each redundant group.
        
        stopifnot(is.numeric(list_of_ranks[[donor]]), !is.null(names(list_of_ranks[[donor]])))
        stopifnot(all(fgsea_res_sig_uncorrected$pathway %in% names(pathway_list)))
        
        # **** **** collapse pathways accross signatures collections ----------
        collapsed_pathways <- collapsePathways( # could be long
            fgseaRes = fgsea_res_sig_uncorrected[padj < 0.05],
            pathways = pathway_list,
            stats    = list_of_ranks[[donor]],
            pval.threshold = 0.01
        )
        fgsea_res_sig <- fgsea_res_sig_uncorrected[fgsea_res_sig_uncorrected$pathway %in% collapsed_pathways$mainPathways, ]
        
        cat(donor, "\t",
            "  before:", dim(fgsea_res_sig_uncorrected[padj < 0.05])[1], "\t",
            "  after:",  dim(fgsea_res_sig)[1], "\n")
        if(dim(fgsea_res_sig)[1] == 0) {
            next
        }
        
        # **** **** emapplots ------------
        emapplot_uncorrected <- fgsea_emap(fgsea_res_sig_uncorrected[padj < 0.05], pathway_list = pathway_list,
                                           edge_source = "gene_sets", color_by = "NES", edge_cutoff = 0.05)
        ggsave(sprintf("%s/data/post-drop/trd-pass/%s_%s_emapplot_all.png",
                       git_folder, exp, donor),
               emapplot_uncorrected, width = 12, height = 9, dpi = 150, bg = "white")
        
        emapplot_collapsed   <- fgsea_emap(fgsea_res_sig, pathway_list = pathway_list,
                                           edge_source = "gene_sets", color_by = "NES", edge_cutoff = 0.05)
        ggsave(sprintf("%s/data/post-drop/trd-pass/%s_%s_emapplot_collapsed.png",
                       git_folder, exp, donor),
               emapplot_collapsed, width = 12, height = 9, dpi = 150, bg = "white")
        
        
        # **** **** tiered genes ----------
        fgsea_res_sig_uncorrected$pathway_genes <- lapply(fgsea_res_sig_uncorrected$pathway, function(x) {
            pathway_list[[x]]
        })
        
        tiered_genes_in_pathways <- list()
        for (tiered_gene in tier_genes) {
            # tiered_gene <- "SLC11A1"
            fgsea_res_with_tiered_genes_in_leading_edge <- fgsea_res_sig_uncorrected[padj < 0.05][sapply(1:dim(fgsea_res_sig_uncorrected[padj < 0.05])[1], function(x) { tiered_gene %in% fgsea_res_sig_uncorrected[padj < 0.05]$leadingEdge[[x]]}), ] %>% 
                mutate(tiered_gene = tiered_gene, evidence = "leading_edge")
            fgsea_res_with_tiered_genes_in_pathway <- fgsea_res_sig_uncorrected[sapply(1:dim(fgsea_res_sig_uncorrected)[1], function(x) { tiered_gene %in% fgsea_res_sig_uncorrected$pathway_genes[[x]]}), ] %>% 
                mutate(tiered_gene = tiered_gene, evidence = "pathway") %>% 
                filter(!(pathway %in% fgsea_res_with_tiered_genes_in_leading_edge$pathway))
            
            cat(tiered_gene, "\t",
                dim(fgsea_res_with_tiered_genes_in_leading_edge)[1], "\t",
                min(fgsea_res_with_tiered_genes_in_leading_edge$padj), "\t",
                dim(fgsea_res_with_tiered_genes_in_pathway)[1], "\t",
                min(fgsea_res_with_tiered_genes_in_pathway$padj),
                "\n")
            
            tiered_genes_in_pathways[[tiered_gene]] <- bind_rows(fgsea_res_with_tiered_genes_in_leading_edge,
                                                                 fgsea_res_with_tiered_genes_in_pathway)
            
        }
        tiered_genes_in_pathways_df <- bind_rows(tiered_genes_in_pathways) %>% 
            dplyr::select(tiered_gene, evidence, everything())
        
        # **** **** save results -----------
        writexl::write_xlsx(
            list(
                "fgsea_res_sig_uncorrected_all" = fgsea_res_sig_uncorrected %>% mutate(across(where(is.list), ~ purrr::map_chr(.x, ~ paste(.x, collapse = ", ")))),
                "fgsea_res_sig" = fgsea_res_sig %>% mutate(across(where(is.list), ~ purrr::map_chr(.x, ~ paste(.x, collapse = ", ")))),
                "tiered_genes_in_pathways" = tiered_genes_in_pathways_df %>% mutate(across(where(is.list), ~ purrr::map_chr(.x, ~ paste(.x, collapse = ", "))))
            ),
            sprintf("%s/data/post-drop/trd-pass/%s_%s_fgsea_results.xlsx",
                    git_folder, exp, donor))
    }
    
    # **** GSEA visualisation ========
    list_of_fgsea_results <- readRDS(sprintf("%s/data/post-drop/trd-pass/%s_list_of_fgsea_results.rds",
                                             git_folder, exp))
    top_gsea_plots <- list() # main top GSEA plots
    top_gsea_plots_tier <- list() # top GSEA plots for signatures with the affected genes in the leading edge list
    for (donor in names(list_of_samples)) {
        fgsea_res <- list_of_fgsea_results[[donor]]
        fgsea_res$signature <- case_when(
            grepl("^WP_", fgsea_res$pathway) ~ "WP",
            grepl("^REACTOME_", fgsea_res$pathway) ~ "REACTOME",
            grepl("^KEGG_", fgsea_res$pathway) ~ "KEGG",
            fgsea_res$pathway %in% hpo$gs_name ~ "HPO",
            grepl("^HALLMARK_", fgsea_res$pathway) ~ "HALLMARK",
            grepl("^GOCC_", fgsea_res$pathway) ~ "GOCC",
            grepl("^GOMF_", fgsea_res$pathway) ~ "GOMF",
            grepl("^GOBP_", fgsea_res$pathway) ~ "GOBP",
            TRUE ~ "OTHER"
        )
        for(sig in names(table(fgsea_res$signature))) {
            cat(donor, "\t", sig, "\n")
            
            # top pathways
            fgsea_res_sig <- fgsea_res[fgsea_res$signature == sig,]
            top_pathways_up   <- fgsea_res_sig[ES > 0][order(padj)][1:8, pathway]
            top_pathways_down <- fgsea_res_sig[ES < 0][order(padj)][1:8, pathway]
            top_pathways <- c(top_pathways_up, rev(top_pathways_down))
            top_pathways <- top_pathways[!is.na(top_pathways)]
            
            if (length(top_pathways) > 0) {
                top_gsea_plots[[donor]][[sig]] <- plotGseaTableHighlighted(
                    pathways = pathway_list[top_pathways],
                    stats = list_of_ranks[[donor]],
                    fgseaRes = fgsea_res,
                    highlight_genes = tier_genes,
                    highlight_color = "#e66c2c",
                    highlight_lwd = 0.8,
                    leading_edge_color = "#2ca6e6",
                    leading_edge_lwd = 0.5,
                    colwidths = c(5, 4, 0.8, 1.2, 2),
                    gseaParam = 0.5,
                    render = FALSE
                )
            }
            
            # top pathways with tier genes
            fgsea_res_tier <- fgsea_res[fgsea_res$signature == sig,]
            fgsea_res_tier <- fgsea_res_tier[sapply(1:dim(fgsea_res_tier)[1], function(x) { any(tier_genes %in% fgsea_res_tier$leadingEdge[[x]])}),]
            top_pathways_up   <- fgsea_res_tier[ES > 0 & padj < 0.05][order(padj)][1:12, pathway]
            top_pathways_down <- fgsea_res_tier[ES < 0 & padj < 0.05][order(padj)][1:12, pathway]
            top_pathways <- c(top_pathways_up, rev(top_pathways_down))
            top_pathways <- top_pathways[!is.na(top_pathways)]
            
            if (length(top_pathways) > 0) {
                top_gsea_plots_tier[[donor]][[sig]] <- plotGseaTableHighlighted(
                    pathways = pathway_list[top_pathways],
                    stats = list_of_ranks[[donor]],
                    fgseaRes = fgsea_res,
                    highlight_genes = tier_genes,
                    highlight_color = "#e66c2c",
                    highlight_lwd = 0.8,
                    leading_edge_color = "#2ca6e6",
                    leading_edge_lwd = 0.5,
                    colwidths = c(5, 4, 0.8, 1.2, 2),
                    gseaParam = 0.5,
                    render = FALSE
                )
            }
        }
        if(purrr::pluck_exists(top_gsea_plots, donor, "HALLMARK")) {
            hallmark_plot <- top_gsea_plots[[donor]][["HALLMARK"]]
            ggsave(
                sprintf("%s/data/post-drop/trd-pass/%s_%s_gsea_hallmark_plot.png",
                        git_folder, exp, donor),
                plot=hallmark_plot,
                width = 2000, height = 800, units = "px", bg = "white", dpi = 165)
        }
        if(purrr::pluck_exists(top_gsea_plots, donor, "GOBP")) {
            gobp_plot <- top_gsea_plots[[donor]][["GOBP"]]
            ggsave(
                sprintf("%s/data/post-drop/trd-pass/%s_%s_gsea_gobp_plot.png",
                        git_folder, exp, donor),
                plot=gobp_plot,
                width = 2000, height = 800, units = "px", bg = "white", dpi = 165)
        }
        if(purrr::pluck_exists(top_gsea_plots, donor, "KEGG")) {
            kegg_plot <- top_gsea_plots[[donor]][["KEGG"]]
            ggsave(
                sprintf("%s/data/post-drop/trd-pass/%s_%s_gsea_kegg_plot.png",
                        git_folder, exp, donor),
                plot=kegg_plot,
                width = 2000, height = 800, units = "px", bg = "white", dpi = 165)
        }
    }
    
    saveRDS(top_gsea_plots,
            sprintf("%s/data/post-drop/trd-pass/%s_list_of_gsea_plots.rds",
                    git_folder, exp))
    saveRDS(top_gsea_plots_tier,
            sprintf("%s/data/post-drop/trd-pass/%s_list_of_gsea_tiered_plots.rds",
                    git_folder, exp))
    
}

# Aggregate ==========
exp <- "ba9_gtex"
list_of_tables <- lapply(names(experiments[[exp]][["samples"]]), function(x) {
    readxl::read_xlsx(sprintf("%s/data/post-drop/trd-pass/ba9_gtex_%s_fgsea_results.xlsx", git_folder, x),
                      sheet = "fgsea_res_sig_uncorrected_all")
})
names(list_of_tables) <- names(experiments[[exp]][["samples"]])
fgsea_res_all <- bind_rows(list_of_tables, .id="donor")

significant_pathways <- readxl::read_xlsx(sprintf("%s/data/post-drop/fst-pass/ba9_gtex_SZ07_fgsea_results.xlsx", git_folder), sheet = "fgsea_res_sig_uncorrected")
significant_pathways$n_donors_p005 <- NA
significant_pathways$n_donors_psmallerSZ07 <- NA
significant_pathways$SZ07_simple_NES <- NA
significant_pathways$other_donors_NES <- NA
significant_pathways$NES_gap <- NA

for (sp in significant_pathways$pathway) {
    # sp <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
    fgsea_res_other_donors_sp <- fgsea_res_all %>% 
        filter(pathway == sp) %>% 
        filter(donor != "SZ07")
    fgsea_res_sz07_sp <- fgsea_res_all %>% 
        filter(pathway == sp) %>% 
        filter(donor == "SZ07")
    
    significant_pathways$n_donors_p005[significant_pathways$pathway == sp] <- sum(fgsea_res_other_donors_sp$padj <= 0.05)
    significant_pathways$n_donors_psmallerSZ07[significant_pathways$pathway == sp] <- sum(fgsea_res_other_donors_sp$padj <= significant_pathways$padj[significant_pathways$pathway == sp])
    significant_pathways$SZ07_simple_NES[significant_pathways$pathway == sp] <- fgsea_res_sz07_sp$NES[fgsea_res_sz07_sp$pathway == sp]  
    significant_pathways$other_donors_NES[significant_pathways$pathway == sp] <- list(fgsea_res_other_donors_sp$NES[fgsea_res_other_donors_sp$pathway == sp])
    significant_pathways$NES_gap[significant_pathways$pathway == sp] <- ifelse(fgsea_res_sz07_sp$NES[fgsea_res_sz07_sp$pathway == sp] > 0,
                                                                               fgsea_res_sz07_sp$NES[fgsea_res_sz07_sp$pathway == sp] - max(fgsea_res_other_donors_sp$NES[fgsea_res_other_donors_sp$pathway == sp]),
                                                                               -(fgsea_res_sz07_sp$NES[fgsea_res_sz07_sp$pathway == sp] - min(fgsea_res_other_donors_sp$NES[fgsea_res_other_donors_sp$pathway == sp])))
}

significant_pathways <- significant_pathways %>% 
    mutate(across(where(is.list), ~ purrr::map_chr(.x, ~ paste(.x, collapse = ", ")))) %>% 
    arrange(desc(NES_gap))

classify <- function(p) {
    P <- toupper(p)
    if (grepl("SYNAP|NEUREXIN|NMDA|POSTSYN|PRESYN|LEARNING|NEUROTRANSMITTER|VESICLE|EXOCYT|IQGAP|PATHWAY_OF_L1|DOPAMINE", P)) return("Synaptic")
    if (grepl("NEURON_PROJECTION|NERVOUS_SYSTEM_DEVELOPMENT|AXON_GUIDANCE|SPINAL", P)) return("Neurodevelopment & Axonal Repair")
    if (grepl("RRNA|PRERIBOSOME|RNA_3_END|MRNA_3_END|RIBOSOM|CLEAVAGE_INVOLVED", P)) return("RNA / Ribosome Biogenesis")
    if (grepl("CALCIUM|ION_TRANSPORT|CHANNEL", P)) return("Calcium")
    if (grepl("TNFA|INFLAMM|TOLL_LIKE|MYD88|INTERLEUKIN", P)) return("TLR & Cytokine Signalling")
    if (grepl("COMPLEMENT|TOXINS|HUMORAL_IMMUNE|OXIDATIVE_DAMAGE", P)) return("Humoral/Complement Immunity")
    "Other"
}
pg <- setNames(vapply(significant_pathways$pathway, classify, character(1)), significant_pathways$pathway)
significant_pathways$pathway_group <- as.vector(pg)
writexl::write_xlsx(significant_pathways %>% 
                        dplyr::select(pathway, pathway_group, everything()) %>% 
                        mutate(pathway_group = factor(pathway_group, levels = c(
                            "Synaptic", "RNA / Ribosome Biogenesis", "Calcium", "TLR & Cytokine Signalling", "Humoral/Complement Immunity", "Neurodevelopment & Axonal Repair", "Other"
                        ))) %>% 
                        arrange(pathway_group, desc(NES_gap)),
                    sprintf("%s/data/post-drop/trd-pass/trd_pass_significant_pathways_fgsea_results.xlsx",
                            git_folder))

out <- plot_donor_specificity_table(
    res             = significant_pathways,
    donor_ids       = names(experiments[[exp]][["samples"]])[names(experiments[[exp]][["samples"]]) != "SZ07"],  # POSITIONAL — see caveat below
    donor_groups    = ifelse(grepl("SZ",names(experiments[[exp]][["samples"]])[names(experiments[[exp]][["samples"]]) != "SZ07"]), "SZ", "HC"), # same order
    pathway_groups  = pg,
    group_order     = c("Synaptic", "RNA / Ribosome Biogenesis", "Calcium", "TLR & Cytokine Signalling", "Humoral/Complement Immunity", "Neurodevelopment & Axonal Repair", "Other"),
    top_n_per_group = c(6, 5, 5, 7, 5, 3, 3),
    rank_by         = "gap", # or "z" / "padj"
    label_fn        = strip_msigdb_prefix,
    label_wrap      = 80, # keeps the WHOLE name
    colwidths       = c(9, 1.2, 1.0, 1.3, 5.5)
)
out$data   # the numbers behind the figure: gap, z, rank, n_more_extreme

# Check the cell-type imbalances ===========

library(dplyr); library(tidyr); library(stringr)
drop_inner <- readxl::read_xlsx(sprintf("%s/data/drop-pipe/drop_summaries.xlsx", git_folder))
m <- read.delim(sprintf("%s/data/rnaseq-pipe/04_qc/00_cell_marker_expression.tsv", git_folder))

# study-specific parsing lives here, not in the workflow
m <- m %>%
    mutate(
        region = str_extract(sample, "BA22p|BA9|BA4"),   # BA22p first: longest wins
        donor  = str_extract(sample, "(HC|SZ)\\d+M?")
    )

scores_all <- m %>%
    filter(sample %in% drop_inner$RNA_ID) %>% 
    filter(!is.na(log2cpm)) %>%
    group_by(region, gene) %>%
    filter(n() >= 3, sd(log2cpm) > 0) %>%            # need spread to z-score
    mutate(z = as.numeric(scale(log2cpm))) %>%       # z per gene WITHIN region
    ungroup() %>%
    group_by(sample, donor, region, cell_type) %>%
    summarise(score = mean(z), n_markers = n(), .groups = "drop") %>%
    pivot_wider(names_from = cell_type,
                values_from = c(score, n_markers))

# composite + the actual question — an analysis choice, so it lives here
scores <- scores_all %>%
    mutate(
        glia = rowMeans(cbind(score_astrocyte, score_oligodendrocyte,
                              score_opc, score_microglia), na.rm = TRUE),
        neu_minus_glia = score_neuron - glia) %>%
    dplyr::select(donor, region, sample, score_oligodendrocyte, score_astrocyte,
           glia, score_neuron, neu_minus_glia) %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    )) %>% 
    mutate(group = factor(group, levels = c("SZ07","SZ","HC"))) %>% 
    dplyr::arrange(group, sample)

scores$sample <- factor(scores$sample, levels = scores$sample)
ggplot(data = scores, aes(x = sample, y = neu_minus_glia, fill = group))+
    geom_col()+
    scale_fill_manual(values = c(
        "SZ07"="#d62728",
        "SZ"  ="#e69f00",
        "HC"  ="#0072b2"
    ))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_grid(cols = vars(region), scales="free_x", space  = "free_x")

scores <- scores_all %>%
    dplyr::select(donor, region, sample,
                  score_bowen19_type2_down,
                  score_bowen19_type2_up,
                  score_cytokine_response_genes,
                  score_fromer16_down,
                  score_fromer16_up,
                  score_lanz19_dlpfc_neuronal_down,
                  score_ruzicka24_down,
                  score_synaptic_readout) %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    )) %>% 
    mutate(group = factor(group, levels = c("SZ07","SZ","HC"))) %>% 
    tidyr::pivot_longer(-c(donor, region, sample, group), names_to = "score", values_to = "z") %>% 
    dplyr::arrange(group, sample) %>% 
    mutate(score = sub("score_", "", score)) %>% 
    mutate(score = case_when(
        score == "bowen19_type2_down" ~ "bowen19\ntype2 (down)",
        score == "bowen19_type2_up" ~ "bowen19\ntype2 (up)",
        score == "cytokine_response_genes" ~ "cytokine\nresponse",
        score == "fromer16_down" ~ "fromer16\ndown",
        score == "fromer16_up" ~ "fromer16\nup",
        score == "lanz19_dlpfc_neuronal_down" ~ "lanz19 DLPFC\nneuronal down",
        score == "ruzicka24_down" ~ "ruzicka24\ndown",
        score == "synaptic_readout" ~ "synaptic\nreadout",
        TRUE ~ NA_character_
    ))

scores$sample <- factor(scores$sample, levels = unique(scores$sample))
ggplot(data = scores, aes(x = sample, y = z, fill = group))+
    geom_col()+
    scale_fill_manual(values = c(
        "SZ07"="#d62728",
        "SZ"  ="#e69f00",
        "HC"  ="#0072b2"
    ))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_grid(rows = vars(score), cols = vars(region), scales="free", space  = "free_x",
               labeller = labeller(
                   .cols = label_wrap_gen(width = 20),
                   .rows = label_wrap_gen(width = 15)
               ))
