rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(readxl)
library(writexl)
library(stringr)

git_folder <- "C:/Users/Nikolay/Dropbox/Git/lpb-ex-ome-presession" # <-- SET THE GIT PATH
dir.create(sprintf("%s/data/post-drop/cell-type-assessment", git_folder), showWarnings = FALSE)

# Check the cell-type imbalances ===========
# **** simple marker based --------
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
    mutate(z = as.numeric(scale(log2cpm))) %>%       # z per gene WITHIN brain region
    ungroup() %>%
    group_by(donor, region, cell_type) %>%
    summarise(score = mean(z), n_markers = n(), .groups = "drop") %>%
    pivot_wider(names_from = cell_type,
                values_from = c(score, n_markers))

# composite + the actual question — an analysis choice, so it lives here
scores <- scores_all %>%
    mutate(
        nonneuronal = rowMeans(cbind(
                              score_astrocyte,
                              score_oligodendroglial,
                              score_microglia,
                              score_mural,
                              score_endothelial,
                              score_fibroblast
                              ), na.rm = TRUE),
        neu_minus_nonneuronal = score_neuron - nonneuronal) %>%
    dplyr::select(donor, region,
                  nonneuronal,
                  score_neuron,
                  neu_minus_nonneuronal) %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    )) %>% 
    mutate(group = factor(group, levels = c("SZ07","SZ","HC"))) %>% 
    dplyr::arrange(group, donor)

scores$donor <- factor(scores$donor, levels = c(
    "SZ07","SZ01","SZ04","SZ05", "SZ06","SZ08","SZ10","SZ11",
    "HC1M","HC24","HC2M","HC31","HC318","HC3M","HC79","HC91"))
scores$region <- factor(scores$region, levels = c("BA9", "BA22p", "BA4"))

markers_plot <- ggplot(data = scores, aes(x = donor, y = neu_minus_nonneuronal, fill = group))+
    geom_col()+
    scale_fill_manual(values = c(
        "SZ07"="#d62728",
        "SZ"  ="#e69f00",
        "HC"  ="#0072b2"
    ))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid=element_blank())+
    facet_grid(cols = vars(region), scales="free_x", space  = "free_x")+
    labs(y = "neuronal-to-non-neuronal\nrepresentation score", fill="", x="")
ggsave(plot= markers_plot,
       sprintf("%s/data/post-drop/cell-type-assessment/markers.png", git_folder),
       width = 8, height = 2.75, dpi = 300, bg = "white")

# alternative plot
ggplot(data = scores, aes(x = donor, y = 1, fill = neu_minus_nonneuronal))+
    geom_tile()+
    scale_fill_continuous(type = "viridis")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid=element_blank(),
          panel.border=element_blank())+
    facet_grid(cols = vars(region), scales="free_x", space  = "free_x")+
    labs(y = "neuronal-to-non-neuronal representation score")

# **** **** markers by individual cell types --------
scores <- scores_all %>%
    dplyr::select(donor, region,
                  score_neuron,
                  score_astrocyte,
                  score_oligodendroglial,
                  score_microglia,
                  score_mural,
                  score_endothelial,
                  score_fibroblast,
                  score_ependymal,
                  score_bam_macrophage,
                  score_choroid_plexus) %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    )) %>% 
    mutate(group = factor(group, levels = c("SZ07","SZ","HC"))) %>% 
    tidyr::pivot_longer(-c(donor, region, group), names_to = "score", values_to = "z") %>% 
    dplyr::arrange(group, donor) %>% 
    mutate(score = sub("score_", "", score)) %>% 
    mutate(score = case_when(
        score == "choroid_plexus" ~ "CP",
        score == "bam_macrophage" ~ "macrophage",
        TRUE ~ score
    ))

scores$donor <- factor(scores$donor, levels = unique(scores$donor))
scores$region <- factor(scores$region, levels = c("BA9", "BA22p", "BA4"))
scores$score <- factor(scores$score, levels = c(
    "neuron",
    "astrocyte",
    "microglia",
    "macrophage",
    "oligodendroglial",
    "mural",
    "fibroblast",
    "endothelial",
    "ependymal",
    "CP"
    ))

ggplot(data = scores, aes(x = donor, y = z, fill = group))+
    geom_col()+
    scale_fill_manual(values = c(
        "SZ07"="#d62728",
        "SZ"  ="#e69f00",
        "HC"  ="#0072b2"
    ))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid=element_blank(),
          panel.background = element_rect(fill = "grey95", colour = "white"),
          strip.text.y.right = element_text(angle = 0))+
    facet_grid(rows = vars(score), cols = vars(region), scales="free", space  = "free_x",
               labeller = labeller(
                   .cols = label_wrap_gen(width = 20),
                   .rows = label_wrap_gen(width = 15)
               ))+
    labs(y = "mean z-score", fill = "")

ggsave(sprintf("%s/data/post-drop/cell-type-assessment/markers_detailed.png", git_folder),
       width = 10, height = 5, dpi = 300, bg = "white")

# **** **** endothelial program in ischemia  --------
scores <- scores_all %>%
    dplyr::select(donor, region,
                  score_bbb_transport,
                  score_endmt,
                  score_endothelial_activation,
                  score_endothelial_core,
                  score_fibroblast_core,
                  score_fibromyocyte,
                  score_fibrotic_ecm,
                  score_sz07_vascular_depleted
                  ) %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    )) %>% 
    mutate(group = factor(group, levels = c("SZ07","SZ","HC"))) %>% 
    tidyr::pivot_longer(-c(donor, region, group), names_to = "score", values_to = "z") %>% 
    dplyr::arrange(group, donor) %>% 
    mutate(score = sub("score_", "", score)) %>% 
    mutate(score = case_when(
        score == "bbb_transport" ~ "BBB transport",
        score == "endmt" ~ "EndMT",
        score == "endothelial_core" ~ "endothelial core",
        score == "endothelial_activation" ~ "endothelial activation",
        score == "fibroblast_core" ~ "fibroblast core",
        score == "fibrotic_ecm" ~ "fibrotic ECM",
        score == "sz07_vascular_depleted" ~ "sz07 vasc. depleted",
        TRUE ~ score
    ))

scores$donor <- factor(scores$donor, levels = unique(scores$donor))
scores$region <- factor(scores$region, levels = c("BA9", "BA22p", "BA4"))
scores$score <- factor(scores$score, levels = c(
    "fibroblast core",
    "endothelial core",
    "endothelial activation",
    "EndMT",
    "BBB transport",
    "fibrotic ECM",
    "fibromyocyte",
    "sz07 vasc. depleted"
    ))

ggplot(data = scores, aes(x = donor, y = z, fill = group))+
    geom_col()+
    scale_fill_manual(values = c(
        "SZ07"="#d62728",
        "SZ"  ="#e69f00",
        "HC"  ="#0072b2"
    ))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid=element_blank(),
          panel.background = element_rect(fill = "grey95", colour = "white"),
          strip.text.y.right = element_text(angle = 0))+
    facet_grid(rows = vars(score), cols = vars(region), scales="free", space  = "free_x",
               labeller = labeller(
                   .cols = label_wrap_gen(width = 20),
                   .rows = label_wrap_gen(width = 15)
               ))+
    labs(y = "mean z-score", fill = "")

ggsave(sprintf("%s/data/post-drop/cell-type-assessment/markers_ischemia_endo_programs.png", git_folder),
       width = 10, height = 4.5, dpi = 300, bg = "white")

# **** **** other markers panels --------
scores <- scores_all %>%
    dplyr::select(donor, region,
                  score_bowen19_type2_down,
                  score_bowen19_type2_up,
                  score_cytokine_response_genes,
                  score_fromer16_down,
                  score_fromer16_up,
                  score_lanz19_dlpfc_neuronal_down,
                  score_ruzicka24_down,
                  score_synaptic_readout,
                  score_microglia_activated
                  ) %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    )) %>% 
    mutate(group = factor(group, levels = c("SZ07","SZ","HC"))) %>% 
    tidyr::pivot_longer(-c(donor, region, group), names_to = "score", values_to = "z") %>% 
    dplyr::arrange(group, donor) %>% 
    mutate(score = sub("score_", "", score)) %>% 
    mutate(score = case_when(
        score == "bowen19_type2_down" ~ "bowen19\ntype2 (down)",
        score == "bowen19_type2_up" ~ "bowen19\ntype2 (up)",
        score == "cytokine_response_genes" ~ "cytokine\nresponse",
        score == "fromer16_down" ~ "fromer16\ndown",
        score == "fromer16_up" ~ "fromer16\nup",
        score == "lanz19_dlpfc_neuronal_down" ~ "lanz19 DLPFC\nneuronal down",
        score == "ruzicka24_down" ~ "ruzicka24\ndown",
        score == "microglia_activated" ~ "microglia\nactivated",
        TRUE ~ score
    ))

scores$donor <- factor(scores$donor, levels = unique(scores$donor))
scores$region <- factor(scores$region, levels = c("BA9", "BA22p", "BA4"))

ggplot(data = scores, aes(x = donor, y = z, fill = group))+
    geom_col()+
    scale_fill_manual(values = c(
        "SZ07"="#d62728",
        "SZ"  ="#e69f00",
        "HC"  ="#0072b2"
    ))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid=element_blank(),
          panel.background = element_rect(fill = "grey95", colour = "white"))+
    facet_grid(rows = vars(score), cols = vars(region), scales="free", space  = "free_x",
               labeller = labeller(
                   .cols = label_wrap_gen(width = 20),
                   .rows = label_wrap_gen(width = 15)
               ))

# **** CIBERSORTx --------
# **** **** neuronal-vs-nonneuronal ----------
cib_siletti <- read.delim(sprintf("%s/data/rnaseq-pipe/09_deconv/_cibersortx/siletti_cortex/results/CIBERSORTx_Adjusted.txt", git_folder))
cib_siletti_export <- cib_siletti %>% 
    filter(Mixture %in% drop_inner$RNA_ID) %>% 
    left_join(m[, c("sample", "donor", "region")] %>% distinct(),
              by = c("Mixture"="sample")) %>% 
    dplyr::select(donor, region, everything()) %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    )) %>% 
    mutate(group = factor(group, levels = c("SZ07", "SZ", "HC"))) %>% 
    arrange(group, region, donor, Mixture) %>% 
    dplyr::rename(Donor=donor, `Brodmann area`=region) %>% 
    dplyr::select(-group)
writexl::write_xlsx(cib_siletti_export, sprintf("%s/data/rnaseq-pipe/09_deconv/_cibersortx/siletti_cortex/results/CIBERSORTx_Adjusted.xlsx", git_folder))

cib_siletti <- cib_siletti %>% 
    filter(Mixture %in% drop_inner$RNA_ID) %>% 
    mutate(donor = sub("^(KH\\d+_|)*([^_]+)_.*$", "\\2", Mixture)) %>% 
    mutate(region = sub("^.*_(.*)_.*", "\\1", Mixture)) %>% 
    group_by(donor, region) %>% 
    summarise(neuronal = mean(neuronal), .groups = "drop") %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    ))
cib_siletti$donor <- factor(cib_siletti$donor, levels = c(
    "SZ07","SZ01","SZ04","SZ05", "SZ06","SZ08","SZ10","SZ11",
    "HC1M","HC24","HC2M","HC31","HC318","HC3M","HC79","HC91"))
cib_siletti$region <- factor(cib_siletti$region, levels = c("BA9", "BA22p", "BA4"))

deconv_plot <- ggplot(data = cib_siletti, aes(x = donor, y = neuronal, fill = group))+
    geom_col()+
    scale_fill_manual(values = c(
        "SZ07"="#d62728",
        "SZ"  ="#e69f00",
        "HC"  ="#0072b2"
    ))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid=element_blank())+
    facet_grid(cols = vars(region), scales="free_x", space  = "free_x")+
    labs(y = "inferred neuronal\nproportion", fill="")
ggsave(plot = deconv_plot,
       sprintf("%s/data/post-drop/cell-type-assessment/deconv.png", git_folder),
       width = 8, height = 2.75, dpi = 300, bg = "white")

# **** **** neuronal and non-neuronal by supercluster----------
cib_siletti_glia_plus <- read.delim(sprintf("%s/data/rnaseq-pipe/09_deconv/_cibersortx/siletti_glia_plus_supercluster/results/CIBERSORTx_Adjusted.txt", git_folder))
cib_siletti_glia_plus_export <- cib_siletti_glia_plus %>% 
    filter(Mixture %in% drop_inner$RNA_ID) %>% 
    left_join(m[, c("sample", "donor", "region")] %>% distinct(),
              by = c("Mixture"="sample")) %>% 
    dplyr::select(donor, region, everything()) %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    )) %>% 
    mutate(group = factor(group, levels = c("SZ07", "SZ", "HC"))) %>% 
    arrange(group, region, donor, Mixture) %>% 
    dplyr::rename(Donor=donor, `Brodmann area`=region) %>% 
    dplyr::select(-group)
writexl::write_xlsx(cib_siletti_glia_plus_export, sprintf("%s/data/rnaseq-pipe/09_deconv/_cibersortx/siletti_glia_plus_supercluster/results/CIBERSORTx_Adjusted.xlsx", git_folder))

cib_siletti_glia_plus <- cib_siletti_glia_plus %>% 
    filter(Mixture %in% drop_inner$RNA_ID) %>% 
    dplyr::rename(
        COP=Committed.oligodendrocyte.precursor,
        `Bergmann glia` = Bergmann.glia,
        OPC = Oligodendrocyte.precursor,
        `Choroid plexus` = Choroid.plexus
    ) %>% 
    mutate(donor = sub("^(KH\\d+_|)*([^_]+)_.*$", "\\2", Mixture)) %>% 
    mutate(region = sub("^.*_(.*)_.*", "\\1", Mixture)) %>% 
    mutate(OPC=OPC+COP) %>% 
    dplyr::select(-COP, -P.value, -Correlation, -RMSE) %>% 
    pivot_longer(-c(Mixture, donor, region), names_to = "cell_type", values_to = "prop") %>% 
    group_by(donor, region, cell_type) %>% 
    summarise(
        prop=mean(prop),
        .groups = "drop") %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("SZ", donor) ~ "SZ",
        grepl("HC", donor) ~ "HC",
        TRUE ~ NA_character_
    ))

cib_siletti_glia_plus$donor <- factor(cib_siletti_glia_plus$donor, levels = c(
    "SZ07","SZ01","SZ04","SZ05", "SZ06","SZ08","SZ10","SZ11",
    "HC1M","HC24","HC2M","HC31","HC318","HC3M","HC79","HC91"))
cib_siletti_glia_plus$region <- factor(cib_siletti_glia_plus$region, levels = c("BA9", "BA22p", "BA4"))
cib_siletti_glia_plus$cell_type <- factor(cib_siletti_glia_plus$cell_type, levels = c(
    "Vascular", "Choroid plexus","Ependymal","Fibroblast","Microglia",
    "Astrocyte","Bergmann glia","Oligodendrocyte","OPC",
    "Neuron"
    ))

deconv_plot2 <- ggplot(data = cib_siletti_glia_plus, aes(x = donor, y = prop, fill = cell_type))+
    geom_col()+
    scale_fill_manual(values = c(
        "Microglia"= "#a41796",
        "Choroid plexus"="#e25f60",
        "Ependymal"="#ed9b9c",
        "Fibroblast"="#aa555e",
        "Vascular"="#d62728",
        "Astrocyte"  = "#e69f00",
        "Bergmann glia"  = "#ffc032",
        "Oligodendrocyte"= "#c2ad0f",
        "OPC"="#f4e46f",
        "Neuron"  ="#0072b2"
    ))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid=element_blank())+
    facet_grid(cols = vars(region), scales="free_x", space  = "free_x")+
    labs(y = "inferred relative cell-type fraction", fill="")
ggsave(plot = deconv_plot2,
       sprintf("%s/data/post-drop/cell-type-assessment/deconv2.png", git_folder),
       width = 10, height = 2.75, dpi = 300, bg = "white")


# Combined plot -------
ggpubr::ggarrange(plotlist = list(markers_plot, grid::nullGrob(), deconv_plot),
                  nrow = 3, ncol = 1,
                  heights = c(1, 0.085, 0.75),
                  common.legend = TRUE, legend = "right")
ggsave(sprintf("%s/data/post-drop/cell-type-assessment/cell_type.png", git_folder),
       width = 7, height = 5, dpi = 300, bg = "white")

# individual markers --------
tpm <- readxl::read_xlsx(sprintf("%s/data/rnaseq-pipe/04_qc/00_expression_summary.xlsx", git_folder))
tpm <- tpm %>% 
    drop_na() %>% 
    pivot_longer(-c(gene_id, gene_name), names_to = "sample", values_to = "TPM") %>% 
    filter(sample %in% drop_inner$RNA_ID) %>% 
    left_join(drop_inner, by = c("sample"="RNA_ID")) %>% 
    mutate(group = case_when(
        donor == "SZ07" ~ "SZ07",
        grepl("^SZ", donor) ~ "SZ",
        TRUE ~ "HC"
    ))
selected_genes <- c("FABP3","RTN4R","CDH5", "ICAM1","OCLN","IFITM1","ABCG2")
selected_genes <- c("DCN","LUM","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL15A1") # every one is increased

tpm$donor <- factor(tpm$donor, levels = c(
    "SZ07","SZ01","SZ04","SZ05", "SZ06","SZ08","SZ10","SZ11",
    "HC1M","HC24","HC2M","HC31","HC318","HC3M","HC79","HC91"))
tpm$group <- factor(tpm$group, levels = c(
    "SZ07","SZ","HC"))

ggplot(data = tpm[tpm$gene_name %in% selected_genes, ] %>% 
           filter(!(sample == "KH28_HC2M_BA22p_S79" & gene_name == "ICAM1")), # massive outlier
       aes(x = group, y = TPM, fill = group))+
    geom_boxplot()+
    scale_fill_manual(values = c(
        "SZ07"="#d62728",
        "SZ"  ="#e69f00",
        "HC"  ="#0072b2"
    ))+
    facet_wrap(brain_region ~ gene_name, scale="free", nrow=3)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
          panel.grid=element_blank())

# check the (turns out to be wrong) discrepancy -----------
cib_siletti_glia_plus <- read.delim(sprintf("%s/data/rnaseq-pipe/09_deconv/_cibersortx/siletti_glia_plus_supercluster/results/CIBERSORTx_Adjusted.txt", git_folder))
combined_mcib <- scores_all %>% 
    select(donor, region, score_fibrotic_ecm) %>% 
    left_join(cib_siletti_glia_plus %>% filter(cell_type == "Fibroblast"), by = c("donor", "region"))
fit <- lm(
    score_fibrotic_ecm ~ prop + region,
    data = subset(combined_mcib, donor != "SZ07")
)
predict_SZ07 <- predict(
    fit,
    newdata = subset(combined_mcib, donor == "SZ07"),
    interval = "prediction"
)

ggplot(combined_mcib, aes(
    x = prop,
    y = score_fibrotic_ecm
))+
    geom_point(aes(color = group))+
    scale_colour_manual(values = c(
        "SZ07"="#d62728",
        "SZ"  ="#e69f00",
        "HC"  ="#0072b2"
    ))+
    facet_wrap(~ region)+
    geom_smooth(
        data = subset(combined_mcib, donor != "SZ07"),
        method = "lm",
        se = TRUE
    )
