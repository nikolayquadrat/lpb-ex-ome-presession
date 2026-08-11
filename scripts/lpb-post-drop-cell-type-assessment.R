rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(writexl)
library(stringr)

git_folder <- "C:/Users/Nikolay/Dropbox/Git/lpb-ex-ome-presession" # <-- SET THE GIT PATH
dir.create(sprintf("%s/data/post-drop/cell-type-assessment", git_folder), showWarnings = TRUE)

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
        glia = rowMeans(cbind(score_astrocyte, score_oligodendrocyte,
                              score_opc, score_microglia), na.rm = TRUE),
        neu_minus_glia = score_neuron - glia) %>%
    dplyr::select(donor, region, score_oligodendrocyte, score_astrocyte,
                  glia, score_neuron, neu_minus_glia) %>% 
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

markers_plot <- ggplot(data = scores, aes(x = donor, y = neu_minus_glia, fill = group))+
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
    labs(y = "neuron-minus-glia z-score", fill="")
ggsave(plot= markers_plot,
       sprintf("%s/data/post-drop/cell-type-assessment/markers.png", git_folder),
       width = 8, height = 2.75, dpi = 300, bg = "white")

ggplot(data = scores, aes(x = donor, y = 1, fill = neu_minus_glia))+
    geom_tile()+
    scale_fill_continuous(type = "viridis")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid=element_blank(),
          panel.border=element_blank())+
    facet_grid(cols = vars(region), scales="free_x", space  = "free_x")+
    labs(y = "neuron-minus-glia z-score")

# **** **** check some specific markers panels --------
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

# **** CIBERSORTx --------
cib_siletti <- read.delim(sprintf("%s/data/rnaseq-pipe/09_deconv/_cibersortx/siletti_cortex/results/CIBERSORTx_Adjusted.txt", git_folder))
writexl::write_xlsx(cib_siletti, sprintf("%s/data/rnaseq-pipe/09_deconv/_cibersortx/siletti_cortex/results/CIBERSORTx_Adjusted.xlsx", git_folder))

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
    labs(y = "inferred neuronal proportion", fill="")
ggsave(plot = deconv_plot,
       sprintf("%s/data/post-drop/cell-type-assessment/deconv.png", git_folder),
       width = 8, height = 2.75, dpi = 300, bg = "white")


ggpubr::ggarrange(plotlist = list(markers_plot, deconv_plot),
                  nrow = 2, ncol = 1,
                  heights = c(1,0.65),
                  common.legend = TRUE, legend = "right")
ggsave(sprintf("%s/data/post-drop/cell-type-assessment/cell_type.png", git_folder),
       width = 8, height = 5, dpi = 300, bg = "white")
