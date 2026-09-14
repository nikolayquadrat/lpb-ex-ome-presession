#' Classify a pathway name into a functional theme (improved recall).
#'
#' Refinement of the original 6-theme classifier: same themes, but each pattern
#' widened to catch the false negatives that were falling into "Other" (e.g.
#' NEUROTROPHIN/TRK, GLIOGENESIS, MYELIN for Neurodevelopment; RNA_POLYMERASE,
#' TRNA, SPLICEOSOME, NUCLEOLAR for RNA; CYTOKINE/INTERFERON/NFKB/T_CELL for
#' TLR; ANTIGEN_PROCESSING/MHC/PHAGOCYTOSIS/FC-receptor for Immune; specific
#' ion channels for Calcium). Order matters -- RNA is tested first because its
#' tokens are least ambiguous; broad immune tokens are kept immune-specific to
#' avoid pulling adhesion/Hippo/organ-development sets.
#'
#' NOTE: even with better recall, most pathways remain "Other" for a
#' GSVA-extreme set, because SZ07's extreme GSVA pathways are dominated by
#' themes NOT among these six (actin/cytoskeleton, carbohydrate/energy
#' metabolism) plus many HPO clinical-phenotype sets. Add themes or drop the
#' HP/HPO collection separately if you want to shrink "Other" further.
#'
#' @param p character scalar pathway name (MSigDB style, underscores).
#' @return theme label (character scalar).

# The old obsolete method:
# classify <- function(p) {
#     P <- toupper(p)
#     if (grepl("SYNAP|NEUREXIN|NMDA|POSTSYN|PRESYN|LEARNING|NEUROTRANSMITTER|VESICLE|EXOCYT|IQGAP|PATHWAY_OF_L1|DOPAMINE", P)) return("Synaptic")
#     if (grepl("NEURON_PROJECTION|NERVOUS_SYSTEM_DEVELOPMENT|AXON_GUIDANCE|SPINAL", P)) return("Neurodevelopment & Axonal Repair")
#     if (grepl("RRNA|PRERIBOSOME|RNA_3_END|MRNA_3_END|RIBOSOM|CLEAVAGE_INVOLVED", P)) return("RNA / Ribosome Biogenesis")
#     if (grepl("CALCIUM|ION_TRANSPORT|CHANNEL", P)) return("Calcium")
#     if (grepl("TNFA|INFLAMM|TOLL_LIKE|MYD88|INTERLEUKIN", P)) return("TLR & Cytokine Signalling")
#     if (grepl("COMPLEMENT|TOXINS|HUMORAL_IMMUNE|OXIDATIVE_DAMAGE", P)) return("Humoral/Complement Immunity")
#     "Other"
# }

classify <- function(p) {
    P <- toupper(p)

    # RNA / ribosome / transcription machinery -- tested first (unambiguous
    # tokens); exclude lncRNA/ncRNA-in-cancer-signalling sets, which are WNT/
    # cancer pathways that merely contain "NCRNA".
    if (grepl(paste0("RRNA|PRERIBOSOME|RIBOSOM|\\bTRNA|SNORNA|SNRNA|SPLICEOSOM|",
                     "PROCESSOME|NUCLEOLUS|NUCLEOLAR|RNA_POLYMERASE|POLYMERASE_I|",
                     "POLYMERASE_II|POLYMERASE_III|TRANSCRIPTION_INITIATION|",
                     "TRANSCRIPTION_ELONGATION|RNA_3_END|MRNA_3_END|MRNA_PROCESS|",
                     "MRNA_CATABOL|RNA_MODIF|RNA_PROCESS|RNA_PHOSPHODIESTER|",
                     "CLEAVAGE_INVOLVED|CAJAL|DEADENYL"), P) &&
        !grepl("LNCRNA|NCRNA", P))
        return("RNA / Ribosome Biogenesis")

    if (grepl(paste0("SYNAP|NEUREXIN|NMDA|POSTSYN|PRESYN|LEARNING|NEUROTRANSMITTER|",
                     "VESICLE|EXOCYT|IQGAP|PATHWAY_OF_L1|DOPAMINE|GLUTAMATE_|GABA|",
                     "DENDRIT|PHOSPHODIESTERASES_IN_NEURONAL"), P))
        return("Synaptic")

    if (grepl(paste0("NEURON_PROJECTION|NERVOUS_SYSTEM_DEVELOPMENT|AXON_GUIDANCE|",
                     "SPINAL|NEUROTROPHIN|_TRK_|NEUROGENESIS|GLIOGENESIS|MYELIN|",
                     "DEMYELINAT|OLIGODENDROCYTE|NEURON_DIFFERENTIATION|NEURON_DEATH|",
                     "NEURON_APOPTO|NEURON_MIGRATION"), P))
        return("Neurodevelopment & Axonal Repair")

    if (grepl(paste0("CALCIUM|ION_TRANSPORT|_CHANNEL|CATION_TRANSPORT|ANION_TRANSPORT|",
                     "VOLTAGE_GATED|MEMBRANE_POTENTIAL|ION_HOMEOSTASIS|",
                     "POTASSIUM_ION_TRANSPORT|SODIUM_ION_TRANSPORT"), P))
        return("Calcium / Ion Transport")

    if (grepl(paste0("TNFA|INFLAMM|TOLL_LIKE|MYD88|INTERLEUKIN|CYTOKINE|CHEMOKINE|",
                     "INTERFERON|NF_KAPPA|NFKB|JAK_STAT|LYMPHOCYTE|\\bT_CELL|\\bB_CELL|",
                     "LEUKOCYTE_MIGRATION|MACROPHAGE|MICROGLI"), P))
        return("TLR & Cytokine Signalling")

    if (grepl(paste0("COMPLEMENT|TOXINS|HUMORAL_IMMUNE|OXIDATIVE_DAMAGE|IMMUNOGLOBULIN|",
                     "ANTIGEN_PROCESS|ANTIGEN_PRESENT|\\bMHC|PHAGOCYT|OPSON|FC_GAMMA|",
                     "FCGAMMA|FC_RECEPTOR"), P))
        return("Humoral/Complement Immunity")

    "Other"
}
