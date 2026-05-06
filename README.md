# Joint Exome/Transcriptome Mutation Prioritisation in Low SZ-PRS Brain
This LPB project, low PRS SZ brain, exome & transcriptome analysis

<! <img src="images/picture.png" alt="Title" width="70%"> >
data in /mnt/data/exome on yc-ncpz

## Reference data acquisition and post-processing (scripts\lpb-exome-priritisation-collect-data.sh)
A bash pipeline (download_references.sh) was developed to assemble all reference and annotation resources required by the variant-calling and annotation workflow. The script is idempotent, supports atomic resumption after interruption, and produces a manifest (TSV format) recording each managed file's path, size, SHA-256 hash, source URL, and validation status.

### Reference genome and known-sites VCFs.
The GRCh38 reference assembly with ALT contigs (Homo_sapiens_assembly38.fasta) and its FAI/dict indices were obtained from the Broad Institute's public Google Cloud bucket (gcp-public-data--broad-references/hg38/v0). Known-sites VCFs for GATK BQSR — HapMap 3.3, the Mills and 1000 Genomes gold-standard indels, and dbSNP 138 — were retrieved from the same source.

### Functional annotation infrastructure.
Ensembl VEP cache release 112 was downloaded from ftp.ensembl.org. VEP plugin source code was cloned from two repositories: the official Ensembl VEP_plugins repository pinned to release/115 (required for dbNSFP v5 column-layout compatibility), and the LOFTEE plugin from the konradjk/loftee repository on the grch38 branch. To accommodate VEP's single-directory plugin-loading model, all .pm files from both repositories were copied into a single canonical directory (plugins/flat/); copies (rather than symlinks) avoid Singularity bind-mount path-resolution failures.

### Plugin data resources.
AlphaMissense scores (AlphaMissense_hg38.tsv.gz) were obtained from the DeepMind public bucket and tabix-indexed. SpliceAI single-nucleotide-variant scores (Ensembl MANE GRCh38 release 110 mirror) were downloaded from the Ensembl FTP; the indel scores were obtained manually via the Illumina BaseSpace CLI (project 66029966) due to licensing constraints. CADD v1.7 SNV and indel scores were retrieved from the University of Washington's CADD distribution. LOFTEE supporting data (human ancestor reference, GERP conservation BigWig, and SQL conservation database) were obtained from personal.broadinstitute.org, with aria2 used preferentially for resilience against intermittent peering issues. dbNSFP v5.3.1a was supplied manually from genos.us and verified against its upstream MD5 checksum.

### REVEL post-processing.
REVEL scores (May 2021 release with Ensembl transcript IDs) require explicit transformation to be readable by VEP's REVEL plugin: (i) the published CSV is converted to TSV, (ii) chromosome names are prefixed with "chr" to match the reference assembly's UCSC-style naming, (iii) the column-header line is prefixed with "#" and the file is indexed via tabix -c '#' rather than tabix -S 1, ensuring tabix -h queries return the header line as expected by the plugin's column-detection routine. Two automated sanity checks validate the resulting file: a BRCA1 lookup at chr17:43106478 and a header-retrieval test.

### Population-frequency annotation.
gnomAD v4.1 exome sites VCFs were downloaded per chromosome (autosomes plus X and Y; ~184 GB total). A helper script (gnomad_strip_concat.sh) was generated and auto-invoked to strip each per-chromosome VCF to relevant AF columns (AF, AF_nfe, AF_afr, AF_amr, AF_eas, AF_sas, AF_fin, AF_asj, nhomalt, AC, AN) using bcftools 1.19 in parallel (8 chromosomes concurrently with GNU parallel), concatenating the stripped files into a single bgzipped VCF (~30 GB final), and removing per-chromosome originals as each was processed to bound peak disk usage. ClinVar GRCh38 was retrieved from a dated NCBI archive snapshot (archive_2.0/<YEAR>/clinvar_<YYYYMMDD>.vcf.gz) for run reproducibility. The gnomAD v4.1 constraint metrics table was downloaded from the gnomAD release directory.

### Gene panels.
 Schizophrenia, bipolar, and autism gene-burden results were obtained as TSV files from the SCHEMA, BipEx, and ASC web applications, respectively. SCHEMA gene results were joined to HGNC symbols via the gnomAD constraint table (Ensembl gene-ID match) for downstream gene-symbol-based filtering. The DDG2P / Genomics England PanelApp panel (ID 484) was retrieved via the panel's TSV download endpoint, with a JSON-API fallback and automated JSON-to-TSV conversion. The capture-kit BED file is supplied externally from the wet-laboratory provider.

### Validation.
Final completeness validation iterates over expected files, confirming non-zero size and presence of required indices (.tbi for VCF/TSV.gz, .fai/.dict for FASTA). Version coherence between VEP cache, plugin branch, and dbNSFP release is enforced at startup. Post-processing failures (REVEL conversion, SCHEMA HGNC join, gnomAD strip+concat, plugin flattening) propagate to the script's exit code.

## Variant-calling and annotation pipeline (scripts\lpb-exome-prioritisation-pipe.smk)
A 27-rule Snakemake workflow processes paired-end whole-exome sequencing data from raw FASTQ to tier-classified candidate variant tables suitable for downstream DROP monoallelic-expression analysis. Each rule executes inside a per-tool Singularity container, ensuring tool-version reproducibility; rules are designed to resume from intermediate outputs after interruption.

### Read processing and alignment.
Adapter and quality trimming was performed with fastp 0.23.4 using default Illumina-adapter detection. Trimmed reads were aligned to GRCh38 (with ALT contigs) using BWA-MEM2 v2.2.1 (rule r02a_bwamem2_align), with read-group tags injected for sample identification. Alignments were coordinate-sorted and indexed with samtools 1.19. PCR/optical duplicates were marked using GATK MarkDuplicatesSpark (GATK 4.5.0.0), which performs sorting and duplicate marking in a single Spark-parallelized pass.

### Base-quality recalibration.
GATK BaseRecalibrator computed per-read-group covariates against HapMap 3.3, Mills/1000G indels, and dbSNP 138 known-sites VCFs (rule r04a_bqsr_table). ApplyBQSR produced recalibrated BAMs.

### Variant calling and joint genotyping.
Per-sample GVCFs were produced by GATK HaplotypeCaller in -ERC GVCF mode. The cohort GVCFs were imported into a per-cohort GenomicsDB datastore (r06_genomicsdbimport), and joint genotyping was performed across all samples by GATK GenotypeGVCFs (r07_genotypegvcfs). Joint calling was chosen over per-sample calling to share evidence at borderline-quality sites; SNV and indel cohort sizes (n=32) below the VQSR threshold made hard filtering preferable to model-based recalibration.

### Hard filtering.
Following GATK best practices for small cohorts, SNVs and indels were extracted and filtered separately. SNV filters: QD < 2.0, FS > 60.0, MQ < 40.0, MQRankSum < -12.5, ReadPosRankSum < -8.0, SOR > 3.0. Indel filters: QD < 2.0, FS > 200.0, ReadPosRankSum < -20.0, SOR > 10.0. Filtered SNVs and indels were merged and PASS-filtered (r08e_merge_and_pass).

### Normalization.
The cohort VCF was left-aligned and multiallelics were split using bcftools norm -m -any --check-ref w against the reference FASTA (rule r09_normalize).

### Internal artifact filtering.
Four rules implement relatedness-aware artifact filtering. PLINK2 (v2.00a5.10) with autosome restriction and MAF/genotype-rate filters (--maf 0.05 --geno 0.05 --hwe 1e-6) produced a BED file from the cohort VCF (r09b_make_plink_bed); KING-format kinship coefficients were computed via --make-king-table (r09c_kinship_table). A custom Python helper (pick_representatives.py) performed connected-component analysis on related-pair edges (KING kinship ≥ 0.0884, the second-degree-relative threshold), retaining one lexicographically-first representative per family (r09d_pick_representatives). The cohort VCF was then subset to representatives and re-tagged with bcftools +fill-tags AC,AN; sites with AC ≥ 3 among representatives that were absent or rare (AF < 0.001) in gnomAD v4.1 were flagged as artifact-suspect tuples (r09e_artifact_blacklist). This procedure removes recurrent batch / mapping / sample-prep artifacts that gnomAD-AF filtering alone cannot detect, while avoiding false positives from variants shared across related samples.

### Per-sample VCF extraction.
The cohort VCF was demultiplexed into per-sample VCFs by bcftools 1.19 with private-variant retention (r10_per_sample_vcf).
Functional annotation. The cohort VCF was annotated by Ensembl VEP 112 (rule r11_vep_annotate_cohort) using the offline cache, MANE-Select / canonical / biotype transcript prioritization (--pick_allele_gene --pick_order mane_select,canonical,biotype), and HGVS notation. Plugins were loaded from a flattened directory and included AlphaMissense, REVEL, LOFTEE (high-confidence loss-of-function flag), CADD v1.7, SpliceAI (max delta scores across acceptor/donor gain/loss), and dbNSFP v5.3.1a fields (MutationTaster_pred, PROVEAN_pred, MetaLR_pred, MetaRNN_pred, M-CAP_pred, PrimateAI_pred, ClinPred_pred, BayesDel_addAF_pred; LRT_pred and FATHMM_pred were retired in dbNSFP v5). VEP --custom annotations attached gnomAD v4.1 exome population-stratified allele frequencies and ClinVar clinical-significance fields.

### Tier classification.
A Python script (scripts\lpb-exome-prioritisation-tier-candidates.py) classifies per-sample variants into three tiers (rule r12_tier_candidates). Common variants (gnomAD popmax AF ≥ 0.001) and internal artifact sites are filtered first. Tier A captures DROP-testable predicted loss-of-function and splice-disrupting variants (LOFTEE high-confidence LoF, or SpliceAI max delta ≥ 0.20). Tier B captures rank-only damaging missense candidates in constrained genes (gnomAD LOEUF < 0.35) by AlphaMissense likely-pathogenic class (score ≥ 0.564) or REVEL ≥ 0.75. Tier C captures any rare protein-altering variant in a curated gene set (SCHEMA at FDR ≤ 0.25, BipEx, ASC, or DDG2P confidence ≥ 2). Per-panel boolean membership flags and panel-specific annotation columns (e.g., SCHEMA Q meta, OR for protein-truncating variants; DDG2P mode-of-inheritance and aggregated phenotypes) are appended to all output tables. Outputs comprise three tier-specific TSVs and a master TSV containing all rare-variant calls with the tier label as the leading column.