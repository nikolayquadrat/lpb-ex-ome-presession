# =============================================================================
# WES variant-calling pipeline for DROP MAE module input
# =============================================================================
# Produces:
#   (1) Joint-called, normalized, per-sample VCF for DROP MAE
#   (2) Tiered annotated mutation tables (Tier A / B / C) for VUS prioritization
#
# Design:
#   * BWA-MEM2 (check the container tag!)
#   * Joint calling via GenomicsDBImport + GenotypeGVCFs across the cohort
#   * Hard filtering split by variant type (SNP / indel) per GATK best practice;
#     VQSR and CNN are not appropriate for this sample size
#   * bcftools norm for left-alignment + multiallelic splitting
#   * Relatedness-aware internal artifact filtering: PLINK2 detects related
#     sample pairs (KING-style kinship), one representative is kept per
#     family, and variants recurrent across the unrelated representatives
#     but absent / ultra-rare in gnomAD are blacklisted before tier
#     classification. Removes batch / mapping / sample-prep artifacts that
#     gnomAD-AF filtering alone misses.
#   * VEP with plugins: SpliceAI, AlphaMissense, LOFTEE, CADD, REVEL, dbNSFP
#   * Per-rule singularity containers preserved (design choice from original)
#
# Reference compatibility:
#   * DNA alignment uses Broad Homo_sapiens_assembly38.fasta (WITH ALT contigs)
#     - standard GATK best-practice reference
#   * Coordinate-compatible with the no-ALT reference used by CMC/GTEx for
#     RNA-seq alignment; variants on primary chromosomes are identical between
#     the two and will pass through DROP MAE cleanly
# =============================================================================

# -----------------------------------------------------------------------------
# Reference paths (keep in sync with the Broad bundle you already use)
# -----------------------------------------------------------------------------
reference_genome_base_wo_path = "resources_broad_hg38_v0_Homo_sapiens_assembly38"
reference_genome_base = "/tmp/genome/genome_for_exome_pipe/" + reference_genome_base_wo_path
reference_genome = reference_genome_base + ".fasta"

# Known-sites for BQSR (from the Broad bundle -- already used in the original)
reference_hapmap    = "/tmp/variation/vcf_for_exome_pipe/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
reference_indels    = "/tmp/variation/vcf_for_exome_pipe/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
reference_dbsnp     = "/tmp/variation/vcf_for_exome_pipe/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"

# Capture target intervals (BED from the kit manufacturer -> list with padding)
capture_bed         = "/tmp/annotation/agilent/S33266436_Regions.bed"
capture_intervals   = "/tmp/annotation/agilent/S33266436_Regions.padded100.interval_list"

# VEP resources
# Cache: Ensembl VEP cache for GRCh38, release aligned with GENCODE v30 (or newer,
# but keep one version for consistency). The "merged" cache (GENCODE + RefSeq)
# is convenient; the "gencode" cache is smaller.
vep_cache_dir       = "/tmp/annotation/vep/cache_grch38"
vep_plugin_dir      = "/tmp/annotation/vep/plugins/flat"

# Plugin data files (download separately, see README)
spliceai_snv        = "/tmp/annotation/vep/plugin_data/spliceai_scores.raw.snv.hg38.vcf.gz"
spliceai_indel      = "/tmp/annotation/vep/plugin_data/spliceai_scores.raw.indel.hg38.vcf.gz" # manually
alphamissense_tsv   = "/tmp/annotation/vep/plugin_data/AlphaMissense_hg38.tsv.gz"
loftee_dir          = "/tmp/annotation/vep/plugin_data/loftee_hg38"
loftee_path_string  = "loftee_path:/tmp/annotation/vep/plugins/loftee_grch38"
cadd_snv            = "/tmp/annotation/vep/plugin_data/whole_genome_SNVs.tsv.gz"
cadd_indel          = "/tmp/annotation/vep/plugin_data/gnomad.genomes.r4.0.indel.tsv.gz"
revel_tsv           = "/tmp/annotation/vep/plugin_data/new_tabbed_revel_grch38.tsv.gz" # manually
dbnsfp_gz           = "/tmp/annotation/vep/plugin_data/dbNSFP5.3.1a_grch38.gz" # manually with registration

# Custom annotation files (not plugins, loaded via --custom)
gnomad_v4_exomes    = "/tmp/annotation/vep/custom/gnomad.exomes.v4.1.sites.stripped.vcf.gz"
clinvar_vcf         = "/tmp/annotation/vep/custom/clinvar.vcf.gz"

# Gene panels for Tier C prioritization
schema_genes        = "/tmp/data/SCHEMA_gene_results_with_hgnc.tsv"
bipex_genes         = "/tmp/data/BipEx_gene_results_with_hgnc.tsv"
asc_genes           = "/tmp/data/ASC_gene_results_with_hgnc.tsv"
ndd_genes           = "/tmp/data/DDG2P_panel.tsv"
loeuf_table         = "/tmp/data/gnomad_v4.1_constraint_metrics.tsv"

# -----------------------------------------------------------------------------
# Internal-recurrence (artifact) filtering -- relatedness-aware
# -----------------------------------------------------------------------------
# Joint-called variants are scanned for sites that recur in multiple unrelated
# samples but are absent or ultra-rare in gnomAD. Such sites are almost always
# batch / mapping / sample-prep artifacts. Relatedness is estimated by PLINK2
# (KING-style kinship), and one representative is kept per family before the
# recurrence count is taken.
# -----------------------------------------------------------------------------

# Kinship coefficient above which two samples are considered "related".
# 0.0884 is KING's documented cutoff for "2nd-degree or closer" relatives
# (grandparent-grandchild, half-siblings, avuncular). 0.177 is "1st-degree".
KINSHIP_THRESHOLD       = 0.0442 # 3rd degree

# A site is flagged as a recurrent artifact when:
#   - it is called in >= ARTIFACT_AC_MIN unrelated representatives, AND
#   - its gnomAD AF is < ARTIFACT_GNOMAD_AF_MAX (i.e. ultra-rare or absent)
# With ~24 representatives, AC>=3 is conservative; tighten to 2 for a stricter
# filter or loosen to 4 if you suspect related-cohort contamination of AC.
ARTIFACT_AC_MIN         = 2
ARTIFACT_GNOMAD_AF_MAX  = 0.001

# Helper script for picking unrelated representatives (greedy graph algorithm).
# Place alongside the tier-candidates script under /tmp/scripts/.
pick_representatives_py = "/tmp/scripts/lpb-exome-prioritisation-pick-family-representatives.py"

# -----------------------------------------------------------------------------
# Sample discovery
# -----------------------------------------------------------------------------
samples ,= glob_wildcards("/tmp/fastq/{sample}_R1_001.fastq.gz") # or whatever fastq is named
# samples = [s for s in samples if s not in ['e1-19-combined']] # testing
# samples = ['e1-19-combined', 'e1-1-combined', 'e1-3-combined'] # testing
print("Samples:", samples)

# -----------------------------------------------------------------------------
# Targets -- only what DROP needs
# -----------------------------------------------------------------------------
rule all:
    input:
        # 
        # (1) Per-sample VCFs for DROP MAE
        expand("/tmp/fastq/10_per_sample_vcf/{sample}.vcf.gz",     sample=samples),
        expand("/tmp/fastq/10_per_sample_vcf/{sample}.vcf.gz.tbi", sample=samples),
        # # (2) Tiered candidate tables for VUS prioritization
        expand("/tmp/fastq/12_tiered/{sample}_tierA_lof_splice.tsv", sample=samples),
        expand("/tmp/fastq/12_tiered/{sample}_tierB_missense.tsv",   sample=samples),
        expand("/tmp/fastq/12_tiered/{sample}_tierC_gene_panels.tsv",sample=samples),
        expand("/tmp/fastq/12_tiered/{sample}_master.tsv",           sample=samples),
        # (3) QC: kinship table + artifact blacklist (built once per cohort)
        # "/tmp/qc/cohort_kinship.kin0",
        # "/tmp/qc/representatives.txt",
        # "/tmp/qc/representatives.log",
        # "/tmp/qc/internal_artifact_sites.tsv",
        # "/tmp/qc/internal_artifact_count.txt",
        "/tmp/qc/inferred_sex.tsv",
        # (4) HLA class I typing with OptiType (best class I accuracy from
        #     exome/DNA data, FASTQ input via razer3 prefilter).
        expand("/tmp/fastq/13_hla/optitype/{sample}/{sample}_result.tsv", sample=samples)

        # "/tmp/fastq/11_vep/cohort.vep.tsv.gz",
        # "/tmp/fastq/11_vep/test_loftee.vep.tsv", # for the testing
    shell: "echo 'DROP-ready outputs produced.'"


# =============================================================================
# Reference preparation
# =============================================================================
rule r000a_samtools_faidx:
    input:  reference_genome
    output: reference_genome + ".fai"
    singularity: "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    shell:  "samtools faidx {input}"

rule r000b_picard_dict:
    """
    Create the FASTA's sequence dictionary (.dict). The GENOME_ASSEMBLY (AS=)
    argument is REQUIRED here -- without it Picard either omits the AS field
    or fills it from the input filename, both of which break downstream
    tools that compare genome metadata across files (notably DROP's
    MonoallelicExpression module, where R's BSgenome compares AS strings
    when merging DNA-VCF and RNA-counts objects, and refuses to merge if
    "hg38" doesn't match the dict's filename-derived label).

    Setting AS=hg38 here ensures consistency across:
      - the dict (this rule)
      - downstream BAM @SQ headers (BWA reads the dict)
      - downstream VCF ##contig lines (GATK reads the dict)
    so DROP's MAE concordance check just works without manual reheadering.
    """
    input:  reference_genome
    output: reference_genome_base + ".dict"
    singularity: "docker://quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
    shell:
        """
        picard CreateSequenceDictionary \
            R={input} \
            O={output} \
            GENOME_ASSEMBLY=hg38 \
            SPECIES=Homo_sapiens

        # Sanity check: every @SQ line should now carry AS:hg38.
        TOTAL=$(grep -c '^@SQ' {output})
        GOOD=$(grep -c '^@SQ.*AS:hg38' {output})
        if [ "$TOTAL" != "$GOOD" ]; then
            echo "ERROR: AS:hg38 not set on all @SQ lines: $GOOD/$TOTAL" >&2
            grep '^@SQ' {output} | head -3 >&2
            exit 1
        fi
        echo "[r000b] AS:hg38 set correctly on all $TOTAL @SQ lines in {output}" >&2
        """

rule r000c_bwamem2_index:
    input:  reference_genome
    output: reference_genome + ".bwt.2bit.64"
    singularity: "docker://quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0"
    shell:  "bwa-mem2 index {input}"

rule r000d_make_intervals:
    # Pad capture BED by 100 bp for boundary variant recovery, convert to
    # interval_list (GATK-native format)
    input:  bed = capture_bed,
            dict_file = reference_genome_base + ".dict"
    output: capture_intervals
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk BedToIntervalList \
            -I {input.bed} \
            -O /tmp/_tmp.interval_list \
            -SD {input.dict_file}
        gatk IntervalListTools \
            -I /tmp/_tmp.interval_list \
            -O {output} \
            --PADDING 100
        rm -f /tmp/_tmp.interval_list
        """

# =============================================================================
# Per-sample: FASTQ -> analysis-ready BAM
# =============================================================================
rule r01_fastp_trim:
    # Adapter/quality trimming + automatic adapter detection + per-sample QC
    input:
        r1 = "/tmp/fastq/{sample}_R1_001.fastq.gz",
        r2 = "/tmp/fastq/{sample}_R2_001.fastq.gz",
    output:
        r1      = temp("/tmp/fastq/01_fastp/{sample}_R1.trim.fastq.gz"),
        r2      = temp("/tmp/fastq/01_fastp/{sample}_R2.trim.fastq.gz"),
        html    = "/tmp/fastq/01_fastp/{sample}.fastp.html",
        json    = "/tmp/fastq/01_fastp/{sample}.fastp.json",
    threads: 2
    resources: disk_mb=10000+55000+10000+10000 
    priority: 10
    singularity: "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 36 \
            --thread {threads} \
            --html {output.html} --json {output.json}
        """

rule r02a_bwamem2_align:
    input:
        r1   = "/tmp/fastq/01_fastp/{sample}_R1.trim.fastq.gz",
        r2   = "/tmp/fastq/01_fastp/{sample}_R2.trim.fastq.gz",
        ref  = reference_genome,
        bwt2 = reference_genome + ".bwt.2bit.64",
    output:
        bam = temp("/tmp/fastq/02_bam/{sample}.unsorted.bam"),
    threads: 6
    resources: disk_mb=10000+55000+10000+10000
    priority: 80
    singularity: "docker://quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0"
    shell:
        """
        bwa-mem2 mem -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA" \
            -K 100000000 -Y \
            {input.ref} {input.r1} {input.r2} \
            > {output.bam}
        """

rule r02b_sort_index:
    input:
        bam = "/tmp/fastq/02_bam/{sample}.unsorted.bam",
    output:
        bam = temp("/tmp/fastq/02_bam/{sample}.sorted.bam"),
        bai = temp("/tmp/fastq/02_bam/{sample}.sorted.bam.bai"),
    threads: 2
    resources: disk_mb=10000+55000+10000+10000
    priority: 90
    singularity: "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule r03_markduplicates_spark:
    # Replaces picard MarkDuplicates + SortSam + samtools index with a single
    # Spark-backed tool. Output is coordinate-sorted and indexed.
    input:
        bam = "/tmp/fastq/02_bam/{sample}.sorted.bam",
        bai = "/tmp/fastq/02_bam/{sample}.sorted.bam.bai",
    output:
        bam     = temp("/tmp/fastq/03_markdup/{sample}.markdup.bam"),
        bai     = temp("/tmp/fastq/03_markdup/{sample}.markdup.bam.bai"),
        metrics = "/tmp/fastq/03_markdup/{sample}.markdup.metrics.txt",
    threads: 6
    resources: disk_mb=10000+55000+10000+10000
    priority: 100
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk MarkDuplicatesSpark \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --conf 'spark.executor.cores={threads}'
        """

# =============================================================================
# HLA class I typing (class I only; class II requires class-II-specific tools)
# =============================================================================
# r03b_optitype : OptiType, the most accurate tool for class I typing from
#                 exome/DNA short reads. Takes FASTQs, prefilters with razer3
#                 against OptiType's bundled IPD-IMGT/HLA reference, then solves
#                 an integer-linear-program for the best allele combination
#                 explaining the read evidence.
#
#
# The trimmed FASTQs from r01 are temp(). Snakemake retains temp files while any
# downstream rule still needs them, so this HLA rule consuming them keeps them
# alive until typing completes.
# =============================================================================
rule r03b_optitype:
    """
    OptiType class I HLA typing from trimmed paired-end FASTQs using the
    bioconda biocontainer (quay.io/biocontainers/optitype:1.3.5--hdfd78af_3).

    The biocontainer ships OptiType under /usr/local/bin/ but does NOT include
    a config.ini -- we generate one at runtime pointing to the bundled
    razers3 and using the GLPK solver (which is what bioconda ships).
    """
    input:
        r1 = "/tmp/fastq/01_fastp/{sample}_R1.trim.fastq.gz",
        r2 = "/tmp/fastq/01_fastp/{sample}_R2.trim.fastq.gz",
    output:
        tsv = "/tmp/fastq/13_hla/optitype/{sample}/{sample}_result.tsv",
        pdf = "/tmp/fastq/13_hla/optitype/{sample}/{sample}_coverage_plot.pdf",
    params:
        outdir = "/tmp/fastq/13_hla/optitype/{sample}",
        # Discovered path inside the bioconda biocontainer (1.3.5--hdfd78af_3):
        hla_ref = "/usr/local/bin/data/hla_reference_dna.fasta",
    threads: 4
    priority: 80
    singularity: "docker://quay.io/biocontainers/optitype:1.3.5--hdfd78af_3"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        TMPDIR=$(mktemp -d -p {params.outdir} optitype.XXXXXX)

        # The bioconda biocontainer doesn't ship a config.ini. Generate one
        # at runtime. razers3 path is fixed (/usr/local/bin/razers3 in the
        # biocontainer); solver is GLPK (the bioconda recipe doesn't include
        # CBC or CPLEX); threads match the rule's allocation.
        cat > $TMPDIR/config.ini <<EOF
[mapping]
razers3=/usr/local/bin/razers3
threads={threads}

[ilp]
solver=glpk
threads=1

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
EOF

        # razer3 prefilter -- map reads against OptiType's bundled HLA ref
        # to drastically reduce the input size for the ILP step.
        razers3 -i 95 -m 1 -dr 0 \
            -tc {threads} \
            -o $TMPDIR/{wildcards.sample}_R1.bam \
            {params.hla_ref} \
            {input.r1}
        razers3 -i 95 -m 1 -dr 0 \
            -tc {threads} \
            -o $TMPDIR/{wildcards.sample}_R2.bam \
            {params.hla_ref} \
            {input.r2}

        # Convert filtered BAMs back to FASTQ for OptiType
        samtools bam2fq $TMPDIR/{wildcards.sample}_R1.bam \
            > $TMPDIR/{wildcards.sample}_R1.fastq
        samtools bam2fq $TMPDIR/{wildcards.sample}_R2.bam \
            > $TMPDIR/{wildcards.sample}_R2.fastq

        # Run OptiType with our generated config
        OptiTypePipeline.py \
            -i $TMPDIR/{wildcards.sample}_R1.fastq $TMPDIR/{wildcards.sample}_R2.fastq \
            --dna \
            -c $TMPDIR/config.ini \
            -o {params.outdir} \
            -p {wildcards.sample}

        rm -rf $TMPDIR
        """

rule r04a_bqsr_table:
    input:
        bam         = "/tmp/fastq/03_markdup/{sample}.markdup.bam",
        bai         = "/tmp/fastq/03_markdup/{sample}.markdup.bam.bai",
        ref         = reference_genome,
        fai         = reference_genome + ".fai",
        dict_file   = reference_genome_base + ".dict",
        hapmap      = reference_hapmap,
        indels      = reference_indels,
        dbsnp       = reference_dbsnp,
        intervals   = capture_intervals,
    output:
        table       = "/tmp/fastq/04_bqsr/{sample}.recal.table",
    threads: 2
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            -L {input.intervals} --interval-padding 100 \
            --known-sites {input.hapmap} \
            --known-sites {input.indels} \
            --known-sites {input.dbsnp} \
            -O {output.table}
        """

rule r04b_apply_bqsr: # probably redundant, present in GATK best practices but absent in DRAGEN
    input:
        bam     = "/tmp/fastq/03_markdup/{sample}.markdup.bam",
        bai     = "/tmp/fastq/03_markdup/{sample}.markdup.bam.bai",
        table   = "/tmp/fastq/04_bqsr/{sample}.recal.table",
        ref     = reference_genome,
    output:
        bam     = "/tmp/fastq/04_bqsr/{sample}.recal.bam",
        bai     = "/tmp/fastq/04_bqsr/{sample}.recal.bai",
    threads: 2
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk ApplyBQSR \
            --java-options "-Xmx8G -Dsamjdk.compression_level=5" \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.table} \
            -O {output.bam}
        """

# =============================================================================
# Per-sample: BAM -> GVCF
# =============================================================================
rule r05_haplotypecaller_gvcf:
    input:
        bam         = "/tmp/fastq/04_bqsr/{sample}.recal.bam",
        bai         = "/tmp/fastq/04_bqsr/{sample}.recal.bai",
        ref         = reference_genome,
        intervals   = capture_intervals,
    output:
        gvcf        = "/tmp/fastq/05_gvcf/{sample}.g.vcf.gz",
        tbi         = "/tmp/fastq/05_gvcf/{sample}.g.vcf.gz.tbi",
    threads: 1
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk HaplotypeCaller \
            --java-options "-Xmx16G" \
            -I {input.bam} \
            -R {input.ref} \
            -L {input.intervals} --interval-padding 100 \
            -ERC GVCF \
            -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
            --native-pair-hmm-threads {threads} \
            -O {output.gvcf}
        """

# =============================================================================
# Cohort-level: GVCFs -> joint-called VCF
# =============================================================================
rule r06_genomicsdbimport:
    # Combines all per-sample GVCFs into a GenomicsDB datastore for joint calling.
    # Even with 8 samples this is beneficial: shared evidence improves rare-variant
    # calibration, and the annotation fields needed by DROP MAE are cleaner.
    input:
        gvcfs       = expand("/tmp/fastq/05_gvcf/{sample}.g.vcf.gz", sample=samples),
        tbis        = expand("/tmp/fastq/05_gvcf/{sample}.g.vcf.gz.tbi", sample=samples),
        intervals   = capture_intervals,
    output:
        gendb       = directory("/tmp/fastq/06_genomicsdb/cohort_db"),
        sentinel    = "/tmp/fastq/06_genomicsdb/import.done",
    params:
        map_string  = lambda wc, input: " ".join(f"-V {g}" for g in input.gvcfs),
    threads: 6
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        rm -rf {output.gendb}
        gatk GenomicsDBImport \
            --java-options "-Xmx16G" \
            {params.map_string} \
            --genomicsdb-workspace-path {output.gendb} \
            -L {input.intervals} --interval-padding 100 \
            --reader-threads {threads} \
            --merge-input-intervals
        touch {output.sentinel}
        """

rule r07_genotypegvcfs:
    input:
        gendb       = "/tmp/fastq/06_genomicsdb/cohort_db",
        sentinel    = "/tmp/fastq/06_genomicsdb/import.done",
        ref         = reference_genome,
        dbsnp       = reference_dbsnp,
    output:
        vcf         = "/tmp/fastq/07_joint/cohort.joint.vcf.gz",
        tbi         = "/tmp/fastq/07_joint/cohort.joint.vcf.gz.tbi",
    threads: 6
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk GenotypeGVCFs \
            --java-options "-Xmx16G" \
            -R {input.ref} \
            -V gendb://{input.gendb} \
            -D {input.dbsnp} \
            -G StandardAnnotation -G AS_StandardAnnotation \
            -O {output.vcf}
        """

# =============================================================================
# Filtering: split by variant type, apply GATK best-practice hard filters
# =============================================================================
# With 8 samples, VQSR is underpowered (needs ~30+ WGS or ~100+ WES for a
# well-calibrated Gaussian mixture). CNN is deprecated in recent GATK. Hard
# filtering by variant type is the appropriate choice here.
# Thresholds from: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112
# -----------------------------------------------------------------------------

rule r08a_select_snps:
    input:  vcf = "/tmp/fastq/07_joint/cohort.joint.vcf.gz",
            ref = reference_genome,
    output: vcf = temp("/tmp/fastq/08_filter/cohort.snps.vcf.gz"),
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk SelectVariants -R {input.ref} -V {input.vcf} \
            --select-type-to-include SNP -O {output.vcf}
        """

rule r08b_select_indels:
    input:  vcf = "/tmp/fastq/07_joint/cohort.joint.vcf.gz",
            ref = reference_genome,
    output: vcf = temp("/tmp/fastq/08_filter/cohort.indels.vcf.gz"),
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk SelectVariants -R {input.ref} -V {input.vcf} \
            --select-type-to-include INDEL --select-type-to-include MIXED \
            -O {output.vcf}
        """

rule r08c_hardfilter_snps:
    input:  vcf = "/tmp/fastq/08_filter/cohort.snps.vcf.gz",
            ref = reference_genome,
    output: vcf = temp("/tmp/fastq/08_filter/cohort.snps.filtered.vcf.gz"),
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk VariantFiltration -R {input.ref} -V {input.vcf} -O {output.vcf} \
            --filter-name "QD2"        --filter-expression "QD < 2.0" \
            --filter-name "FS60"       --filter-expression "FS > 60.0" \
            --filter-name "MQ40"       --filter-expression "MQ < 40.0" \
            --filter-name "MQRS-12.5"  --filter-expression "MQRankSum < -12.5" \
            --filter-name "RPRS-8"     --filter-expression "ReadPosRankSum < -8.0" \
            --filter-name "SOR3"       --filter-expression "SOR > 3.0"
        """

rule r08d_hardfilter_indels:
    input:  vcf = "/tmp/fastq/08_filter/cohort.indels.vcf.gz",
            ref = reference_genome,
    output: vcf = temp("/tmp/fastq/08_filter/cohort.indels.filtered.vcf.gz"),
    singularity: "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        """
        gatk VariantFiltration -R {input.ref} -V {input.vcf} -O {output.vcf} \
            --filter-name "QD2"        --filter-expression "QD < 2.0" \
            --filter-name "FS200"      --filter-expression "FS > 200.0" \
            --filter-name "RPRS-20"    --filter-expression "ReadPosRankSum < -20.0" \
            --filter-name "SOR10"      --filter-expression "SOR > 10.0"
        """

rule r08e_merge_and_pass:
    # Merge SNP + indel filtered VCFs, drop non-PASS records, drop sites with no
    # ALT genotype called in any sample (happens after joint calling).
    input:
        snps    = "/tmp/fastq/08_filter/cohort.snps.filtered.vcf.gz",
        indels  = "/tmp/fastq/08_filter/cohort.indels.filtered.vcf.gz",
    output:
        vcf     = "/tmp/fastq/08_filter/cohort.filtered.vcf.gz",
        tbi     = "/tmp/fastq/08_filter/cohort.filtered.vcf.gz.tbi",
    singularity: "docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    shell:
        """
        bcftools concat -a -Oz -o /tmp/_merged.vcf.gz {input.snps} {input.indels}
        bcftools sort -Oz -o /tmp/_sorted.vcf.gz /tmp/_merged.vcf.gz
        bcftools view -f PASS -Oz -o {output.vcf} /tmp/_sorted.vcf.gz
        bcftools index -t {output.vcf}
        rm -f /tmp/_merged.vcf.gz /tmp/_sorted.vcf.gz
        """

# =============================================================================
# Normalization -- CRITICAL step before annotation / DROP MAE
# -----------------------------------------------------------------------------
# Splits multiallelic sites into biallelic records and left-aligns indels.
# Required for correct allele matching in VEP (plugins like SpliceAI and
# AlphaMissense only match left-aligned biallelic records) and for correct
# allele counting in the DROP MAE module.
# =============================================================================
rule r09_normalize:
    input:
        vcf = "/tmp/fastq/08_filter/cohort.filtered.vcf.gz",
        ref = reference_genome,
    output:
        vcf = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz",
        tbi = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz.tbi",
    singularity: "docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    shell:
        """
        bcftools norm -m -any -f {input.ref} --check-ref w \
            -Oz -o {output.vcf} {input.vcf}
        bcftools index -t {output.vcf}
        """

# =============================================================================
# Kinship analysis and internal-artifact blacklist
# -----------------------------------------------------------------------------
# Four rules running between joint-call normalization and per-sample VCF
# extraction. The job: detect which samples are related, pick one
# representative per family, then flag (CHROM, POS, REF, ALT) tuples that
# recur in unrelated representatives but are absent / rare in gnomAD --
# almost certainly batch / sample-prep / mapping artifacts. The tier
# script consumes this blacklist and drops matched variants before tier
# classification.
# =============================================================================

rule r09b_make_plink_bed:
    # PLINK2 BED file used as input for kinship analysis.
    # Restrict to autosomes; apply MAF/genotype-rate filters because kinship
    # is best estimated from common, confidently-genotyped variants.
    input:
        vcf = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz",
        tbi = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz.tbi",
    output:
        bed = "/tmp/qc/cohort_for_king.bed",
        bim = "/tmp/qc/cohort_for_king.bim",
        fam = "/tmp/qc/cohort_for_king.fam",
    threads: 2
    singularity: "docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
    shell:
        """
        mkdir -p /tmp/qc
        plink2 --vcf {input.vcf} \
               --autosome \
               --allow-extra-chr \
               --max-alleles 2 \
               --maf 0.05 \
               --geno 0.05 \
               --hwe 1e-6 \
               --threads {threads} \
               --make-bed \
               --out /tmp/qc/cohort_for_king
        """

rule r09c_kinship_table:
    # Compute KING-format kinship coefficients via plink2's native
    # --make-king-table (produces .kin0 in the same format as the
    # standalone KING tool, no separate binary needed).
    # If no related pairs exist, plink2 may not produce .kin0 at all --
    # write an empty header-only placeholder so downstream rules don't
    # break.
    input:
        bed = "/tmp/qc/cohort_for_king.bed",
        bim = "/tmp/qc/cohort_for_king.bim",
        fam = "/tmp/qc/cohort_for_king.fam",
    output:
        kin0 = "/tmp/qc/cohort_kinship.kin0",
    threads: 2
    singularity: "docker://quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
    shell:
        """
        plink2 --bfile /tmp/qc/cohort_for_king \
               --make-king-table \
               --king-table-filter 0 \
               --threads {threads} \
               --out /tmp/qc/cohort_kinship
        if [ ! -s {output.kin0} ]; then
            printf "#FID1\\tIID1\\tFID2\\tIID2\\tNSNP\\tHETHET\\tIBS0\\tKINSHIP\\n" \
                > {output.kin0}
        fi
        """

rule r09d_pick_representatives:
    # Greedy graph algorithm: connected-component analysis on the
    # related-pair graph, one (lexicographically first) representative per
    # family. See pick_representatives.py for details.
    input:
        kin0   = "/tmp/qc/cohort_kinship.kin0",
        fam    = "/tmp/qc/cohort_for_king.fam",
        script = pick_representatives_py,
    output:
        reps = "/tmp/qc/representatives.txt",
        log  = "/tmp/qc/representatives.log",
    params:
        kinship_threshold = KINSHIP_THRESHOLD,
    singularity: "docker://quay.io/biocontainers/pandas:1.5.2"
    shell:
        """
        python {input.script} \
            --kin0 {input.kin0} \
            --fam  {input.fam} \
            --kinship_threshold {params.kinship_threshold} \
            --out_representatives {output.reps} \
            --out_log {output.log}
        """

rule r09e_artifact_blacklist:
    # Compute the artifact blacklist as a TSV of (CHROM, POS, REF, ALT)
    # tuples. Steps:
    #   1. subset the cohort VCF to representatives only
    #   2. recompute INFO/AC across just those representatives
    #   3. keep sites where AC >= ARTIFACT_AC_MIN
    #   4. annotate with gnomAD AF
    #   5. keep only sites with AF < ARTIFACT_GNOMAD_AF_MAX or AF missing
    # The output is the blacklist; intermediate files are removed.
    #
    # NOTES on past bugs avoided here:
    #   * Use bgzipped VCF (not BCF) for the gnomAD reference, because
    #     bcftools annotate --output-type Oz handles VCF.gz cleanly and
    #     it's consistent with the format used downstream by VEP --custom.
    #   * Don't try to rely on bcftools' implicit popmax field; we use AF
    #     directly because the stripped gnomAD VCF was reduced to overall
    #     AF + per-population AFs in the upstream stripping step.
    #   * Single-quote the bcftools -i expressions so shell doesn't
    #     interpolate, and double the curly braces if used (none here).
    input:
        vcf      = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz",
        tbi      = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz.tbi",
        reps     = "/tmp/qc/representatives.txt",
        gnomad   = gnomad_v4_exomes,
    output:
        sites    = "/tmp/qc/internal_artifact_sites.tsv",
        n_sites  = "/tmp/qc/internal_artifact_count.txt",
    params:
        ac_min   = ARTIFACT_AC_MIN,
        af_max   = ARTIFACT_GNOMAD_AF_MAX,
    threads: 2
    singularity: "docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    shell:
        r"""
        # Step 1+2+3: subset to representatives, recompute AC, keep recurrent
        bcftools view --threads {threads} \
                      -S {input.reps} \
                      --force-samples \
                      -Ou {input.vcf} \
        | bcftools +fill-tags -Ou - -- -t AC,AN \
        | bcftools view -i 'INFO/AC >= {params.ac_min}' \
                      -Oz -o /tmp/qc/_recurrent.vcf.gz
        bcftools index -t /tmp/qc/_recurrent.vcf.gz

        # Step 4: annotate with gnomAD AF.
        # The gnomAD reference here must be a bgzipped VCF (not BCF) so it
        # can be read both by bcftools and downstream by VEP --custom.
        bcftools annotate --threads {threads} \
            -a {input.gnomad} \
            -c INFO/gnomAD_AF:=AF \
            -Oz -o /tmp/qc/_recurrent_with_gnomad.vcf.gz \
            /tmp/qc/_recurrent.vcf.gz
        bcftools index -t /tmp/qc/_recurrent_with_gnomad.vcf.gz

        # Step 5: keep only sites where gnomAD AF is rare or absent, and
        # write the blacklist as plain tab-separated tuples.
        bcftools view -i 'INFO/gnomAD_AF < {params.af_max} || INFO/gnomAD_AF = "."' \
            -Ou /tmp/qc/_recurrent_with_gnomad.vcf.gz \
        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
            > {output.sites}

        n=$(wc -l < {output.sites})
        echo "$n" > {output.n_sites}
        echo "Internal artifact sites detected: $n" >&2

        rm -f /tmp/qc/_recurrent.vcf.gz /tmp/qc/_recurrent.vcf.gz.tbi \
              /tmp/qc/_recurrent_with_gnomad.vcf.gz /tmp/qc/_recurrent_with_gnomad.vcf.gz.tbi
        """

rule r09f_sex_coverage_per_sample:
    """Per-sample chrX/chrY/autosome mean depth over capture targets."""
    input:
        bam = "/tmp/fastq/04_bqsr/{sample}.recal.bam",
        bai = "/tmp/fastq/04_bqsr/{sample}.recal.bai",
        bed = capture_bed,
    output:
        tsv = "/tmp/qc/sex/{sample}.sex_depth.tsv",
    threads: 2
    singularity: "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv})"
        WORK=$(mktemp -d); trap 'rm -rf "${{WORK:-}}"' EXIT INT TERM

        awk 'BEGIN{{OFS="\t"}} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {{print $1,$2,$3}}' {input.bed} > "$WORK/auto.bed"
        awk 'BEGIN{{OFS="\t"}} $1 == "chrX" {{print $1,$2,$3}}' {input.bed} > "$WORK/chrX.bed"
        awk 'BEGIN{{OFS="\t"}} $1 == "chrY" {{print $1,$2,$3}}' {input.bed} > "$WORK/chrY.bed"

        md() {{
            [ -s "$2" ] || {{ echo 0; return; }}
            samtools bedcov -Q 20 "$2" "$1" \
                | awk '{{cov+=$NF; len+=($3-$2)}} END{{if(len>0) printf "%.6f",cov/len; else print 0}}'
        }}
        printf '%s\t%s\t%s\t%s\n' "{wildcards.sample}" "$(md {input.bam} $WORK/auto.bed)" \
            "$(md {input.bam} $WORK/chrX.bed)" "$(md {input.bam} $WORK/chrY.bed)" > {output.tsv}
        """

rule r09g_infer_sex:
    """Aggregate per-sample depths into the cohort sex-inference table."""
    input:
        depths = expand("/tmp/qc/sex/{sample}.sex_depth.tsv", sample=samples),
    output:
        tsv = "/tmp/qc/inferred_sex.tsv",
    run:
        with open(output.tsv, "w") as out:
            out.write("sample\tmean_depth_auto\tmean_depth_chrX\tmean_depth_chrY\t"
                      "x_ratio\ty_ratio\tinferred_sex\n")
            for path in input.depths:
                s, da, dx, dy = open(path).read().split()
                da, dx, dy = float(da), float(dx), float(dy)
                if da <= 0:
                    out.write(f"{s}\t{da}\t{dx}\t{dy}\tNA\tNA\tambiguous\n"); continue
                xr, yr = dx/da, dy/da
                sex = "XX" if (xr >= 0.80 and yr < 0.15) else \
                      "XY" if (xr < 0.65 and yr >= 0.15) else "ambiguous"
                out.write(f"{s}\t{da}\t{dx}\t{dy}\t{xr:.4f}\t{yr:.4f}\t{sex}\n")

# =============================================================================
# Per-sample VCFs for DROP MAE
# -----------------------------------------------------------------------------
# DROP's MAE module expects one VCF per RNA sample (pointed at via DNA_VCF_FILE
# in the sample annotation sheet). Produce per-sample VCFs that drop sites
# where that sample is 0/0 or ./. -- keeps DROP's input tidy.
# =============================================================================
rule r10_per_sample_vcf:
    input:
        vcf = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz",
        tbi = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz.tbi",
    output:
        vcf = "/tmp/fastq/10_per_sample_vcf/{sample}.vcf.gz",
        tbi = "/tmp/fastq/10_per_sample_vcf/{sample}.vcf.gz.tbi",
    params:
        unsorted_vcf = lambda wc: "/tmp/fastq/10_per_sample_vcf/" + wc.sample + ".unsorted.vcf.gz",
    singularity: "docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    shell:
        """
        bcftools view -s {wildcards.sample} -Ou {input.vcf} \
        | bcftools view -e 'GT="0/0" || GT="./." || GT="0|0" || GT=".|."' \
            -Oz -o {params.unsorted_vcf}
        bcftools sort -Oz -o {output.vcf} {params.unsorted_vcf}
        bcftools index -t {output.vcf}
        """

# =============================================================================
# Annotation -- VEP with plugins
# -----------------------------------------------------------------------------
# Big upgrade over the original:
#   * SpliceAI  -- essential for splice-VUS interpretation (Tier A candidates).
#                  Uses split_output=1 so the 4 delta-scores (DS_AG/AL/DG/DL)
#                  and 4 delta-positions (DP_AG/AL/DG/DL) are emitted as
#                  separate columns, making it easier to identify the splicing
#                  mechanism (e.g., donor loss vs acceptor gain) rather than
#                  collapsing into a single SpliceAI_pred string.
#   * AlphaMissense -- state-of-the-art missense pathogenicity (Tier B)
#   * LOFTEE    -- high-confidence LoF classification (Tier A)
#   * CADD      -- general deleteriousness score
#   * REVEL     -- ensemble missense score (complements AlphaMissense)
#   * dbNSFP    -- pulls 8 in-silico predictors
#   * gnomAD v4 -- 807k samples, much better ancestry resolution than v2.0.1
#   * ClinVar   -- curated clinical assertions
# Annotation done on the cohort VCF (single VEP run, faster than per-sample)
# and split per sample afterwards.
# =============================================================================
rule r11_vep_annotate_cohort:
    input:
        vcf             = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz",
        tbi             = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz.tbi",
        ref             = reference_genome,
        spliceai_snv    = spliceai_snv,
        spliceai_indel  = spliceai_indel,
        alphamissense   = alphamissense_tsv,
        cadd_snv        = cadd_snv,
        cadd_indel      = cadd_indel,
        revel           = revel_tsv,
        dbnsfp          = dbnsfp_gz,
        gnomad          = gnomad_v4_exomes,
        clinvar         = clinvar_vcf,
        lof_pm       = "/tmp/annotation/vep/plugins/flat/LoF.pm",
        lof_src_pm   = "/tmp/annotation/vep/plugins/loftee_grch38/LoF.pm",
        lof_ancestor = "/tmp/annotation/vep/plugin_data/loftee_hg38/human_ancestor.fa.gz",
        lof_sql      = "/tmp/annotation/vep/plugin_data/loftee_hg38/loftee.sql",
        lof_gerp     = "/tmp/annotation/vep/plugin_data/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
    output:
        tsv     = "/tmp/fastq/11_vep/cohort.vep.tsv.gz",
        stats   = "/tmp/fastq/11_vep/cohort.vep.stats.html",
    params:
        cache_dir   = vep_cache_dir,
        plugin_dir  = vep_plugin_dir,
        loftee_dir  = loftee_dir,
        loftee_path_string  = loftee_path_string,
    threads: 6
    # singularity: "docker://ensemblorg/ensembl-vep:release_112.0"
    singularity: "/tmp/repo/sing/ensembl-vep-loftee-112.simg"
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --output_file STDOUT \
            --stats_file {output.stats} \
            --fasta {input.ref} \
            --cache --offline --dir_cache {params.cache_dir} \
            --dir_plugins {params.plugin_dir} \
            --assembly GRCh38 \
            --fork {threads} \
            --buffer_size 20000 \
            --force_overwrite \
            --tab --compress_output gzip \
            --symbol --canonical --biotype --hgvs --numbers \
            --mane --pick_allele_gene --pick_order mane_select,canonical,biotype \
            --sift b --polyphen b \
            --check_existing --no_check_alleles \
            --plugin SpliceAI,snv={input.spliceai_snv},indel={input.spliceai_indel},split_output=1 \
            --plugin AlphaMissense,file={input.alphamissense} \
            --plugin "LoF,{params.loftee_path_string},human_ancestor_fa:{input.lof_ancestor},conservation_file:{input.lof_sql},gerp_bigwig:{input.lof_gerp}" \
            --plugin CADD,{input.cadd_snv},{input.cadd_indel} \
            --plugin REVEL,{input.revel} \
            --plugin dbNSFP,{input.dbnsfp},MutationTaster_pred,PROVEAN_pred,\
MetaLR_pred,MetaRNN_pred,M-CAP_pred,PrimateAI_pred,\
ClinPred_pred,BayesDel_addAF_pred,\
            --custom file={input.gnomad},short_name=gnomADv4,format=vcf,type=exact,\
fields=AF%AF_nfe%AF_afr%AF_amr%AF_eas%AF_sas%AF_fin%AF_asj%nhomalt \
            --custom file={input.clinvar},short_name=ClinVar,format=vcf,type=exact,\
fields=CLNSIG%CLNREVSTAT%CLNDN%CLNDISDB \
        > {output.tsv}

        chmod a+wx /tmp/fastq/11_vep
        """

# rule r11_vep_annotate_cohort_test_bloody_loftee_plugin_arrgh:
#     input:
#         vcf          = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz",
#         tbi          = "/tmp/fastq/09_normalized/cohort.norm.vcf.gz.tbi",
#         ref          = reference_genome,

#         lof_pm       = "/tmp/annotation/vep/plugins/flat/LoF.pm",
#         lof_src_pm   = "/tmp/annotation/vep/plugins/loftee_grch38/LoF.pm",
#         lof_ancestor = "/tmp/annotation/vep/plugin_data/loftee_hg38/human_ancestor.fa.gz",
#         lof_sql      = "/tmp/annotation/vep/plugin_data/loftee_hg38/loftee.sql",
#         lof_gerp     = "/tmp/annotation/vep/plugin_data/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
#     output:
#         tsv = "/tmp/fastq/11_vep/test_loftee.vep.tsv",
#         ok  = "/tmp/fastq/11_vep/test_loftee.ok",
#     params:
#         cache_dir  = vep_cache_dir,
#         plugin_dir = vep_plugin_dir,
#         loftee_dir = loftee_dir,
#         loftee_path_string = loftee_path_string,
#     threads: 1
#     # singularity: "docker://ensemblorg/ensembl-vep:release_112.0"
#     singularity: "/tmp/repo/sing/ensembl-vep-loftee-112.simg"
#     shell:
#         r"""
#         set -euo pipefail

#         mkdir -p /tmp/fastq/11_vep

#         vep \
#             --input_file {input.vcf} \
#             --output_file {output.tsv} \
#             --fasta {input.ref} \
#             --cache --offline --dir_cache {params.cache_dir} \
#             --dir_plugins {params.plugin_dir} \
#             --assembly GRCh38 \
#             --force_overwrite \
#             --tab \
#             --no_stats \
#             --safe \
#             --chr chr21 \
#             --symbol --canonical --biotype --numbers \
#             --no_check_alleles \
#             --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,SYMBOL,LoF,LoF_filter,LoF_flags,LoF_info" \
#             --plugin "LoF,{params.loftee_path_string},human_ancestor_fa:{input.lof_ancestor},conservation_file:{input.lof_sql},gerp_bigwig:{input.lof_gerp}"

#         echo "Checking for LOFTEE columns..." >&2
#         grep -m1 '^#Uploaded_variation' {output.tsv} \
#             | tr '\t' '\n' \
#             | grep -qx 'LoF' || {{
#                 echo "ERROR: VEP output is missing LOFTEE column 'LoF'" >&2
#                 exit 1
#             }}
#         grep -m1 '^#Uploaded_variation' {output.tsv} \
#             | tr '\t' '\n' \
#             | grep -qx 'LoF_filter' || {{
#                 echo "ERROR: VEP output is missing LOFTEE column 'LoF_filter'" >&2
#                 exit 1
#             }}
#         grep -m1 '^#Uploaded_variation' {output.tsv} \
#             | tr '\t' '\n' \
#             | grep -qx 'LoF_flags' || {{
#                 echo "ERROR: VEP output is missing LOFTEE column 'LoF_flags'" >&2
#                 exit 1
#             }}
#         touch {output.ok}
#         chmod a+wx /tmp/fastq/11_vep
#         """

# =============================================================================
# r12 -- THE CRITICAL RULE: tier-based filtering of VEP output
# -----------------------------------------------------------------------------
# Replaces the original r14_filter_annotations.R script. Key changes:
#   * Variant-class-aware: LoF, splice, and missense handled separately
#   * Uses current-gen predictors (SpliceAI, AlphaMissense, LOFTEE)
#                 instead of SIFT/PolyPhen-2 from the original
#   * Uses gnomAD v4 AF with ancestry stratification
#   * Incorporates gene-level constraint (LOEUF) and disease-gene panels
#   * Produces three tiers aligned with what RNA-seq can empirically validate:
#       Tier A: DROP-testable (LoF + splice)       -- RNA confirmation expected
#       Tier B: Missense in constrained genes      -- RNA usually silent, rank only
#       Tier C: Any rare variant in SZ/NDD panels  -- hand-curate regardless
#   * Consequence deduplication: one row per variant (canonical / MANE-Select
#     transcript), not one row per VEP consequence
# =============================================================================
rule r12_tier_candidates:
    input:
        vep_tsv         = "/tmp/fastq/11_vep/cohort.vep.tsv.gz",
        per_sample_vcf  = "/tmp/fastq/10_per_sample_vcf/{sample}.vcf.gz",
        schema_genes    = schema_genes,
        bipex_genes     = bipex_genes,
        asc_genes       = asc_genes,
        ndd_genes       = ndd_genes,
        loeuf_table     = loeuf_table,
        artifact_sites  = "/tmp/qc/internal_artifact_sites.tsv",
    output:
        tier_a          = "/tmp/fastq/12_tiered/{sample}_tierA_lof_splice.tsv",
        tier_b          = "/tmp/fastq/12_tiered/{sample}_tierB_missense.tsv",
        tier_c          = "/tmp/fastq/12_tiered/{sample}_tierC_gene_panels.tsv",
        master          = "/tmp/fastq/12_tiered/{sample}_master.tsv",
    params:
        script          = "/tmp/scripts/lpb-exome-prioritisation-tier-candidates.py",
    singularity: "docker://quay.io/biocontainers/pandas:1.5.2"
    shell:
        """
        python {params.script} \
            --vep_tsv {input.vep_tsv} \
            --sample_vcf {input.per_sample_vcf} \
            --sample_id {wildcards.sample} \
            --schema_genes {input.schema_genes} \
            --bipex_genes {input.bipex_genes} \
            --asc_genes {input.asc_genes} \
            --ndd_genes {input.ndd_genes} \
            --loeuf_table {input.loeuf_table} \
            --artifact_sites {input.artifact_sites} \
            --out_tier_a {output.tier_a} \
            --out_tier_b {output.tier_b} \
            --out_tier_c {output.tier_c} \
            --out_master {output.master}

        chmod a+wx /tmp/fastq/12_tiered
        """
