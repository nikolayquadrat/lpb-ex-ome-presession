# =============================================================================
# RNA-seq alignment pipeline, GTEx V11-compatible
# =============================================================================
# Purpose: produce STAR-aligned BAM files for the 8-donor / ~36-sample brain
# RNA-seq cohort, using exactly the alignment parameters and reference files
# of the GTEx V11 release (https://gtexportal.org/home/downloads/adult-gtex/
# bulk_tissue_expression). This makes the BAMs and gene-level read counts
# directly mergeable with the GTEx V11 expression matrices for DROP/OUTRIDER
# reference panel expansion.
#
# Reference setup (GTEx V11):
#   - genome:    Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta (TOPMed ref)
#   - annotation: GENCODE v47 (gencode.v47.GRCh38.annotation.gtf)
#   - STAR:      v2.7.11b (use this exact version for full reproducibility)
#   - RSEM:      v1.3.3
#   - sjdbOverhang: depends on YOUR read length (NOT 75, unless your reads are
#                   2x76bp like GTEx). Set READ_LENGTH below to match your data.
#
# Outputs per sample:
#   - 03_bam_star/{sample}/Aligned.sortedByCoord.out.patched.md.bam
#       (the canonical GTEx-style output: STAR-aligned, coordinate-sorted,
#        duplicate-marked; this is the one to feed to DROP)
#   - 03_bam_star/{sample}/Aligned.toTranscriptome.out.bam
#       (transcriptome-coordinate BAM, for downstream RSEM if desired)
#   - 03_bam_star/{sample}/{sample}.ReadsPerGene.out.tab
#       (STAR's gene-level counts, mergeable with GTEx V11 count matrix)
#   - 03_bam_star/{sample}/{sample}.SJ.out.tab
#       (splice junction file, used by FRASER)
#
# Differences from your previous pipeline:
#   - STAR parameters now match GTEx V11 exactly (twopassMode, chim*, etc.)
#   - Genome FASTA is the GTEx-specific one (no ALT/HLA/decoy contigs) - this
#     matters because STAR auto-prefers ALT contigs when present, which would
#     produce different alignments than GTEx
#   - Reference annotation is GENCODE v47, not Ensembl 109
#   - STAR index uses --sjdbOverhang appropriate for YOUR reads
#   - Picard MarkDuplicates is run after sorting (GTEx pattern)
#   - Salmon/HISAT2/contamination paths from your earlier pipeline are removed
#     (this Snakefile is DROP-focused; keep the other one for general work)
# =============================================================================

import pandas as pd
import os

# -----------------------------------------------------------------------------
# User-tunable parameters
# -----------------------------------------------------------------------------

# Read length in your FASTQs. CRITICAL: STAR's --sjdbOverhang must be
# (read_length - 1). GTEx itself uses 75 because their reads are 2x76bp.
# To check: zcat your_R1.fq.gz | head -2 | tail -1 | wc -c   (subtract 1 for newline)
READ_LENGTH = 150   # <-- ADJUST TO ACTUAL READ LENGTH
SJDB_OVERHANG = READ_LENGTH - 1

# Sample table: same format as your previous pipeline, with
# columns: name, path, dataset1, dataset2, extension, include_in_analysis
SAMPLE_FILE = "/tmp/data/00_additional_files/sample_data.txt"
sample_data = pd.read_csv(SAMPLE_FILE, sep="\t")
sample_data = sample_data[sample_data["include_in_analysis"] == 1]

# Filter out samples whose forward FASTQ is not present on disk.
# This makes the pipeline robust to a master sample table that lists more
# samples than are currently available (e.g., when sequencing is rolling in
# in batches, or some samples have been moved/archived). The warning lists
# every dropped sample by name so a typo'd filename still surfaces clearly.
def _fwd_fastq_path(row):
    return f"/tmp/data/{row['path']}/{row['dataset1']}.{row['extension']}"

_present_mask = sample_data.apply(
    lambda r: os.path.exists(_fwd_fastq_path(r)), axis=1)

if (~_present_mask).any():
    missing = sample_data.loc[~_present_mask, "name"].tolist()
    print("=" * 72)
    print(f"WARNING: {len(missing)} sample(s) listed in {SAMPLE_FILE}")
    print("with include_in_analysis=1 but their forward FASTQ file is missing")
    print("on disk. These samples will be SILENTLY SKIPPED in this run.")
    print("If this is unintended (e.g., typo in 'path' or 'dataset1'),")
    print("fix the sample table or restore the missing files before re-running.")
    print()
    print("Skipped samples (and their expected R1 paths):")
    for name in missing:
        row = sample_data[sample_data["name"] == name].iloc[0]
        print(f"  - {name}: {_fwd_fastq_path(row)}")
    print("=" * 72)
    sample_data = sample_data.loc[_present_mask].reset_index(drop=True)

samples = sample_data["name"].tolist()
print(f"Aligning {len(samples)} samples")


# Reference paths.
# IMPORTANT: download these files to the indicated locations before running.
# See "Reference download instructions" at the bottom of this file.
REF_DIR     = "/tmp/data/00_additional_files/gtex_v11_refs"
GENOME_FA   = REF_DIR + "/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"
GTF         = REF_DIR + "/gencode.v47.GRCh38.annotation.gtf"
STAR_INDEX  = REF_DIR + f"/star_index_oh{SJDB_OVERHANG}"
RSEM_INDEX  = REF_DIR + "/rsem_reference"

# Container versions (matching GTEx V11 release where possible)
STAR_CONTAINER   = "docker://quay.io/biocontainers/star:2.7.11b--h43eeafb_0"
SAMTOOLS_CONTAINER = "docker://biocontainers/samtools:v1.9-4-deb_cv1"
PICARD_CONTAINER = "docker://broadinstitute/picard:2.27.5"
TRIM_CONTAINER   = "docker://mskaccess/trim_galore" 

def get_fastq_path(wildcards):
    """Same logic as your existing get_fastq_path(): paired or single-end."""
    row = sample_data[sample_data.name == wildcards.sample].iloc[0]
    fwd = f"/tmp/data/{row['path']}/{row['dataset1']}.{row['extension']}"
    if str(row['dataset2']).strip().lower() == "nan":
        return [fwd, fwd]   # double up for SE
    rev = f"/tmp/data/{row['path']}/{row['dataset2']}.{row['extension']}"
    return [fwd, rev]


# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------

rule all:
    input:
        # STAR BAMs ready for DROP, plus splice-junction files
        expand("/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam", sample=samples),
        # expand("/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai", sample=samples),
        expand("/tmp/data/03_bam_star/{sample}/{sample}.SJ.out.tab", sample=samples),
        expand("/tmp/data/03_bam_star/{sample}/{sample}.ReadsPerGene.out.tab", sample=samples),
        # Mapping stat summary
        "/tmp/data/03_bam_star/00_mapping_stat/mapping_stat.txt",
        # QC: per-sample Picard CollectRnaSeqMetrics + cohort summary table
        # expand("/tmp/data/04_qc/{sample}/{sample}.rnaseq_metrics.txt", sample=samples),
        "/tmp/data/04_qc/00_qc_summary.tsv",
    shell: "echo 'GTEx-V11-compatible alignment + QC complete.'"

# -----------------------------------------------------------------------------
# Trimming (kept identical to your existing pipeline for consistency)
# -----------------------------------------------------------------------------
# Note: GTEx itself does NOT trim. They feed raw FASTQ to STAR. Trimming is
# generally OK for STAR but produces slightly different alignments than GTEx.
# If you want maximum GTEx fidelity, set USE_TRIMMING = False below and edit
# rule r03d_STAR_mapping_GTEx to read directly from the un-trimmed FASTQ.
# For DROP use, trimmed-vs-untrimmed does not significantly affect outlier
# detection, so trimming is left enabled by default.

USE_TRIMMING = True

rule r02a_trim_galore:
    input:
        get_fastq_path,
    output:
        R1 = "/tmp/data/02_trimmed_fastq/{sample}_R1_val_1.fq.gz",
    singularity: TRIM_CONTAINER
    shell:
        r"""
        mkdir -p /tmp/data/02_trimmed_fastq/temp/{wildcards.sample}
        n=$(perl -e 'if (qq|{input}| =~ /^(\S+)\s\1$/) {{print 1}} else {{print 0}}')
        if (($n==1)); then
            trim_galore {input[0]} -o /tmp/data/02_trimmed_fastq/temp/{wildcards.sample}
            sleep 5
            file_to_move=$(find /tmp/data/02_trimmed_fastq/temp/{wildcards.sample} \
                             -type f -name "*trimmed\.fq\.gz*")
            mv "$file_to_move" /tmp/data/02_trimmed_fastq/{wildcards.sample}_R1_val_1.fq.gz
        else
            trim_galore --paired {input[0]} {input[1]} \
                -o /tmp/data/02_trimmed_fastq/temp/{wildcards.sample}
            sleep 5
            file_to_move1=$(find /tmp/data/02_trimmed_fastq/temp/{wildcards.sample} \
                              -type f -name "*_val_1\.fq\.gz*")
            mv "$file_to_move1" /tmp/data/02_trimmed_fastq/{wildcards.sample}_R1_val_1.fq.gz
            file_to_move2=$(find /tmp/data/02_trimmed_fastq/temp/{wildcards.sample} \
                              -type f -name "*_val_2\.fq\.gz*")
            mv "$file_to_move2" /tmp/data/02_trimmed_fastq/{wildcards.sample}_R2_val_2.fq.gz
        fi
        rm -rf /tmp/data/02_trimmed_fastq/temp/{wildcards.sample}
        """


# -----------------------------------------------------------------------------
# STAR index, GTEx V11 parameters
# -----------------------------------------------------------------------------

rule r03a_STAR_index:
    """
    Build the STAR index using GENCODE v47 annotation and sjdbOverhang
    matching YOUR read length (read_length - 1).
    GTEx itself uses sjdbOverhang=75 because their reads are 2x76bp.
    """
    input:
        genome = GENOME_FA,
        gtf    = GTF,
    output:
        checkpoint = STAR_INDEX + "/SA",
    params:
        index_dir = STAR_INDEX,
        sjdb_overhang = SJDB_OVERHANG,
    threads: 16
    singularity: STAR_CONTAINER
    shell:
        """
        mkdir -p {params.index_dir}
        STAR \
            --runMode genomeGenerate \
            --genomeDir {params.index_dir} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdb_overhang} \
            --runThreadN {threads}
        """


# -----------------------------------------------------------------------------
# STAR alignment with GTEx V11 parameters
# -----------------------------------------------------------------------------

rule r03d_STAR_mapping_GTEx:
    """
    STAR alignment using exactly the parameter set from the GTEx V11 pipeline
    (Broad's run_STAR.py defaults). The flags below come verbatim from
    https://github.com/broadinstitute/gtex-pipeline/blob/master/rnaseq/src/run_STAR.py
    Together with --twopassMode Basic, these are what GTEx V11 used to
    produce the count matrices on the GTEx portal. Reproducing them exactly
    means our gene-level counts can be merged into a unified expression
    matrix with GTEx samples for DROP/OUTRIDER reference panel expansion.

    Notable parameters:
      --twopassMode Basic              : two-pass alignment for better novel-junction detection
      --outFilterMultimapNmax 20       : up to 20 multi-mapping locations per read
      --outFilterMismatchNoverLmax 0.1 : <=10% mismatches per read length
      --outFilterType BySJout          : only retain junctions supported in 2nd pass
      --quantMode TranscriptomeSAM GeneCounts : output transcriptome BAM + STAR gene counts
      --chimSegmentMin 15              : enable chimeric detection (TOPMed/GTEx convention)
      --outSAMtype BAM SortedByCoordinate : GTEx writes a coordinate-sorted BAM directly
      --outSAMunmapped Within          : keep unmapped reads in the BAM (GTEx convention)
    """
    input:
        R1 = "/tmp/data/02_trimmed_fastq/{sample}_R1_val_1.fq.gz" if USE_TRIMMING
              else lambda wc: get_fastq_path(wc)[0],
        index_checkpoint = ancient(STAR_INDEX + "/SA"),
    output:
        bam_sorted = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bam_transcripts = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.toTranscriptome.out.bam",
        gene_counts = "/tmp/data/03_bam_star/{sample}/{sample}.ReadsPerGene.out.tab",
        sj_tab = "/tmp/data/03_bam_star/{sample}/{sample}.SJ.out.tab",
        chimeric = "/tmp/data/03_bam_star/{sample}/{sample}.Chimeric.out.junction",
        log_final = "/tmp/data/03_bam_star/{sample}/{sample}.Log.final.out",
    params:
        out_prefix = lambda wc: f"/tmp/data/03_bam_star/{wc.sample}/{wc.sample}.",
        out_dir    = lambda wc: f"/tmp/data/03_bam_star/{wc.sample}/",
        tmp_dir    = lambda wc: f"/tmp/data/03_bam_star/{wc.sample}/_STAR_tmp",
        star_index = STAR_INDEX,
        # Read group line for downstream tools. Replace SM tag with sample id;
        # ID is set to the sample id as well for simplicity. PL=ILLUMINA, LB=lib1.
        rg_line = lambda wc: f"ID:{wc.sample} SM:{wc.sample} PL:ILLUMINA LB:lib1",
        r2_path = lambda wc: f"/tmp/data/02_trimmed_fastq/{wc.sample}_R2_val_2.fq.gz",
    threads: 8
    singularity: STAR_CONTAINER
    shell:
        r"""
        # STAR refuses to write into a non-empty tmp dir; clean before running
        rm -rf {params.tmp_dir} || true
        mkdir -p {params.out_dir}

        # Decide PE vs SE based on whether the R2 file exists
        if [ -f "{params.r2_path}" ]; then
            R2_ARG="{params.r2_path}"
        else
            R2_ARG=""
        fi

        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.star_index} \
            --twopassMode Basic \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.1 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFilterType BySJout \
            --outFilterScoreMinOverLread 0.33 \
            --outFilterMatchNminOverLread 0.33 \
            --limitSjdbInsertNsj 1200000 \
            --readFilesIn {input.R1} $R2_ARG \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.out_prefix} \
            --outTmpDir {params.tmp_dir} \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs None \
            --alignSoftClipAtReferenceEnds Yes \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --genomeLoad NoSharedMemory \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --chimOutType Junctions WithinBAM SoftClip \
            --chimMainSegmentMultNmax 1 \
            --outSAMattributes NH HI AS nM NM ch \
            --outSAMattrRGline {params.rg_line}

        rm -rf {params.tmp_dir} || true
        """


# -----------------------------------------------------------------------------
# Index BAM (samtools)
# -----------------------------------------------------------------------------

rule r03e_index_bam:
    input:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bai = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
    threads: 4
    singularity: SAMTOOLS_CONTAINER
    shell:
        "samtools index -@ {threads} {input.bam} {output.bai}"


# -----------------------------------------------------------------------------
# Picard MarkDuplicates
# -----------------------------------------------------------------------------
# GTEx applies MarkDuplicates AFTER coordinate-sorting from STAR. The output
# is the canonical "<sample>.Aligned.sortedByCoord.out.patched.md.bam" -
# this is the file that gets fed to RNA-SeQC and RSEM in the GTEx pipeline,
# and it is also the file DROP MAE / OUTRIDER / FRASER expect as input.
#
# Note: bamsync (which carries QC flags from the original GTEx CRAM) is NOT
# applied here because we are starting from FASTQ, not from existing BAM.
# The "patched" suffix is kept in the filename only for naming consistency
# with the GTEx convention.

rule r03f_markduplicates:
    input:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
    output:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        metrics = "/tmp/data/03_bam_star/{sample}/{sample}.markdup_metrics.txt",
    threads: 4
    singularity: PICARD_CONTAINER
    shell:
        """
        java -jar /usr/picard/picard.jar MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            ASSUME_SORT_ORDER=coordinate \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
            VALIDATION_STRINGENCY=LENIENT
        """


rule r03g_index_md_bam:
    input:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
    output:
        bai = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai",
    threads: 4
    singularity: SAMTOOLS_CONTAINER
    shell:
        "samtools index -@ {threads} {input.bam} {output.bai}"


# -----------------------------------------------------------------------------
# Aggregate mapping stats
# -----------------------------------------------------------------------------

rule r03h_mapping_stat:
    input:
        log_files = expand("/tmp/data/03_bam_star/{sample}/{sample}.Log.final.out",
                           sample=samples),
    output:
        "/tmp/data/03_bam_star/00_mapping_stat/mapping_stat.txt",
    run:
        with open(output[0], "w") as fout:
            fout.write("Sample\tInput_reads\tUniquely_mapped\tUniquely_mapped_pct\t"
                       "Multi_mapped\tMulti_mapped_pct\tToo_many_loci\tUnmapped\n")
            for log_file in input.log_files:
                sample = os.path.basename(log_file).replace(".Log.final.out", "")
                metrics = {"input_reads": "", "unique": "", "unique_pct": "",
                           "multi": "", "multi_pct": "", "too_many": "", "unmapped": ""}
                with open(log_file) as f:
                    for line in f:
                        line = line.strip()
                        if "Number of input reads" in line:
                            metrics["input_reads"] = line.split("|")[1].strip()
                        elif "Uniquely mapped reads number" in line:
                            metrics["unique"] = line.split("|")[1].strip()
                        elif "Uniquely mapped reads %" in line:
                            metrics["unique_pct"] = line.split("|")[1].strip()
                        elif "Number of reads mapped to multiple loci" in line:
                            metrics["multi"] = line.split("|")[1].strip()
                        elif "% of reads mapped to multiple loci" in line:
                            metrics["multi_pct"] = line.split("|")[1].strip()
                        elif "% of reads mapped to too many loci" in line:
                            metrics["too_many"] = line.split("|")[1].strip()
                        elif "% of reads unmapped: too short" in line:
                            metrics["unmapped"] = line.split("|")[1].strip()
                fout.write(f"{sample}\t{metrics['input_reads']}\t"
                           f"{metrics['unique']}\t{metrics['unique_pct']}\t"
                           f"{metrics['multi']}\t{metrics['multi_pct']}\t"
                           f"{metrics['too_many']}\t{metrics['unmapped']}\n")


# -----------------------------------------------------------------------------
# Picard CollectRnaSeqMetrics: RNA-seq QC including 3' bias
# -----------------------------------------------------------------------------
# Picard CollectRnaSeqMetrics reports several metrics critical for catching
# RNA degradation and library-prep issues:
#
#   - MEDIAN_3PRIME_BIAS / MEDIAN_5PRIME_BIAS / MEDIAN_5PRIME_TO_3PRIME_BIAS:
#     coverage uniformity along transcripts. A high 3' bias (>0.5 typically,
#     or much higher than the cohort median) indicates RNA degradation -
#     polyA-selected libraries lose 5' coverage as RNA degrades because
#     fragments without a polyA tail are lost preferentially.
#   - PCT_RIBOSOMAL_BASES: rRNA contamination. Should be <5% in a well-
#     prepared polyA-selected library; anything >10% indicates poor depletion.
#   - PCT_CODING_BASES / PCT_UTR_BASES / PCT_INTRONIC_BASES / PCT_INTERGENIC_BASES:
#     where reads come from. Low coding+UTR fraction (<60% combined) suggests
#     genomic DNA contamination or unspliced pre-mRNA dominance, both of which
#     can confound expression analysis.
#   - PCT_MRNA_BASES: percent of reads in coding+UTR (the "useful" fraction
#     for differential expression).
#
# Why this matters for cell-type composition: RNA degradation can selectively
# affect long transcripts (e.g., GFAP at ~10 kb) more than short ones,
# producing apparent "depletion" of cell types whose markers happen to be
# long transcripts. A high 3' bias on a sample with apparently missing
# astrocyte signal (GFAP, AQP4 are both long) would explain the pattern
# without needing a real biological loss-of-cell-type story.
# -----------------------------------------------------------------------------

rule r04a_make_refflat:
    """
    Convert GENCODE v47 GTF to UCSC refFlat format, which is what Picard
    CollectRnaSeqMetrics requires. We use ucsc-gtfToGenePred (with
    -genePredExt to retain gene_id) and then reformat columns to refFlat.

    refFlat columns: geneName, name (transcript), chrom, strand, txStart,
    txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds.
    genePredExt has these in different order; the awk fixes that.
    """
    input:
        gtf = GTF,
    output:
        refflat = REF_DIR + "/gencode.v47.GRCh38.refFlat.txt",
    singularity: "docker://quay.io/biocontainers/ucsc-gtftogenepred:482--h0b57e2e_0"
    shell:
        r"""
        # gtfToGenePred -genePredExt produces 15-column genePredExt; refFlat
        # wants 11 columns with geneName as the FIRST column. We swap and trim.
        gtfToGenePred -genePredExt -geneNameAsName2 {input.gtf} /dev/stdout \
            | awk -v OFS='\t' '{{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}' \
            > {output.refflat}
        """


rule r04b_make_rrna_intervals:
    """
    Build a Picard interval-list of ribosomal RNA loci from the GENCODE GTF.
    Picard CollectRnaSeqMetrics uses this to compute PCT_RIBOSOMAL_BASES.

    Selects gene_type "rRNA" and "Mt_rRNA" entries. The interval-list format
    is a SAM-like header (from the BAM) plus one BED-like line per region.
    We extract the header from any sample BAM (they all share the same
    reference, so any one works).
    """
    input:
        gtf = GTF,
        # Borrow the header from the first sample's STAR BAM. The header
        # contains @SQ lines that identify the reference contigs - required
        # for a valid interval-list file.
        ref_bam = expand(
            "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
            sample=samples[:1]),
    output:
        rrna_list = REF_DIR + "/gencode.v47.GRCh38.rRNA.interval_list",
    singularity: SAMTOOLS_CONTAINER
    shell:
        r"""
        # 1. Header from the BAM (only @HD and @SQ lines; @PG lines from STAR
        #    confuse Picard's interval-list parser)
        samtools view -H {input.ref_bam[0]} \
            | grep -E '^@HD|^@SQ' > {output.rrna_list}

        # 2. rRNA gene loci: from the GTF, gene_type "rRNA" or "Mt_rRNA",
        #    feature == "gene", emit chrom, start, end, strand, gene_id.
        #    Picard interval-list uses 1-based coordinates so no -1 on start.
        awk 'BEGIN{{OFS="\t"}} \
            $3=="gene" && /gene_type "rRNA"|gene_type "Mt_rRNA"/ {{ \
                match($0, /gene_id "[^"]+"/); \
                gid=substr($0, RSTART+9, RLENGTH-10); \
                print $1, $4, $5, $7, gid \
            }}' {input.gtf} >> {output.rrna_list}
        """


rule r04c_picard_rnaseq_metrics:
    """
    Per-sample CollectRnaSeqMetrics. Strand-specificity is set to NONE
    because we don't know your library protocol - if your library is
    strand-specific (most modern protocols are), update this to
    SECOND_READ_TRANSCRIPTION_STRAND (TruSeq dUTP/Stranded mRNA) or
    FIRST_READ_TRANSCRIPTION_STRAND. Wrong strand setting inflates the
    PCT_CORRECT_STRAND_READS metric but does NOT affect 3'/5' bias - so
    NONE is safe for our primary QC purpose.
    """
    input:
        bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        bai = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai",
        refflat = REF_DIR + "/gencode.v47.GRCh38.refFlat.txt",
        rrna_list = REF_DIR + "/gencode.v47.GRCh38.rRNA.interval_list",
    output:
        metrics = "/tmp/data/04_qc/{sample}/{sample}.rnaseq_metrics.txt",
    threads: 2
    singularity: PICARD_CONTAINER
    shell:
        """
        mkdir -p /tmp/data/04_qc/{wildcards.sample}
        java -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
            I={input.bam} \
            O={output.metrics} \
            REF_FLAT={input.refflat} \
            RIBOSOMAL_INTERVALS={input.rrna_list} \
            STRAND_SPECIFICITY=NONE \
            VALIDATION_STRINGENCY=LENIENT \
            ASSUME_SORTED=true
        """


rule r04d_qc_summary:
    """
    Aggregate per-sample CollectRnaSeqMetrics outputs into one table, flagging
    outliers. Picard's metrics file has two main sections:
      - METRICS section: one line of values per sample
      - HISTOGRAM section: normalized coverage by percentile of transcript
        position (0=5'-end, 100=3'-end); we don't aggregate this here but
        leave the per-sample file available for plotting.

    The output table makes it easy to spot:
      - Samples with high MEDIAN_3PRIME_BIAS (RNA degradation candidates)
      - Samples with high PCT_RIBOSOMAL_BASES (rRNA depletion failure)
      - Samples with low PCT_MRNA_BASES (DNA contamination or pre-mRNA
        dominance)
      - Samples with low PCT_USABLE_BASES (general QC poor)
    """
    input:
        metrics_files = expand("/tmp/data/04_qc/{sample}/{sample}.rnaseq_metrics.txt",
                                sample=samples),
    output:
        summary = "/tmp/data/04_qc/00_qc_summary.tsv",
    run:
        import os, statistics

        # Columns of interest from Picard's METRICS section
        cols_of_interest = [
            "PF_BASES", "PF_ALIGNED_BASES",
            "PCT_RIBOSOMAL_BASES",
            "PCT_CODING_BASES", "PCT_UTR_BASES",
            "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES",
            "PCT_MRNA_BASES", "PCT_USABLE_BASES",
            "MEDIAN_CV_COVERAGE",
            "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS",
            "MEDIAN_5PRIME_TO_3PRIME_BIAS",
            "PCT_CORRECT_STRAND_READS",
        ]

        rows = {}
        for f in input.metrics_files:
            sample = os.path.basename(f).replace(".rnaseq_metrics.txt", "")
            with open(f) as fh:
                lines = fh.readlines()
            # Find the METRICS header line (starts with "PF_BASES")
            header_idx = None
            for i, line in enumerate(lines):
                if line.startswith("PF_BASES"):
                    header_idx = i
                    break
            if header_idx is None:
                continue
            header = lines[header_idx].strip().split("\t")
            values = lines[header_idx + 1].strip().split("\t")
            row = dict(zip(header, values))
            rows[sample] = {c: row.get(c, "") for c in cols_of_interest}

        # Compute cohort medians for the bias metrics so we can flag outliers
        def floats(col):
            return [float(rows[s][col]) for s in rows
                    if rows[s].get(col, "") not in ("", "?", "NA")]
        bias_3p = floats("MEDIAN_3PRIME_BIAS")
        rrna    = floats("PCT_RIBOSOMAL_BASES")
        med_3p_bias = statistics.median(bias_3p) if bias_3p else 0
        mad_3p_bias = statistics.median([abs(x - med_3p_bias) for x in bias_3p]) \
                      if bias_3p else 0
        med_rrna    = statistics.median(rrna) if rrna else 0

        # Write summary
        with open(output.summary, "w") as fout:
            fout.write("sample\t" + "\t".join(cols_of_interest) +
                       "\tflag_3prime_bias_high\tflag_rrna_high\tflag_mrna_low\n")
            for sample in sorted(rows):
                row = rows[sample]
                # Flagging logic (heuristic; tune to your cohort):
                #   3' bias > median + 3*MAD -> degradation candidate
                #   rRNA > 10% -> rRNA depletion failure
                #   PCT_MRNA_BASES < 0.6 -> DNA / pre-mRNA contamination
                try:
                    b3 = float(row.get("MEDIAN_3PRIME_BIAS", "0") or 0)
                    flag_3p = "YES" if (mad_3p_bias > 0 and
                                        b3 > med_3p_bias + 3*mad_3p_bias) else ""
                except ValueError:
                    flag_3p = "?"
                try:
                    rb = float(row.get("PCT_RIBOSOMAL_BASES", "0") or 0)
                    flag_rrna = "YES" if rb > 0.10 else ""
                except ValueError:
                    flag_rrna = "?"
                try:
                    pm = float(row.get("PCT_MRNA_BASES", "0") or 0)
                    flag_mrna = "YES" if pm < 0.60 else ""
                except ValueError:
                    flag_mrna = "?"

                vals = [row.get(c, "") for c in cols_of_interest]
                fout.write(f"{sample}\t" + "\t".join(vals) +
                           f"\t{flag_3p}\t{flag_rrna}\t{flag_mrna}\n")


# =============================================================================
# Reference download instructions (run ONCE before the pipeline)
# =============================================================================
# 1. GENOME FASTA (GTEx-specific: no ALT, no HLA, no decoy contigs)
#    This file is on Google Cloud Storage:
#      gsutil cp gs://gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
#                /tmp/data/00_additional_files/gtex_v11_refs/
#    (gsutil is from the gcloud SDK; alternatively use rclone or the public
#    URL if you don't have a GCP account.)
#
#    If you cannot use gsutil, the same FASTA is served by the Broad's TOPMed
#    public bucket:
#      wget https://storage.googleapis.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
#
# 2. GENCODE v47 GTF
#      cd /tmp/data/00_additional_files/gtex_v11_refs/
#      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
#      gunzip gencode.v47.annotation.gtf.gz
#      mv gencode.v47.annotation.gtf gencode.v47.GRCh38.annotation.gtf
#
# 3. (Optional) Pre-built STAR index from GTEx
#    Building from scratch takes ~30 min on 16 cores. If you want exact
#    parity with GTEx, you can also download their pre-built index:
#      gsutil cp -r gs://gtex-resources/references/star_index_oh75 \
#                /tmp/data/00_additional_files/gtex_v11_refs/
#    BUT note that this index is for sjdbOverhang=75 only, so it only matches
#    if your reads are 2x76bp. For other read lengths, build your own.
#
# 4. (Optional, for DROP-OUTRIDER reference panel) GTEx V11 gene counts
#    Used downstream to merge with your samples for OUTRIDER:
#      Download from https://gtexportal.org/home/downloads/adult-gtex/
#      bulk_tissue_expression -- pick the
#      "GTEx_Analysis_v11_RNASeQC_2.4.3_RNAseQC_*_reads.gct.gz" file.
#      Filter to the brain-cortex-related tissues (Brain - Frontal Cortex BA9,
#      Brain - Anterior cingulate cortex BA24, etc.) for OUTRIDER reference.
#
# =============================================================================
# Sample table format:
# =============================================================================
# Tab-separated, columns: name path dataset1 dataset2 extension include_in_analysis
#   name                = unique sample id
#   path                = subdirectory under /tmp/data/ where FASTQs live
#   dataset1            = R1 file basename (without extension)
#   dataset2            = R2 file basename (or "NaN" for SE)
#   extension           = "fq.gz" or "fastq.gz"
#   include_in_analysis = 1 (process) or 0 (skip)
#
# =============================================================================
