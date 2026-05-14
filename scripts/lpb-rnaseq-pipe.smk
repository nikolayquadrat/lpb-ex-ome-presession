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
import sys

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
# bwa-mem2 + samtools 1.13 mulled image. Both binaries are on PATH inside,
# so r05b can stream from samtools fastq -> bwa-mem2 mem -> samtools view
# in a single container invocation.
BWA_CONTAINER = (
    "docker://quay.io/biocontainers/"
    "mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:"
    "c5b8c4b7735290369693e2b63cfc1ea0732fde07-0"
)

# BBTools (bbduk.sh in particular). Used by r05c_entropy_filter to drop
# low-complexity / low-entropy reads (poly-A tails, simple repeats,
# AT-rich runs) before they're aligned to the contamination reference.
# Removing these reads at the FASTQ level is the most effective single
# intervention for reducing cross-mapping false positives, since
# low-complexity reads align cleanly to many references regardless of
# bwa-mem2 stringency settings.
BBMAP_CONTAINER = "docker://quay.io/biocontainers/bbmap:39.06--h92535d8_0"

# bedtools, used by r05a0_mask_conserved_regions to hard-mask conserved
# rRNA/tRNA regions in the contamination reference. These regions are
# extremely sequence-conserved across species (bacterial 16S/23S has
# stretches 100% identical to human rRNA, mitochondrial rRNA reads
# cross-map to apicomplexan plastid rRNA, etc.), so masking them at the
# reference level prevents host-derived rRNA reads from generating
# spurious contamination hits regardless of alignment stringency.
BEDTOOLS_CONTAINER = "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"

# NCBI datasets CLI container, used by r05_prep_contamination_refs when
# the user supplies species.txt and asks the pipeline to build the
# contamination reference on the fly. Ships `datasets` and `unzip`.
# DATASETS_CONTAINER = "docker://quay.io/biocontainers/ncbi-datasets-cli:14.26.0"
DATASETS_CONTAINER = "docker://biocontainers/ncbi-datasets-cli:16.22.1_cv1"

# -----------------------------------------------------------------------------
# Optional contamination screen
# -----------------------------------------------------------------------------
# Three modes, chosen by which files are present in /tmp/data/00_additional_files/contamination/:
#
#   (A) PREBUILT: if contamination.fa + contamination.gtf + gtf_seqnames.tsv
#       are all present and non-empty, they're used as-is.
#
#   (B) ON-DEMAND: if (A)'s files are not all present but species.txt IS
#       present, the pipeline downloads each species' genome+GTF from NCBI
#       using the `datasets` CLI tool (mirroring the user-supplied
#       download_species.sh pattern), concatenates the FASTAs into
#       contamination.fa, the GTFs into contamination.gtf, and emits a
#       gtf_seqnames.tsv that maps each FASTA seqname (NCBI accession) to
#       a species_id of the form <name>_<assembly_accession>.
#
#   (C) DISABLED: neither (A) nor (B) is satisfied. r05* rules are skipped
#       and the QC summary omits contamination columns.
#
# species.txt format (tab-separated, no header):
#   species_name<TAB>assembly_accession
# e.g.:
#   mesomycoplasma_hyorhinis<TAB>GCF_001705605.1
#
# The contig naming inside the auto-built contamination.fa preserves the
# native NCBI seqnames (e.g. ">NC_010163.1 Acholeplasma laidlawii PG-8A").
# r05b's awk uses the gtf_seqnames.tsv to attribute reads to species.
#
# Alignment parameters (r05b): bwa-mem2 mem is run with conservative
# settings to reduce false positives from rRNA homologs etc:
#   -T 60  : minimum alignment score (default 30) -- requires ~60bp of
#            near-perfect identity for a read to be reported
#   -k 25  : minimum seed length (default 19)    -- longer seeds reject
#            short noisy matches
#   -r 2.0 : re-seed threshold (default 1.5)     -- less aggressive
#            re-seeding; fewer marginal alignments
# Combined with the MAPQ>=20 post-filter, the result is much closer to
# "reads that genuinely come from this organism" than the defaults.
CONTAM_DIR          = "/tmp/data/00_additional_files/contamination"
CONTAM_FA           = f"{CONTAM_DIR}/contamination.fa"
CONTAM_GTF          = f"{CONTAM_DIR}/contamination.gtf"
CONTAM_SEQNAMES     = f"{CONTAM_DIR}/gtf_seqnames.tsv"
CONTAM_SPECIES_LIST = f"{CONTAM_DIR}/species.txt"

_prebuilt = all(
    os.path.exists(p) and os.path.getsize(p) > 0
    for p in (CONTAM_FA, CONTAM_GTF, CONTAM_SEQNAMES)
)
_have_species_list = (
    os.path.exists(CONTAM_SPECIES_LIST)
    and os.path.getsize(CONTAM_SPECIES_LIST) > 0
)
CONTAMINATION_ENABLED  = _prebuilt or _have_species_list
CONTAMINATION_AUTOBUILD = (not _prebuilt) and _have_species_list

if CONTAMINATION_ENABLED:
    if CONTAMINATION_AUTOBUILD:
        print(f"[contamination screen] enabled (auto-build mode)",
              file=sys.stderr)
        print(f"[contamination screen] species list: {CONTAM_SPECIES_LIST}",
              file=sys.stderr)
        print(f"[contamination screen] will download genomes from NCBI and "
              f"build {CONTAM_FA} / {CONTAM_GTF} / {CONTAM_SEQNAMES}",
              file=sys.stderr)
    else:
        print(f"[contamination screen] enabled (pre-built files)",
              file=sys.stderr)
        print(f"[contamination screen] FASTA: {CONTAM_FA}", file=sys.stderr)
        print(f"[contamination screen] seqname map: {CONTAM_SEQNAMES}",
              file=sys.stderr)
else:
    print(f"[contamination screen] disabled (need either pre-built "
          f"contamination.fa+gtf+gtf_seqnames.tsv or a species.txt in "
          f"{CONTAM_DIR})", file=sys.stderr)

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
# Trimming
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
        R1 = temp("/tmp/data/02_trimmed_fastq/{sample}_R1_val_1.fq.gz"),
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
    Per-sample CollectRnaSeqMetrics. Strand-specificity could be set to NONE
    update this to SECOND_READ_TRANSCRIPTION_STRAND (TruSeq dUTP/Stranded mRNA) or
    FIRST_READ_TRANSCRIPTION_STRAND. Wrong strand setting inflates the
    PCT_CORRECT_STRAND_READS metric but does NOT affect 3'/5' bias - so
    NONE is safe for primary QC purpose.
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
            STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
            VALIDATION_STRINGENCY=LENIENT \
            ASSUME_SORTED=true
        """

# =============================================================================
# Optional contamination screen
# =============================================================================
# Rules:
#   r05_prep_contamination_refs : (only when CONTAMINATION_AUTOBUILD)
#       Downloads each species' genome+GTF from NCBI using `datasets` and
#       concatenates them into contamination.fa / contamination.gtf /
#       gtf_seqnames.tsv. Mirrors the structure of the user's
#       download_species.sh helper.
#   r05a0_mask_conserved_regions: hard-mask rRNA/tRNA regions in
#                                 contamination.fa using bedtools
#                                 maskfasta. These regions cross-map to
#                                 host rRNA/mitochondrial RNA at 100%
#                                 identity regardless of alignment
#                                 stringency, so masking them is the only
#                                 way to eliminate that noise source.
#   r05a_contamination_index    : builds the bwa-mem2 index (run once)
#   r05b_extract_unmapped       : per-sample BAM -> unmapped.fq.gz
#   r05c_entropy_filter         : per-sample bbduk entropy filter to drop
#                                 low-complexity reads (poly-A, simple
#                                 repeats, etc.) before alignment. This is
#                                 the strongest single intervention
#                                 against the cross-mapping false-positive
#                                 floor.
#   r05d_contamination_screen   : per-sample bwa-mem2 alignment of the
#                                 entropy-filtered reads to the
#                                 contamination index, plus per-species
#                                 read attribution via the seqname map.
#
# r04d_qc_summary depends on r05d's per-sample counts via the conditional
# `unpack(...)` input, and folds the contamination metrics into the
# single cohort QC TSV (no separate summary file).
#
# Contig naming: contamination.fa uses the native NCBI seqnames
# (e.g. ">NC_010163.1 Acholeplasma laidlawii PG-8A"). The seqname->species
# mapping comes from gtf_seqnames.tsv (header: species_id<TAB>seqname);
# r05_prep_contamination_refs generates this file in auto-build mode by
# reading FASTA headers of each downloaded genome.

if CONTAMINATION_ENABLED:

    if CONTAMINATION_AUTOBUILD:

        rule r05_prep_contamination_refs:
            """
            On-demand contamination-reference builder. Reads species.txt
            (tab-separated: species_name<TAB>assembly_accession), downloads
            each genome+GTF from NCBI using the `datasets` CLI, then
            concatenates them into the three files r05a/r05b/r04d expect.

            Some NCBI assemblies, especially fungi/environmental genomes, provide
            genome FASTA but no GTF. For contamination screening we only need a
            syntactically valid annotation file, so generate one mock whole-contig
            gene/transcript/exon per contig. Do not interpret these as real genes.

            Internet access is required at run time (NCBI download). The
            species.txt file should be small (a few dozen species at
            most); the per-species download is ~3-50 MB and serialized,
            so the total wall-clock for a 10-species panel is ~2-5 min.

            Idempotent re-runs: if all three output files exist and are
            non-empty, Snakemake will not re-fire this rule. To force a
            rebuild after editing species.txt, delete the three output
            files.
            """
            input:
                species_list = CONTAM_SPECIES_LIST,
            output:
                fa       = CONTAM_FA,
                gtf      = CONTAM_GTF,
                seqnames = CONTAM_SEQNAMES,
            threads: 2
            singularity: DATASETS_CONTAINER
            shell:
                r"""
                set -euo pipefail
                CONTAM_DIR=$(dirname {output.fa})
                WORK="$CONTAM_DIR/_dl"
                mkdir -p "$WORK"

                # Initialize outputs to empty so any failure mid-way leaves
                # zero-byte files (which CONTAMINATION_ENABLED treats as
                # absent on next run).
                : > {output.fa}
                : > {output.gtf}
                printf "species_id\tseqname\n" > {output.seqnames}

                # Iterate species.txt. Format mirrors download_species.sh:
                #   species_name<TAB>assembly_accession   (no header)
                while IFS=$'\t' read -r species accession || \
                    [[ -n "${{species:-}}${{accession:-}}" ]]; do
                    [[ -z "${{species// }}" ]] && continue
                    [[ "${{species:0:1}}" == "#" ]] && continue

                    # Strip CRLF and leading/trailing whitespace
                    species=$(printf '%s' "$species" | sed 's/\r$//; s/^[[:space:]]*//; s/[[:space:]]*$//')
                    accession=$(printf '%s' "$accession" | sed 's/\r$//; s/^[[:space:]]*//; s/[[:space:]]*$//')

                    [[ -z "$species" ]] && continue
                    [[ "${{species:0:1}}" == "#" ]] && continue
                    [[ -z "$accession" ]] && continue

                    if [[ ! "$accession" =~ ^GC[AF]_[0-9]+\.[0-9]+$ ]]; then
                        echo "[r05_prep] ERROR: malformed accession for $species: [$accession]" >&2
                        printf '[r05_prep] shell-escaped accession: %q\n' "$accession" >&2
                        exit 1
                    fi

                    safe_species=$(printf '%s' "$species" \
                        | sed 's/[[:space:]]\+/_/g; s#[/\\]#_#g')
                    species_id="${{safe_species}}_${{accession}}"

                    zipfile="$WORK/${{species_id}}.zip"
                    destdir="$WORK/${{species_id}}"

                    printf '[r05_prep] accession shell-escaped: %q\n' "$accession" >&2
                    datasets version >&2 || true

                    echo "[r05_prep] Downloading: $species ($accession)" >&2
                    rm -f "$zipfile"

                    if ! datasets download genome accession "$accession" \
                            --include genome,gtf \
                            --filename "$zipfile" \
                            --no-progressbar; then
                        echo "[r05_prep] WARNING: genome,gtf download failed for $species ($accession); trying genome,gff3" >&2
                        rm -f "$zipfile"
                        datasets download genome accession "$accession" \
                            --include genome,gff3 \
                            --filename "$zipfile" \
                            --no-progressbar
                    fi

                    mkdir -p "$destdir"
                    unzip -oq "$zipfile" -d "$destdir"
                    rm -f "$zipfile"

                    # Locate the genome FASTA and GTF inside the unzipped
                    # tree (datasets puts them under
                    # ncbi_dataset/data/<accession>/).
                    asm_dir="$destdir/ncbi_dataset/data/${{accession}}"
                    if [[ ! -d "$asm_dir" ]]; then
                        echo "[r05_prep] ERROR: assembly dir missing for $species" >&2
                        exit 1
                    fi

                    # FASTA: usually *_genomic.fna (no fixed prefix).
                    # Use the first .fna found that is not the genomic_gaps
                    # variant.
                    asm_fa=$(find "$asm_dir" -maxdepth 2 -name "*.fna" \
                            ! -name "*_genomic_gaps.fna" | head -1)
                    asm_gtf=$(find "$asm_dir" -maxdepth 2 -name "*.gtf" | head -1)

                    if [[ -z "$asm_fa" ]]; then
                        echo "[r05_prep] ERROR: missing genome FASTA for $species" >&2
                        echo "[r05_prep]   asm_fa=$asm_fa" >&2
                        exit 1
                    fi

                    if [[ -z "$asm_gtf" ]]; then
                        echo "[r05_prep] WARNING: missing GTF for $species ($accession); generating mock whole-contig GTF" >&2

                        asm_gtf="$destdir/${{species_id}}.mock.gtf"

                        awk -v sp="$species_id" '
                            BEGIN {{
                                OFS = "\t"
                                seq = ""
                                len = 0
                            }}

                            /^>/ {{
                                if (seq != "" && len > 0) {{
                                    gene_id = sp "__" seq
                                    tx_id   = gene_id ".t1"

                                    print seq, "mock_gtf", "gene",       1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; gene_name \"" gene_id "\";"

                                    print seq, "mock_gtf", "transcript", 1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; transcript_id \"" tx_id "\"; gene_name \"" gene_id "\";"

                                    print seq, "mock_gtf", "exon",       1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; transcript_id \"" tx_id "\"; exon_number \"1\";"
                                }}

                                seq = substr($1, 2)
                                len = 0
                                next
                            }}

                            {{
                                gsub(/[[:space:]]/, "", $0)
                                len += length($0)
                            }}

                            END {{
                                if (seq != "" && len > 0) {{
                                    gene_id = sp "__" seq
                                    tx_id   = gene_id ".t1"

                                    print seq, "mock_gtf", "gene",       1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; gene_name \"" gene_id "\";"

                                    print seq, "mock_gtf", "transcript", 1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; transcript_id \"" tx_id "\"; gene_name \"" gene_id "\";"

                                    print seq, "mock_gtf", "exon",       1, len, ".", "+", ".", \
                                        "gene_id \"" gene_id "\"; transcript_id \"" tx_id "\"; exon_number \"1\";"
                                }}
                            }}
                        ' "$asm_fa" > "$asm_gtf"
                    fi

                    # Concatenate FASTA + GTF onto the cohort files (we do
                    # NOT rename seqnames; native NCBI accessions are used
                    # throughout, with gtf_seqnames.tsv carrying the
                    # accession-to-species mapping).
                    cat "$asm_fa"  >> {output.fa}
                    cat "$asm_gtf" >> {output.gtf}

                    # Extract seqnames from the FASTA headers and emit
                    # one row per seqname in gtf_seqnames.tsv.
                    grep '^>' "$asm_fa" | awk -v sp="$species_id" '
                        BEGIN {{ OFS="\t" }}
                        {{
                            seqname = substr($1, 2)
                            print sp, seqname
                        }}
                    ' >> {output.seqnames}

                    rm -rf "$destdir"
                done < {input.species_list}

                rm -rf "$WORK"

                # Sanity: ensure each output is non-empty
                for f in {output.fa} {output.gtf} {output.seqnames}; do
                    if [[ ! -s "$f" ]]; then
                        echo "[r05_prep] ERROR: output is empty: $f" >&2
                        exit 1
                    fi
                done

                n_species=$(awk 'NR>1 {{print $1}}' {output.seqnames} \
                            | sort -u | wc -l)
                n_seqs=$(awk 'NR>1' {output.seqnames} | wc -l)
                echo "[r05_prep] Built $n_species species, $n_seqs contigs" >&2
                """

    rule r05a0_mask_conserved_regions:
        """
        Hard-mask conserved rRNA/tRNA regions in the contamination FASTA
        so that host-derived rRNA reads don't generate cross-mapping
        contamination hits.

        Mechanism: bacterial 16S/23S rRNA has stretches 100% identical
        to human cytosolic rRNA; mitochondrial rRNA reads (very abundant
        in any RNA-seq library that didn't deplete mt-RNA) cross-map to
        apicomplexan apicoplast rRNA and to bacterial rRNA broadly;
        fungal small-subunit rRNA shares conserved cores with all of the
        above. No alignment-parameter tightening can fix this -- the
        reads ARE identical sequences. The only effective fix is to
        prevent bwa-mem2 from seeing those regions in the reference at
        all, by hard-masking them with N's. bwa-mem2 won't seed into N
        runs, so reads hitting masked regions fail to align entirely.

        Feature types masked (extracted from contamination.gtf):
          - rRNA      : the primary problem; all bacterial / eukaryotic
                        species have annotated rRNA features
          - tRNA      : also evolutionarily conserved; secondary noise
                        source
          - ncRNA     : misc conserved ncRNAs (riboswitches, etc.)
          - misc_RNA  : catch-all for additional conserved RNA features

        Predicted effect: drops Plasmodium/Toxoplasma/Penicillium/
        Cladosporium/Cutibacterium "contamination" counts to near-zero
        in samples that don't actually carry those organisms, since
        their current signal is dominated by rRNA cross-mapping. Real
        contamination (which has reads across the whole genome, not
        just at rRNA loci) is preserved.
        """
        input:
            fa = CONTAM_FA,
            gtf = CONTAM_GTF,
        output:
            masked_fa = "/tmp/data/05_contamination/index/contamination.masked.fa",
            bed = "/tmp/data/05_contamination/index/conserved_regions.bed",
        threads: 2
        singularity: BEDTOOLS_CONTAINER
        shell:
            r"""
            set -euo pipefail
            OUTDIR=$(dirname {output.masked_fa})
            mkdir -p "$OUTDIR"

            # Extract conserved-RNA features from the GTF/GFF. Handles
            # both formats because GTF column 3 (feature type) is the
            # same in GFF3; we just look for feature-type values
            # matching rRNA / tRNA / ncRNA / misc_RNA. Note: column 4 is
            # 1-based inclusive (GTF/GFF), but BED is 0-based half-open,
            # so we subtract 1 from start. We also drop features with
            # malformed coordinates (start > end), which can happen for
            # malformed annotations.
            awk -F'\t' 'BEGIN {{ OFS="\t" }}
                /^#/ {{ next }}
                NF >= 8 && ($3 == "rRNA" || $3 == "tRNA" ||
                            $3 == "ncRNA" || $3 == "misc_RNA") {{
                    if ($4 + 0 > 0 && $5 + 0 >= $4 + 0) {{
                        print $1, $4 - 1, $5
                    }}
                }}
            ' {input.gtf} | sort -k1,1 -k2,2n -u > {output.bed}

            N_FEATURES=$(wc -l < {output.bed})
            echo "[r05a0_mask] extracted $N_FEATURES conserved-RNA features to mask" >&2

            # Even if no features were found (e.g., user supplied a
            # FASTA-only reference or all-virus species.txt), continue
            # gracefully: just copy the FASTA through unchanged.
            if [[ "$N_FEATURES" -eq 0 ]]; then
                echo "[r05a0_mask] WARNING: no rRNA/tRNA features in GTF; passing FASTA through unmasked" >&2
                cp {input.fa} {output.masked_fa}
                exit 0
            fi

            # Hard-mask with N's. bwa-mem2 cannot seed into N runs, so
            # reads that would have mapped into these regions simply
            # fail to align in r05d.
            bedtools maskfasta \
                -fi {input.fa} \
                -bed {output.bed} \
                -fo {output.masked_fa} \
                -mc N

            # Sanity check: report how many bases were masked
            N_BASES_ORIG=$(grep -v '^>' {input.fa}        | tr -d '\n' | tr -cd 'Nn' | wc -c)
            N_BASES_MASK=$(grep -v '^>' {output.masked_fa} | tr -d '\n' | tr -cd 'Nn' | wc -c)
            N_BASES_ADDED=$((N_BASES_MASK - N_BASES_ORIG))
            echo "[r05a0_mask] masked $N_BASES_ADDED bp of conserved-RNA across $N_FEATURES features" >&2
            """

    rule r05a_contamination_index:
        """
        Build a bwa-mem2 index for the contamination reference. Run once,
        shared across all samples. Consumes the masked FASTA from
        r05a0_mask_conserved_regions so that conserved rRNA/tRNA
        regions are N-masked before indexing. The seqname->species
        mapping lives in gtf_seqnames.tsv (either supplied by the user
        or auto-generated by r05_prep_contamination_refs); we do NOT
        re-derive it from FASTA headers because the contamination FASTA
        uses bare NCBI accessions (e.g. NC_010163.1) which don't carry
        species info on their own.

        The index is built under a name that points at the masked
        FASTA, so r05b/r05d (which derive the FASTA path from the
        index_done file's directory) automatically pick up the masked
        version too.
        """
        input:
            fa = "/tmp/data/05_contamination/index/contamination.masked.fa",
        output:
            done = "/tmp/data/05_contamination/index/index.done",
        threads: 4
        singularity: BWA_CONTAINER
        shell:
            r"""
            set -euo pipefail
            INDEX_DIR=$(dirname {output.done})
            mkdir -p "$INDEX_DIR"
            # r05d's $INDEX_FA expects the file named contamination.fa
            # in the index dir. Since the masked FASTA already lives
            # there, just symlink the conventional name.
            ln -sf $(basename {input.fa}) "$INDEX_DIR/contamination.fa"
            cd "$INDEX_DIR"
            bwa-mem2 index contamination.fa
            touch {output.done}
            """

    rule r05b_extract_unmapped:
        """
        Extract a sample's unmapped reads from the STAR BAM as a single
        gzipped interleaved FASTQ. Separated from the contamination
        screen so the entropy-filter step (r05c) can operate on the
        FASTQ before alignment.

        Filter flags:
          -f 4    : keep only unmapped reads
          -F 256  : drop secondary alignments
          -N      : retain /1, /2 read-pair suffix in QNAME so the
                    downstream interleaved-stream consumers can pair
                    reads correctly
        """
        input:
            bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        output:
            fq = temp("/tmp/data/05_contamination/{sample}/{sample}.unmapped.fq.gz"),
        threads: 2
        singularity: BWA_CONTAINER
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.fq})
            samtools fastq -f 4 -F 256 -N {input.bam} 2>/dev/null \
                | gzip > {output.fq}
            """

    rule r05c_entropy_filter:
        """
        Drop low-complexity / low-entropy reads using bbduk.sh before the
        contamination alignment. This is the strongest single intervention
        against the cross-mapping false-positive floor: reads like
        poly-A tails, ATx50 simple repeats, and other low-complexity
        sequences align cleanly to many references (especially AT-rich
        ones like Plasmodium falciparum or AT-rich bacterial regions)
        regardless of bwa-mem2 score / seed-length stringency settings.

        Parameters:
          entropy=0.7        Shannon entropy threshold (1.0 = max entropy).
                             Typical RNA-seq reads sit at 0.85-0.95; poly-A
                             at ~0; simple repeats at ~0.5. 0.7 drops the
                             low-complexity tail without eating into real
                             reads.
          entropywindow=20   Window size for entropy calculation (default 50;
                             20 catches shorter low-complexity stretches).
          entropyk=5         k-mer size for entropy calculation (default 5).
          int=t              Treat stream as interleaved paired-end (matches
                             samtools fastq -N output for paired data; safe
                             for single-end since singletons are emitted to
                             the same stream and processed independently).
        """
        input:
            fq = "/tmp/data/05_contamination/{sample}/{sample}.unmapped.fq.gz",
        output:
            fq = temp("/tmp/data/05_contamination/{sample}/{sample}.unmapped.filtered.fq.gz"),
        threads: 2
        singularity: BBMAP_CONTAINER
        shell:
            r"""
            set -euo pipefail
            bbduk.sh \
                in={input.fq} \
                out={output.fq} \
                entropy=0.7 \
                entropywindow=20 \
                entropyk=5 \
                int=t \
                threads={threads} \
                overwrite=t 2>&1 | grep -E "Input|Result|entropy|Time" || true
            """

    rule r05d_contamination_screen:
        """
        Align the entropy-filtered unmapped reads to the contamination
        reference with conservative bwa-mem2 parameters, then attribute
        per-species read counts via the seqname->species map.

        bwa-mem2 parameters:
          -T 95   : minimum alignment score (default 30). A 95bp read
                    needs near-perfect identity to be reported. Pairs
                    with the entropy filter to push the false-positive
                    rate as low as practical.
          -k 35   : minimum seed length (default 19). Longer seeds
                    reject short noisy matches before extension.
          -r 2.0  : re-seed threshold multiplier (default 1.5). Less
                    aggressive re-seeding, fewer marginal alignments.
        Plus the post-filter samtools view -q 30 (MAPQ confidence >= 30)
        to drop ambiguous multi-mappers across organisms.

        Per-species attribution uses the user-supplied gtf_seqnames.tsv
        (header: species_id<TAB>seqname). Reads whose accession is
        absent from the map are bucketed into 'unknown' rather than
        dropped, so additions to the FASTA without matching seqname-map
        entries are still counted and visibly attributed.
        """
        input:
            fq = "/tmp/data/05_contamination/{sample}/{sample}.unmapped.filtered.fq.gz",
            bam = "/tmp/data/03_bam_star/{sample}/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
            index_done = "/tmp/data/05_contamination/index/index.done",
            species_map = CONTAM_SEQNAMES,
        output:
            counts = "/tmp/data/05_contamination/{sample}/{sample}.contamination_counts.tsv",
            bam = temp("/tmp/data/05_contamination/{sample}/{sample}.contamination.bam"),
        threads: 4
        singularity: BWA_CONTAINER
        shell:
            r"""
            set -euo pipefail
            OUTDIR=$(dirname {output.counts})
            mkdir -p "$OUTDIR"

            INDEX_FA=$(dirname {input.index_done})/contamination.fa

            # bwa-mem2 with -p reads interleaved paired-end from stdin
            # (matches bbduk's int=t output and samtools fastq -N).
            # Works correctly for SE input too since singletons stream
            # through the same path.
            zcat {input.fq} \
                | bwa-mem2 mem -t {threads} -p -T 95 -k 35 -r 2.0 -B 6 -O 8 -L 10 \
                      "$INDEX_FA" - 2>/dev/null \
                | samtools view -b -F 4 - > {output.bam}

            # Totals: total unmapped reads in the host BAM (denominator
            # context only; the QC summary uses Input_reads as the
            # percentage denominator), and total confidently-mapped to
            # any contaminant genome.
            TOTAL_UNMAPPED=$(samtools view -c -f 4 -F 256 {input.bam})
            TOTAL_MAPPED=$(samtools view -c -F 260 -q 30 {output.bam})

            # Per-species counts via samtools view + awk join on the
            # gtf_seqnames.tsv map. All species from the map are emitted
            # (with count 0 if absent in this sample) so the cohort
            # summary has consistent columns across samples. NOTE: the
            # seqname-map header is "species_id<TAB>seqname" so we key
            # the lookup on the SECOND column (seqname / NCBI accession)
            # -> FIRST column (species_id).
            samtools view -F 260 -q 30 {output.bam} \
                | awk -v map="{input.species_map}" \
                      -v sample="{wildcards.sample}" \
                      -v total_unmapped="$TOTAL_UNMAPPED" \
                      -v total_mapped="$TOTAL_MAPPED" \
                      -v out_counts="{output.counts}" '
                    BEGIN {{
                        OFS="\t"
                        # load seqname -> species_id; build species set
                        while ((getline line < map) > 0) {{
                            n = split(line, a, "\t")
                            if (n == 2 && a[1] != "species_id") {{
                                # a[1] = species_id, a[2] = seqname
                                seqname2species[a[2]] = a[1]
                                all_species[a[1]] = 1
                            }}
                        }}
                        close(map)
                    }}
                    {{
                        species = (($3 in seqname2species) ? seqname2species[$3] : "unknown")
                        counts[species]++
                        all_species[species] = 1
                    }}
                    END {{
                        print "metric", "value" > out_counts
                        print "sample", sample > out_counts
                        print "total_unmapped_in_star_bam", total_unmapped > out_counts
                        print "total_mapped_to_contamination", total_mapped > out_counts
                        # emit a deterministic order of species rows
                        n = 0
                        for (s in all_species) sp_arr[++n] = s
                        # naive insertion sort
                        for (i = 2; i <= n; i++) {{
                            key = sp_arr[i]; j = i - 1
                            while (j > 0 && sp_arr[j] > key) {{
                                sp_arr[j + 1] = sp_arr[j]; j--
                            }}
                            sp_arr[j + 1] = key
                        }}
                        for (i = 1; i <= n; i++) {{
                            print "species:" sp_arr[i], (sp_arr[i] in counts ? counts[sp_arr[i]] : 0) > out_counts
                        }}
                    }}
                '
            """


rule r04d_qc_summary:
    """
    Aggregate per-sample CollectRnaSeqMetrics outputs (Picard) and STAR
    mapping statistics into one table, flagging outliers. Picard's metrics
    file has two main sections:
      - METRICS section: one line of values per sample
      - HISTOGRAM section: normalized coverage by percentile of transcript
        position (0=5'-end, 100=3'-end); we don't aggregate this here but
        leave the per-sample file available for plotting.

    The STAR Log.final.out values come from r03h_mapping_stat's aggregated
    output (mapping_stat.txt) so we don't re-parse the per-sample logs here.
    The uniquely-mapped read count drives the library-complexity flag.

    If the optional contamination screen (r05a/r05b) is enabled, the
    per-sample contamination counts are folded into the same cohort table.
    For each species in the contamination reference we add two data
    columns (<species>_reads and <species>_rpm) plus a per-species flag.
    A combined flag_contamination fires if total contaminant reads exceed
    100 reads per million unmapped (the Olarerin-George & Hogenesch 2015
    threshold).

    Thresholds chosen for DROP / OUTRIDER analysis:
      - 3' bias high: cohort-relative outlier (median + 3*MAD on MEDIAN_3PRIME_BIAS)
      - rRNA high:    > 10%  (PCT_RIBOSOMAL_BASES > 0.10)
      - mRNA low:     < 60%  (PCT_MRNA_BASES < 0.60)
      - low complexity: < 10 M uniquely-mapped reads (OUTRIDER/FRASER lose
                        statistical power below this; ENCODE's full-quality
                        bar is 30 M, but 10 M is the working floor)
      - contamination: combined contaminant reads > 0.01% of total input reads
      - per-species:   single-species contaminant reads > 0.001% of total input reads

    Contamination percentages use TOTAL input reads (STAR's "Number of
    input reads") as denominator, so percentages are directly comparable
    across samples regardless of host-alignment rate. The 0.01% threshold
    corresponds to 100 contaminant reads per million total reads -- the
    same Olarerin-George & Hogenesch (2015, NAR) bar as the original
    rpm-of-unmapped formulation, just rescaled to a denominator that
    doesn't change with sample alignment quality.
    """
    input:
        unpack(lambda wc: {
            "metrics_files": expand(
                "/tmp/data/04_qc/{sample}/{sample}.rnaseq_metrics.txt",
                sample=samples),
            "mapping_stat":
                "/tmp/data/03_bam_star/00_mapping_stat/mapping_stat.txt",
            **(
                {
                    "contam_counts": expand(
                        "/tmp/data/05_contamination/{sample}/{sample}.contamination_counts.tsv",
                        sample=samples),
                    "contam_species_map": CONTAM_SEQNAMES,
                }
                if CONTAMINATION_ENABLED else {}
            ),
        }),
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

        # Extra columns we surface from STAR's mapping_stat.txt. Input_reads
        # is needed as the denominator for contamination percentages so that
        # the value is comparable across samples (independent of host-alignment
        # rate, library size, etc).
        mapping_cols = ["Input_reads", "Uniquely_mapped", "Uniquely_mapped_pct"]

        # Outlier flags from RNA-quality / complexity metrics
        flag_cols = [
            "flag_3p_bias_high",
            "flag_rrna_high",
            "flag_mrna_low",
            "flag_low_complexity",
        ]

        # -- Parse Picard metrics files --
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

        # -- Parse STAR mapping_stat.txt for read-count columns --
        # Format (from r03h_mapping_stat):
        #   Sample  Input_reads  Uniquely_mapped  Uniquely_mapped_pct
        #           Multi_mapped Multi_mapped_pct Too_many_loci Unmapped
        with open(input.mapping_stat) as f:
            ms_lines = f.readlines()
        if ms_lines:
            ms_header = ms_lines[0].rstrip("\n").split("\t")
            for line in ms_lines[1:]:
                parts = line.rstrip("\n").split("\t")
                if not parts or not parts[0]:
                    continue
                sample = parts[0]
                if sample not in rows:
                    # STAR-only sample; will appear as Picard-NaN row below
                    rows[sample] = {c: "" for c in cols_of_interest}
                ms_row = dict(zip(ms_header, parts))
                for c in mapping_cols:
                    rows[sample][c] = ms_row.get(c, "")

        # -- Parse contamination counts (when CONTAMINATION_ENABLED) --
        # Format per-sample (long-form key-value, from r05b):
        #   metric                              value
        #   sample                              S001
        #   total_unmapped_in_star_bam          12345678
        #   total_mapped_to_contamination       42
        #   species:mycoplasma_hyorhinis_GCF... 37
        #   species:candida_albicans_GCF...     5
        #   ...
        #
        # Percentages use the sample's TOTAL input reads (STAR's
        # "Number of input reads", surfaced as Input_reads) as denominator,
        # so they are directly comparable across samples regardless of how
        # much of the library aligned to the host genome.
        contam_cols = []           # Data columns appended to row table
        contam_flag_cols = []      # Flag columns appended after the RNA flags
        species_sorted = []        # Stable ordering for column emission
        if CONTAMINATION_ENABLED:
            # Build species list from the user-supplied seqname map. The
            # map's first column is species_id (species name + GCF
            # accession joined by underscores), second column is seqname
            # (NCBI accession of one contig from that species).
            species_set = set()
            with open(input.contam_species_map) as f:
                next(f)  # skip header (species_id<TAB>seqname)
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == 2:
                        species_set.add(parts[0])
            species_sorted = sorted(species_set)

            contam_cols = (
                ["total_mapped_to_contamination",
                 "contamination_pct_of_input"]
                + [f"{sp}_reads" for sp in species_sorted]
                + [f"{sp}_pct"   for sp in species_sorted]
            )
            contam_flag_cols = (
                ["flag_contamination"]
                + [f"flag_{sp}" for sp in species_sorted]
            )

            for f_path in input.contam_counts:
                data = {}
                with open(f_path) as fh:
                    next(fh)  # skip header
                    for line in fh:
                        parts = line.rstrip("\n").split("\t", 1)
                        if len(parts) == 2:
                            data[parts[0]] = parts[1]

                sample = data.get("sample", "")
                if not sample:
                    continue
                if sample not in rows:
                    rows[sample] = {c: "" for c in cols_of_interest}

                # Denominator: TOTAL input reads (from mapping_stat, already
                # in the row dict via the mapping_cols pass above).
                try:
                    total_input = int(rows[sample].get("Input_reads", 0) or 0)
                except (ValueError, TypeError):
                    total_input = 0
                total_contam = int(data.get("total_mapped_to_contamination", 0) or 0)
                pct = (total_contam / total_input * 100) if total_input > 0 else 0.0

                rows[sample]["total_mapped_to_contamination"] = total_contam
                rows[sample]["contamination_pct_of_input"] = f"{pct:.6f}"

                for sp in species_sorted:
                    sp_reads = int(data.get(f"species:{sp}", 0) or 0)
                    sp_pct = (sp_reads / total_input * 100) if total_input > 0 else 0.0
                    rows[sample][f"{sp}_reads"] = sp_reads
                    rows[sample][f"{sp}_pct"]   = f"{sp_pct:.6f}"

        # -- Compute cohort statistics for outlier detection --
        def floats(col):
            return [float(rows[s][col]) for s in rows
                    if rows[s].get(col, "") not in ("", "?", "NA")]
        bias_3p = floats("MEDIAN_3PRIME_BIAS")
        med_3p_bias = statistics.median(bias_3p) if bias_3p else 0
        mad_3p_bias = statistics.median([abs(x - med_3p_bias) for x in bias_3p]) \
                      if bias_3p else 0
        bias_3p_threshold = med_3p_bias + 3 * mad_3p_bias

        # Library-complexity threshold: 10M uniquely-mapped reads. ENCODE's
        # gold-standard bar is 30M, but OUTRIDER/FRASER autoencoders can still
        # work down to ~10M with reduced power; below that, statistical
        # outlier detection becomes unreliable.
        complexity_threshold = 10_000_000

        # Contamination thresholds (only used when CONTAMINATION_ENABLED).
        # Both are expressed as percent of total input reads to match the
        # contamination column units. 100 reads per million total reads =
        # 0.01% combined; per-species 10 reads per million = 0.001%.
        contam_total_pct_threshold   = 0.01
        contam_species_pct_threshold = 0.001

        def is_outlier(sample):
            r = rows[sample]
            flags = {}
            # 3' bias: cohort-relative MAD outlier
            try:
                flags["flag_3p_bias_high"] = float(r["MEDIAN_3PRIME_BIAS"]) > bias_3p_threshold
            except (ValueError, KeyError):
                flags["flag_3p_bias_high"] = ""
            # rRNA: absolute threshold
            try:
                flags["flag_rrna_high"] = float(r["PCT_RIBOSOMAL_BASES"]) > 0.10
            except (ValueError, KeyError):
                flags["flag_rrna_high"] = ""
            # mRNA: absolute threshold
            try:
                flags["flag_mrna_low"] = float(r["PCT_MRNA_BASES"]) < 0.60
            except (ValueError, KeyError):
                flags["flag_mrna_low"] = ""
            # Library complexity: absolute threshold on uniquely-mapped reads
            try:
                flags["flag_low_complexity"] = int(r["Uniquely_mapped"]) < complexity_threshold
            except (ValueError, KeyError):
                flags["flag_low_complexity"] = ""

            # Contamination flags (only when CONTAMINATION_ENABLED)
            if CONTAMINATION_ENABLED:
                try:
                    flags["flag_contamination"] = (
                        float(r["contamination_pct_of_input"]) > contam_total_pct_threshold
                    )
                except (ValueError, KeyError):
                    flags["flag_contamination"] = ""
                for sp in species_sorted:
                    try:
                        flags[f"flag_{sp}"] = (
                            float(r[f"{sp}_pct"]) > contam_species_pct_threshold
                        )
                    except (ValueError, KeyError):
                        flags[f"flag_{sp}"] = ""
            return flags

        # -- Write summary --
        all_data_cols = cols_of_interest + mapping_cols + contam_cols
        all_flag_cols = flag_cols + contam_flag_cols
        with open(output.summary, "w") as fout:
            fout.write("sample\t" + "\t".join(all_data_cols)
                       + "\t" + "\t".join(all_flag_cols) + "\n")
            for sample in sorted(rows):
                vals = [str(rows[sample].get(c, "")) for c in all_data_cols]
                flags = is_outlier(sample)
                flag_vals = [str(flags.get(c, "")) for c in all_flag_cols]
                fout.write(sample + "\t" + "\t".join(vals)
                           + "\t" + "\t".join(flag_vals) + "\n")



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
#    if reads are 2x76bp. For other read lengths, build your own.
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
