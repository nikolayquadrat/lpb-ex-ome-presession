#!/usr/bin/env bash
# =============================================================================
# lpb-rnaseq-set-up-arcashla.sh
# =============================================================================
# Build a complete arcasHLA reference for the RNA-seq HLA-typing rules in
# lpb-rnaseq-pipe.smk.
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# The bioconda biocontainer quay.io/biocontainers/arcas-hla:0.6.0--hdfd78af_2
# ships only a PARTIAL reference. To run arcasHLA you need, under dat/:
#   - IMGTHLA/        a clone of the ANHIG/IMGTHLA database (with hla.dat
#                     UNZIPPED -- it ships compressed and arcasHLA needs the
#                     uncompressed file plus the wmda/ nomenclature tables)
#   - info/           small JSON/TSV files the container bundles (parameters.json,
#                     hla_freq.tsv, ...)
#   - ref/            the derived files: hla.fasta, hla.idx, hla.p.json,
#                     hla.convert.json, hla_partial.fasta, hla_partial.idx,
#                     hla_partial.p.json  -- built by `arcasHLA reference
#                     --rebuild`, which parses IMGTHLA and builds the Kallisto
#                     indices.
#
# This script does all four steps, in order, idempotently, with verification
# after each. It is SELF-CONTAINED: it does not assume any file already exists
# and is safe to run on a fresh machine.
#
# QUIRKS HANDLED (learned the hard way)
# -----------------------------------------------------------
#  1. hla.dat ships as hla.dat.zip inside the IMGTHLA repo -- must be unzipped
#     (plus the other *.zip archives). Uncompressed hla.dat is ~323 MB.
#  2. The container's bundled info/ and ref/ must be copied to the HOST dat/
#     dir. When you bind the host dat/ over the container's dat/ at run time,
#     the bind HIDES the container's bundled files, so the host copy must
#     contain them. We use `cp -rT` (NOT `cp -r`) to avoid the
#     "destination dir already exists -> nested copy" footgun.
#  3. `arcasHLA reference --rebuild` must be run with the host dat/ bound OVER
#     the container's internal dat/ path
#     (/usr/local/share/arcas-hla-0.6.0-2/dat), so the rebuild writes into the
#     host directory.
#  4. apptainer/singularity is required. Container operations may need to run
#     under sudo depending on how apptainer is configured on the host (suid
#     install vs unprivileged userns). This script auto-detects and can be
#     told to use sudo via ARCASHLA_USE_SUDO=1.
#
# USAGE
# -----
#   bash install-arcashla-ref.sh
#
# Environment overrides (all optional):
#   ARCASHLA_REF_DIR     install target (default below)
#   ARCASHLA_CONTAINER   biocontainer URI
#   ARCASHLA_USE_SUDO    set to 1 to prefix container calls with sudo
#   SING_BIN             force "apptainer" or "singularity"
#
# After it completes, run the RNA-seq pipeline with the host dat/ bound over
# the container's dat/. The exact bind line for this machine is printed at the
# end.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
ARCASHLA_REF_DIR="${ARCASHLA_REF_DIR:-/mnt/data/rnaseq/rnaseq-drop/00_additional_files/arcashla_ref}"
ARCASHLA_CONTAINER="${ARCASHLA_CONTAINER:-docker://quay.io/biocontainers/arcas-hla:0.6.0--hdfd78af_2}"
CONTAINER_DAT_PATH="/usr/local/share/arcas-hla-0.6.0-2/dat"
IMGTHLA_REPO="https://github.com/ANHIG/IMGTHLA.git"

HOST_DAT="$ARCASHLA_REF_DIR/dat"
IMGTHLA_DIR="$HOST_DAT/IMGTHLA"

# -----------------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------------
log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
err()  { printf '[%s] ERROR: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die()  { err "$*"; exit 1; }

# -----------------------------------------------------------------------------
# 0. Tool detection
# -----------------------------------------------------------------------------
log "==== arcasHLA reference installer ===="
log "install dir : $ARCASHLA_REF_DIR"
log "container   : $ARCASHLA_CONTAINER"

# git + unzip are needed on the host
command -v git   >/dev/null 2>&1 || die "git not found on PATH (needed to clone IMGTHLA)"
command -v unzip >/dev/null 2>&1 || die "unzip not found on PATH (needed to decompress hla.dat.zip)"

# apptainer or singularity
SING_BIN="${SING_BIN:-}"
if [[ -z "$SING_BIN" ]]; then
    if command -v apptainer >/dev/null 2>&1; then
        SING_BIN="apptainer"
    elif command -v singularity >/dev/null 2>&1; then
        SING_BIN="singularity"
    else
        die "neither apptainer nor singularity found on PATH. Install one, e.g.:
       sudo apt-get install -y apptainer
     or see https://apptainer.org/docs/admin/main/installation.html"
    fi
fi
log "container runtime: $SING_BIN"

# sudo prefix for container calls (some apptainer installs need it)
SUDO=""
if [[ "${ARCASHLA_USE_SUDO:-0}" == "1" ]]; then
    SUDO="sudo"
    log "ARCASHLA_USE_SUDO=1 -> container calls will use sudo"
fi

# Convenience wrapper for "run a command inside the arcasHLA container"
sing_exec() {
    # usage: sing_exec <bind_spec_or_empty> <command...>
    local bind="$1"; shift
    if [[ -n "$bind" ]]; then
        $SUDO "$SING_BIN" exec --bind "$bind" "$ARCASHLA_CONTAINER" "$@"
    else
        $SUDO "$SING_BIN" exec "$ARCASHLA_CONTAINER" "$@"
    fi
}

# -----------------------------------------------------------------------------
# Idempotency: are we already done?
# -----------------------------------------------------------------------------
KEY_FILES=(
    "$HOST_DAT/IMGTHLA/hla.dat"
    "$HOST_DAT/IMGTHLA/wmda/hla_nom_p.txt"
    "$HOST_DAT/info/parameters.json"
    "$HOST_DAT/ref/hla.p.json"
    "$HOST_DAT/ref/hla.convert.json"
    "$HOST_DAT/ref/hla.idx"
    "$HOST_DAT/ref/hla_partial.idx"
)

all_present() {
    local f
    for f in "${KEY_FILES[@]}"; do
        [[ -s "$f" ]] || return 1
    done
    return 0
}

if all_present; then
    log "reference already fully built at $HOST_DAT -- nothing to do."
    log "key files present:"
    for f in "${KEY_FILES[@]}"; do
        log "  $(du -h "$f" | awk '{printf "%-8s", $1}')  $f"
    done
    # Everything below is guarded by SKIP_BUILD; verification + bind-line
    # printout still run so the user gets the run instructions.
    SKIP_BUILD=1
else
    SKIP_BUILD=0
fi

mkdir -p "$HOST_DAT"

# =============================================================================
# Step 1: clone ANHIG/IMGTHLA
# =============================================================================
if [[ "$SKIP_BUILD" == "0" ]]; then
    if [[ -s "$IMGTHLA_DIR/hla.dat" || -s "$IMGTHLA_DIR/hla.dat.zip" ]]; then
        log "[1/4] IMGTHLA already cloned at $IMGTHLA_DIR -- skipping clone"
    else
        log "[1/4] cloning ANHIG/IMGTHLA (~1.2 GB, depth 1) ..."
        # Remove a partial/empty dir if a previous run died mid-clone.
        if [[ -d "$IMGTHLA_DIR" ]] && [[ -z "$(ls -A "$IMGTHLA_DIR" 2>/dev/null)" ]]; then
            rmdir "$IMGTHLA_DIR"
        fi
        git clone --depth 1 "$IMGTHLA_REPO" "$IMGTHLA_DIR" \
            || die "git clone of IMGTHLA failed"
        log "[1/4] IMGTHLA cloned"
    fi
fi

# =============================================================================
# Step 2: unzip hla.dat.zip (and other archives) inside IMGTHLA
# =============================================================================
if [[ "$SKIP_BUILD" == "0" ]]; then
    if [[ -s "$IMGTHLA_DIR/hla.dat" ]]; then
        log "[2/4] hla.dat already unzipped -- skipping"
    else
        log "[2/4] unzipping archives in $IMGTHLA_DIR ..."
        shopt -s nullglob
        zips=("$IMGTHLA_DIR"/*.zip)
        shopt -u nullglob
        if (( ${#zips[@]} == 0 )); then
            # Some IMGTHLA layouts keep hla.dat uncompressed already; only fail
            # if hla.dat is genuinely absent.
            [[ -s "$IMGTHLA_DIR/hla.dat" ]] \
                || die "no *.zip archives and no hla.dat in $IMGTHLA_DIR"
        fi
        for z in "${zips[@]}"; do
            log "  unzip $(basename "$z")"
            unzip -o -q "$z" -d "$IMGTHLA_DIR" || die "unzip failed: $z"
        done
        [[ -s "$IMGTHLA_DIR/hla.dat" ]] \
            || die "hla.dat still missing after unzip (expected ~323 MB)"
        log "[2/4] hla.dat unzipped ($(du -h "$IMGTHLA_DIR/hla.dat" | cut -f1))"
    fi

    # Sanity: the nomenclature table arcasHLA reads must be present.
    [[ -s "$IMGTHLA_DIR/wmda/hla_nom_p.txt" ]] \
        || die "missing $IMGTHLA_DIR/wmda/hla_nom_p.txt (IMGTHLA clone incomplete)"
fi

# =============================================================================
# Step 3: seed dat/info and dat/ref from the biocontainer
# =============================================================================
# The container bundles small files under its internal dat/info and dat/ref.
# We must copy them into the HOST dat/ so they survive the run-time bind that
# overlays the host dat/ on top of the container's dat/.
#
# Critical: use `cp -rT SRC DST` so DST is treated as the target directory
# itself (not nested into DST/SRC) regardless of whether DST already exists.
if [[ "$SKIP_BUILD" == "0" ]]; then
    seed_done=1
    for f in info/parameters.json ref/cDNA.json; do
        [[ -s "$HOST_DAT/$f" ]] || seed_done=0
    done

    if (( seed_done == 1 )); then
        log "[3/4] dat/info and dat/ref already seeded -- skipping"
    else
        log "[3/4] seeding dat/info + dat/ref from container ..."

        # Remove pre-existing (possibly partial) dirs so a re-seed is clean.
        if [[ -d "$HOST_DAT/info" ]]; then
            log "  removing existing $HOST_DAT/info before re-seed"
            rm -rf "$HOST_DAT/info"
        fi
        if [[ -d "$HOST_DAT/ref" ]]; then
            log "  removing existing $HOST_DAT/ref before re-seed"
            rm -rf "$HOST_DAT/ref"
        fi

        # Bind the host dat/ at /host_dat and copy the container's bundled
        # info/ and ref/ into it with cp -rT.
        sing_exec "$HOST_DAT:/host_dat" bash -c '
            set -e
            src=/usr/local/share/arcas-hla-0.6.0-2/dat
            [ -d "$src" ] || { echo "container dat dir missing: $src" >&2; exit 1; }
            mkdir -p /host_dat
            cp -rT "$src/info" /host_dat/info
            cp -rT "$src/ref"  /host_dat/ref
            echo "--- seeded /host_dat/info ---"; ls -la /host_dat/info | head
            echo "--- seeded /host_dat/ref ----"; ls -la /host_dat/ref  | head
        ' || die "seeding dat/info + dat/ref from container failed"

        # Verify the seed produced the files we expect on the host.
        for f in info/parameters.json info/hla_freq.tsv ref/cDNA.json; do
            [[ -s "$HOST_DAT/$f" ]] \
                || die "seed verification failed: missing $HOST_DAT/$f"
        done
        log "[3/4] dat/info ($(ls "$HOST_DAT/info" | wc -l) files) + dat/ref seeded"
    fi
fi

# =============================================================================
# Step 4: build the derived reference (hla.idx, hla.p.json, ...)
# =============================================================================
# `arcasHLA reference --rebuild` parses IMGTHLA and builds the Kallisto
# indices. It must run with the host dat/ bound OVER the container's internal
# dat/ so it writes the derived files into the host directory.
if [[ "$SKIP_BUILD" == "0" ]]; then
    rebuild_done=1
    for f in ref/hla.p.json ref/hla.convert.json ref/hla.idx ref/hla_partial.idx; do
        [[ -s "$HOST_DAT/$f" ]] || rebuild_done=0
    done

    if (( rebuild_done == 1 )); then
        log "[4/4] derived reference already built -- skipping rebuild"
    else
        log "[4/4] running 'arcasHLA reference --rebuild' (5-15 min) ..."
        sing_exec "$HOST_DAT:$CONTAINER_DAT_PATH" \
            arcasHLA reference --rebuild -v \
            || die "arcasHLA reference --rebuild failed"

        for f in ref/hla.p.json ref/hla.convert.json ref/hla.idx ref/hla_partial.idx; do
            [[ -s "$HOST_DAT/$f" ]] \
                || die "rebuild verification failed: missing $HOST_DAT/$f"
        done
        log "[4/4] rebuild complete:"
        log "  hla.idx         $(du -h "$HOST_DAT/ref/hla.idx"         | cut -f1)"
        log "  hla_partial.idx $(du -h "$HOST_DAT/ref/hla_partial.idx" | cut -f1)"
    fi
fi

# =============================================================================
# Final verification
# =============================================================================
if ! all_present; then
    err "reference build INCOMPLETE -- the following key files are missing:"
    for f in "${KEY_FILES[@]}"; do
        if [[ -s "$f" ]]; then
            log "  present: $f"
        else
            err "  MISSING: $f"
        fi
    done
    die "see messages above; re-run after resolving the issue."
fi

log "================================================================"
log "SUCCESS: arcasHLA reference fully built at"
log "  $HOST_DAT"
log ""
log "Key files:"
for f in "${KEY_FILES[@]}"; do
    log "  $(du -h "$f" | awk '{printf "%-8s", $1}')  $f"
done
log "================================================================"

# -----------------------------------------------------------------------------
# Print the bind line for the snakemake run
# -----------------------------------------------------------------------------
# The pipeline expects the host arcashla_ref dir at /tmp/data/.../arcashla_ref
# (i.e. your data root bound to /tmp/data), AND the dat/ subdir additionally
# bound over the container's internal dat/ path.
log ""
log "NEXT STEP -- run the RNA-seq pipeline with these binds. With apptainer:"
log ""
log "  snakemake -s lpb-rnaseq-pipe.smk \\"
log "    --use-singularity \\"
log "    --singularity-args \"--bind <DATA_ROOT>:/tmp/data --bind $HOST_DAT:$CONTAINER_DAT_PATH\" \\"
log "    [other options]"
log ""
log "where <DATA_ROOT> is the host directory that should appear as /tmp/data"
log "inside the container (the one containing 00_additional_files, 03_bam_star, ...)."
log "For this install that is likely:"
log "  $(dirname "$(dirname "$ARCASHLA_REF_DIR")")"
log ""
log "The pipeline auto-detects the reference at"
log "  /tmp/data/00_additional_files/arcashla_ref/dat"
log "and enables the r06a/r06b/r06c HLA rules only when it is present."
