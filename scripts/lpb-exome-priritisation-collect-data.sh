#!/usr/bin/env bash
# =============================================================================
# Reference / annotation downloader for the DROP-directed WES pipeline (v4)
# -----------------------------------------------------------------------------
# Designed to run on a clean Ubuntu instance with Docker installed.
#
# Improvements over v3:
#  1. Checksums and a per-file manifest (TSV at $MANIFEST_FILE).
#     Files with an upstream .md5/.sha256 are verified after download;
#     all files (downloaded and post-processed) get a locally-computed
#     SHA-256 entry in the manifest. Re-running detects content drift.
#  2. TLS verification ON by default. A small allowlist of hosts where
#     cert validation is known-broken can opt-out via $TLS_INSECURE_HOSTS.
#  3. Post-processing failures (REVEL conversion, SCHEMA HGNC join, etc.)
#     count as real failures - they increment N_FAIL and the script
#     exits non-zero.
#  4. Pinned, dated/versioned URLs. ClinVar uses a dated archive snapshot
#     instead of the rolling latest. Other resources already had versions
#     in their paths.
#  5. Index-file completeness checks. After downloads, every VCF/TSV.gz
#     is required to have a non-empty .tbi alongside it; FASTA must have
#     .fai and .dict. Missing indexes are listed as failures.
#  6. The gnomAD strip+concat helper now exports KEEP into parallel jobs,
#     uses 'parallel --halt' to surface the first failure, and validates
#     the final concatenated VCF + index before declaring success.
#  7. Version coherence is enforced explicitly: VEP_CACHE_RELEASE,
#     VEP_PLUGINS_RELEASE, DBNSFP_VERSION are the single source of truth,
#     printed in a banner at the start, and a sanity-check function
#     warns if you've manually broken alignment (e.g., changed the
#     plugin branch but left the dbNSFP version on v4.x).
#  8. Final completeness validation walks the full expected-files list
#     and exits non-zero if anything is missing or malformed.
#  9. Plugin .pm files from both git repos (Ensembl VEP_plugins + LOFTEE)
#     are flattened into a single directory after cloning, so VEP's
#     single --dir_plugins flag can find every plugin. The .pm files are
#     COPIED (not symlinked) to avoid Singularity bind-mount breakage:
#     absolute symlinks point at host paths invisible inside containers,
#     while copies are visible everywhere. Critical plugins are verified
#     present after flattening; missing .pm files are flagged as failures
#     rather than silently dropped at runtime.
#     IMPORTANT: in your Snakefile, point VEP at the flat directory:
#       --dir_plugins /tmp/annotation/vep/plugins/flat
#     and set --plugin LoF,loftee_path:/tmp/annotation/vep/plugins/loftee_grch38
#     (the loftee_path is the source-code directory, NOT the plugin_data
#     directory that holds the supporting .fa.gz / .sql / .bw files).
# 10. REVEL conversion produces a bgzipped TSV that VEP's REVEL plugin can
#     actually read. Three transformations are applied to the source CSV:
#       (a) comma-separated -> tab-separated
#       (b) 'chr' prefix added to data rows (REVEL ships bare numeric
#           chromosomes, but our pipeline uses 'chr1', 'chr2', ...)
#       (c) header line prefixed with '#' AND tabix indexed with '-c #'
#           (instead of '-S 1') so 'tabix -h' returns the header line
#     Without (b) tabix queries return nothing for chr-prefixed inputs.
#     Without (c) VEP fails with "Could not read headers" at runtime even
#     though regular queries work fine; the REVEL plugin reads its column
#     names via 'tabix -h', which only returns lines registered as
#     comments.
#     After conversion, two sanity checks run: (i) BRCA1 site returns a
#     record, (ii) 'tabix -h' returns a '#'-prefixed header line. Failure
#     of either is treated as a real error. On re-runs, an existing REVEL
#     file is verified via the same checks before being trusted; old
#     broken files are auto-rebuilt.
# 11. The gnomAD strip+concat helper auto-runs by default after the
#     per-chromosome downloads complete. The helper deletes per-chromosome
#     ~184 GB of source files as it processes them, leaving only the
#     ~30 GB merged stripped VCF that VEP --custom uses. Set
#     SKIP_GNOMAD_PROCESSING=1 in the env to opt out (e.g. when you want
#     to do the heavy processing on a different machine). Failures count
#     as real failures.
# 12. gnomAD per-chromosome downloads are now gated on the absence of the
#     merged stripped file. A re-run with the merged file already present
#     skips the entire 184 GB download cycle (previously it would have
#     re-downloaded every per-chromosome file just to be re-deleted by the
#     helper). Also: any per-chromosome files left over from an interrupted
#     previous run are auto-cleaned when the merged file is detected on
#     re-run, reclaiming disk space without manual intervention.
#
# === ITEMS THAT REQUIRE HUMAN ACTION ===
#   A. Capture-kit BED file - from your wet-lab supplier
#   B. SpliceAI INDEL VCF - via Illumina BaseSpace CLI (free account)
#   C. dbNSFP v5.3.1a - via genos.us (license form per session)
#   D. REVEL - via Google Sites (license form). Auto-converted on re-run.
#   E. SCHEMA / BipEx / ASC TSVs - via web-app "Download" buttons
#   F. LOFTEE supporting data - if peering to personal.broadinstitute.org
#      fails, see end-of-run notes for AWS-staging workaround
#
# Idempotent: non-empty existing files are SKIPped. Atomic: writes to
# .partial and renames on success. Safe to re-run after crashes.
# =============================================================================

set -uo pipefail

# -----------------------------------------------------------------------------
# Single source of truth: pinned versions
# -----------------------------------------------------------------------------
VEP_CACHE_RELEASE="112"           # VEP cache used at runtime; matches container
VEP_PLUGINS_RELEASE="release/115" # plugin code; intentionally newer than cache
                                  # so dbNSFP v5.x parsing logic is available
DBNSFP_VERSION="5.3.1a"           # current academic dbNSFP release
GNOMAD_RELEASE="4.1"              # gnomAD VCF release line
CADD_VERSION="v1.7"               # CADD score version
CLINVAR_DATE="20260420"           # dated ClinVar archive snapshot (see notes)

# -----------------------------------------------------------------------------
# Layout - keep in sync with paths in the Snakefile
# -----------------------------------------------------------------------------
GENOME_DIR=/mnt/data/exome/genome/genome_for_exome_pipe
VARIATION_DIR=/mnt/data/exome/variation/vcf_for_exome_pipe
VEP_CACHE=/mnt/data/exome/annotation/vep/cache_grch38
VEP_PLUGINS=/mnt/data/exome/annotation/vep/plugins
VEP_PLUGIN_DATA=/mnt/data/exome/annotation/vep/plugin_data
VEP_CUSTOM=/mnt/data/exome/annotation/vep/custom
DATA_DIR=/mnt/data/exome/data
CAPTURE_DIR=/mnt/data/exome/annotation/capture
CAPTURE_BED="$CAPTURE_DIR/capture.bed"
GNOMAD_HELPER=/mnt/data/exome/scripts/gnomad_strip_concat.sh

mkdir -p "$GENOME_DIR" "$VARIATION_DIR" \
         "$VEP_CACHE" "$VEP_PLUGINS" "$VEP_PLUGIN_DATA" "$VEP_CUSTOM" \
         "$DATA_DIR" "$CAPTURE_DIR" "$(dirname "$GNOMAD_HELPER")"

# Manifest file: TSV with one row per managed file
MANIFEST_FILE="${MANIFEST_FILE:-/mnt/data/exome/MANIFEST.tsv}"

# TLS opt-out allowlist: hosts where cert validation can be skipped if needed.
# Default: empty (TLS validated for all hosts). To allow a host explicitly,
# set TLS_INSECURE_HOSTS as a colon-separated list of host substrings.
TLS_INSECURE_HOSTS="${TLS_INSECURE_HOSTS:-}"

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
LOG_FILE="${DOWNLOAD_LOG:-$(dirname "$(realpath "$0")")/download_references_$(date +%Y%m%d_%H%M%S).log}"
: > "$LOG_FILE"

log() {
    local level="$1"; shift
    local ts
    ts=$(date '+%Y-%m-%d %H:%M:%S')
    printf '[%s] %-5s %s\n' "$ts" "$level" "$*" | tee -a "$LOG_FILE"
}

section() {
    printf '\n' | tee -a "$LOG_FILE"
    log INFO  "==== $1 ===="
}

declare -i N_SKIP=0 N_OK=0 N_FAIL=0 N_MANUAL=0
FAILED_ITEMS=()
MANUAL_ITEMS=()

# -----------------------------------------------------------------------------
# Manifest helpers
# -----------------------------------------------------------------------------
# Record one file in the manifest. Atomic via write-then-rename so a
# concurrent run can't see a half-written state.
manifest_record() {
    local path="$1"; local status="$2"; local source_url="${3:-}"
    if [[ ! -f "$path" || ! -s "$path" ]]; then
        return 0
    fi
    local size sha256 timestamp
    size=$(stat -c%s "$path" 2>/dev/null || echo 0)
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    # Re-use existing SHA if path+size+mtime unchanged from last run
    local mtime cached
    mtime=$(stat -c%Y "$path" 2>/dev/null || echo 0)
    if [[ -f "$MANIFEST_FILE" ]]; then
        cached=$(awk -F'\t' -v p="$path" -v s="$size" -v m="$mtime" \
            '$1==p && $4==s && $5==m {print $6; exit}' "$MANIFEST_FILE")
    else
        cached=""
    fi
    if [[ -n "$cached" ]]; then
        sha256="$cached"
    else
        sha256=$(sha256sum "$path" 2>/dev/null | awk '{print $1}')
    fi

    # Build new manifest atomically
    local tmp="${MANIFEST_FILE}.new"
    {
        if [[ -f "$MANIFEST_FILE" ]]; then
            head -1 "$MANIFEST_FILE" 2>/dev/null
            grep -v "^${path}	" "$MANIFEST_FILE" 2>/dev/null | tail -n +2
        else
            printf 'path\tstatus\tsource\tsize\tmtime\tsha256\ttimestamp\n'
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$path" "$status" "$source_url" "$size" "$mtime" "$sha256" "$timestamp"
    } > "$tmp"
    mv "$tmp" "$MANIFEST_FILE"
}

# Verify against an upstream checksum file. Returns 0 on success, 1 on mismatch.
# Supported algorithms: md5, sha256.
verify_checksum() {
    local file="$1"; local sumfile="$2"; local algo="${3:-md5}"
    if [[ ! -f "$sumfile" || ! -s "$sumfile" ]]; then
        return 2  # no checksum file - caller decides what to do
    fi

    # Source files often look like:
    #   "<sum>  <filename>" (standard) or just "<sum>" (bare hash)
    local expected=""
    expected=$(awk '/^[a-fA-F0-9]+/ {print $1; exit}' "$sumfile")
    if [[ -z "$expected" ]]; then
        log WARN "  could not parse checksum from $sumfile"
        return 2
    fi

    local actual
    case "$algo" in
        md5)    actual=$(md5sum    "$file" | awk '{print $1}') ;;
        sha256) actual=$(sha256sum "$file" | awk '{print $1}') ;;
        *)      log WARN "  unknown checksum algorithm $algo"; return 2 ;;
    esac

    if [[ "$expected" == "$actual" ]]; then
        log INFO "  $algo OK ($actual)"
        return 0
    else
        log ERROR "  $algo MISMATCH for $file"
        log ERROR "    expected: $expected"
        log ERROR "    actual:   $actual"
        return 1
    fi
}

# -----------------------------------------------------------------------------
# Step 0: tool prerequisites
# -----------------------------------------------------------------------------
section "0. Tool prerequisites"
log INFO "version pins: VEP cache=$VEP_CACHE_RELEASE, VEP plugins=$VEP_PLUGINS_RELEASE,"
log INFO "              dbNSFP=$DBNSFP_VERSION, gnomAD=$GNOMAD_RELEASE,"
log INFO "              CADD=$CADD_VERSION, ClinVar snapshot=$CLINVAR_DATE"
log INFO "manifest:     $MANIFEST_FILE"

NO_INSTALL="${NO_INSTALL:-0}"

# List of tools the script actually needs
REQUIRED_TOOLS="wget curl aria2c git unzip tabix bcftools python3 sha256sum"
MISSING_TOOLS=()
for t in $REQUIRED_TOOLS; do
    command -v "$t" >/dev/null 2>&1 || MISSING_TOOLS+=("$t")
done

if (( ${#MISSING_TOOLS[@]} == 0 )); then
    log INFO "all required tools present, skipping apt-get"
elif [[ "$NO_INSTALL" == "1" ]]; then
    log WARN "missing tools (${MISSING_TOOLS[*]}) but NO_INSTALL=1 set; continuing"
else
    if command -v apt-get >/dev/null 2>&1; then
        log INFO "missing tools: ${MISSING_TOOLS[*]} -- attempting apt install"
        SUDO=""
        [[ $EUID -ne 0 ]] && SUDO="sudo"

        # Only run update if we actually need to install something.
        # Wrap in `timeout` so a hung mirror doesn't stall forever.
        timeout 120 $SUDO apt-get update -qq 2>&1 | tee -a "$LOG_FILE" \
            || log WARN "apt-get update timed out or failed (continuing with cached lists)"

        timeout 300 $SUDO DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
            "${MISSING_TOOLS[@]}" 2>&1 | tee -a "$LOG_FILE" \
            || log WARN "apt-get install reported errors (continuing)"
    else
        log WARN "missing tools (${MISSING_TOOLS[*]}) but apt-get unavailable"
    fi
fi

# Tool detection
HAVE_WGET=0; HAVE_CURL=0; HAVE_ARIA2=0; HAVE_TABIX=0; HAVE_BCFTOOLS=0
HAVE_BGZIP=0; HAVE_GIT=0; HAVE_UNZIP=0; HAVE_PYTHON=0; HAVE_SHA256=0
command -v wget       >/dev/null 2>&1 && HAVE_WGET=1
command -v curl       >/dev/null 2>&1 && HAVE_CURL=1
command -v aria2c     >/dev/null 2>&1 && HAVE_ARIA2=1
command -v tabix      >/dev/null 2>&1 && HAVE_TABIX=1
command -v bcftools   >/dev/null 2>&1 && HAVE_BCFTOOLS=1
command -v bgzip      >/dev/null 2>&1 && HAVE_BGZIP=1
command -v git        >/dev/null 2>&1 && HAVE_GIT=1
command -v unzip      >/dev/null 2>&1 && HAVE_UNZIP=1
command -v python3    >/dev/null 2>&1 && HAVE_PYTHON=1
command -v sha256sum  >/dev/null 2>&1 && HAVE_SHA256=1

log INFO "tool detection:"
for kv in "wget=$HAVE_WGET" "curl=$HAVE_CURL" "aria2c=$HAVE_ARIA2" \
          "tabix=$HAVE_TABIX" "bcftools=$HAVE_BCFTOOLS" "bgzip=$HAVE_BGZIP" \
          "git=$HAVE_GIT" "unzip=$HAVE_UNZIP" "python3=$HAVE_PYTHON" \
          "sha256sum=$HAVE_SHA256"; do
    log INFO "  $kv"
done

if (( HAVE_WGET == 0 && HAVE_CURL == 0 )); then
    log ERROR "neither wget nor curl is available; cannot proceed"
    exit 1
fi
if (( HAVE_PYTHON == 0 )); then
    log WARN "python3 unavailable - post-processing helpers will fail"
fi
if (( HAVE_SHA256 == 0 )); then
    log WARN "sha256sum unavailable - manifest will lack hash entries"
fi

# Set TMPDIR onto the data volume (cloud VMs commonly have tiny /tmp)
DEFAULT_TMP=/mnt/data/exome_tmp
if [[ "${TMPDIR:-}" == "" ]]; then
    mkdir -p "$DEFAULT_TMP"
    export TMPDIR="$DEFAULT_TMP"
    log INFO "set TMPDIR=$TMPDIR (on the big disk)"
else
    log INFO "TMPDIR=$TMPDIR (already set)"
fi

# Version-coherence sanity check. Flags common misalignments at startup
# rather than letting them surface as a runtime failure later.
section "Version coherence sanity check"
coherence_ok=1
case "$DBNSFP_VERSION" in
    5.*)
        if [[ "$VEP_PLUGINS_RELEASE" =~ ^release/11[0-2]$ ]]; then
            log ERROR "version mismatch: dbNSFP $DBNSFP_VERSION (v5+) needs"
            log ERROR "  VEP_PLUGINS_RELEASE >= release/115 (you have $VEP_PLUGINS_RELEASE)"
            coherence_ok=0
        fi
        ;;
    4.*)
        log INFO "dbNSFP v4 series; older plugin branches OK but v115 still works"
        ;;
esac
if (( coherence_ok == 1 )); then
    log INFO "version pins are mutually consistent"
else
    log ERROR "edit the version constants at the top of this script and re-run"
    exit 1
fi

# -----------------------------------------------------------------------------
# Wget / aria2 with TLS opt-in policy
# -----------------------------------------------------------------------------
# By default, TLS certs are validated. To allow a specific host whose cert is
# known-broken, add a substring of its hostname to TLS_INSECURE_HOSTS.
url_is_tls_insecure() {
    local url="$1"
    if [[ -z "$TLS_INSECURE_HOSTS" ]]; then return 1; fi
    local IFS=':'
    for h in $TLS_INSECURE_HOSTS; do
        [[ -z "$h" ]] && continue
        [[ "$url" == *"$h"* ]] && return 0
    done
    return 1
}

wget_rc_meaning() {
    case "$1" in
        1) echo "(generic error)";;
        2) echo "(parse error)";;
        3) echo "(file I/O error)";;
        4) echo "(network failure)";;
        5) echo "(SSL verification failure)";;
        6) echo "(auth failure)";;
        7) echo "(protocol error)";;
        8) echo "(server error response - 4xx/5xx)";;
        *) echo "";;
    esac
}

get_single_wget() {
    local url="$1"; local partial="$2"
    local tls_args=()
    if url_is_tls_insecure "$url"; then
        tls_args+=(--no-check-certificate)
        log WARN "  TLS validation DISABLED for $url (per TLS_INSECURE_HOSTS)"
    fi
    wget -c "${tls_args[@]}" \
         --timeout=300 --read-timeout=300 \
         --tries=3 --waitretry=30 \
         --show-progress --progress=dot:giga \
         -a "$LOG_FILE" \
         "$url" -O "$partial"
}

get_single_aria2() {
    local url="$1"; local partial="$2"
    local dir; dir=$(dirname "$partial")
    local file; file=$(basename "$partial")
    local tls_args=(--check-certificate=true)
    if url_is_tls_insecure "$url"; then
        tls_args=(--check-certificate=false)
        log WARN "  TLS validation DISABLED for $url (per TLS_INSECURE_HOSTS)"
    fi
    aria2c --continue=true --max-tries=5 --retry-wait=30 \
           --max-connection-per-server=4 --split=4 \
           --timeout=300 --connect-timeout=60 \
           --lowest-speed-limit=10K \
           "${tls_args[@]}" \
           --console-log-level=warn --summary-interval=30 \
           --allow-overwrite=true \
           --dir="$dir" --out="$file" \
           "$url" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

get_single() {
    local url="$1"; local target="$2"
    local engine="wget"
    if (( HAVE_ARIA2 == 1 )) && [[ "$url" == *personal.broadinstitute.org* ]]; then
        engine="aria2"
    fi

    local partial="${target}.partial"
    local start_ts elapsed bytes http_status rc
    start_ts=$(date +%s)

    http_status=""
    if (( HAVE_CURL == 1 )); then
        local curl_tls=""
        url_is_tls_insecure "$url" && curl_tls="-k"
        http_status=$(curl -sIL $curl_tls --max-time 30 -o /dev/null -w '%{http_code}' "$url" 2>/dev/null || echo "???")
    fi

    log INFO "  engine=$engine  HTTP=$http_status"
    case "$engine" in
        aria2) get_single_aria2 "$url" "$partial"; rc=$? ;;
        wget)  get_single_wget  "$url" "$partial"; rc=$? ;;
    esac

    elapsed=$(( $(date +%s) - start_ts ))
    bytes=0
    [[ -f "$partial" ]] && bytes=$(stat -c%s "$partial" 2>/dev/null || echo 0)

    if [[ $rc -eq 0 && -s "$partial" ]]; then
        mv "$partial" "$target"
        log INFO "  got $(du -h "$target" | cut -f1) in ${elapsed}s -> $target"
        return 0
    fi

    log ERROR "  download FAILED"
    log ERROR "    URL:          $url"
    log ERROR "    engine:       $engine"
    log ERROR "    HTTP status:  ${http_status:-not probed}"
    if [[ "$engine" == "wget" ]]; then
        log ERROR "    wget exit:    $rc $(wget_rc_meaning $rc)"
    else
        log ERROR "    aria2 exit:   $rc"
    fi
    log ERROR "    elapsed:      ${elapsed}s"
    log ERROR "    bytes saved:  $bytes"
    [[ $bytes -gt 0 ]] && log ERROR "    .partial kept at $partial for resume"
    return 1
}

# -----------------------------------------------------------------------------
# get URL [URL_ALT...] [-- TARGET]
# -----------------------------------------------------------------------------
get() {
    local urls=() target=""
    local saw_sep=0
    for arg in "$@"; do
        if [[ "$arg" == "--" ]]; then saw_sep=1; continue; fi
        if (( saw_sep )); then target="$arg"; else urls+=("$arg"); fi
    done
    if [[ -z "$target" ]]; then
        if (( ${#urls[@]} == 2 )) && [[ ! "${urls[1]}" =~ ^https?:// ]]; then
            target="${urls[1]}"; urls=("${urls[0]}")
        else
            target="$(basename "${urls[0]}")"
        fi
    fi

    if [[ -f "$target" && -s "$target" ]]; then
        log SKIP  "$(printf '%-68s (%s)' "$target" "$(du -h "$target" | cut -f1)")"
        N_SKIP+=1
        manifest_record "$(realpath "$target")" "skipped" "${urls[0]:-}"
        return 0
    fi

    local i=0 n=${#urls[@]} ok=0 used_url=""
    for url in "${urls[@]}"; do
        i=$((i+1))
        if (( n > 1 )); then log INFO "attempt $i/$n for $target"; fi
        if get_single "$url" "$target"; then ok=1; used_url="$url"; break; fi
        if (( i < n )); then log WARN "  trying next mirror..."; fi
    done

    if (( ok == 1 )); then
        N_OK+=1
        manifest_record "$(realpath "$target")" "downloaded" "$used_url"
        return 0
    fi
    log ERROR "  ALL ${n} source(s) failed for $target"
    FAILED_ITEMS+=("$target  (from: ${urls[*]})")
    N_FAIL+=1
    return 1
}

manual_needed() {
    local target="$1"; local instr="$2"
    if [[ -f "$target" && -s "$target" ]]; then
        log SKIP  "$(printf '%-68s (%s)' "$target" "$(du -h "$target" | cut -f1)")"
        N_SKIP+=1
        manifest_record "$(realpath "$target")" "manual_provided" ""
        return 0
    fi
    log WARN  "MANUAL  $target"
    log WARN  "        $instr"
    MANUAL_ITEMS+=("$target")
    N_MANUAL+=1
}

# -----------------------------------------------------------------------------
# Post-processing helpers (counted as real failures on error)
# -----------------------------------------------------------------------------
post_process_step() {
    # post_process_step DESCRIPTION COMMAND [ARGS...]
    # Run the command; if it exits non-zero, count as a failure.
    local desc="$1"; shift
    log INFO "post-processing: $desc"
    if "$@"; then
        log INFO "  OK"
        return 0
    else
        local rc=$?
        log ERROR "  $desc FAILED (exit $rc)"
        FAILED_ITEMS+=("post:$desc")
        N_FAIL+=1
        return 1
    fi
}

extract_if_missing() {
    local tarball="$1"; local marker="$2"
    if [[ -d "$marker" ]]; then
        log SKIP "extraction of $tarball ($marker exists)"
        return 0
    fi
    if [[ ! -f "$tarball" ]]; then
        log ERROR "cannot extract $tarball -- file does not exist"
        FAILED_ITEMS+=("extract:$tarball")
        N_FAIL+=1
        return 1
    fi
    log INFO  "extracting $tarball -> $marker"
    if tar -xzf "$tarball" 2>>"$LOG_FILE"; then
        log INFO "  extraction complete"
        return 0
    else
        log ERROR "  extraction FAILED (see log)"
        FAILED_ITEMS+=("extract:$tarball")
        N_FAIL+=1
        return 1
    fi
}

clone_if_missing() {
    local url="$1"; local dir="$2"; local branch="${3:-}"
    if [[ -d "$dir/.git" ]]; then
        log SKIP "clone $url ($dir exists)"
        return 0
    fi
    log INFO "cloning $url -> $dir${branch:+ (branch: $branch)}"
    local git_cmd=(git clone --depth 1 "$url" "$dir")
    [[ -n "$branch" ]] && git_cmd=(git clone --branch "$branch" --depth 1 "$url" "$dir")
    if "${git_cmd[@]}" 2>>"$LOG_FILE"; then
        log INFO "  clone complete"
        return 0
    else
        log ERROR "  clone FAILED"
        FAILED_ITEMS+=("clone:$dir")
        N_FAIL+=1
        return 1
    fi
}

# Track expected-files for the final completeness check
declare -a EXPECTED_FILES=()
expect_file() { EXPECTED_FILES+=("$1"); }

# =============================================================================
# 1. GENOME
# =============================================================================
section "1. Broad hg38 reference FASTA"
BROAD_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
cd "$GENOME_DIR"
get "$BROAD_URL/Homo_sapiens_assembly38.fasta"     -- resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
get "$BROAD_URL/Homo_sapiens_assembly38.fasta.fai" -- resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai
get "$BROAD_URL/Homo_sapiens_assembly38.dict"      -- resources_broad_hg38_v0_Homo_sapiens_assembly38.dict
expect_file "$GENOME_DIR/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
expect_file "$GENOME_DIR/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai"
expect_file "$GENOME_DIR/resources_broad_hg38_v0_Homo_sapiens_assembly38.dict"

# =============================================================================
# 2. VARIATION
# =============================================================================
section "2. Known-sites VCFs for BQSR"
cd "$VARIATION_DIR"
get "$BROAD_URL/hapmap_3.3.hg38.vcf.gz"     -- resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz
get "$BROAD_URL/hapmap_3.3.hg38.vcf.gz.tbi" -- resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz.tbi
get "$BROAD_URL/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"     -- resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
get "$BROAD_URL/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" -- resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
get "$BROAD_URL/Homo_sapiens_assembly38.dbsnp138.vcf"     -- resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
get "$BROAD_URL/Homo_sapiens_assembly38.dbsnp138.vcf.idx" -- resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.idx
expect_file "$VARIATION_DIR/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
expect_file "$VARIATION_DIR/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz.tbi"
expect_file "$VARIATION_DIR/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
expect_file "$VARIATION_DIR/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
expect_file "$VARIATION_DIR/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
expect_file "$VARIATION_DIR/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.idx"

# =============================================================================
# 3. VEP CACHE
# =============================================================================
section "3. VEP cache (Ensembl release $VEP_CACHE_RELEASE)"
cd "$VEP_CACHE"
get "https://ftp.ensembl.org/pub/release-${VEP_CACHE_RELEASE}/variation/indexed_vep_cache/homo_sapiens_vep_${VEP_CACHE_RELEASE}_GRCh38.tar.gz"
post_process_step "extract VEP cache" \
    extract_if_missing "homo_sapiens_vep_${VEP_CACHE_RELEASE}_GRCh38.tar.gz" "homo_sapiens"
expect_file "$VEP_CACHE/homo_sapiens"

# =============================================================================
# 4. VEP PLUGIN CODE
# -----------------------------------------------------------------------------
# We clone two separate git repos:
#   - Ensembl VEP_plugins (AlphaMissense, REVEL, CADD, dbNSFP, SpliceAI, ...)
#   - LOFTEE (LoF.pm + helper modules)
#
# VEP's --dir_plugins flag accepts only ONE directory, so the .pm files from
# both repos must be visible in a single place. After both clones succeed,
# we symlink all plugin .pm files into a "flat" directory and point the
# Snakefile's --dir_plugins at THAT directory rather than the parent.
#
# Symlinks (not copies) are used so that `git pull` on either source repo
# automatically propagates to the flat directory.
# =============================================================================
section "4. VEP plugin code"
cd "$VEP_PLUGINS"
post_process_step "clone VEP_plugins" \
    clone_if_missing "https://github.com/Ensembl/VEP_plugins.git" ensembl_plugins "$VEP_PLUGINS_RELEASE"
post_process_step "clone loftee" \
    clone_if_missing "https://github.com/konradjk/loftee.git" loftee_grch38 grch38

# Flatten plugin .pm files into one directory that --dir_plugins can read.
# CRITICAL: the Snakefile must point --dir_plugins at $VEP_PLUGINS/flat
# (i.e. /tmp/annotation/vep/plugins/flat by default).
#
# We COPY the .pm files (not symlink) because Singularity bind-mounts
# break absolute symlinks: a symlink whose target uses the host path
# (e.g. /mnt/data/exome/...) will dangle inside the container, where
# only the bind-mounted view (/tmp/annotation/...) is visible. Copies
# avoid this. The .pm files are tiny (~10 KB each), so the disk cost is
# negligible. Re-runs of this script overwrite cleanly.
flatten_plugins() {
    local flat="$VEP_PLUGINS/flat"
    mkdir -p "$flat"

    # Copy Ensembl plugins
    local n_ens=0
    if [[ -d "$VEP_PLUGINS/ensembl_plugins" ]]; then
        for pm in "$VEP_PLUGINS/ensembl_plugins"/*.pm; do
            [[ -e "$pm" ]] || continue
            cp -f "$pm" "$flat/$(basename "$pm")"
            n_ens=$((n_ens + 1))
        done
    fi

    # Copy LOFTEE plugin and compile-time helper scripts.
    # LoF.pm uses Perl require() for helper .pl files such as utr_splice.pl.
    # These must be visible in VEP's @INC, which includes --dir_plugins.
    local n_loftee=0
    if [[ -d "$VEP_PLUGINS/loftee_grch38" ]]; then
        for f in "$VEP_PLUGINS/loftee_grch38"/*.pm \
                "$VEP_PLUGINS/loftee_grch38"/*.pl; do
            [[ -e "$f" ]] || continue
            cp -f "$f" "$flat/$(basename "$f")"
            n_loftee=$((n_loftee + 1))
        done

        # Subdirectories required by LoF.pm at compile/runtime.
        for d in maxEntScan; do
            if [[ -d "$VEP_PLUGINS/loftee_grch38/$d" ]]; then
                cp -a "$VEP_PLUGINS/loftee_grch38/$d" "$flat/"
                n_loftee=$((n_loftee + 1))
            fi
        done
    fi

    log INFO "  flattened $n_ens Ensembl + $n_loftee LOFTEE .pm files -> $flat"

    # Verify the critical plugins are present. A missing .pm file would
    # silently fail at VEP runtime; catch it here instead.
    local missing=()
    for required in AlphaMissense.pm REVEL.pm LoF.pm SpliceAI.pm CADD.pm dbNSFP.pm; do
        if [[ ! -f "$flat/$required" ]]; then
            missing+=("$required")
        fi
    done

    if (( ${#missing[@]} > 0 )); then
        log ERROR "  flat plugin dir missing required .pm files: ${missing[*]}"
        return 1
    fi
    return 0
}

post_process_step "flatten plugin .pm files into single directory" flatten_plugins
expect_file "$VEP_PLUGINS/flat/AlphaMissense.pm"
expect_file "$VEP_PLUGINS/flat/REVEL.pm"
expect_file "$VEP_PLUGINS/flat/LoF.pm"
expect_file "$VEP_PLUGINS/flat/SpliceAI.pm"
expect_file "$VEP_PLUGINS/flat/CADD.pm"
expect_file "$VEP_PLUGINS/flat/dbNSFP.pm"

# =============================================================================
# 5. VEP PLUGIN DATA
# =============================================================================
section "5. VEP plugin data"
cd "$VEP_PLUGIN_DATA"

# --- 5a. AlphaMissense -----------------------------------------------------
get "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
if [[ -f AlphaMissense_hg38.tsv.gz && ! -f AlphaMissense_hg38.tsv.gz.tbi ]]; then
    if (( HAVE_TABIX == 1 )); then
        post_process_step "tabix-index AlphaMissense" \
            tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz
    else
        log ERROR "tabix not installed - AlphaMissense will not be indexed"
        FAILED_ITEMS+=("post:tabix-AlphaMissense (no tabix tool)")
        N_FAIL+=1
    fi
fi
expect_file "$VEP_PLUGIN_DATA/AlphaMissense_hg38.tsv.gz"
expect_file "$VEP_PLUGIN_DATA/AlphaMissense_hg38.tsv.gz.tbi"

# --- 5b. SpliceAI SNV (Ensembl mirror) -------------------------------------
ENSEMBL_SPLICEAI_URL="https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins"
get "$ENSEMBL_SPLICEAI_URL/spliceai_scores.raw.snv.ensembl_mane.grch38.110.vcf.gz" \
    -- spliceai_scores.raw.snv.hg38.vcf.gz
get "$ENSEMBL_SPLICEAI_URL/spliceai_scores.raw.snv.ensembl_mane.grch38.110.vcf.gz.tbi" \
    -- spliceai_scores.raw.snv.hg38.vcf.gz.tbi
expect_file "$VEP_PLUGIN_DATA/spliceai_scores.raw.snv.hg38.vcf.gz"
expect_file "$VEP_PLUGIN_DATA/spliceai_scores.raw.snv.hg38.vcf.gz.tbi"

# --- 5c. SpliceAI indel -- manual ------------------------------------------
manual_needed "spliceai_scores.raw.indel.hg38.vcf.gz" \
    "Use Illumina BaseSpace CLI ('bs'). Accept project 66029966 once via the browser, then: bs file download --id=16534036123 -o $VEP_PLUGIN_DATA"
manual_needed "spliceai_scores.raw.indel.hg38.vcf.gz.tbi" \
    "(same project, file id 16534036125)"
expect_file "$VEP_PLUGIN_DATA/spliceai_scores.raw.indel.hg38.vcf.gz"
expect_file "$VEP_PLUGIN_DATA/spliceai_scores.raw.indel.hg38.vcf.gz.tbi"

# --- 5d. LOFTEE supporting data --------------------------------------------
LOFTEE_URL="https://personal.broadinstitute.org/konradk/loftee_data/GRCh38"
mkdir -p loftee_hg38
cd loftee_hg38
get "$LOFTEE_URL/human_ancestor.fa.gz"      -- human_ancestor.fa.gz
get "$LOFTEE_URL/human_ancestor.fa.gz.fai"  -- human_ancestor.fa.gz.fai
get "$LOFTEE_URL/human_ancestor.fa.gz.gzi"  -- human_ancestor.fa.gz.gzi

if [[ ! -f loftee.sql || ! -s loftee.sql ]]; then
    if get "$LOFTEE_URL/loftee.sql.gz" -- loftee.sql.gz; then
        post_process_step "decompress loftee.sql.gz" \
            gunzip -k loftee.sql.gz
    else
        get "$LOFTEE_URL/loftee.sql" -- loftee.sql
    fi
else
    log SKIP "loftee.sql (already present)"
    N_SKIP+=1
fi

get "$LOFTEE_URL/gerp_conservation_scores.homo_sapiens.GRCh38.bw" \
    -- gerp_conservation_scores.homo_sapiens.GRCh38.bw
cd "$VEP_PLUGIN_DATA"

expect_file "$VEP_PLUGIN_DATA/loftee_hg38/human_ancestor.fa.gz"
expect_file "$VEP_PLUGIN_DATA/loftee_hg38/human_ancestor.fa.gz.fai"
expect_file "$VEP_PLUGIN_DATA/loftee_hg38/human_ancestor.fa.gz.gzi"
expect_file "$VEP_PLUGIN_DATA/loftee_hg38/loftee.sql"
expect_file "$VEP_PLUGIN_DATA/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw"

# --- 5e. CADD --------------------------------------------------------------
CADD_BASE="https://krishna.gs.washington.edu/download/CADD/${CADD_VERSION}/GRCh38"
get "$CADD_BASE/whole_genome_SNVs.tsv.gz"
get "$CADD_BASE/whole_genome_SNVs.tsv.gz.tbi"
get "$CADD_BASE/gnomad.genomes.r4.0.indel.tsv.gz"
get "$CADD_BASE/gnomad.genomes.r4.0.indel.tsv.gz.tbi"
expect_file "$VEP_PLUGIN_DATA/whole_genome_SNVs.tsv.gz"
expect_file "$VEP_PLUGIN_DATA/whole_genome_SNVs.tsv.gz.tbi"
expect_file "$VEP_PLUGIN_DATA/gnomad.genomes.r4.0.indel.tsv.gz"
expect_file "$VEP_PLUGIN_DATA/gnomad.genomes.r4.0.indel.tsv.gz.tbi"

# --- 5f. REVEL: manual download, automated post-processing -----------------
revel_convert() {
    if [[ -f revel_with_transcript_ids.zip && ! -f revel_with_transcript_ids ]]; then
        if (( HAVE_UNZIP == 0 )); then
            log ERROR "  unzip not installed"
            return 1
        fi
        unzip -o revel_with_transcript_ids.zip 2>&1 | tee -a "$LOG_FILE" || return 1
    fi
    if [[ ! -f revel_with_transcript_ids ]]; then
        log ERROR "  revel_with_transcript_ids not present after unzip"
        return 1
    fi
    if (( HAVE_BGZIP == 0 || HAVE_TABIX == 0 )); then
        log ERROR "  bgzip/tabix required for REVEL conversion"
        return 1
    fi
    mkdir -p sort_tmp

    # The REVEL data file needs THREE specific transformations to be usable
    # by VEP's REVEL plugin:
    #
    #   (a) Comma-separated -> tab-separated. The REVEL distribution is a CSV;
    #       VEP plugins read TSV.
    #
    #   (b) Chromosome column gets a 'chr' prefix. REVEL ships chromosomes as
    #       bare numbers ('1', '2', ..., 'X', 'Y'), but our pipeline (matching
    #       GRCh38 reference, gnomAD, VEP cache) uses 'chr1', 'chr2', etc.
    #       Without the prefix, tabix lookups for 'chr17:...' return nothing.
    #
    #   (c) Header line gets a '#' prefix, AND we index with '-c #' instead of
    #       '-S 1'. The REVEL plugin's source code (see ensembl_plugins/REVEL.pm)
    #       reads its column-name header via 'tabix -h ...', which only returns
    #       lines that were registered as comments in the index. Indexing with
    #       '-S 1' (skip first line) makes tabix queries skip the header but
    #       does NOT register it as a comment, so 'tabix -h' returns no header
    #       and the plugin dies with "Could not read headers". Indexing with
    #       '-c #' tells tabix to treat '#'-prefixed lines as comments, which
    #       both skips them in regular queries AND surfaces them in 'tabix -h'.
    #
    # Without (b) the file looks fine to tabix queries but is unusable; without
    # (c) the file looks fine to tabix queries AND to 'tabix -l', but VEP
    # rejects it with a misleading error.
    {
        head -n1 revel_with_transcript_ids | tr "," "\t" | sed 's/^/#/'
        tail -n +2 revel_with_transcript_ids \
            | tr "," "\t" \
            | awk -F'\t' -v OFS='\t' '$3 ~ /^[0-9]+$/ {$1="chr"$1; print}' \
            | sort -k1,1 -k3,3n \
                   --temporary-directory=./sort_tmp \
                   --buffer-size=4G \
                   --parallel=4
    } | bgzip -c -@ 4 > new_tabbed_revel_grch38.tsv.gz \
      && tabix -f -s 1 -b 3 -e 3 -c '#' new_tabbed_revel_grch38.tsv.gz
    local rc=$?
    rm -rf sort_tmp
    if (( rc != 0 )); then
        return $rc
    fi

    # Sanity-check 1: a known REVEL site should now return a record. BRCA1
    # c.181T>G (chr17:43106478) is well-covered. If the lookup is empty,
    # something went wrong with the chromosome prefix or tabix indexing
    # and the file is broken (VEP would fail at runtime anyway).
    if ! tabix new_tabbed_revel_grch38.tsv.gz chr17:43106478-43106478 \
            2>/dev/null | grep -q '^chr17'; then
        log ERROR "  REVEL sanity check FAILED: tabix lookup at known site"
        log ERROR "  (chr17:43106478 BRCA1) returned no record. The conversion"
        log ERROR "  produced a file VEP cannot use. Inspect the source"
        log ERROR "  revel_with_transcript_ids file's column layout and adjust"
        log ERROR "  the conversion accordingly."
        return 1
    fi

    # Sanity-check 2: 'tabix -h' should return the '#'-prefixed header line.
    # This is what the REVEL plugin actually uses to read column names; if
    # this returns empty, VEP will fail with "Could not read headers" at
    # runtime even though the data lookup works.
    if ! tabix -h new_tabbed_revel_grch38.tsv.gz chr17:43106478-43106478 \
            2>/dev/null | head -1 | grep -q '^#'; then
        log ERROR "  REVEL sanity check FAILED: tabix -h returns no header line"
        log ERROR "  The header is not registered as a comment in the tabix index."
        log ERROR "  This means VEP's REVEL plugin will fail with 'Could not read"
        log ERROR "  headers' even though tabix queries work. Conversion bug."
        return 1
    fi
    log INFO "  REVEL sanity checks OK (BRCA1 lookup + tabix -h header)"
    return 0
}

if [[ -f new_tabbed_revel_grch38.tsv.gz && -f new_tabbed_revel_grch38.tsv.gz.tbi ]]; then
    # Existing file: verify it actually works before trusting it. Two failure
    # modes have been seen with files produced by earlier script revisions:
    #
    #   1. No 'chr' prefix on data rows -> tabix lookups for chr17:... return
    #      nothing. Detect via lookup at known BRCA1 site.
    #
    #   2. Header line not '#'-prefixed and tabix indexed with -S 1 instead of
    #      -c '#' -> 'tabix -h' returns no header line, so the REVEL plugin
    #      fails with "Could not read headers" at VEP runtime even though
    #      regular tabix queries work fine. Detect via 'tabix -h ...' check.
    #
    # If either check fails, back up the broken file and rebuild from source.
    revel_ok=1
    if (( HAVE_TABIX == 0 )); then
        revel_ok=0
        revel_skip_reason="tabix not installed; cannot validate"
    elif ! tabix new_tabbed_revel_grch38.tsv.gz chr17:43106478-43106478 \
                2>/dev/null | grep -q '^chr17'; then
        revel_ok=0
        revel_skip_reason="tabix lookup at BRCA1 returns no record (likely no 'chr' prefix on data)"
    elif ! tabix -h new_tabbed_revel_grch38.tsv.gz chr17:43106478-43106478 \
                2>/dev/null | head -1 | grep -q '^#'; then
        revel_ok=0
        revel_skip_reason="tabix -h returns no header line (likely indexed with -S 1 instead of -c '#'); VEP plugin will fail"
    fi

    if (( revel_ok == 1 )); then
        log SKIP "REVEL: already converted and valid ($(du -h new_tabbed_revel_grch38.tsv.gz | cut -f1))"
        N_SKIP+=1
    else
        log WARN "REVEL: existing file failed validation:"
        log WARN "       $revel_skip_reason"
        log WARN "       Backing up and rebuilding."
        mv new_tabbed_revel_grch38.tsv.gz     new_tabbed_revel_grch38.tsv.gz.broken
        mv new_tabbed_revel_grch38.tsv.gz.tbi new_tabbed_revel_grch38.tsv.gz.tbi.broken 2>/dev/null
        if [[ -f revel_with_transcript_ids.zip ]] || [[ -f revel_with_transcript_ids ]]; then
            post_process_step "rebuild REVEL (chr prefix + #header)" revel_convert
        else
            log ERROR "  source revel_with_transcript_ids[.zip] missing; cannot rebuild"
            FAILED_ITEMS+=("post:REVEL rebuild (source ZIP missing)")
            N_FAIL+=1
        fi
    fi
elif [[ -f revel_with_transcript_ids.zip ]] || [[ -f revel_with_transcript_ids ]]; then
    post_process_step "convert REVEL to bgzipped TSV with tabix index" revel_convert
else
    manual_needed "new_tabbed_revel_grch38.tsv.gz" \
        "Get revel_with_transcript_ids.zip from https://sites.google.com/site/revelgenomics/downloads (license form), drop into $VEP_PLUGIN_DATA, re-run this script (auto-converts)."
fi
expect_file "$VEP_PLUGIN_DATA/new_tabbed_revel_grch38.tsv.gz"
expect_file "$VEP_PLUGIN_DATA/new_tabbed_revel_grch38.tsv.gz.tbi"

# --- 5g. dbNSFP -- manual download with auto-verification ------------------
DBNSFP_FILE="dbNSFP${DBNSFP_VERSION}_grch38.gz"
if [[ -f "$DBNSFP_FILE" && -s "$DBNSFP_FILE" ]]; then
    log SKIP "$(printf '%-68s (%s)' "$DBNSFP_FILE" "$(du -h "$DBNSFP_FILE" | cut -f1)")"
    N_SKIP+=1
    if [[ -f "${DBNSFP_FILE}.md5" ]]; then
        if ! verify_checksum "$DBNSFP_FILE" "${DBNSFP_FILE}.md5" md5; then
            FAILED_ITEMS+=("checksum:$DBNSFP_FILE")
            N_FAIL+=1
        fi
    else
        log WARN "  no .md5 file - cannot verify upstream checksum"
    fi
    manifest_record "$(realpath "$DBNSFP_FILE")" "manual_provided" ""
else
    manual_needed "$DBNSFP_FILE" \
        "Fill the form at https://www.dbnsfp.org/download (pick version $DBNSFP_VERSION). Download ${DBNSFP_FILE}, ${DBNSFP_FILE}.tbi, ${DBNSFP_FILE}.md5 into $VEP_PLUGIN_DATA. (v5 retired LRT_pred and FATHMM_pred -- the Snakefile already accounts for this.)"
fi
expect_file "$VEP_PLUGIN_DATA/$DBNSFP_FILE"
expect_file "$VEP_PLUGIN_DATA/${DBNSFP_FILE}.tbi"

# =============================================================================
# 6. VEP CUSTOM
# =============================================================================
section "6. VEP --custom annotations"
cd "$VEP_CUSTOM"

# --- 6a. gnomAD per-chromosome VCFs ----------------------------------------
# Two-stage flow:
#   Stage 1: if the merged stripped VCF is already present, we have nothing
#            to do. Skip the per-chromosome downloads entirely (they would
#            be redundant: ~184 GB of network transfer just to be deleted
#            again by the helper). Also clean up any per-chromosome leftovers
#            that may have survived an interrupted previous run.
#   Stage 2: otherwise, download the 24 per-chromosome files and run the
#            strip+concat helper (which deletes per-chromosome originals as
#            it processes them).
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/${GNOMAD_RELEASE}/vcf/exomes"
GNOMAD_MERGED="$VEP_CUSTOM/gnomad.exomes.v${GNOMAD_RELEASE}.sites.stripped.vcf.gz"

if [[ -f "$GNOMAD_MERGED" && -f "${GNOMAD_MERGED}.tbi" ]]; then
    # Merged file already present: skip downloads and post-processing.
    log SKIP "gnomAD: merged stripped file already present ($(du -h "$GNOMAD_MERGED" | cut -f1))"
    N_SKIP+=1

    # Opportunistic cleanup: if a previous (interrupted) run left per-chromosome
    # files lying around, remove them now. They're useless once the merged file
    # exists and they take up ~184 GB.
    leftover_count=0
    for chr in {1..22} X Y; do
        leftover="$VEP_CUSTOM/gnomad.exomes.v${GNOMAD_RELEASE}.sites.chr${chr}.vcf.bgz"
        if [[ -f "$leftover" ]]; then
            leftover_count=$((leftover_count + 1))
        fi
    done
    if (( leftover_count > 0 )); then
        log INFO "  found $leftover_count per-chromosome leftover(s); reclaiming disk space"
        for chr in {1..22} X Y; do
            rm -f "$VEP_CUSTOM/gnomad.exomes.v${GNOMAD_RELEASE}.sites.chr${chr}.vcf.bgz" \
                  "$VEP_CUSTOM/gnomad.exomes.v${GNOMAD_RELEASE}.sites.chr${chr}.vcf.bgz.tbi"
        done
        log INFO "  removed leftover per-chromosome files"
    fi
else
    # Merged file not present: download per-chromosome files.
    log INFO "downloading gnomAD per-chromosome files (~184 GB total)..."
    for chr in {1..22} X Y; do
        get "$GNOMAD_BASE/gnomad.exomes.v${GNOMAD_RELEASE}.sites.chr${chr}.vcf.bgz"
        get "$GNOMAD_BASE/gnomad.exomes.v${GNOMAD_RELEASE}.sites.chr${chr}.vcf.bgz.tbi"
    done
fi

# Helper script for the strip+concat (downloaded but not auto-run)
if [[ -f "$GNOMAD_MERGED" && -f "${GNOMAD_MERGED}.tbi" ]]; then
    : # merged file already present, skip helper-write step entirely
elif [[ ! -f "$GNOMAD_HELPER" ]]; then
    log INFO "writing gnomAD strip+concat helper to $GNOMAD_HELPER"
    cat > "$GNOMAD_HELPER" <<GNOMAD_EOF
#!/usr/bin/env bash
# =============================================================================
# gnomAD v${GNOMAD_RELEASE} exomes: per-chromosome -> stripped+concatenated bgzipped VCF
# -----------------------------------------------------------------------------
# Why: per-chromosome files total ~184 GB; stripping to AF columns drops to
# ~25-30 GB. VEP --custom needs bgzipped VCF (NOT BCF), so we keep VCF.gz.
# Disk: ~45 GB peak; final output ~30 GB. Time: ~30-60 min on 8-core SSD.
# =============================================================================

set -euo pipefail

CUSTOM_DIR=/mnt/data/exome/annotation/vep/custom
cd "\$CUSTOM_DIR"

# AF columns the pipeline uses (matches the VEP --custom field list)
KEEP="INFO/AF,INFO/AF_nfe,INFO/AF_afr,INFO/AF_amr,INFO/AF_eas,INFO/AF_sas,INFO/AF_fin,INFO/AF_asj,INFO/nhomalt,INFO/AC,INFO/AN"
export KEEP   # CRITICAL: parallel subshells won't see KEEP without export

mkdir -p stripped

strip_one() {
    local chr=\$1
    local in=gnomad.exomes.v${GNOMAD_RELEASE}.sites.chr\${chr}.vcf.bgz
    local out=stripped/gnomad.chr\${chr}.stripped.bcf
    [[ ! -f "\$in" ]] && { echo "[chr\${chr}] MISSING: \$in" >&2; return 1; }
    [[ -s "\$out" ]]  && { echo "[chr\${chr}] already stripped"; return 0; }

    echo "[chr\${chr}] stripping..."
    bcftools annotate -x "^\${KEEP}" --threads 2 -Ob -o "\$out" "\$in" || return 1
    bcftools index "\$out" || return 1

    local in_size out_size
    in_size=\$(stat -c%s "\$in")
    out_size=\$(stat -c%s "\$out")
    if [[ \$out_size -lt \$((in_size / 2)) ]]; then
        echo "  [chr\${chr}] \${in_size}B -> \${out_size}B; deleting original"
        rm -f "\$in" "\${in}.tbi"
    else
        echo "  [chr\${chr}] WARNING: stripped file is suspiciously large; keeping original" >&2
        return 1
    fi
}
export -f strip_one

# Run 8 chromosomes in parallel. --halt soon,fail=1 makes the first failure
# stop the whole batch so we don't keep working on a broken pipeline.
if command -v parallel >/dev/null 2>&1; then
    parallel --halt soon,fail=1 -j 8 strip_one ::: {1..22} X Y
else
    # Manual fallback: track per-chromosome exit codes
    declare -A child_status
    for chr in {1..22} X Y; do
        strip_one \$chr &
        child_status[\$!]=\$chr
        if [[ \$(jobs -r | wc -l) -ge 8 ]]; then wait -n; fi
    done
    failed=0
    for pid in "\${!child_status[@]}"; do
        if ! wait \$pid; then
            echo "ERROR: chromosome \${child_status[\$pid]} failed" >&2
            failed=1
        fi
    done
    [[ \$failed -ne 0 ]] && exit 1
fi

# Concatenate. Brace expansion can't mix numeric range and explicit values,
# so spell out 1..22 and X,Y as two separate brace expressions.
echo "concatenating into bgzipped VCF..."
bcftools concat --threads 4 -Oz \\
    -o gnomad.exomes.v${GNOMAD_RELEASE}.sites.stripped.vcf.gz \\
    stripped/gnomad.chr{1..22}.stripped.bcf \\
    stripped/gnomad.chr{X,Y}.stripped.bcf
bcftools index -t gnomad.exomes.v${GNOMAD_RELEASE}.sites.stripped.vcf.gz

# Validate the concatenated output before declaring success
if [[ ! -s gnomad.exomes.v${GNOMAD_RELEASE}.sites.stripped.vcf.gz \\
   || ! -s gnomad.exomes.v${GNOMAD_RELEASE}.sites.stripped.vcf.gz.tbi ]]; then
    echo "ERROR: final file or index is missing/empty" >&2
    exit 1
fi
record_count=\$(bcftools view -H gnomad.exomes.v${GNOMAD_RELEASE}.sites.stripped.vcf.gz | head -100 | wc -l)
if [[ \$record_count -lt 50 ]]; then
    echo "ERROR: final file has suspiciously few records (\$record_count in first 100 read)" >&2
    exit 1
fi

echo
echo "result:"
ls -lh gnomad.exomes.v${GNOMAD_RELEASE}.sites.stripped.vcf.gz*
echo
echo "header sanity:"
bcftools view -h gnomad.exomes.v${GNOMAD_RELEASE}.sites.stripped.vcf.gz | tail -5

rm -rf stripped/
echo
echo "DONE. Update wes_for_drop.smk:"
echo "  gnomad_v4_exomes = \"\$CUSTOM_DIR/gnomad.exomes.v${GNOMAD_RELEASE}.sites.stripped.vcf.gz\""
GNOMAD_EOF
    chmod +x "$GNOMAD_HELPER"
    log INFO "  helper written: $GNOMAD_HELPER"
fi

# Auto-run the gnomAD strip+concat helper unless the user has opted out.
# The helper produces the single merged stripped VCF that VEP --custom
# expects, AND deletes the per-chromosome ~184 GB of source files as it
# processes each chromosome (peak working disk ~45 GB; final merged file
# ~30 GB). Without this step you'd be sitting on 184 GB of unprocessed
# input until you remembered to run the helper later.
#
# To skip auto-run, set SKIP_GNOMAD_PROCESSING=1 in the environment.
# Useful when:
#   - you want to download files now and post-process during off-hours
#   - the host is short on CPU and you'll move data to a beefier machine
#   - you want to verify per-chromosome integrity manually before stripping
SKIP_GNOMAD_PROCESSING="${SKIP_GNOMAD_PROCESSING:-0}"

if [[ -f "$GNOMAD_MERGED" && -f "${GNOMAD_MERGED}.tbi" ]]; then
    # Either the upstream stage-1 block already SKIPped (merged file pre-existing
    # at script start), or the auto-run we're about to do has already happened
    # in this same run (impossible by control flow, but defensive). Either way,
    # nothing to do here.
    :
elif [[ "$SKIP_GNOMAD_PROCESSING" == "1" ]]; then
    log WARN "SKIP_GNOMAD_PROCESSING=1 set; not auto-running the helper."
    log WARN "Run it manually when ready: bash $GNOMAD_HELPER"
    log WARN "Per-chromosome files (~184 GB) are still in $VEP_CUSTOM"
    MANUAL_ITEMS+=("gnomAD strip+concat (auto-run skipped via env var)")
    N_MANUAL+=1
else
    # Confirm at least chr1's source file exists before invoking the helper -
    # if downloads failed, the helper would loop through every chromosome
    # generating identical "MISSING" errors. One quick check up front is
    # clearer.
    if [[ ! -s "$VEP_CUSTOM/gnomad.exomes.v${GNOMAD_RELEASE}.sites.chr1.vcf.bgz" ]]; then
        log ERROR "gnomAD chr1 source missing; cannot run strip+concat helper"
        log ERROR "  expected: $VEP_CUSTOM/gnomad.exomes.v${GNOMAD_RELEASE}.sites.chr1.vcf.bgz"
        FAILED_ITEMS+=("gnomAD strip+concat (source files missing)")
        N_FAIL+=1
    else
        log INFO "auto-running gnomAD strip+concat helper..."
        log INFO "  expected runtime: ~30-60 minutes on 8-core SSD"
        log INFO "  per-chromosome originals will be deleted as they are processed"
        log INFO "  set SKIP_GNOMAD_PROCESSING=1 next time to skip auto-run"
        post_process_step "gnomAD strip+concat" bash "$GNOMAD_HELPER"
    fi
fi
expect_file "$GNOMAD_MERGED"
expect_file "${GNOMAD_MERGED}.tbi"

# --- 6b. ClinVar GRCh38 (dated archive snapshot) ---------------------------
# NCBI's archive directory holds dated snapshots. Pinning to a specific date
# gives reproducible runs; the previous "latest" symlink could change between
# re-runs. ClinVar archive layout: /pub/clinvar/vcf_GRCh38/archive_2.0/<YYYY>/clinvar_<YYYYMMDD>.vcf.gz
CLINVAR_YEAR="${CLINVAR_DATE:0:4}"
get "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/${CLINVAR_YEAR}/clinvar_${CLINVAR_DATE}.vcf.gz" \
    -- clinvar.vcf.gz
get "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/${CLINVAR_YEAR}/clinvar_${CLINVAR_DATE}.vcf.gz.tbi" \
    -- clinvar.vcf.gz.tbi
expect_file "$VEP_CUSTOM/clinvar.vcf.gz"
expect_file "$VEP_CUSTOM/clinvar.vcf.gz.tbi"

# =============================================================================
# 7. GENE PANELS & CONSTRAINT
# =============================================================================
section "7. Gene panels & gnomAD constraint"
cd "$DATA_DIR"

manual_needed "SCHEMA_gene_results.tsv" \
    "https://schema.broadinstitute.org/results -> 'Download' button -> save into $DATA_DIR. Re-run THIS script - it auto-joins HGNC symbols."
manual_needed "BipEx_gene_results.tsv" \
    "https://bipex.broadinstitute.org/downloads -> 'BipEx gene results (TSV)' -> save into $DATA_DIR. Re-run THIS script - it auto-joins HGNC symbols and filters to the canonical group (BIPEX_GROUP, default 'Bipolar Disorder')."
manual_needed "ASC_gene_results.tsv" \
    "https://asc.broadinstitute.org/downloads -> 'ASC gene results (TSV)' -> save into $DATA_DIR. Re-run THIS script - it auto-joins HGNC symbols."

# DDG2P via PanelApp
PANELAPP_TSV="https://panelapp.genomicsengland.co.uk/panels/484/download/All/"
PANELAPP_JSON="https://panelapp.genomicsengland.co.uk/api/v1/panels/484/?format=json"

ddg2p_json_to_tsv() {
    python3 - <<'PYEOF'
import json, csv, sys
with open("DDG2P_panel.json", encoding="utf-8", errors="replace") as f:
    panel = json.load(f)
with open("DDG2P_panel.tsv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene_symbol", "confidence_level", "mode_of_inheritance", "phenotypes"])
    n = 0
    for g in panel.get("genes", []):
        gd = g.get("gene_data", {}) or {}
        symbol = gd.get("hgnc_symbol") or gd.get("gene_symbol") or ""
        conf = g.get("confidence_level", "")
        moi  = g.get("mode_of_inheritance", "")
        phen = "; ".join(g.get("phenotypes", []) or [])
        if symbol:
            w.writerow([symbol, conf, moi, phen]); n += 1
print(f"wrote {n} genes")
sys.exit(0 if n > 0 else 1)
PYEOF
}

if [[ -f DDG2P_panel.tsv && -s DDG2P_panel.tsv ]]; then
    log SKIP "DDG2P_panel.tsv (already present)"
    N_SKIP+=1
    manifest_record "$(realpath DDG2P_panel.tsv)" "skipped" "$PANELAPP_TSV"
elif get "$PANELAPP_TSV" -- DDG2P_panel.tsv; then
    : # success
else
    log WARN "PanelApp TSV download failed; trying JSON API"
    if get "$PANELAPP_JSON" -- DDG2P_panel.json; then
        if (( HAVE_PYTHON == 1 )); then
            post_process_step "convert DDG2P JSON to TSV" ddg2p_json_to_tsv
        else
            log ERROR "python3 unavailable; cannot convert DDG2P JSON"
            FAILED_ITEMS+=("post:DDG2P JSON->TSV (no python)")
            N_FAIL+=1
        fi
    fi
fi
expect_file "$DATA_DIR/DDG2P_panel.tsv"

# gnomAD constraint metrics
get "https://storage.googleapis.com/gcp-public-data--gnomad/release/${GNOMAD_RELEASE}/constraint/gnomad.v${GNOMAD_RELEASE}.constraint_metrics.tsv" \
    -- gnomad_v${GNOMAD_RELEASE}_constraint_metrics.tsv
expect_file "$DATA_DIR/gnomad_v${GNOMAD_RELEASE}_constraint_metrics.tsv"

# SCHEMA HGNC join
schema_hgnc_join() {
    python3 - <<PYEOF
import csv, sys
ensg_to_hgnc = {}
with open("gnomad_v${GNOMAD_RELEASE}_constraint_metrics.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    fields = r.fieldnames or []
    ensg_col = next((c for c in fields if c in ("gene_id", "ensembl_gene_id")), None)
    hgnc_col = next((c for c in fields if c in ("gene", "symbol", "gene_symbol")), None)
    if not ensg_col or not hgnc_col:
        print(f"missing ENSG/HGNC columns: {fields}", file=sys.stderr); sys.exit(1)
    for row in r:
        ensg = row[ensg_col].split(".")[0]
        hgnc = row[hgnc_col]
        if ensg and ensg not in ensg_to_hgnc:
            ensg_to_hgnc[ensg] = hgnc
n_total = n_matched = 0
with open("SCHEMA_gene_results.tsv") as fin, \
     open("SCHEMA_gene_results_with_hgnc.tsv", "w", newline="") as fout:
    r = csv.reader(fin, delimiter="\t")
    w = csv.writer(fout, delimiter="\t")
    header = next(r)
    w.writerow([header[0], "hgnc_symbol"] + header[1:])
    for row in r:
        n_total += 1
        ensg = row[0].split(".")[0]
        hgnc = ensg_to_hgnc.get(ensg, "")
        if hgnc: n_matched += 1
        w.writerow([row[0], hgnc] + row[1:])
print(f"matched {n_matched}/{n_total} ({100*n_matched/n_total:.1f}%)" if n_total > 0 else "no rows")
sys.exit(0)
PYEOF
}

if [[ -f SCHEMA_gene_results_with_hgnc.tsv && -s SCHEMA_gene_results_with_hgnc.tsv ]]; then
    log SKIP "SCHEMA_gene_results_with_hgnc.tsv (already present)"
    N_SKIP+=1
elif [[ -f SCHEMA_gene_results.tsv && -f gnomad_v${GNOMAD_RELEASE}_constraint_metrics.tsv ]]; then
    if (( HAVE_PYTHON == 1 )); then
        post_process_step "join SCHEMA Ensembl IDs to HGNC symbols" schema_hgnc_join
    else
        log ERROR "python3 unavailable; cannot run SCHEMA HGNC join"
        FAILED_ITEMS+=("post:SCHEMA HGNC join (no python)")
        N_FAIL+=1
    fi
fi
expect_file "$DATA_DIR/SCHEMA_gene_results_with_hgnc.tsv"

# -----------------------------------------------------------------------------
# BipEx and ASC HGNC join + group filter
# -----------------------------------------------------------------------------
# BipEx_gene_results.tsv and ASC_gene_results.tsv from the Broad exome
# browsers are keyed only on Ensembl gene_id (no HGNC symbol column) and
# contain multiple rows per gene -- one per analysis stratum (BipEx: BD /
# BD1 / BD2; ASC: typically just "All", but other strata may exist in
# future releases). tier_candidates.py needs a single row per gene with an
# hgnc_symbol column matching VEP's SYMBOL output. This post-processor:
#
#   - drops summary / non-ENSG rows
#   - filters to the canonical group per panel (configurable below)
#   - joins gene_id -> hgnc_symbol via the gnomAD constraint table
#   - keeps only the ID/group + analysis columns we use downstream
#   - de-duplicates so each gene appears at most once
#
# The two output files (BipEx_gene_results_with_hgnc.tsv and
# ASC_gene_results_with_hgnc.tsv) are what tier_candidates.py reads via
# the --bipex_genes / --asc_genes CLI arguments.
#
# To change which group rows are kept (e.g. BD1-only for a bipolar-I-
# specific gene set), edit the BIPEX_GROUP / ASC_GROUP arrays below.
# The values here are the canonical / paper-cited groups.
BIPEX_GROUP=("Bipolar Disorder")
ASC_GROUP=()  # empty = no filter; ASC currently has only "All"
 
panel_hgnc_join() {
    local input="$1"; local output="$2"; local panel_name="$3"
    local group_filter="$4"  # space-separated literal group names, may be empty
    shift 4
    local report_cols=("$@")
 
    # Export to subprocess env -- bash arrays don't traverse into python heredocs
    INPUT="$input" OUTPUT="$output" PANEL="$panel_name" \
    GROUP_FILTER="$group_filter" \
    REPORT_COLS="${report_cols[*]}" \
    CONSTRAINT_TSV="gnomad_v${GNOMAD_RELEASE}_constraint_metrics.tsv" \
    python3 - <<'PYEOF'
import csv, os, sys
 
panel       = os.environ["PANEL"]
input_path  = os.environ["INPUT"]
output_path = os.environ["OUTPUT"]
constraint  = os.environ["CONSTRAINT_TSV"]
group_filt  = [g for g in os.environ["GROUP_FILTER"].split("\t") if g]
report_cols = [c for c in os.environ["REPORT_COLS"].split() if c]
 
# Build gene_id -> hgnc_symbol mapping from the constraint table
ensg_to_hgnc = {}
with open(constraint) as f:
    r = csv.DictReader(f, delimiter="\t")
    fields = r.fieldnames or []
    ensg_col = next((c for c in fields if c in ("gene_id", "ensembl_gene_id")), None)
    hgnc_col = next((c for c in fields if c in ("gene", "symbol", "gene_symbol")), None)
    if not ensg_col or not hgnc_col:
        print(f"  {panel}: missing ENSG/HGNC columns in constraint: {fields}",
              file=sys.stderr)
        sys.exit(1)
    for row in r:
        ensg = row[ensg_col].split(".")[0]
        hgnc = row[hgnc_col]
        if ensg and ensg not in ensg_to_hgnc and hgnc:
            ensg_to_hgnc[ensg] = hgnc
 
# Read panel file
with open(input_path) as fin:
    r = csv.DictReader(fin, delimiter="\t")
    src_cols = r.fieldnames or []
    if "gene_id" not in src_cols:
        print(f"  {panel}: input has no gene_id column: {src_cols[:6]}",
              file=sys.stderr)
        sys.exit(1)
 
    # Resolve which report_cols are actually present
    keep_report = [c for c in report_cols if c in src_cols]
    missing     = [c for c in report_cols if c not in src_cols]
    if missing:
        print(f"  {panel}: requested columns absent from source "
              f"(will be skipped): {missing}", file=sys.stderr)
 
    # Extra columns that always go through if present
    extra_id_cols = [c for c in ("group", "n_cases", "n_controls") if c in src_cols]
 
    out_cols = ["hgnc_symbol", "gene_id"] + extra_id_cols + keep_report
 
    rows_out = []
    seen_genes = set()
    n_raw = n_after_ensg = n_after_group = n_with_hgnc = n_unique = 0
 
    for row in r:
        n_raw += 1
        ensg = row["gene_id"].split(".")[0]
        if not ensg.startswith("ENSG"):
            continue
        n_after_ensg += 1
        if group_filt:
            if row.get("group", "") not in group_filt:
                continue
        n_after_group += 1
        hgnc = ensg_to_hgnc.get(ensg, "")
        if not hgnc:
            continue
        n_with_hgnc += 1
        if hgnc in seen_genes:
            continue
        seen_genes.add(hgnc)
        n_unique += 1
        rows_out.append(
            [hgnc, ensg] + [row.get(c, "") for c in extra_id_cols]
            + [row.get(c, "") for c in keep_report]
        )
 
with open(output_path, "w", newline="") as fout:
    w = csv.writer(fout, delimiter="\t")
    w.writerow(out_cols)
    w.writerows(rows_out)
 
print(f"  {panel}: raw={n_raw} -> ensg-only={n_after_ensg} -> "
      f"group-filt={n_after_group} -> hgnc-joined={n_with_hgnc} -> "
      f"unique-genes={n_unique}", file=sys.stderr)
print(f"  {panel}: wrote {output_path} ({len(out_cols)} cols)",
      file=sys.stderr)
sys.exit(0 if n_unique > 0 else 1)
PYEOF
}
 
bipex_hgnc_join() {
    # Delimiter for GROUP_FILTER is a tab so multi-word group names like
    # "Bipolar Disorder" survive. Multiple groups would join with tabs.
    local groups; groups=$(IFS=$'\t'; printf '%s' "${BIPEX_GROUP[*]}")
    panel_hgnc_join \
        "BipEx_gene_results.tsv" \
        "BipEx_gene_results_with_hgnc.tsv" \
        "BipEx" \
        "$groups" \
        ptv_case_count ptv_control_count \
        ptv_fisher_gnom_non_psych_pval ptv_fisher_gnom_non_psych_OR \
        damaging_missense_case_count damaging_missense_control_count \
        damaging_missense_fisher_gnom_non_psych_pval \
        damaging_missense_fisher_gnom_non_psych_OR
}
 
asc_hgnc_join() {
    local groups; groups=$(IFS=$'\t'; printf '%s' "${ASC_GROUP[*]}")
    panel_hgnc_join \
        "ASC_gene_results.tsv" \
        "ASC_gene_results_with_hgnc.tsv" \
        "ASC" \
        "$groups" \
        xcase_dn_ptv xcont_dn_ptv \
        xcase_dn_misb xcont_dn_misb \
        xcase_dn_misa xcont_dn_misa \
        qval
}
 
# Run BipEx join (gated on python availability and source files present)
if [[ -f BipEx_gene_results_with_hgnc.tsv && -s BipEx_gene_results_with_hgnc.tsv ]]; then
    log SKIP "BipEx_gene_results_with_hgnc.tsv (already present)"
    N_SKIP+=1
elif [[ -f BipEx_gene_results.tsv && -f gnomad_v${GNOMAD_RELEASE}_constraint_metrics.tsv ]]; then
    if (( HAVE_PYTHON == 1 )); then
        post_process_step "join BipEx Ensembl IDs to HGNC symbols (group=${BIPEX_GROUP[*]})" \
            bipex_hgnc_join
    else
        log ERROR "python3 unavailable; cannot run BipEx HGNC join"
        FAILED_ITEMS+=("post:BipEx HGNC join (no python)")
        N_FAIL+=1
    fi
fi
expect_file "$DATA_DIR/BipEx_gene_results_with_hgnc.tsv"
 
# Run ASC join
if [[ -f ASC_gene_results_with_hgnc.tsv && -s ASC_gene_results_with_hgnc.tsv ]]; then
    log SKIP "ASC_gene_results_with_hgnc.tsv (already present)"
    N_SKIP+=1
elif [[ -f ASC_gene_results.tsv && -f gnomad_v${GNOMAD_RELEASE}_constraint_metrics.tsv ]]; then
    if (( HAVE_PYTHON == 1 )); then
        _asc_msg="join ASC Ensembl IDs to HGNC symbols"
        if (( ${#ASC_GROUP[@]} > 0 )); then
            _asc_msg="$_asc_msg (group=${ASC_GROUP[*]})"
        fi
        post_process_step "$_asc_msg" asc_hgnc_join
        unset _asc_msg
    else
        log ERROR "python3 unavailable; cannot run ASC HGNC join"
        FAILED_ITEMS+=("post:ASC HGNC join (no python)")
        N_FAIL+=1
    fi
fi
expect_file "$DATA_DIR/ASC_gene_results_with_hgnc.tsv"

# =============================================================================
# 8. CAPTURE BED -- cannot be automated
# =============================================================================
section "8. Capture-kit BED"
manual_needed "$CAPTURE_BED" \
    "Get GRCh38 target-region BED from your wet-lab supplier or kit vendor (Agilent SureSelect: register at https://earray.chem.agilent.com/suredesign/). Place at $CAPTURE_BED"
expect_file "$CAPTURE_BED"

# =============================================================================
# Final completeness check + manifest update
# =============================================================================
section "Completeness validation"
declare -a MISSING=()
declare -a INDEX_PROBLEMS=()

for f in "${EXPECTED_FILES[@]}"; do
    if [[ ! -e "$f" ]]; then
        MISSING+=("$f")
        continue
    fi
    if [[ -d "$f" ]]; then
        continue
    fi
    if [[ ! -s "$f" ]]; then
        INDEX_PROBLEMS+=("$f (zero-byte)")
        continue
    fi
    # Refresh manifest entry (records current SHA-256 if file changed)
    manifest_record "$f" "validated" ""
done

# Index-pair check: any *.vcf.gz / *.tsv.gz must have a *.tbi alongside.
# *.vcf must have *.idx. *.fasta must have *.fai and *.dict.
for f in "${EXPECTED_FILES[@]}"; do
    [[ -d "$f" ]] && continue
    [[ ! -s "$f" ]] && continue
    case "$f" in
        *.vcf.gz|*.tsv.gz|*.vcf.bgz)
            if [[ ! -s "${f}.tbi" ]]; then
                INDEX_PROBLEMS+=("missing index: ${f}.tbi")
            fi
            ;;
        *.vcf)
            if [[ ! -s "${f}.idx" ]]; then
                INDEX_PROBLEMS+=("missing index: ${f}.idx")
            fi
            ;;
        *.fasta|*.fa)
            if [[ ! -s "${f}.fai" ]]; then
                INDEX_PROBLEMS+=("missing index: ${f}.fai")
            fi
            ;;
    esac
done

if (( ${#MISSING[@]} > 0 )); then
    log ERROR "${#MISSING[@]} expected file(s) missing:"
    for m in "${MISSING[@]}"; do log ERROR "  - $m"; done
    N_FAIL+=$((${#MISSING[@]}))
    FAILED_ITEMS+=("missing: ${#MISSING[@]} files (see log)")
fi

if (( ${#INDEX_PROBLEMS[@]} > 0 )); then
    log ERROR "${#INDEX_PROBLEMS[@]} index/empty-file problem(s):"
    for ip in "${INDEX_PROBLEMS[@]}"; do log ERROR "  - $ip"; done
    N_FAIL+=$((${#INDEX_PROBLEMS[@]}))
    FAILED_ITEMS+=("indexes: ${#INDEX_PROBLEMS[@]} problems (see log)")
fi

if (( ${#MISSING[@]} == 0 && ${#INDEX_PROBLEMS[@]} == 0 )); then
    log INFO "all ${#EXPECTED_FILES[@]} expected files present and indexed"
fi

# =============================================================================
# Summary
# =============================================================================
section "Summary"
log INFO "results: $N_OK downloaded, $N_SKIP skipped, $N_FAIL failed, $N_MANUAL need manual action"

if (( N_FAIL > 0 )); then
    log ERROR ""
    log ERROR "Failed items (${#FAILED_ITEMS[@]}):"
    for item in "${FAILED_ITEMS[@]}"; do log ERROR "  - $item"; done
fi

if (( N_MANUAL > 0 )); then
    log WARN ""
    log WARN "Files requiring manual action (${#MANUAL_ITEMS[@]}):"
    for item in "${MANUAL_ITEMS[@]}"; do log WARN "  - $item"; done
fi

cat <<EOF | tee -a "$LOG_FILE"

==============================================================================
NEXT STEPS - things requiring human action
==============================================================================

A. CAPTURE-KIT BED (BLOCKER for variant-calling rules)
   Path: $CAPTURE_BED
   From your wet-lab supplier; vendor URL printed above.

B. SpliceAI INDEL VCF
   Path:  $VEP_PLUGIN_DATA/spliceai_scores.raw.indel.hg38.vcf.gz
   1. Visit https://basespace.illumina.com/, accept project 66029966
   2. Install bs CLI; bs auth
   3. cd $VEP_PLUGIN_DATA
      bs file download --id=16534036123 -o .
      bs file download --id=16534036125 -o .

C. dbNSFP v$DBNSFP_VERSION
   Path: $VEP_PLUGIN_DATA/dbNSFP${DBNSFP_VERSION}_grch38.gz (+ .tbi + .md5)
   1. https://www.dbnsfp.org/download (license form)
   2. Pick version $DBNSFP_VERSION, download three files into $VEP_PLUGIN_DATA
   IMPORTANT: dbNSFP v5 retired LRT_pred and FATHMM_pred. The pipeline
   already accounts for this.

D. REVEL
   Path: $VEP_PLUGIN_DATA/new_tabbed_revel_grch38.tsv.gz (+ .tbi)
   1. https://sites.google.com/site/revelgenomics/downloads (license form)
   2. Drop revel_with_transcript_ids.zip into $VEP_PLUGIN_DATA
   3. Re-run THIS script - it auto-converts.

E. SCHEMA / BipEx / ASC TSVs
   Path: $DATA_DIR/{SCHEMA,BipEx,ASC}_gene_results.tsv
   Web-app "Download" buttons:
     SCHEMA: https://schema.broadinstitute.org/results
     BipEx:  https://bipex.broadinstitute.org/results
     ASC:    https://asc.broadinstitute.org/results
   The script auto-creates SCHEMA_gene_results_with_hgnc.tsv on next run.

F. gnomAD strip + concat (auto-runs by default)
   The strip+concat helper auto-runs at the end of section 6 unless you
   set SKIP_GNOMAD_PROCESSING=1 in the environment. It deletes the
   per-chromosome ~184 GB of source files as it processes each chromosome
   and produces the merged ~30 GB stripped VCF that VEP --custom uses.
   Runtime: ~30-60 min on 8-core SSD. To run manually if you skipped
   auto-run:
     bash $GNOMAD_HELPER

G. LOFTEE network workaround (if section 5d failed)
   personal.broadinstitute.org has poor peering with some networks.
   If aria2 couldn't fetch, stage via an AWS instance:
     # On AWS:
     mkdir -p loftee_hg38 && cd loftee_hg38
     for f in human_ancestor.fa.gz human_ancestor.fa.gz.fai \\
              human_ancestor.fa.gz.gzi loftee.sql.gz \\
              gerp_conservation_scores.homo_sapiens.GRCh38.bw; do
         aria2c -x 4 -s 4 \\
             https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/\$f
     done
     gunzip -k loftee.sql.gz && tar czf loftee_hg38.tar.gz *
     # rsync / s3 cp the tarball back to the target host

==============================================================================
ALSO AUTOMATED
==============================================================================
  - apt-package installation (when sudo+apt available)
  - VEP plugin .pm flattening into single directory at:
      $VEP_PLUGINS/flat
    Use --dir_plugins $VEP_PLUGINS/flat in the Snakefile (NOT the parent
    $VEP_PLUGINS, which contains no .pm files directly).
  - AlphaMissense tabix index
  - SpliceAI SNV: Ensembl mirror, no auth
  - LOFTEE supporting data: aria2 from Broad
  - DDG2P: TSV with JSON->TSV fallback
  - SCHEMA HGNC join: auto-runs once SCHEMA TSV is in place
  - REVEL conversion: auto-runs once you place the ZIP
  - dbNSFP md5 verification: auto-runs once .md5 is in place
  - gnomAD strip+concat: auto-runs after per-chromosome downloads complete;
      deletes ~184 GB of per-chromosome originals as it processes them and
      produces the ~30 GB merged stripped VCF. Set SKIP_GNOMAD_PROCESSING=1
      to skip the auto-run (in which case run "bash $GNOMAD_HELPER" later).
  - Manifest at $MANIFEST_FILE (paths, sizes, sha256, source URLs, status)
  - Index-completeness validation
  - Version-coherence sanity check at startup

==============================================================================
SNAKEFILE WIRING NOTE - LOFTEE has TWO directories
==============================================================================
  Source code (.pm files):    $VEP_PLUGINS/loftee_grch38/
  Supporting data (.fa.gz/.sql/.bw):  $VEP_PLUGIN_DATA/loftee_hg38/

  These are easy to confuse. The VEP --plugin LoF arguments must match:
    --plugin LoF,loftee_path:$VEP_PLUGINS/loftee_grch38,human_ancestor_fa:$VEP_PLUGIN_DATA/loftee_hg38/human_ancestor.fa.gz,conservation_file:$VEP_PLUGIN_DATA/loftee_hg38/loftee.sql,gerp_bigwig:$VEP_PLUGIN_DATA/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw

  loftee_path: points at the source-code directory (so LOFTEE can find its
  helper Perl modules). human_ancestor_fa, conservation_file, gerp_bigwig
  point at the supporting-data files. Mixing these up makes LOFTEE silently
  fail to register, which then makes Tier A's HC-LoF filter unreachable.

==============================================================================
Logs and manifest
==============================================================================
  Log:      $LOG_FILE
  Manifest: $MANIFEST_FILE

  Re-running this script is idempotent - it retries only failures and
  picks up post-processing for any newly-supplied manual files. The
  manifest gives you a per-file audit trail (size, sha256, source URL,
  validation status) so you can confirm content stability across runs.

EOF

exit $(( N_FAIL > 0 ? 1 : 0 ))
