#!/usr/bin/env bash
set -uo pipefail

# ============================
# Config
# ============================
AGAT_ENV="${AGAT_ENV:-/home/yangy/anaconda3/envs/agat_env}"
AGAT_PERL="${AGAT_PERL:-${AGAT_ENV}/bin/perl}"

AGAT_CONVERT="${AGAT_CONVERT:-${AGAT_ENV}/bin/agat_convert_sp_gxf2gxf.pl}"
AGAT_FIX_OVERLAP="${AGAT_FIX_OVERLAP:-${AGAT_ENV}/bin/agat_sp_fix_overlaping_genes.pl}"
AGAT_KEEP_LONGEST="${AGAT_KEEP_LONGEST:-${AGAT_ENV}/bin/agat_sp_keep_longest_isoform.pl}"
AGAT_FIX_PHASE="${AGAT_FIX_PHASE:-${AGAT_ENV}/bin/agat_sp_fix_cds_phases.pl}"

# Inputs/Outputs
GFF_DIR="."
GFF_PATTERN="*.gff"    # change if needed
OUTDIR="out"
DRYRUN=0

# Optional phase fix (needs genome fasta). If you don't want phase fixing at all: --no_phase
ENABLE_PHASE=1
GENOME_DIR=""  # empty means same dir as GFF
GENOME_SUFFIXES=("fa" "fasta" "fna" "fa.gz" "fasta.gz" "fna.gz")

# ============================
# Helpers
# ============================
die() { echo "[ERROR] $*" >&2; exit 1; }

count_features() {
  local gff="$1"
  local type="$2"
  awk -F'\t' -v t="$type" 'BEGIN{c=0} $0 !~ /^#/ && $3==t {c++} END{print c+0}' "$gff" 2>/dev/null || echo "0"
}

resolve_genome() {
  local sample="$1"
  local gff_path="$2"
  local gff_dir; gff_dir="$(dirname "$gff_path")"
  local search_dir="$gff_dir"
  [[ -n "$GENOME_DIR" ]] && search_dir="$GENOME_DIR"

  local suf cand
  for suf in "${GENOME_SUFFIXES[@]}"; do
    cand="${search_dir}/${sample}.${suf}"
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  echo ""
}

log_msg() {
  # log_msg <message>
  local msg="$1"
  echo "$msg" | tee -a "$CURRENT_LOG" >/dev/null
}

run_cmd() {
  # run_cmd <step_name> <cmd_string>
  local step="$1"
  local cmd="$2"

  CURRENT_STEP="$step"
  CURRENT_CMD="$cmd"

  log_msg "+ [$CURRENT_SAMPLE][$step] $cmd"
  if [[ "$DRYRUN" -eq 1 ]]; then
    return 0
  fi

  eval "$cmd"
  local rc=$?
  if [[ $rc -ne 0 ]]; then
    log_msg "[FAIL] [$CURRENT_SAMPLE][$step] rc=$rc"
    return $rc
  fi
  return 0
}

run_if_missing() {
  # run_if_missing <step_name> <outfile> <cmd_string>
  local step="$1"
  local outfile="$2"
  local cmd="$3"

  if [[ -s "$outfile" ]]; then
    log_msg "[SKIP] [$CURRENT_SAMPLE][$step] exists: $outfile"
    return 0
  fi
  run_cmd "$step" "$cmd"
}

write_summary_line() {
  # write_summary_line <status> <failed_step> <phase_used> <genome>
  local status="$1"
  local failed_step="$2"
  local phase_used="$3"
  local genome="$4"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$CURRENT_SAMPLE" "$status" "$failed_step" \
    "${genes_in:-NA}" "${mrna_in:-NA}" \
    "${genes_norm:-NA}" "${mrna_norm:-NA}" \
    "${genes_fixov:-NA}" "${mrna_fixov:-NA}" \
    "${genes_longest:-NA}" "${mrna_longest:-NA}" \
    "${genes_final:-NA}" "${mrna_final:-NA}" \
    "${phase_used}\t${genome}" >> "$SUMMARY"
}

process_one_sample() {
  local gff="$1"

  local bn sample sample_out logfile genome phase_used final_gff
  bn="$(basename "$gff")"
  sample="${bn%.*}"
  sample_out="$OUTDIR/$sample"
  mkdir -p "$sample_out"

  logfile="$sample_out/run.log"
  : > "$logfile"

  # expose globals for logging
  CURRENT_SAMPLE="$sample"
  CURRENT_LOG="$logfile"
  CURRENT_STEP=""
  CURRENT_CMD=""

  # init counters
  genes_in="NA"; mrna_in="NA"
  genes_norm="NA"; mrna_norm="NA"
  genes_fixov="NA"; mrna_fixov="NA"
  genes_longest="NA"; mrna_longest="NA"
  genes_final="NA"; mrna_final="NA"
  phase_used="0"

  log_msg "=============================="
  log_msg "Sample: $sample"
  log_msg "GFF: $gff"
  log_msg "Start: $(date +'%F %T')"

  genome="$(resolve_genome "$sample" "$gff")"
  if [[ "$ENABLE_PHASE" -eq 1 && -n "$genome" ]]; then
    log_msg "Genome: $genome"
  else
    log_msg "Genome: NOT USED (phase fixing disabled or genome not found)"
    genome="NA"
  fi

  genes_in="$(count_features "$gff" "gene")"
  mrna_in="$(count_features "$gff" "mRNA")"
  log_msg "Counts input: gene=$genes_in mRNA=$mrna_in"

  # Step 0: convert/normalize
  local cmd0
  cmd0="\"$AGAT_PERL\" \"$AGAT_CONVERT\" -g \"$gff\" -o \"$sample_out/00.norm.gff3\""
  if ! run_if_missing "convert" "$sample_out/00.norm.gff3" "$cmd0"; then
    return 10
  fi
  genes_norm="$(count_features "$sample_out/00.norm.gff3" "gene")"
  mrna_norm="$(count_features "$sample_out/00.norm.gff3" "mRNA")"
  log_msg "Counts 00.norm: gene=$genes_norm mRNA=$mrna_norm"

  # Step 1: merge CDS-overlap genes
  local cmd1
  cmd1="\"$AGAT_FIX_OVERLAP\" -f \"$sample_out/00.norm.gff3\" -o \"$sample_out/01.fix_cds_overlap.gff3\""
  if ! run_if_missing "fix_overlap" "$sample_out/01.fix_cds_overlap.gff3" "$cmd1"; then
    return 11
  fi
  genes_fixov="$(count_features "$sample_out/01.fix_cds_overlap.gff3" "gene")"
  mrna_fixov="$(count_features "$sample_out/01.fix_cds_overlap.gff3" "mRNA")"
  log_msg "Counts 01.fix: gene=$genes_fixov mRNA=$mrna_fixov"

  # Step 2: keep longest isoform
  local cmd2
  cmd2="\"$AGAT_KEEP_LONGEST\" --gff \"$sample_out/01.fix_cds_overlap.gff3\" -o \"$sample_out/02.longest_isoform.gff3\""
  if ! run_if_missing "keep_longest" "$sample_out/02.longest_isoform.gff3" "$cmd2"; then
    return 12
  fi
  genes_longest="$(count_features "$sample_out/02.longest_isoform.gff3" "gene")"
  mrna_longest="$(count_features "$sample_out/02.longest_isoform.gff3" "mRNA")"
  log_msg "Counts 02.longest: gene=$genes_longest mRNA=$mrna_longest"

  final_gff="$sample_out/02.longest_isoform.gff3"

  # Step 3 (optional): fix CDS phases
  if [[ "$ENABLE_PHASE" -eq 1 && "$genome" != "NA" ]]; then
    local cmd3
    cmd3="\"$AGAT_FIX_PHASE\" --gff \"$sample_out/02.longest_isoform.gff3\" --fa \"$genome\" -o \"$sample_out/03.phase_fixed.gff3\""
    if ! run_if_missing "fix_phase" "$sample_out/03.phase_fixed.gff3" "$cmd3"; then
      return 13
    fi
    final_gff="$sample_out/03.phase_fixed.gff3"
    phase_used="1"
  fi

  final_name="$sample_out/${sample}.final.gff3"
  cp "$final_gff" "$final_name"

  genes_final="$(count_features "$final_name" "gene")"
  mrna_final="$(count_features "$final_name" "mRNA")"
  log_msg "Counts ${sample}.final.gff3: gene=$genes_final mRNA=$mrna_final"

  log_msg "End: $(date +'%F %T')"
  return 0
}

# ============================
# Args
# ============================
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--gff_dir) GFF_DIR="$2"; shift 2;;
    -p|--pattern) GFF_PATTERN="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -g|--genome_dir) GENOME_DIR="$2"; shift 2;;
    --agat_env)
      AGAT_ENV="$2"
      AGAT_PERL="${AGAT_ENV}/bin/perl"
      AGAT_CONVERT="${AGAT_ENV}/bin/agat_convert_sp_gxf2gxf.pl"
      AGAT_FIX_OVERLAP="${AGAT_ENV}/bin/agat_sp_fix_overlaping_genes.pl"
      AGAT_KEEP_LONGEST="${AGAT_ENV}/bin/agat_sp_keep_longest_isoform.pl"
      AGAT_FIX_PHASE="${AGAT_ENV}/bin/agat_sp_fix_cds_phases.pl"
      shift 2;;
    --no_phase) ENABLE_PHASE=0; shift 1;;
    --dryrun) DRYRUN=1; shift 1;;
    -h|--help)
      cat <<EOF
Usage:
  $(basename "$0") [-i DIR] [-p PATTERN] [-o OUTDIR] [-g GENOME_DIR] [--agat_env DIR] [--no_phase] [--dryrun]

Behavior:
  - Run AGAT only: convert -> fix_overlap -> keep_longest -> (optional) fix_cds_phases
  - Continue to next sample on failure
  - Per-sample log: out/<sample>/run.log
  - Failure list: out/failed.list
  - Summary: out/summary_agat.tsv (always has lines, even for FAIL)
EOF
      exit 0;;
    *) die "Unknown argument: $1";;
  esac
done

# ============================
# Sanity
# ============================
[[ -x "$AGAT_PERL" ]] || die "AGAT perl not found/executable: $AGAT_PERL"
[[ -x "$AGAT_CONVERT" ]] || die "AGAT convert script not found: $AGAT_CONVERT"
[[ -x "$AGAT_FIX_OVERLAP" ]] || die "AGAT overlap script not found: $AGAT_FIX_OVERLAP"
[[ -x "$AGAT_KEEP_LONGEST" ]] || die "AGAT keep-longest script not found: $AGAT_KEEP_LONGEST"
[[ "$ENABLE_PHASE" -eq 0 || -x "$AGAT_FIX_PHASE" ]] || die "AGAT fix-phase script not found: $AGAT_FIX_PHASE"

mkdir -p "$OUTDIR"

shopt -s nullglob
GFFS=( "$GFF_DIR"/$GFF_PATTERN )
shopt -u nullglob
[[ ${#GFFS[@]} -gt 0 ]] || die "No GFF files found: ${GFF_DIR}/${GFF_PATTERN}"

FAILED_LIST="$OUTDIR/failed.list"
SUMMARY="$OUTDIR/summary_agat.tsv"

: > "$FAILED_LIST"
printf "sample\tstatus\tfailed_step\tgenes_in\tmrna_in\tgenes_norm\tmrna_norm\tgenes_fixov\tmrna_fixov\tgenes_longest\tmrna_longest\tgenes_final\tmrna_final\tphase_used\tgenome\n" > "$SUMMARY"

# ============================
# Main loop
# ============================
for gff in "${GFFS[@]}"; do
  # pre-set globals for summary fallback
  CURRENT_SAMPLE=""
  CURRENT_LOG=""
  CURRENT_STEP=""
  CURRENT_CMD=""

  if process_one_sample "$gff"; then
    # OK summary
    # phase_used/genome are captured in process_one_sample scope; we re-read from log variables:
    # easiest: infer from final files existence
    phase_used="0"
    [[ -s "$OUTDIR/$CURRENT_SAMPLE/03.phase_fixed.gff3" ]] && phase_used="1"
    genome_note="NA"
    # if phase_used==1, genome was used; read from log (best-effort)
    # keep genome in summary as NA if not present
    write_summary_line "OK" "-" "$phase_used" "$genome_note"
  else
    rc=$?
    # Determine failed step from CURRENT_STEP
    failed_step="${CURRENT_STEP:-unknown}"
    echo -e "${CURRENT_SAMPLE:-NA}\t${failed_step}\trc=${rc}\tcmd=${CURRENT_CMD:-NA}\tlog=${CURRENT_LOG:-NA}" >> "$FAILED_LIST"

    # Always append a FAIL summary line (counts might be NA or partially filled)
    phase_used="0"
    write_summary_line "FAIL" "$failed_step" "$phase_used" "NA"

    # Also print a concise message to stderr
    echo "[WARN] sample=${CURRENT_SAMPLE:-NA} failed_step=${failed_step} rc=${rc}. See ${CURRENT_LOG:-NA}" >&2
    continue
  fi
done

echo "All samples completed."
echo "Summary: $SUMMARY"
echo "Failed list: $FAILED_LIST"
