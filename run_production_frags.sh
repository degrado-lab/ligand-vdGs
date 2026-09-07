#!/bin/bash
# Generate and submit the full production vdG-library fleet, longest job first.
#
# MODE IS REQUIRED, and deliberately has no default. The two modes differ in
# what they do to an existing library, and the destructive reading is the one a
# forgotten flag would silently select:
#
#   --mode threshold-plus-include
#       Builds every fragment over $MIN_INSTANCES, PLUS anything in $INCLUDE.
#       Demands an empty library: this is a full build, not an addition.
#
#   --mode include-only
#       Builds ONLY the fragments named in $INCLUDE, ignoring the threshold.
#       Demands an existing library, and leaves what is already in it alone.
#       This is the top-up: use it to add a fragment you realised you needed
#       after the main build.
#
# The names say what gets built rather than when you run it, because "full" and
# "top-up" do not by themselves tell you whether $INCLUDE is also honoured. It
# is, in both -- the difference is only whether the threshold pass runs too.
#
# Usage:
#   ./run_production_frags.sh --mode threshold-plus-include --dry-run
#   ./run_production_frags.sh --mode threshold-plus-include
#   MAX_H_RT=24:00:00 ./run_production_frags.sh --mode threshold-plus-include
#
#   # main build, plus fragments you need whatever the sampling estimate said
#   INCLUDE='cnnnn' ./run_production_frags.sh --mode threshold-plus-include
#
#   # later: add fragments to the existing library, no threshold pass
#   INCLUDE='cnnnn' ./run_production_frags.sh --mode include-only

set -euo pipefail

PDB_DIR="${PDB_DIR:-/wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/}"
PROBE_DIR="${PROBE_DIR:-/wynton/group/degradolab/skt/docking/databases/probe_output/}"
VDG_LIB_DIR="${VDG_LIB_DIR:-/wynton/home/degradolab/skt/docking/frag_lib}"
LOG_DIR="${LOG_DIR:-/wynton/home/degradolab/skt/docking/frag_sge_logs}"
ESTIMATE="${ESTIMATE:-resources/frag_cost_estimate.tsv}"

MAX_H_RT="${MAX_H_RT:-36:00:00}"
# 250 CG occurrences, the generator default. The threshold cannot be made
# reliable by tuning: over five seeds of the 3000-structure sample the selected
# set holds ~720 fragments every time but ~40% of the union changes membership,
# so anything needed specifically belongs in INCLUDE rather than under a lower
# cutoff. 100 would build 1103 fragments with no less seed noise.
MIN_INSTANCES="${MIN_INSTANCES:-250}"
MEM_FREE="${MEM_FREE:-4G}"
SCRATCH="${SCRATCH:-20G}"

# Fragments to build regardless of their estimated count, space-separated. The
# threshold is a sampling estimate: across five seeds of the 3000-structure
# sample ~40% of the selected set changes membership while its size barely
# moves, so a fragment you actually need can miss the cut on the draw you took.
INCLUDE="${INCLUDE:-}"

usage() {
    cat >&2 <<'USAGE'
Usage: run_production_frags.sh --mode <threshold-plus-include|include-only> [--dry-run]

  --mode threshold-plus-include   Every fragment over $MIN_INSTANCES, plus $INCLUDE.
                                  Requires an EMPTY library: this rebuilds everything.
  --mode include-only             Only the fragments in $INCLUDE, threshold skipped.
                                  Requires an EXISTING library; adds to it in place.
  --dry-run                       Generate scripts and print the order; submit nothing.

$INCLUDE is honoured in BOTH modes; the mode decides only whether the count
threshold runs as well. There is no default mode -- a forgotten flag would
otherwise pick the full rebuild silently.
USAGE
    exit 2
}

DRY_RUN=0
MODE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --mode)    MODE="${2:-}"; shift 2 || usage ;;
        --mode=*)  MODE="${1#*=}"; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage ;;
        *) echo "ERROR: unknown argument '$1'." >&2; usage ;;
    esac
done

case "$MODE" in
    threshold-plus-include) TOP_UP=0 ;;
    include-only)           TOP_UP=1 ;;
    "") echo "ERROR: --mode is required." >&2; usage ;;
    *)  echo "ERROR: unknown --mode '$MODE'." >&2; usage ;;
esac

# Defaulted only now: include-only must not land in the full build's script dir,
# because the generate step wipes whatever directory it is pointed at.
if [ "$TOP_UP" -eq 1 ]; then
    # Timestamped per run, because the generate step rm -rf's whatever it is
    # pointed at: a fixed name would make a second top-up delete the first's
    # scripts. SGE spools the script at submit time so queued jobs survive that,
    # but the scripts are then gone for inspection and resubmission.
    SCRIPT_DIR="${SCRIPT_DIR:-ligand_vdgs/generate_vdgs/frag_submit_scripts_topup_$(date +%Y%m%d_%H%M%S)}"
    [ -n "$INCLUDE" ] || { echo "ERROR: --mode include-only needs INCLUDE='<smarts> [...]'." >&2; exit 1; }
else
    SCRIPT_DIR="${SCRIPT_DIR:-ligand_vdgs/generate_vdgs/frag_submit_scripts}"
fi

# --- preflight -------------------------------------------------------------
if [ "$TOP_UP" -eq 1 ]; then
    # The inverse of the fresh-build check: a top-up is meaningless without a
    # library to add to, and set_up_outdir refuses any fragment already built,
    # so an accidental repeat fails per-fragment rather than corrupting anything.
    [ -d "$VDG_LIB_DIR" ] || { echo "ERROR: --mode include-only needs an existing $VDG_LIB_DIR." >&2; exit 1; }
else
    # The library must hold no fragment directories. A run that writes new-key
    # fragment dirs alongside pre-ring-encoding ones produces a library no consumer
    # can tell apart. Only subdirectories count: the generator itself drops
    # fragment_aliases.tsv in this root, so a re-run must not trip over its own file.
    if [ -d "$VDG_LIB_DIR" ] && [ -n "$(find "$VDG_LIB_DIR" -mindepth 1 -maxdepth 1 -type d -print -quit)" ]; then
        echo "ERROR: $VDG_LIB_DIR already holds fragment directories, and" >&2
        echo "  --mode threshold-plus-include rebuilds the whole library." >&2
        echo "  To ADD fragments to what is already there, use:" >&2
        echo "    INCLUDE='<smarts> [...]' $0 --mode include-only" >&2
        echo "  To rebuild from scratch, move $VDG_LIB_DIR aside first." >&2
        exit 1
    fi
fi
[ -f "$ESTIMATE" ] || { echo "ERROR: missing $ESTIMATE (run estimate_frag_cost.py)." >&2; exit 1; }

mkdir -p "$LOG_DIR"

# --- generate --------------------------------------------------------------
# Regenerated from scratch every run: the generator refuses a non-empty output
# directory, and the scripts are derived output (gitignored, reproducible from
# the dict + estimate + flags). Reusing them would silently keep a stale
# MAX_H_RT, which is the one value most likely to change between runs.
rm -rf "$SCRIPT_DIR"

# --include takes the fragments verbatim; word-splitting INCLUDE is intended, so
# each SMARTS becomes its own argument.
INCLUDE_ARGS=()
if [ -n "$INCLUDE" ]; then
    INCLUDE_ARGS=(--include $INCLUDE)
    [ "$TOP_UP" -eq 1 ] && INCLUDE_ARGS+=(--include-only)
fi

# --frag-cost-estimate rather than the inline default so the submission order
# below is derived from exactly the counts the tiers were assigned from.
python ligand_vdgs/generate_vdgs/make_sge_scripts_for_frags.py \
    --frags-dict resources/database_frags_dict.pkl \
    --frag-cost-estimate "$ESTIMATE" \
    --min-instances "$MIN_INSTANCES" \
    ${INCLUDE_ARGS[@]+"${INCLUDE_ARGS[@]}"} \
    --pdb-dir "$PDB_DIR" \
    --probe-dir "$PROBE_DIR" \
    --vdg-lib-dir "$VDG_LIB_DIR" \
    --log-dir "$LOG_DIR" \
    --sge-out-dir "$SCRIPT_DIR" \
    --max-h-rt "$MAX_H_RT" \
    --mem-free "$MEM_FREE" \
    --scratch "$SCRATCH" \
    --subset-sizes 1 2

# --- order -----------------------------------------------------------------
# smiles_to_filename is the same encoder the generator used, so every line here
# names a script that exists; a missing one is a real mismatch, not a skip.
ORDER_FILE=$(mktemp)
trap 'rm -f "$ORDER_FILE"' EXIT
python - "$ESTIMATE" "$SCRIPT_DIR" > "$ORDER_FILE" <<'PY'
import os, sys
from ligand_vdgs.functions.utils import smiles_to_filename
from ligand_vdgs.generate_vdgs.estimate_frag_cost import read_estimate_tsv

estimate_path, script_dir = sys.argv[1:]
est, _ = read_estimate_tsv(estimate_path)

built = {f[:-3] for f in os.listdir(script_dir) if f.endswith('.sh')}
by_cg = {}
for smiles, count in est.items():
    cg = smiles_to_filename(smiles)
    if cg in built:
        by_cg[cg] = max(count, by_cg.get(cg, 0))

missing = built - set(by_cg)
if missing:
    sys.exit(f'{len(missing)} generated script(s) absent from the estimate, '
             f'e.g. {sorted(missing)[:3]}')

# Descending estimated cost: longest-running job first, so it goes in while the
# maintenance window is still wide enough to schedule it.
ordered = sorted(by_cg, key=lambda cg: -by_cg[cg])
for cg in ordered:
    print(f'{by_cg[cg]}\t{os.path.join(script_dir, cg + ".sh")}')
PY

# In include-only mode, refuse fragments the library already holds. Nothing
# downstream checks this before submission: the job would qsub cleanly, sit in
# the queue, start, and only then die in set_up_outdir with "output directory is
# not empty" -- so without this the run reports "1 job submitted" and the
# fragment silently never appears. An EMPTY directory is fine (set_up_outdir
# accepts it), which is what a fragment killed before it wrote anything leaves.
if [ "$TOP_UP" -eq 1 ]; then
    KEPT_FILE=$(mktemp)
    already=0
    partial=0
    while IFS=$'\t' read -r est_count script; do
        cg=$(basename "$script" .sh)
        d="$VDG_LIB_DIR/$cg"
        if [ ! -d "$d" ] || [ -z "$(ls -A "$d" 2>/dev/null)" ]; then
            # Absent, or present but empty -- set_up_outdir accepts an empty dir,
            # which is what a job killed before it wrote anything leaves.
            printf '%s\t%s\n' "$est_count" "$script" >> "$KEPT_FILE"
        elif grep -q '^Job completed\.' "$d/${cg}_log" 2>/dev/null; then
            echo "ALREADY BUILT, skipping: $cg" >&2
            already=$((already + 1))
        else
            # Non-empty but no completion marker: a crashed or clamped job. Never
            # silently skipped -- that would strand the fragment permanently,
            # since every later top-up would see the same non-empty directory and
            # skip it again. set_up_outdir will not overwrite it either, so the
            # only way forward is an explicit removal by the operator.
            echo "PARTIAL (no 'Job completed.'), NOT submitted: $cg" >&2
            echo "    rm -rf $d   # then re-run this top-up" >&2
            partial=$((partial + 1))
        fi
    done < "$ORDER_FILE"
    mv "$KEPT_FILE" "$ORDER_FILE"
    if [ ! -s "$ORDER_FILE" ]; then
        echo "Nothing to submit: $already already built, $partial partial." >&2
        [ "$partial" -gt 0 ] && exit 1
        exit 0
    fi
    [ "$already" -gt 0 ] && echo "Skipped $already already-built fragment(s)."
    [ "$partial" -gt 0 ] && echo "WARNING: $partial partial fragment dir(s) left alone; see above."
fi

TOTAL=$(wc -l < "$ORDER_FILE")
echo
echo "Submission order: $TOTAL job(s), longest estimated runtime first."
echo "Maintenance boundary: 2026-09-07 15:00 PDT. h_rt ceiling: $MAX_H_RT."
echo

if [ "$DRY_RUN" -eq 1 ]; then
    echo "--- first 15 ---"; head -15 "$ORDER_FILE"
    echo "--- last 10 ---";  tail -10 "$ORDER_FILE"
    echo
    echo "Dry run: nothing submitted. Scripts are in $SCRIPT_DIR."
    exit 0
fi

# qsub failures are collected rather than fatal. Under `set -e` a single
# transient rejection -- the per-user pending-job cap is easy to hit at this
# scale -- would abort the loop and leave every remaining fragment unsubmitted
# with nothing recorded about it.
FAILED_FILE="qsub_failures.txt"
: > "$FAILED_FILE"
n=0
ok=0
while IFS=$'\t' read -r est_count script; do
    n=$((n + 1))
    printf '[%4d/%4d] %7s structures  ' "$n" "$TOTAL" "$est_count"
    if qsub "$script"; then
        ok=$((ok + 1))
    else
        echo "QSUB FAILED: $script" >&2
        echo "$script" >> "$FAILED_FILE"
    fi
done < "$ORDER_FILE"

echo
echo "Submitted $ok of $n job(s). Watch with: qstat -u \$USER"
if [ -s "$FAILED_FILE" ]; then
    echo "WARNING: $(wc -l < "$FAILED_FILE") job(s) were rejected; see $FAILED_FILE."
    echo "  Resubmit with: while read -r s; do qsub \"\$s\"; done < $FAILED_FILE"
else
    rm -f "$FAILED_FILE"
fi
# A production run writes no marker file -- _finish() appends 'Job completed.' to
# <cg>_log. That line is the only way to tell a finished fragment from one the
# maintenance window killed, since both leave a populated nr_vdgs/.
# Driven off the generated scripts, not off $VDG_LIB_DIR/*/: a fragment that was
# never scheduled, was rejected at qsub, or died before set_up_outdir ran has no
# library directory at all, so globbing the library cannot see the failures that
# matter most.
echo "Afterwards, find fragments that did not finish with:"
echo "  for s in $SCRIPT_DIR/*.sh; do cg=\$(basename \"\$s\" .sh); \\"
echo "    grep -q '^Job completed\.' \"$VDG_LIB_DIR/\$cg/\${cg}_log\" 2>/dev/null \\"
echo "      || echo \"INCOMPLETE \$cg\"; done"
