#!/usr/bin/env bash
# Reproduce the G1 attribution (revision plan, 17 November afternoon entry) and the
# stall check. Usage, from the repository root:
#   bash "2025 sdf paper/jfinecon_referees/archive/g1_attribution/reproduce.sh" /path/to/empty/outdir
# Needs the local Preqin pools (data_prepared_2026_b and _c, VYP files) and the G1
# results in r_project/results/data_out_2026-emp-c. About 20 minutes (the old code
# is serial). Creates two temporary git worktrees inside OUT and removes them.
set -euo pipefail
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
REPO=$(git rev-parse --show-toplevel)
HERE="$REPO/2025 sdf paper/jfinecon_referees/archive/g1_attribution"
OUT=${1:?output directory}; mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
LOG="$OUT/cells.log"
# data_prepared_2026 holds q_factors.csv (not tracked at d2bdf6e)
link_data() { for d in data_prepared_2026 data_prepared_2026_b data_prepared_2026_c; do [ -e "$1/r_project/empirical/$d" ] || ln -s "$REPO/r_project/empirical/$d" "$1/r_project/empirical/$d"; done; }
cleanup() { git -C "$REPO" worktree remove --force "$OUT/oldcode" 2>/dev/null || true; git -C "$REPO" worktree remove --force "$OUT/ownend" 2>/dev/null || true; }
trap cleanup EXIT

# 1) code of the paper's tables
git -C "$REPO" worktree add --detach -q "$OUT/oldcode" d2bdf6e; link_data "$OUT/oldcode"
# 2) current code with the old window endpoint (one line in helper/objective.R)
git -C "$REPO" worktree add --detach -q "$OUT/ownend" HEAD; link_data "$OUT/ownend"
python3 - "$OUT/ownend/r_project/helper/objective.R" <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
a = "pad(d, t0[[paste(d$type[1], d$Vintage[1])]], t.max))"
assert s.count(a) == 1, "patch target not found"
open(p, "w").write(s.replace(a, "pad(d, t0[[paste(d$type[1], d$Vintage[1])]], max(d$Date)))"))
PY

run() { (cd "$1/r_project" && Rscript "$HERE/run_cell.R" "$2" "$OUT/$3" "$LOG" > "$OUT/$3.out" 2>&1); tail -2 "$LOG"; }
run "$OUT/oldcode" data_prepared_2026_b old_b
run "$OUT/oldcode" data_prepared_2026_c old_c
run "$REPO"        data_prepared_2026_b new_b
run "$OUT/ownend"  data_prepared_2026_b ownend_b
run "$OUT/ownend"  data_prepared_2026_c ownend_c
rm -f "$OUT/oldcode/r_project/Rplots.pdf" "$OUT/ownend/r_project/Rplots.pdf" "$REPO/r_project/Rplots.pdf"

cd "$REPO/r_project"
Rscript "$HERE/compare.R" "$OUT"
if cmp -s "$OUT/g1_attribution_north_america.csv" "$HERE/../g1_attribution_north_america.csv"; then
  echo "attribution table identical to archive/g1_attribution_north_america.csv"
else
  echo "attribution table DIFFERS from archive/g1_attribution_north_america.csv"
fi
Rscript "$HERE/check_stall.R" "$OUT"
