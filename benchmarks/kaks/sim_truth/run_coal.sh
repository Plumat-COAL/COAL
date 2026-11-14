#!/bin/bash
# Accuracy of COAL's Ka/Ks against simulated truth (docs/kaks.md, "Accuracy against simulated truth"). From the
# repository root, with numpy, scipy, pandas (and matplotlib for the figures):
#   export COAL_SIM_DIR=work/sim_truth COAL_CPP_BIN=build/coal-cpp
#   python benchmarks/kaks/sim_truth/validate_sim.py     # the simulator against its model (optional)
#   python benchmarks/kaks/sim_truth/gen.py              # the grid: 120 AXT files, 240,000 pairs, ~500 MB
#   bash benchmarks/kaks/sim_truth/run_coal.sh           # this script: every method and option run (hours, 1 core)
#   python benchmarks/kaks/sim_truth/analyze.py && python benchmarks/kaks/sim_truth/perpair.py
#   python benchmarks/kaks/sim_truth/make_tables.py      # -> $COAL_SIM_DIR/res/tables.md
#   python benchmarks/kaks/sim_truth/make_figs.py        # -> $COAL_SIM_DIR/plots/
# Resumable COAL estimation over the grid. One gzipped TSV per (input file, run); skipped if present.
set -u
W=${COAL_SIM_DIR:-work/sim_truth}
C=$(realpath "${COAL_CPP_BIN:-build/coal-cpp}")
cd "$W"
mkdir -p res/coal
run() {  # name, input, args...
  local name=$1 in=$2; shift 2
  local out=res/coal/$(basename $in .axt).$name.tsv.gz
  [ -s $out ] && return 0
  nice -n 10 $C "$@" -j 1 -q -o res/coal/tmp.$name.tsv -i $in 2>>res/coal/stderr.log
  gzip -c res/coal/tmp.$name.tsv > $out.tmp && mv $out.tmp $out && rm -f res/coal/tmp.$name.tsv
}
for f in data/f*.axt; do
  run main $f dnds -m YN -m MYN -m NG -m LWL -m MLWL -m LPB -m MLPB -m Comeron -m Li
  run cgest $f dnds -m Comeron --comeron-flavour gestimator
  run cms0 $f dnds -m Comeron --comeron-min-sites 0
  run cms10 $f dnds -m Comeron --comeron-min-sites 10
  run liseq $f dnds -m Li --li-flavour seqinr
  run kaks $f kaks -m YN -m MYN
  case $f in *_k2_*) run ynpaml $f dnds -m YN --yn-flavour paml; run ynw0 $f dnds -m YN --yn-weighting 0;; esac
  echo "$(date +%H:%M) $f" >> res/coal/done.log
done
echo ALLDONE >> res/coal/done.log
