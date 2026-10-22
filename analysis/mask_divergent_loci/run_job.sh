#!/bin/bash
# Run one mask-divergent-loci job for the sensitivity analysis (see MASK_DIVERGENT_LOCI.md).
#   run_job.sh <dataset> <k> <window size> <outlier threshold>
# <dataset> is "orig" (the original MSAs in DATA) or the name of a directory $WORK/<dataset>_msas
# (e.g. planted, planted_nodonor, control_nodonor).
# Environment:
#   WORK     working directory containing planted_msas/ (required)
#   DATA     original MSA directory (default: example fungal marker set)
#   THREADS  threads per run (default: 24)
#   KEEP=1   keep the comparisons tables; results go to $WORK/details/ instead of $WORK/runs/
#   UPPER=1  also apply the high-identity test, with the same k, window size and threshold (by default
#            it is disabled with --no-upper-mask so the sweeps evaluate the low-identity test alone)
#   UPPER_SWEEP=1  sweep the high-identity test: k, window size and threshold apply to it, the low-identity
#            test uses its defaults (3, 120, 25), and results go to <dataset>_uk<k>_un<N>[_ut<T>]
set -euo pipefail
: "${WORK:?Set WORK to the analysis working directory}"
DATA=${DATA:-/srv/projects/donovan/gtrim/data/gtdb_r232_fun77_markers}
THREADS=${THREADS:-24}
BIN=$(cd "$(dirname "$0")/../.." && pwd)/target/release/gtrim

ds=$1; k=$2; n=$3; t=$4
if [ "$ds" = orig ]; then IN=$DATA; else IN=$WORK/${ds}_msas; fi
UPPER_ARG=--no-upper-mask
[ "${UPPER:-0}" = 1 ] && UPPER_ARG="--upper-neighbours $k --upper-window-size $n --upper-outlier-threshold $t"
LOWER_ARGS="--neighbours $k --window-size $n --outlier-threshold $t"
NAME=${ds}_k${k}_n${n}
SUFFIX=_t${t}
if [ "${UPPER_SWEEP:-0}" = 1 ]; then
  UPPER_ARG="--upper-neighbours $k --upper-window-size $n --upper-outlier-threshold $t"
  LOWER_ARGS=""
  NAME=${ds}_uk${k}_un${n}
  SUFFIX=_ut${t}
fi

if [ "${KEEP:-0}" = 1 ]; then
  d=$WORK/details/$NAME
else
  d=$WORK/runs/$NAME
fi
[ "$t" != 25 ] && d=${d}${SUFFIX}
[ -f "$d/divergent_loci_masked.tsv" ] && exit 0

if ! "$BIN" mask-divergent-loci -i "$IN" -o "$d" -t "$THREADS" $LOWER_ARGS $UPPER_ARG > /dev/null 2>&1; then
  echo "FAILED $d"
  exit 1
fi

rm -rf "$d/masked_msas"
[ "${KEEP:-0}" = 1 ] || rm -f "$d/divergent_loci_comparisons.tsv" "$d/divergent_loci_upper_comparisons.tsv"
echo "$(date +%H:%M:%S) done $ds k=$k n=$n t=$t"
