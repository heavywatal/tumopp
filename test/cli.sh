set -eux
TMP_OUT=tumopp_cli.sh

./tumopp -C hex -L step -P mindrag -N 255 -o $TMP_OUT
rm -r $TMP_OUT

./tumopp -Chex -Lstep -Pmindrag -N255 -o$TMP_OUT
rm -r $TMP_OUT
