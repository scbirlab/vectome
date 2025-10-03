#!/usr/bin/env python

set -euox pipefail

OUTDIR=test/outputs
CACHE=$OUTDIR/.cache
TEST_LIST=$OUTDIR/test-inputs.txt

mkdir -p $OUTDIR

vectome build 1 --cache "$CACHE"

queries=("Escherichia coli" 83332 83333 "Klebsiella pneumoniae")
for q in "${queries[@]}"
do
    echo "$q" >> $TEST_LIST
done

vectome embed "$TEST_LIST" --cache "$CACHE" > $OUTDIR/test1.tsv
vectome embed "$TEST_LIST" --projection 3 --cache "$CACHE" > $OUTDIR/test2.tsv
vectome embed "$TEST_LIST" --projection 3 --cache "$CACHE" > $OUTDIR/test3.tsv
vectome embed "$TEST_LIST" --projection 3 --seed 0 --cache "$CACHE" > $OUTDIR/test3.tsv
vectome embed "$TEST_LIST" --cache "$CACHE" > $OUTDIR/test3.tsv
vectome embed "$TEST_LIST" --method landmarks --projection 3 --cache "$CACHE" > $OUTDIR/test3.tsv

>&2 echo "[$(date)] Done!"
