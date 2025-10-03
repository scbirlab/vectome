#!/usr/bin/env python

set -euox pipefail

OUTDIR=test/outputs
CACHE=$OUTDIR/.cache
TEST_LIST=$OUTDIR/test-inputs.txt
GROUP=1
PROJ=8

mkdir -p $OUTDIR

vectome build $GROUP --force --cache "$CACHE"
vectome info
vectome info --cache "$CACHE"

queries=("Escherichia coli" 83332 83333 "Klebsiella pneumoniae")
for q in "${queries[@]}"
do
    echo "$q" >> $TEST_LIST
done

vectome embed "$TEST_LIST" --cache "$CACHE" > $OUTDIR/test1.tsv
vectome embed "$TEST_LIST" --projection $PROJ --seed 0 --cache "$CACHE" > $OUTDIR/test2.tsv
vectome embed "$TEST_LIST" --projection $PROJ --seed 42 --cache "$CACHE" > $OUTDIR/test3.tsv

if [ "$(diff $OUTDIR/test2.tsv $OUTDIR/test3.tsv | wc -l)" -eq "0" ]
then
    >&2 echo "Seed did not change output!"
    exit 1
fi

vectome embed "$TEST_LIST" --method landmark --group $GROUP --cache "$CACHE" > $OUTDIR/test4.tsv
vectome embed "$TEST_LIST" --method landmark --group $GROUP --projection $PROJ --cache "$CACHE" > $OUTDIR/test5.tsv

for f in $OUTDIR/test?.tsv
do 
    >&2 echo "::: $f"
    >&2 head "$f"
done

>&2 echo "[$(date)] Done!"
