for file in ~/inst/cbcbench/*.mps.gz;
do
    ./preprocess $file -probeAndCut
done
for file in ~/inst/bpconf/lp/*.lp;
do
    ./preprocess $file -probeAndCut
done

