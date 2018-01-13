for file in ~/inst/cbcbench/*.mps.gz;
do
    ./preprocess $file
done
for file in ~/inst/bpconf/lp/*.lp;
do
    ./preprocess $file
done

