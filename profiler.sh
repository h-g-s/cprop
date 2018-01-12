clang++ preprocess.c cprop.c lp.cpp containers.c memory.c strutils.c cp_cuts.c \
     -O1 -g3 -std=gnu++11 -Wno-deprecated \
     -DCBC -I/usr/include/coin/ \
     -L/usr//lib -lCbcSolver -lCbc -lpthread -lrt -lreadline -lncurses -lCgl -lOsiClp -lClpSolver -lClp -lreadline -lncurses \
     -lOsi -lCoinUtils -lreadline -lncurses -lbz2 -lz -llapack -lblas -lm -lprofiler -fopenmp -fopenmp=libiomp5 \
     -o preprocess
CPUPROFILE=cpuprofile.log ./preprocess ~/inst/cbcbench/23844-T1_0.mps.gz
pprof --gv preprocess cpuprofile.log
