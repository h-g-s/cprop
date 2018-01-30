g++ `pkg-config --cflags cbc` \
    containers.c cp_cuts.c cprop.c memory.c strutils.c \
    preprocess.c cprop_lp.c lp.cpp -DCBC -g -O2 \
    mincut.c \
    `pkg-config --libs cbc` \
    -o preprocess 
