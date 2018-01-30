g++ `pkg-config --cflags cbc`  \
    containers.c cp_cuts.c cprop.c memory.c strutils.c \
    genbip.c cprop_lp.c lp.cpp -DCBC -Ofast \
    `pkg-config --libs cbc` \
    -o genbip
