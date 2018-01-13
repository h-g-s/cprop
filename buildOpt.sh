g++ `pkg-config --cflags cbc` -Ofast -fexpensive-optimizations -flto -DNDEBUG \
    containers.c cp_cuts.c cprop.c memory.c strutils.c \
    preprocess.c cprop_lp.c lp.cpp -DCBC \
    `pkg-config --libs cbc` \
    -o preprocess 
