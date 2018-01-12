g++ `pkg-config --cflags cbc` -fsanitize=address \
    containers.c cp_cuts.c cprop.c memory.c strutils.c \
    preprocess.c cprop_lp.c lp.cpp -DCBC -O0 -g3 \
    `pkg-config --libs cbc` \
    -o preprocess 
