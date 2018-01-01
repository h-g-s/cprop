add_executable(preprocess preprocess.c cprop.c containers.c memory.c strutils.c lp.cpp cp_cuts.c)
add_executable(cprop main.c cprop.c containers.c memory.c strutils.c lp.cpp cp_cuts.c)
#add_executable(relax relax.c cprop.c containers.c memory.c strutils.c lp.cpp cp_cuts.c)
#add_executable(genbip genbip.c lp.cpp)

SET(CMAKE_C_COMPILER clang)
SET(CMAKE_CXX_COMPILER clang++)

#SET( CPLEXDIR "/opt/ibm/ILOG/CPLEX_Studio1271/cplex/")

SET( CMAKE_CXX_FLAGS "-Wall -std=c++11 -DCBC  -fopenmp -I/usr/include/coin/ -I/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.1/include/ -g" )
    
SET( CMAKE_C_FLAGS "-Wall -std=c99  -fopenmp -I/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.1/include/ -g")
    
SET( CMAKE_EXE_LINKER_FLAGS "-fopenmp=libiomp5" )

target_link_libraries(preprocess -L/usr//lib -lCbcSolver -lCbc -lpthread -lrt -lreadline -lncurses -lCgl -lOsiClp -lClpSolver -lClp -lreadline -lncurses -lOsi -lCoinUtils -lreadline -lncurses -lbz2 -lz -llapack -lblas -lm -lprofiler )

target_link_libraries(cprop -L/usr//lib -lCbcSolver -lCbc -lpthread -lrt -lreadline -lncurses -lCgl -lOsiClp -lClpSolver -lClp -lreadline -lncurses -lOsi -lCoinUtils -lreadline -lncurses -lbz2 -lz -llapack -lblas -lm)
#target_link_libraries(relax cplex m pthread)
#target_link_libraries(genbip -L/usr//lib -lCbcSolver -lCbc -lpthread -lrt -lreadline -lncurses -lCgl -lOsiClp -lClpSolver -lClp -lreadline -lncurses -lOsi -lCoinUtils -lreadline -lncurses -lbz2 -lz -llapack -lblas -lm)
