add_executable(cprop main.c cprop.c containers.c memory.c strutils.c lp.cpp)

SET( CPLEXDIR "/opt/ibm/ILOG/CPLEX_Studio1271/cplex/")

SET( CMAKE_CXX_FLAGS "-Wall -std=c++11 -DCPX -I${CPLEXDIR}/include/ilcplex/ -fopenmp -DDEBUG -O0 -g3" )
    
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -L${CPLEXDIR}lib/x86-64_linux/static_pic/" )

target_link_libraries(cprop cplex m pthread)
