# cmake version:
cmake_minimum_required( VERSION 3.0 )

#project name:
project( fdtdWidget )

#flags
set(CMAKE_CXX_FLAGS)

#files

#include
include_directories( ./include/fdtd ./src )

SET (COMMON_SRCS
          ./src/field_updates.c
          ./src/grid_initialization.c
          ./src/memory_allocation.c
          ./src/array_functions.c
          ./src/ricker.c
          ./src/grid_free.c
          ./src/structure_functions.c
     )


SET (EMPTY_TARGET
          ./src/empty_fdtd.c
     )

SET (STRUCTURE_TARGET
          ./src/structure_fdtd.c
     )

SET (TEST_TARGET
          ./src/test_memory.c
     )


#target
add_executable( empty ${COMMON_SRCS} ${EMPTY_TARGET} )
add_executable( structure ${COMMON_SRCS} ${STRUCTURE_TARGET} )
add_executable( test_mem ${COMMON_SRCS} ${TEST_TARGET} )

#link libraries
target_link_libraries(empty m)
target_link_libraries(structure m)
target_link_libraries(test_mem m)
