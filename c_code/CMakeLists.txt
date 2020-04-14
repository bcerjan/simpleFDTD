###############################################################################
#    Copyright (c) 2020 Ben Cerjan
#
#    This file is part of simpleFDTD.
#
#    simpleFDTD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    simpleFDTD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with simpleFDTD.  If not, see <https://www.gnu.org/licenses/>.
###############################################################################

# cmake version:
cmake_minimum_required( VERSION 3.0 )

#project name:
project( fdtdWidget )

#flags
set(CMAKE_CXX_FLAGS)

# Build-type Optimizations
SET (DEBUG_OPTIONS "-Wall -Wextra -g")
SET (RELEASE_OPTIONS "-O2")

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

SET (MULTIPLE_TARGET
         ./src/multiple_empty.c
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
add_executable( multi ${COMMON_SRCS} ${MULTIPLE_TARGET} )

#link libraries
target_link_libraries(empty m)
target_link_libraries(structure m)
target_link_libraries(test_mem m)
target_link_libraries(multi m)


# Finally apply our additional changes:
target_compile_options(multi PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
target_compile_options(multi PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")
