project(ilp_part)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH               "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -DIL_STD -m64 -Wall -ldl -fPIC")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX   "-isystem ")
set(LINKER_OPTIONS                  "-Wl,--export-dynamic")

find_package(Cplex)
find_package(Threads REQUIRED)
find_package(ortools REQUIRED)

add_executable(ilp_part Main.cpp IlpReadHypergraph.cpp)
target_link_libraries(ilp_part
    PUBLIC
    ${CMAKE_THREAD_LIBS_INIT}
    ortools::ortools)

# CPLEX is optional. When found (e.g. commercial install on the build host) the
# ILP uses the CPLEX backend, matching the original behavior. When absent (e.g.
# the shippable OR-Tools-only container) the open-source OR-Tools backend is
# compiled in instead -- no source edits needed.
if(CPLEX_FOUND)
    message(STATUS "ilp_part: CPLEX found -> using CPLEX backend")
    target_compile_definitions(ilp_part PRIVATE USE_CPLEX)
    target_include_directories(ilp_part SYSTEM PRIVATE ${CPLEX_INCLUDE_DIRS})
    target_link_libraries(ilp_part PUBLIC ${CPLEX_LIBRARIES})
else()
    message(STATUS "ilp_part: CPLEX not found -> using OR-Tools backend")
endif()