project(ilp_part)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH               "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -DIL_STD -m64 -Wall -ldl -fPIC")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX   "-isystem ")
set(LINKER_OPTIONS                  "-Wl,--export-dynamic")

find_package(Cplex)
find_package(Threads REQUIRED)
include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS})
add_executable(ilp_part Main.cpp IlpReadHypergraph.cpp)
target_link_libraries(ilp_part
    PUBLIC
    ${CPLEX_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT})