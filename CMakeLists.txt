cmake_minimum_required(VERSION 3.16)
project(hdr_generator)

set(CMAKE_CXX_STANDARD 11)

#set(OPTIMISATIONS_OPTIONS -Ofast -march=native -ftree-vectorize -ffast-math -fprofile-generate -ftree-parallelize-loops=2 -fno-inline -v)  #pour la generation du profiling
set(OPTIMISATIONS_OPTIONS -Ofast -march=native -ftree-vectorize -ffast-math  -fno-inline -ftree-parallelize-loops=2 -fprofile-use -v)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(hdr_generator main.cpp)
target_link_libraries(hdr_generator Threads::Threads)
target_link_libraries(hdr_generator gcov)  #pour la generation du profiling
target_link_libraries(hdr_generator OpenMP::OpenMP_CXX)

target_compile_options(hdr_generator PUBLIC "$<$<CONFIG:DEBUG>:${OPTIMISATIONS_OPTIONS}>")
