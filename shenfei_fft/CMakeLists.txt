cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)
project(KISS_FFT)

include_directories(./src/)

link_directories(/usr/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

find_library(mLib m) # libm.so

set(targetLib kissfft)
file (GLOB lib_srcs src/*.c)
add_library(${targetLib} SHARED ${lib_srcs})
#set(target_name "kiss_fft_main")
add_executable(kissfft2d_example ./kissfft2d_example.cpp)
target_link_libraries(kissfft2d_example ${targetLib} ${mLib})