cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)
project(KISS_FFT)


include_directories(./)
include_directories(./tools/)

link_directories(./)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

set(targetLib kissfft)
file (GLOB lib_srcs tools/*.c ./*.c)
add_library(${targetLib} SHARED ${lib_srcs})
set(target_name "kiss_fft_main")
add_executable(${target_name} test/test_real.c)
target_link_libraries(${target_name} ${targetLib})