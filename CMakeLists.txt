cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)
project(KISS_FFT)


include_directories(./)
include_directories(./tools/)

link_directories(/usr/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

find_library(mLib m) # libm.so

set(targetLib kissfft)
file (GLOB lib_srcs tools/*.c ./*.c)
add_library(${targetLib} SHARED ${lib_srcs})
#set(target_name "kiss_fft_main")
add_executable(test_real test/test_real.c)
target_link_libraries(test_real ${targetLib} ${mLib})