cmake_minimum_required(VERSION 2.8)
project(libKrylov)

set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_C_COMPILER clang)

set(SHARED_FLAGS "-DDEBUG -Wall -Wextra -O0 -g -fcolor-diagnostics")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SHARED_FLAGS} -stdlib=libc++ -std=c++11 -fpic")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SHARED_FLAGS} -std=c11")

include_directories(
  include
  gtest/include
)

add_subdirectory(gtest)

add_library(kry_core SHARED 
  src/Math.cxx src/Runtime.cxx src/Utility.cxx
  src/Arnoldi.cxx
)

add_library(kry_input SHARED
  src/PSSE.cxx
  src/Common.cc
)

add_executable(BasicMathTests test/BasicMathTests.cxx)
target_link_libraries(BasicMathTests kry_core gtest gtest_main pthread)

add_executable(ArnoldiTests test/ArnoldiTests.cxx)
target_link_libraries(ArnoldiTests kry_core gtest gtest_main pthread)

add_executable(ieee14Test test/ieee14.cxx)
target_link_libraries(ieee14Test kry_core kry_input gtest gtest_main pthread)
