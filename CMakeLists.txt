cmake_minimum_required(VERSION 3.6)
project(ScientificComputingProj)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp Prob1.cpp Prob1.h AbstractProb.cpp AbstractProb.h Prob2.cpp Prob2.h Prob3.cpp Prob3.h Prob4.cpp Prob4.h Prob5.cpp Prob5.h Prob6.cpp Prob6.h)
add_executable(ScientificComputingProj ${SOURCE_FILES})