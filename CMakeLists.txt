cmake_minimum_required(VERSION 3.15)

project(
    IncBeta
    VERSION 0.0.1
    LANGUAGES CXX
)

#================================================================
# Default build type
#================================================================
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message("-- Using default value of CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
else()
    message("-- Obtained CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

#================================================================
# Add executable
#================================================================
add_executable(
    IncBeta
    incbeta.cpp
)
target_compile_features(IncBeta PUBLIC cxx_std_11)