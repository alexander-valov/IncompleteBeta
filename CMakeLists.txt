cmake_minimum_required(VERSION 3.10)

project(
    IncompleteBeta
    VERSION 0.0.1
    LANGUAGES CXX
)

#================================================================
# Check if IncompleteBeta is built as subproject or it's main project
#================================================================
set(INC_BETA_MAIN_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(INC_BETA_MAIN_PROJECT ON)
endif()

#================================================================
# Define options
#================================================================
if (POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()

option(INC_BETA_BUILD_TESTS "Build the unit tests" ${INC_BETA_MAIN_PROJECT})

#================================================================
# Define configuration
#================================================================
set(INC_BETA_TARGET_NAME ${PROJECT_NAME})
set(INC_BETA_INCLUDE_BUILD_DIR "${PROJECT_SOURCE_DIR}/include")

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
# Add INTERFACE library
#================================================================
add_library(${INC_BETA_TARGET_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${INC_BETA_TARGET_NAME} ALIAS ${INC_BETA_TARGET_NAME})
target_compile_features(${INC_BETA_TARGET_NAME} INTERFACE cxx_std_11)
target_include_directories(
    ${INC_BETA_TARGET_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${INC_BETA_INCLUDE_BUILD_DIR}>
    $<INSTALL_INTERFACE:include>
)

#================================================================
# Buld unit tests
#================================================================
# if (INC_BETA_BUILD_TESTS)
#     enable_testing()
#     add_subdirectory(test)
# endif()