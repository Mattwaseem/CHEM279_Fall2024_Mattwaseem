cmake_minimum_required(VERSION 3.10)

# Project Name
project(LJ_Project)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(include)

# Set the output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Add the executable with source files
add_executable(ljpotenial src/main.cpp src/file_io.cpp src/LJ_PE.cpp)

# Ensure the input directory exists
add_custom_command(
    OUTPUT sample_input/Energy
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/sample_input/Energy
    COMMENT "Ensuring sample_input/Energy directory exists"
)

# Run ljpotenial to generate the truncation_error_Problem2.txt
# Set the WORKING_DIRECTORY to the root project folder to access the correct input files.
add_custom_command(
    OUTPUT truncation_error_Problem2.txt
    COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ljpotenial
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}  # Ensure the correct working directory
    DEPENDS sample_input/Energy
    COMMENT "Running ljpotenial to generate truncation_error_Problem2.txt"
)

# Add a custom target to run ljpotenial
add_custom_target(
    run_ljpotenial ALL
    DEPENDS truncation_error_Problem2.txt
)

# Custom target to clean the bin directory
add_custom_target(
    clean-bin
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*
    COMMENT "Cleaning up bin directory"
)

# Custom target to clean all build files (including bin)
add_custom_target(
    clean-all
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/* ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*
    COMMENT "Cleaning all build and bin directories"
)
