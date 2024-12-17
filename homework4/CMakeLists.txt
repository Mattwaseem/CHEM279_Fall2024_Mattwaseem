cmake_minimum_required(VERSION 3.10)
project(CNDO2)

set(CMAKE_CXX_STANDARD 17)

# Find Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Set output directories to the `bin` directory one level up from the `build` directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Add include directory
include_directories(include)

# Add source files
add_executable(CNDO2
    src/main.cpp
    src/DensityMatrix.cpp
    src/FockMatrix.cpp
    src/GammaCalculator.cpp
    src/Molecule.cpp
    src/FileInputParser.cpp
    src/Constants.cpp
    src/OverlapMatrix.cpp
    src/CartesianGaussian.cpp
)

# Link Armadillo
target_link_libraries(CNDO2 ${ARMADILLO_LIBRARIES})

# Run command to execute binary from bin directory
add_custom_target(run
    COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/CNDO2 ${CMAKE_SOURCE_DIR}/sample_input/H2.txt ${CMAKE_SOURCE_DIR}/sample_input/HF.txt ${CMAKE_SOURCE_DIR}/sample_input/HO.txt ${CMAKE_SOURCE_DIR}/sample_input/N2.txt ${CMAKE_SOURCE_DIR}/sample_input/O2.txt
    DEPENDS CNDO2
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running CNDO2 program"
)

# Custom target to clean build and bin directories
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/Makefile
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/cmake_install.cmake
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/* # Clean contents within build
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_SOURCE_DIR}/bin/* # Clean contents within bin
    COMMENT "Cleaning all build and binary files"
)
