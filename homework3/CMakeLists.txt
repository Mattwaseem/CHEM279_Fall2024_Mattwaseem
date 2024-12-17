cmake_minimum_required(VERSION 3.10)
project(MoleculeSetup)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_OSX_DEPLOYMENT_TARGET "14.0")

include_directories(include)

set(ARMADILLO_INCLUDE_DIR /opt/homebrew/opt/armadillo/include)
set(ARMADILLO_LIBRARY /opt/homebrew/opt/armadillo/lib/libarmadillo.dylib)

include_directories(${ARMADILLO_INCLUDE_DIR})

add_executable(MoleculeSetup
    src/main.cpp
    src/molecule.cpp
    src/input_parser.cpp
    src/CartesianGaussian.cpp
    src/OverlapMatrix.cpp
    src/Hamiltonian.cpp
    src/basis_function.cpp
)

target_link_libraries(MoleculeSetup ${ARMADILLO_LIBRARY})

target_compile_options(MoleculeSetup PRIVATE -g -Wall -Wextra -Wpedantic)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

add_custom_target(run
    COMMAND ${CMAKE_SOURCE_DIR}/bin/MoleculeSetup ${CMAKE_SOURCE_DIR}/sample_input/C2H2.txt > ${CMAKE_SOURCE_DIR}/calculated_output/C2H2_output.txt
    COMMAND ${CMAKE_SOURCE_DIR}/bin/MoleculeSetup ${CMAKE_SOURCE_DIR}/sample_input/C2H4.txt > ${CMAKE_SOURCE_DIR}/calculated_output/C2H4_output.txt
    COMMAND ${CMAKE_SOURCE_DIR}/bin/MoleculeSetup ${CMAKE_SOURCE_DIR}/sample_input/H2.txt > ${CMAKE_SOURCE_DIR}/calculated_output/H2_output.txt
    DEPENDS MoleculeSetup
    COMMENT "Running MoleculeSetup with C2H2.txt, C2H4.txt, and H2.txt and generating output for each"
    VERBATIM
)

# Custom clean-all target to remove build files, executables, and output
add_custom_target(clean-all
    COMMAND rm -rf ${CMAKE_BINARY_DIR}/*
    COMMAND rm -rf ${CMAKE_SOURCE_DIR}/bin/*
    COMMAND rm -rf ${CMAKE_SOURCE_DIR}/calculated_output/*
    COMMENT "Removing all build files, executables, and output files"
)
