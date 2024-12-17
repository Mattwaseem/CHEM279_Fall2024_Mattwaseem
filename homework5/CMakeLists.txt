# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)

project(CNDO_2_Gradient)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${CMAKE_SOURCE_DIR}/include)

set(SOURCES
    ${CMAKE_SOURCE_DIR}/src/main.cpp
    ${CMAKE_SOURCE_DIR}/src/Molecule.cpp
    ${CMAKE_SOURCE_DIR}/src/CNDO.cpp
    ${CMAKE_SOURCE_DIR}/src/factorial.cpp
    ${CMAKE_SOURCE_DIR}/src/FileInputParser.cpp
    ${CMAKE_SOURCE_DIR}/src/overlapIntegrals.cpp
)

add_executable(main ${SOURCES})

find_package(Armadillo REQUIRED)

if(ARMADILLO_FOUND)
    include_directories(${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(main PRIVATE ${ARMADILLO_LIBRARIES})
else()
    message(FATAL_ERROR "Armadillo not found. Please install the Armadillo library.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2")

# New code to run main on each input file
# Get the list of input files
file(GLOB SAMPLE_INPUT_FILES "${CMAKE_SOURCE_DIR}/sample_input/*")

set(CALCULATED_OUTPUTS_DIR "${CMAKE_SOURCE_DIR}/calculated_outputs")

foreach(input_file ${SAMPLE_INPUT_FILES})
    get_filename_component(input_filename ${input_file} NAME)
    set(output_file "${CALCULATED_OUTPUTS_DIR}/${input_filename}")

    add_custom_command(
        OUTPUT "${output_file}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CALCULATED_OUTPUTS_DIR}"
        COMMAND "${CMAKE_CURRENT_BINARY_DIR}/main" "${input_file}" > "${output_file}"
        DEPENDS main "${input_file}"
        COMMENT "Generating output for ${input_file}"
    )
    list(APPEND OUTPUT_FILES "${output_file}")
endforeach()

add_custom_target(run_all_outputs ALL DEPENDS ${OUTPUT_FILES})
