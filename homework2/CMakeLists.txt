cmake_minimum_required(VERSION 3.10)
project(homework2)

set(CMAKE_CXX_STANDARD 17)

# Include header files
include_directories(include)

# Find Armadillo library - to help me solving matrices easier per TA and Daktoas rec.
find_package(Armadillo REQUIRED)

# Include Armadillo's directories
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Set the output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Add source files for the project
add_executable(homework2_exec
    src/main.cpp
    src/Gaussian1D.cpp
    src/InputParser.cpp
    src/Integrator1D.cpp
    src/CartesianGaussian.cpp
    src/OverlapIntegral.cpp
)

# Link Armadillo library AFTER defining the executable
target_link_libraries(homework2_exec PRIVATE ${ARMADILLO_LIBRARIES})

# List of input files for numerical and analytical cases
set(NUMERICAL_INPUT_FILES
    ${CMAKE_SOURCE_DIR}/sample_input/numerical/1.txt
    ${CMAKE_SOURCE_DIR}/sample_input/numerical/2.txt
    ${CMAKE_SOURCE_DIR}/sample_input/numerical/3.txt
)

set(ANALYTICAL_INPUT_FILES
    ${CMAKE_SOURCE_DIR}/sample_input/analytical/1.txt
    ${CMAKE_SOURCE_DIR}/sample_input/analytical/2.txt
    ${CMAKE_SOURCE_DIR}/sample_input/analytical/3.txt
)

# Handle numerical inputs (Problem 1)
foreach(INPUT_FILE ${NUMERICAL_INPUT_FILES})
    # Get the name of the input file without extension
    get_filename_component(FILE_NAME ${INPUT_FILE} NAME_WE)

    # Define the output file path in calculated_output_p1
    set(OUTPUT_FILE ${CMAKE_SOURCE_DIR}/calculated_output_p1/${FILE_NAME}-output.txt)

    # Add a custom target to run the executable and redirect output
    add_custom_target(run_numerical_${FILE_NAME}
        COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/homework2_exec ${INPUT_FILE} > ${OUTPUT_FILE}
        DEPENDS homework2_exec
        COMMENT "Running homework2_exec with numerical input file ${FILE_NAME}.txt"
    )
endforeach()

# Handle analytical inputs (Problem 2)
foreach(INPUT_FILE ${ANALYTICAL_INPUT_FILES})
    # Get the name of the input file without extension
    get_filename_component(FILE_NAME ${INPUT_FILE} NAME_WE)

    # Define the output file path in calculated_output_p2
    set(OUTPUT_FILE ${CMAKE_SOURCE_DIR}/calculated_output_p2/${FILE_NAME}-output.txt)

    # Add a custom target to run the executable and redirect output
    add_custom_target(run_analytical_${FILE_NAME}
        COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/homework2_exec ${INPUT_FILE} > ${OUTPUT_FILE}
        DEPENDS homework2_exec
        COMMENT "Running homework2_exec with analytical input file ${FILE_NAME}.txt"
    )
endforeach()

# Add a target that runs all numerical and analytical input files
add_custom_target(run-all
    COMMAND ${CMAKE_COMMAND} --build . --target run_numerical_1
    COMMAND ${CMAKE_COMMAND} --build . --target run_numerical_2
    COMMAND ${CMAKE_COMMAND} --build . --target run_numerical_3
    COMMAND ${CMAKE_COMMAND} --build . --target run_analytical_1
    COMMAND ${CMAKE_COMMAND} --build . --target run_analytical_2
    COMMAND ${CMAKE_COMMAND} --build . --target run_analytical_3
    COMMENT "Running all numerical and analytical input files..."
)

# Custom target to clean all build artifacts and output files
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove -f
        ${CMAKE_SOURCE_DIR}/calculated_output_p1/*.txt
        ${CMAKE_SOURCE_DIR}/calculated_output_p2/*.txt
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*
        ${CMAKE_BINARY_DIR}/*
    COMMENT "Cleaning all build files, executables, and output text files"
)
