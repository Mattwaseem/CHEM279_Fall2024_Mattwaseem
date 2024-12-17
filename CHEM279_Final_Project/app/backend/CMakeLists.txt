cmake_minimum_required(VERSION 3.14)

# Project name and version
project(WaterSim VERSION 0.1 DESCRIPTION "Simple C++ NodeJS Water Simulation" LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Manually specify the Node.js and node-addon-api include paths. These are used for compiling our C++
# code into a node js module.
set(NODE_JS_INCLUDE_DIRS "/usr/include/nodejs/src")
set(NODE_ADDON_API_INCLUDE_DIRS "/repo/app/frontend/node_modules/node-addon-api/")

# Include directories
include_directories(
    "./include/"
    ${NODE_JS_INCLUDE_DIRS} 
    ${NODE_ADDON_API_INCLUDE_DIRS}
)

# Add the shared library target (Node.js addon)
add_library(${PROJECT_NAME} SHARED src/water_sim.cpp)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(OpenMP REQUIRED)

# Set output file to be a .node (Node.js addon format)
set_target_properties(${PROJECT_NAME} PROPERTIES
    PREFIX ""
    SUFFIX ".node"
)

# Specify the path to the Node.js shared library (libnode.so)
set(NODE_JS_LIBRARY "/usr/lib/aarch64-linux-gnu/libnode.so")

# Link Node.js library
target_link_libraries(${PROJECT_NAME}   PRIVATE ${NODE_JS_LIBRARY}
                                        PRIVATE ${ARMADILLO_LIBRARIES}
                                        PRIVATE OpenMP::OpenMP_CXX)
target_compile_options(${PROJECT_NAME} PRIVATE -fopenmp)
