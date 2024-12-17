# Numerical Algorithms Applied to Computational Quantum Chemistry
## Homework 2: Analytical and Numerical Integration for Overlap Integrals

**Link to GitHub Repo:** [https://github.com/Mattwaseem/chem279-Waseem-hw](https://github.com/Mattwaseem/chem279-Waseem-hw)
**This is my personal Repo and it is made:** Publicly available.
**Username:** MattWaseem

---

## Overview

This homework focused on numerical concepts covered in first five weeks of lecture on approaches to estimate energy of MO through a summation of atomic orbitals by calculating the overlap integrals between Gaussian functions using both numerical and analytical approaches. This was divided into two problems:

- **Problem 1:** Numerical 1D Overlap Integral
- **Problem 2:** Analytical 3D Overlap Integral

Each problem is approached by implementing specific algorithms and executed code for integration in modular method, and the results are output into calculated output directory after processing the input files.

---

## Directory Structure

```bash
.
├── CMakeLists.txt                  # CMake build configuration file
├── Plots                           # Directory to store any generated plots (currently empty)
├── bin                             # Directory for compiled executable 'homework2_exec'
│   └── homework2_exec
├── build                           # Directory containing build-related files
├── calculated_output_p1            # Directory for output files of Problem 1 (Numerical)
│   ├── 1-output.txt
│   ├── 2-output.txt
│   └── 3-output.txt
├── calculated_output_p2            # Directory for output files of Problem 2 (Analytical)
│   ├── 1-output.txt
│   ├── 2-output.txt
│   └── 3-output.txt
├── include                         # Directory containing header files for the project
│   ├── CartesianGaussian.hpp
│   ├── Gaussian1D.hpp
│   ├── InputParser.hpp
│   ├── Integrator1D.hpp
│   └── OverlapIntegral.hpp
├── sample_input                    # Directory containing sample input files
│   ├── analytical                  # Input files for Problem 2 (Analytical 3D)
│   ├── numerical                   # Input files for Problem 1 (Numerical 1D)
├── sample_output                   # Sample output files for comparison
├── src                             # Directory containing source code files
│   ├── CartesianGaussian.cpp
│   ├── Gaussian1D.cpp
│   ├── InputParser.cpp
│   ├── Integrator1D.cpp
│   ├── OverlapIntegral.cpp
│   └── main.cpp


## How to Build and Run the Code
1. Go to the build directory
2. Run the following commands to utilize the Cmake build system set up to execute this code
        cd build
        cmake ..
        make
3. To run the executable
    make run-all
4. This command will execute the program for all input files (numerical and analytical).
    The output will be saved in the calculated_output_p1 and calculated_output_p2 directories, corresponding to Problem 1 and Problem 2, respectively.
5. TO delete output files and restart you can do the following:
    cd build
    make clean-all

For any errors or issues please refer to PDF. Thank you!
