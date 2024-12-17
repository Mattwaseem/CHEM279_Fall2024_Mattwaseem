# Numerical Algorithms Applied to Computational Quantum Chemistry
## Homework 3: Implementation of Extended Huckel Theory (EHT)

### Link to Github Repo: https://github.com/Mattwaseem/Chem279-Numerical-Algos-HW3
### Repo: publicly available.
### Username: MattWaseem


—

#### Assignment Overview:
<p> The goal here was to extend the concepts of overlap integrals, focusing on the next steps in molecular orbital theory calculations using Gaussian basis sets by:
Parse hydrocarbon molecule input and set up basis functions.
Compute the overlap matrix using contracted Gaussian-type orbitals (cGTOs).
Assemble and diagonalize the Hamiltonian matrix to solve for molecular orbitals.
In similar manner the project was organized in three sections with a supplementary section.
Problem 1: focused on molecule input and basis function set up
Problem 2: focused on Matrix calculation
Problem 3: focused on Hamiltonian and Diagonalization
The best way to approach these problems was to modularize the code, so that subsequent tasks can reuse code from previous tasks. To begin below is the directory structure for the assignment. For the purpose of this discussion files within each directory will not be listed, please refer to my github for file specific details. This is to get a broad overview picture of the directory structure and ease of navigation. </p>

***
\documentclass{article}
\usepackage{verbatim}

\begin{document}

\begin{verbatim}
.
├── CMakeLists.txt
├── Chem279_HW3__pdf.pdf
├── README.md
├── bin
│   └── MoleculeSetup
├── build
├── Makefile
├── cmake_install.cmake
├── calculated_output
│   ├── C2H2_output.txt
│   ├── C2H4_output.txt
│   └── H2_output.txt
├── homework3.pdf
├── include
│   ├── CartesianGaussian.hpp
│   ├── Hamiltonian.hpp
│   ├── OverlapMatrix.hpp
│   ├── basis_function.hpp
│   ├── input_parser.h
│   └── molecule.hpp
├── sample_input
│   ├── C2H2.txt
│   ├── C2H4.txt
│   └── H2.txt
├── sample_output
│   ├── C2H2.out
│   ├── C2H4.out
│   └── H2.out
└── src
    ├── CartesianGaussian.cpp
    ├── Hamiltonian.cpp
    ├── OverlapMatrix.cpp
    ├── basis_function.cpp
    ├── input_parser.cpp
    ├── main.cpp
    └── molecule.cpp
\end{verbatim}

\end{document}
***
How to Compile and Run:
Ensure you have armadillo library installed in your system
Also make sure you have the input files under the input directory.
CD to build directory (if there is not one, mkdir build, please).
Once in build directory, run the following command to use CMakeList.txt to compile code:
Cmake ..
Make
Make run
Results from the terminal should be generated in the calculated_output directory for each input file.
Once, satisfied with response, you can run the following to clean all files and artifacts:
Make clean-all



Problem 1: Molecule Input and Basis Function Setup
First task was to read in molecular coordinates and use the formula to derive the number of basis functions and electrons. A hydrocarbon molecule was specified by the input files in the format E X Y Z, where E is the element, and X, Y, Z are the atomic coordinates. The number of basis functions (N) was calculated using the relation N = 4a + b, where a is the number of carbon atoms and b is the number of hydrogen atoms.
Next the basis functions were constructed using contracted Gaussians, where each basis function consists of three primitive Gaussians, the sum of those primitive Gaussians made up the contracted Gaussian.. The contracted Gaussian functions are organized by shells (S and P shells) and centered on the atomic coordinates. For each basis function, the was stored:
Center: Atomic position (X, Y, Z)
Quantum Numbers: Angular momentum (l, m, n)
Exponents (α): Radial decay constants for the Gaussian
Contraction Coefficients (d): Coefficients for combining primitive Gaussians
Normalization Constants (N): Ensuring the overlap of each primitive with itself is one.
This set up the atomic orbital basis and was helpful in  allowing for further overlap and Hamiltonian matrix calculations.
Problem 2: Overlap Matrix Calculation
I used my homework on two primitive Gaussian overlap integrals implementations for this problem to calculate the overlap matrix between contracted Gaussian atomic orbitals. Since contract Guassian is a sum of primitive Guassains the only modification needed was to implement the normalization. Indeed the normalization aspect did introduce difficulty as it was hard to determine how to implement a sum of normalized values for each Guassian. The overlap matrix elements  Sμν were calculated by summing over the overlap integrals of the primitive Gaussians comprising the contracted functions:

Where Skl is the overlap integral between two unnormalized primitive Gaussians. Using the normalization constants calculated earlier, the program computed the full overlap matrix for the molecule.


Problem 3: Hamiltonian Assembly and Diagonalization
Once the molecular input and overlap matrix was computed, the next task was to assemble the Hamiltonian matrix. The diagonal elements of the Hamiltonian correspond to the atomic orbital energies, using standard parameters provided for hydrogen and carbon atoms ( -13.6 eV for H, -21.4 eV for C(2s), -11.4 eV for C(2p)). The off-diagonal elements were computed as:

Where K = 1.75 is a scaling factor, and A(μ) and L(μ) are the atom and angular momentum correlated with the basis function ωμ.
The next step was solving the generalized eigenvalue problem to find the molecular orbital coefficients and eigenvalues:

Using Armadillo, the program diagonalizes the Hamiltonian matrix, obtaining the molecular orbital coefficients (C) and orbital energies (ε).
Finally, the total energy of the molecule was computed by summing over the occupied molecular orbitals:

For the H2 molecule, the calculated energy was approximately -35 eV, consistent with reference values. The program was also tested on C2H2 and C2H4, however, I was not able to get an accurate energy value for these two more complex molecules. I assume this has to do with my energy calculation, my scaling factor, and the bond between C-C, C-H which is more complex compared to just an H-H bond. Hence, my energy calculation for H2 being close to the reference, but my total energy not getting anything near what the reference energy was for C2H2 and C2H4.


Problem 4: GOING FURTHER Resonance Energy of Benzene:
