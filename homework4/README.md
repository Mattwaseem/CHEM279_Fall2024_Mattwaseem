<!DOCTYPE html>
# Implementation of CNDO/2 Method for Molecular Calculations
## Numerical Algorithms Applied to Computational Quantum Chemistry
### Homework 4: Implementation of Complete Neglect of Differential Overlap (CNDO/2) Method for Molecular Calculations
** ----------------------------------------------------------------**
#### Problem 1: Building the CNDO/2 Fock Matrix for Simple Molecules
#### Problem 2: Evaluating the Open-Shell CNDO/2 Self-Consistent Field Equations

<h2>Assignment Overview</h2>

<p>The goal of this assignment was to extend concepts from previous homework (EHT, Matrix Overlap calculations, etc.) to implement the CNDO/2 method, which is a semi-empirical Hartree-Fock calculation.</p>

<p>The key steps for this assignment included:</p>

<ol>
    <li><strong>Parsing molecule input files</strong> (code recycled from previous homework)</li>
    <li><strong>Building the CNDO/2 Fock Matrix</strong></li>
    <li><strong>Solving the Self Consistent Field (SCF) equations</strong></li>
    <li><strong>Calculating the total energy</strong></li>
</ol>

<p>The assignment was organized into two main problems with a supplementary section bonus section :) :</p>

<ul>
    <li><strong>Problem 1:</strong> Correctly Building the CNDO/2 Fock Matrix</li>
    <li><strong>Problem 2:</strong> Calculating the SCF equations with the open-shell CNDO/2 method</li>
    <li><strong>Problem 3:</strong> Applying the CNDO/2 method to more complex molecules such as N₂ and O₂</li>
</ul>

<hr>

<h2>How to Compile and Run</h2>

<h3>Dependencies:</h3>

<p>Make sure you have the <strong>Armadillo</strong> library installed on your system.</p>

<h3>Preparation:</h3>

<p>Make sure you have the input files under the <code>sample_input</code> directory.</p>

<h3>Compilation:</h3>

<ol>
    <li>Navigate to the <code>build</code> directory (create it if it doesn't exist):

    <pre><code>mkdir build
cd build</code></pre>
    </li>
    <li>Run CMake to configure the build:

    <pre><code>cmake ..</code></pre>
    </li>
    <li>Compile the code using Make:

    <pre><code>make</code></pre>
    </li>
</ol>

<h3>Execution:</h3>

<ol start="4">
    <li>Run the program:

    <pre><code>make run</code></pre>
    </li>
</ol>

<p>Results from the terminal will be generated in the <code>calculated_outputs</code> directory for each input file. Each input file will have a corresponding output file with a similar name (e.g., <code>H2.txt</code> → <code>H2_calculated_output.txt</code>).</p>

<h3>Cleaning Up:</h3>

<ol start="5">
    <li>To clean all files and artifacts:

    <pre><code>make clean-all</code></pre>
    </li>
</ol>

<hr>

<h2>Directory Structure</h2>

<pre><code>.
├── CMakeLists.txt
├── GammaPseudoCode.cpp
├── README.md
├── STO3G_info
│   ├── C_STO3G.txt
│   ├── F_STO3G.txt
│   ├── H_STO3G.txt
│   ├── N_STO3G.txt
│   └── O_STO3G.txt
├── bin
│   └── CNDO2
├── build
├── calculated_outputs
│   ├── H2_calculated_output.txt
│   ├── HF_calculated_output.txt
│   ├── HO_calculated_output.txt
│   ├── N2_calculated_output.txt
│   └── O2_calculated_output.txt
├── homework4.pdf
├── include
│   ├── CartesianGaussian.h
│   ├── Constants.h
│   ├── DensityMatrix.h
│   ├── FileInputParser.h
│   ├── FockMatrix.h
│   ├── GammaCalculator.h
│   ├── Molecule.h
│   ├── OverlapMatrix.h
│   └── STO3GBasis.h
├── sample_input
│   ├── H2.txt
│   ├── HF.txt
│   ├── HO.txt
│   ├── N2.txt
│   └── O2.txt
├── sample_output
│   ├── H2.txt
│   ├── HF.txt
│   └── HO.txt
└── src
    ├── CartesianGaussian.cpp
    ├── Constants.cpp
    ├── DensityMatrix.cpp
    ├── FileInputParser.cpp
    ├── FockMatrix.cpp
    ├── GammaCalculator.cpp
    ├── Molecule.cpp
    ├── OverlapMatrix.cpp
    ├── STO3GBasis.cpp
    └── main.cpp
</code></pre>

<hr>

<h2>Problem 1: Building the CNDO/2 Fock Matrix for Simple Molecules</h2>

<p>To implement the CNDO/2 method, the first step was to build the Fock matrix by correctly parsing the input files to obtain the molecular coordinates. The input file parser was designed to handle both closed-shell and open-shell molecules. After parsing the files, the correct number of atomic orbitals was determined based on how each atom contributes a specific number of atomic orbitals corresponding to their valence electrons.</p>

<p>In the CNDO/2 method, a minimal basis set consisting only of valence electrons is used. The total number of basis functions (<strong>N</strong>) was calculated by summing the contributions from each atom. The atomic orbital information was stored to facilitate the calculation of overlap integrals (using the overlap matrix code from previous homework) and Fock matrix elements.</p>

<p>The overlap matrix elements between atomic orbitals were computed. Calculating the overlap integrals was essential for solving the generalized eigenvalue problems in the SCF procedure in the second part of the assignment. The density matrices for alpha and beta spins were also initialized and used as the starting point for the iterative SCF process.</p>

<h3>Fock Matrix Construction</h3>

<p>Lastly, the CNDO/2 Fock matrix was built using several empirical parameters and approximations. The diagonal elements were calculated using the following formula:</p>

<h4>Diagonal Elements:</h4>

<p>
\[ F^\alpha_{\mu\mu} = -\frac{1}{2} (I_\mu + A_\mu) + \left[ \left( P^\text{tot}_{AA} - Z_A \right) - \left( P^\alpha_{\mu\mu} - \frac{1}{2} \right) \right] \gamma_{AA} + \sum_{B \ne A} \left( P^\text{tot}_{BB} - Z_B \right) \gamma_{AB} \]
</p>

<p>where:</p>

<ul>
    <li><code>F<sup>α</sup><sub>μμ</sub></code> is the diagonal element of the alpha-spin Fock matrix.</li>
    <li><code>I<sub>μ</sub></code> and <code>A<sub>μ</sub></code> are the ionization potential and electron affinity of orbital <code>μ</code>.</li>
    <li><code>P<sup>α</sup><sub>μμ</sub></code> is the diagonal element of the alpha-spin density matrix.</li>
    <li><code>P<sup>tot</sup><sub>AA</sub></code> is the total electron density on atom <code>A</code>.</li>
    <li><code>Z<sub>A</sub></code> is the valence atomic number of atom <code>A</code>.</li>
    <li><code>γ<sub>AA</sub></code> is the two-electron repulsion integral for atom <code>A</code>.</li>
</ul>

<h4>Off-Diagonal Elements:</h4>

<p>
\[ F^\alpha_{\mu\nu} = \frac{1}{2} (\beta_A + \beta_B) S_{\mu\nu} - P^\alpha_{\mu\nu} \gamma_{AB} \]
</p>

<p>where:</p>

<ul>
    <li><code>β<sub>A</sub></code> and <code>β<sub>B</sub></code> are the bonding parameters for atoms <code>A</code> and <code>B</code>.</li>
    <li><code>S<sub>μν</sub></code> is the overlap matrix element between atomic orbitals <code>μ</code> and <code>ν</code>.</li>
    <li><code>P<sup>α</sup><sub>μν</sub></code> is the off-diagonal element of the alpha-spin density matrix.</li>
    <li><code>γ<sub>AB</sub></code> is the two-electron repulsion integral between atoms <code>A</code> and <code>B</code>.</li>
</ul>

<h3>Gamma Values:</h3>

<p>The <code>γ<sub>AB</sub></code> values were determined based on:</p>

<p>
\[ \gamma_{AB} = \frac{14.397}{R_{AB} + \frac{1}{2} \left( \frac{1}{U_A} + \frac{1}{U_B} \right)} \]
</p>

<p>Finally, unit conversions were performed to calc. all energies in electron volts (eV) from atomic units. The Fock matrix was built, and the calculated outputs were written to corresponding output files with names similar to the input files (e.g., <code>H2.txt</code> → <code>H2_calculated_output.txt</code>).</p>

<p>Overall, this process built the foundational steps fo`r a semi-empirical Hartree-Fock calculations, bridging the molecular input data with the SCF process to calculate the total energies of the system.</p>

<hr>

<h2>Problem 2: Evaluating the Open-Shell CNDO/2 Self-Consistent Field Equations</h2>

<p>The Fock matrix from Problem 1 was used to implement the code for <code>f<sup>α</sup></code> and <code>f<sup>β</sup></code> to develop the self-consistent field (SCF) algorithm for the CNDO/2 method. This involved solving the SCF equations in the atomic orbital (AO) basis for both the alpha and beta spins:</p>

<p>
\[ \begin{align} f^\alpha \mathbf{c}^\alpha &= \mathbf{c}^\alpha \boldsymbol{\epsilon}^\alpha \\ f^\beta \mathbf{c}^\beta &= \mathbf{c}^\beta \boldsymbol{\epsilon}^\beta \end{align} \]
</p>

<p>where:</p>

<ul>
    <li><code>f<sup>α</sup></code> and <code>f<sup>β</sup></code> are the Fock matrices for alpha and beta spins, respectively.</li>
    <li><code>c<sup>α</sup></code> and <code>c<sup>β</sup></code> are the molecular orbital (MO) coefficient matrices.</li>
    <li><code>ε<sup>α</sup></code> and <code>ε<sup>β</sup></code> are the diagonal matrices of orbital energies (eigenvalues).</li>
</ul>

<h3>Implementation Steps:</h3>

<ol>
    <li><strong>Initial Guess:</strong>
        <ul>
            <li>The density matrices for alpha and beta spins were initialized to zero matrices.</li>
        </ul>
    </li>
    <li><strong>Fock Matrix Construction:</strong>
        <ul>
            <li>The Fock matrices <code>f<sup>α</sup></code> and <code>f<sup>β</sup></code> were constructed using the current density matrices.</li>
        </ul>
    </li>
    <li><strong>Solving the Eigenvalue Problems:</strong>
        <ul>
            <li>The eigenvalue problems were solved to calculate the molecular orbital coefficients and eigenvalues using the equations above.</li>
            <li>The Armadillo library was used for efficient matrix operations.</li>
        </ul>
    </li>
    <li><strong>Updating Density Matrices:</strong>
        <ul>
            <li>The old density matrices were stored.</li>
            <li>New density matrices were assembled by occupying the <code>p</code> lowest energy alpha MOs and the <code>q</code> lowest energy beta MOs:</li>
        </ul>

        <p>
        \[ P^\alpha_{\mu\nu} = \sum_{i=1}^{p} c^\alpha_{\mu i} c^\alpha_{\nu i} \]
        \[ P^\beta_{\mu\nu} = \sum_{i=1}^{q} c^\beta_{\mu i} c^\beta_{\nu i} \]
        </p>
    </li>
    <li><strong>Convergence Check:</strong>
        <ul>
            <li>The maximum change in the density matrices was calculated.</li>
            <li>If the change was less than the specified tolerance (e.g., \( 1 \times 10^{-6} \)), the SCF process was considered converged.</li>
        </ul>
    </li>
    <li><strong>Total Energy Calculation:</strong>
        <ul>
            <li>Upon convergence, the total energy was calculated using:</li>
        </ul>

        <p>
        \[ E_{\text{CNDO/2}} = \frac{1}{2} \sum_{\mu\nu} P^\alpha_{\mu\nu} \left( h_{\mu\nu} + f^\alpha_{\mu\nu} \right) + \frac{1}{2} \sum_{\mu\nu} P^\beta_{\mu\nu} \left( h_{\mu\nu} + f^\beta_{\mu\nu} \right) + \sum_{A} \sum_{B<A} \frac{Z_A Z_B}{R_{AB}} \]
        </p>

        <p>where <code>h<sub>μν</sub></code> is the core Hamiltonian matrix.</p>
    </li>
</ol>

<h3>Results:</h3>

<p>The implementation of the open-shell CNDO/2 was successful. The SCF algorithm converged, and the total energies were calculated. The eigenvalues and molecular orbital coefficients were obtained after convergence. While the output total energies did not exactly match the expected energies from sample output files, the overall implementation was successful, and the total energy calculations could have been improved if I had more time.</p>

<hr>

<h2>Problem 3: Going Further – Building CNDO/2 for N₂ and O₂</h2>

<p>The extension of the CNDO/2 code to N₂ and O₂ molecules provides insights into the chemical properties and validated the code's capability to handle different electron configurations. Although the total energies did not match experimental results, the ability to build the CNDO/2 method for other systems shows the flexibility of the current implementation.</p>

<h3>N₂ Simulation:</h3>

<ul>
    <li><strong>Results:</strong>
        <ul>
            <li>A strong triple bond and high stability were reflected in the results, with a higher energy output in the <code>N2_calculated_output.txt</code> file.</li>
            <li>The code accurately represented a closed-shell diatomic molecule.</li>
        </ul>
    </li>
</ul>

<h3>O₂ Simulation:</h3>

<ul>
    <li><strong>Results:</strong>
        <ul>
            <li>Open-shell configuration modeled.</li>
            <li>The simulation demonstrated the code's ability to handle unpaired electrons.</li>
        </ul>
    </li>
</ul>

<p><strong>Overall, the CNDO/2 code proved to be a useful for exploring molecular properties and behaviors.</strong></p>

<hr>

<h2>Conclusion</h2>

<p>This project successfully implemented the CNDO/2 method for calculating molecular energies and properties. By building the Fock matrices, solving the SCF equations, and calculating total energies, the code can handle both simple molecules and more complex systems like N₂ and O₂. While there is room for improvement in matching experimental energies, the current implementation shows the potential of semi-empirical methods in computational quantum chemistry.</p>

<hr>

<p><em>Note: For detailed code and implementation specifics, please refer to the source files in the repository.</em></p>

<hr>

</body>
</html>
