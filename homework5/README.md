<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CHEM 279 - Homework 5</title>
</head>
<body>
    <header>
        <h1>Chemistry</h1>
        <h2>Homework 5: Implementation of Complete Neglect of Differential Overlap 2 (CNDO/2) Method for Molecular Calculations</h2>
        <p>
            <strong>Link to Github Repo:</strong>
            <a href="https://github.com/Mattwaseem/Chem279-Numerical-Algos-HW5" target="_blank">https://github.com/Mattwaseem/Chem279-Numerical-Algos-HW5</a><br>
            Repo: Publicly Available & TA invited.<br>
            Username: MattWaseem
        </p>
    </header>
    <hr>
    <section>
        <h3>Assignment Overview</h3>
        <p>
            The goal here was to extend concepts from previous homework (EHT, Matrix Overlap calculations, Fock matrix, Density Matrix, CNDO/2, etc.) to first derive and then implement the full ab initio SCF energy for a molecular system. This is a semi-empirical calculation that builds on the concepts from previous homework. Rather than determining the energy, the objective is to obtain the gradient of it with respect to the position of the nuclei of the atoms. The key steps for this assignment included:
        </p>
        <ul>
            <li>Deriving the theory first.</li>
            <li>Parsing molecule input file (code recycled from previous homework).</li>
            <li>Building the CNDO/2 Fock Matrix from Homework 4.</li>
            <li>Solving the Self-Consistent Field equation.</li>
            <li>Taking derivative of SCF with respect to positions RA and RB.</li>
            <li>Calculating the gradient of energy.</li>
        </ul>
        <p>
            In a similar manner to previous homeworks, the assignment was organized in two sections with a supplementary section:
        </p>
        <ul>
            <li><strong>Problem 1:</strong> Correctly derive the energy gradient and use the CNDO/2 method.</li>
            <li><strong>Problem 2:</strong> Calculating the SCF equation with the open-shell CNDO/2.</li>
        </ul>
    </section>
    <section>
        <h3>Problem 1: Deriving The Theory Energy Gradient</h3>
        <p>
            CNDO/2 or Complete Neglect of Differential Overlap two is the second method of the CNDO approach that splits the total energy of a molecule into parts focused on:
        </p>
        <ul>
            <li><strong>Overlap:</strong> Between atomic orbitals.</li>
            <li><strong>Electron Repulsion:</strong></li>
        </ul>
        <p>
            The total energy represented by the CNDO/2 method is as follows:
        </p>
        <p>
            <strong>Where the symbols mean:</strong>
        </p>
        <ul>
            <li>&psi;<sub>i</sub> and &psi;<sub>j</sub> refer to different atomic orbitals.</li>
            <li>P represents the density matrix - which is focused on electron density around a molecule.</li>
            <li>h is the integral determining the attraction between electrons and nuclei.</li>
            <li>f is a representation of the Fock Matrix, which is a matrix of electron repulsions.</li>
            <li>Z<sub>A</sub> and Z<sub>B</sub> are the nuclear charges for corresponding atoms A and B.</li>
            <li>R<sub>AB</sub> is the distance between the nuclei of such atoms.</li>
        </ul>
        <p>
            The above terms represent a simplified summation of interactions between electrons and nuclei within a molecule that make up the total energy. The difference in energy with respect to distance between atoms of a molecule must be calculated to approximate gradients by taking the derivative of the CNDO/2 (total energy) equation with respect to the position of the nuclei of two atoms RA and RB.(please refer to pdf my hand written derivation)
        </p>
        <p>
            <strong>Derivative:</strong> Refer to handwritten derivations for detailed calculations.
        </p>
    </section>
    <footer>
        <p><em>Handwritten derivations are available PDF notes.</em></p>
    </footer>
</body>
</html>
