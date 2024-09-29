# Hexagonal Layout Mapping

## For Details

Refer to our publish work:
    Abhoy Kole, Kamalika Datta, Indranil Sengupta, and Rolf Drechsler. 2024. Exploiting the Extended Neighborhood of Hexagonal Qubit Architecture for Mapping Quantum Circuits. J. Emerg. Technol. Comput. Syst. Just Accepted (August 2024). https://doi.org/10.1145/3688391


## For Installation

Download (clone) the source to a local directory and run the following commands:

### For source build:

    1. chmod +x install.sh

    2. ./install.sh cmake


### For executing the mapping tool:

    1. Compiling for a fixed lookahead window:

          <<Build Path>>/HEX_QxMapping <path to source qasm file> <<path to target qasm file>> <<layout width>> <<layout height>> <<cost metric: (RCNOT = 0 and Coupling Cost = 1)>> <<approach: (RCNOT = 0, SWAP = 1 and SWAP + RCNOT = 2)>>

    2. Compiling for a variable lookahead window: 

          <<Build Path>>/HEX_QxMapping_II <path to source qasm file> <<path to target location>> <<layout width>> <<layout height>> <<cost metric: (RCNOT = 0 and Coupling Cost = 1)>>

   
    Option 2 require subfolders M0A1 and M0A2 in target location for getting compiled circuit.    