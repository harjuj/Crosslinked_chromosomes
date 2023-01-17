# Crosslinked chromosome

## Introduction
This repository is based on [the MaxEnt model for a Caulobacter crescentus chromosome](https://doi.org/10.1038/s41467-021-22189-x). It has been edited by:

1. Removing the Hi-C map constraints for the polymer
2. Adding crosslinks that can bind two monomers that are on the same lattice site.

Crosslinks can bind together two sites on the polymer. They have a chemical potential and a site-associated binding energy. A monomer with a crosslink cannot move until it is unbound.

## Installation and dependencies
The code can simply be installed using git. In addition, two open-source C++ libraries are required:

1. [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
2. [Boost](https://www.boost.org/users/download/)

Once these are installed, you should edit the 'CMakeLists' file to include the right paths, or link the libraries using your compiler of choice.

The program uses 'pthread', so a compatible compiler is required.

## Set-up and running
The 'main.cpp' file has two directory paths that we need; 'forward_output_directory' and 'dir'. These default to 'out' and 'initial_configs' respectively, but can be edited in main.cpp.

'forward_output_directory' is the absolute path to the folder where the simulation output directory (default name date and potential details) is created.

'dir' is the absolute path to the folder where optional previous configurations that can be used as input are stored.

Once these paths are set, the program is ready to be used.

    cmake .
    make
    ./Forward_crosslinks [0=random potential/1=cosine potential] [chemical potential] [energy amplitude] ([input_directory])

'input_directory' is an optional argument; if provided, the program will load initial configurations from './initial_configs/[input_directory] instead of starting with a generated coiled polymer configuration. 

## Help and contact details
For help or questions with the repository, start an issue or email j.k.harju[at]vu.nl.
