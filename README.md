# SMMB documentation
## Stochastic Multiple Markov Blankets (SMMB)
SMMB is a cross-platform software package, as soon as compilation is available on desired workstation.

## Needs for compilation
* g++ version compatible with C++11 functionalities
* boost library

## Compilation instructions
### 1. Change the line n°13 in the Makefile : `BOOST_FOLDER=value`.
It is needed to replace `value` with the path to the installed boost library

### 2. Execute following command lines:
    cd \<path to smmb project directory>
    make`
    
## Clean the project for a complete re-compilation
    cd <path to smmb project directory>
    make remove

## Execution with command line
    ./bin/executable <path_to_genotypes> <path_to_phenotypes>

## Parameters
The complete list of parameters can be found in the file `PARAMETERS.txt`.
Each parameter is explained and can be tuned by the user.
