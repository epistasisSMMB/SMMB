# SMMB documentation
## Stochastic Multiple Markov Blankets (SMMB)
SMMB is a cross-platform software package, as soon as compilation is available on desired workstation.

## Needs for compilation
* g++ version compatible with c++11 functionalities
* boost library

## Compilation instructions
* Change the line n°13 in the Makefile :
    BOOST_FOLDER=value
It is needed to replace "value" with the path to the installed boost library

* Execute following command lines :
    cd <path to smmb project directory>
    make


## Clean the project for a complete re-compilation
    cd <path to smmb project directory>
    make remove


## Execution with command line
    ./bin/executable <path_to_genotypes> <path_to_phenotypes>
