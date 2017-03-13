# SMMB documentation
## Stochastic Multiple Markov Blankets (SMMB)
SMMB is a cross-platform software package, provided compilation is available on desired workstation.

## Needs for compilation
* g++ version compatible with C++11 functionalities
* boost library. The version 1.61.0 was used by the authors. It can be found at http://www.boost.org/users/history/version_1_61_0.html

## Compilation instructions
### 1. Change the first line in the Makefile : `BOOST_FOLDER=<value>`.
Replace `<value>` with the path to the installed boost library.

### 2. Execute following command lines:
    cd <path to SMMB project directory>
    make
    
## Clean the project for a complete re-compilation
    cd <path to SMMB project directory>
    make remove

## Execute command line
    ./SMMB <path_to_genotypes> <path_to_phenotypes>

## Parameters
The complete list of parameters can be found in the file `PARAMETERS.txt`.

Each parameter is explained in `PARAMETERS.txt` and can be tuned by the user.

## Launch the analysis of the toy dataset
A Bash script is provided to launch the analysis of the toy dataset in a simple way on Linux plateforms.
The command line to execute is:

    ./launch_toy_dataset_analysis.sh

## Download or clone the SMMB source code
### 1. ZIP download 
An ZIP compressed archive is downloadable in the main web page of the SMMB project : https://github.com/epistasisSMMB/SMMB.
From this web page, the user finds a "Clone or download" button. After clicking on it, the user is invited to download a ZIP version of the project.

### 2. Clone with Git
If Git is installed on the user's workstation, SMMB can also be cloned with the following git command line:

    git clone https://github.com/epistasisSMMB/SMMB.git

Then, a new repository is created locally, which contains the SMMB project (source code, Makefile, parameters file...).
