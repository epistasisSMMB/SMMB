#!/bin/bash

# Compilation if needed
make

# SMMB command line : 
# ./SMMB <path_to_genotypes> <path_to_phenotypes>
./SMMB toy_dataset/genotypes_toy_dataset.txt toy_dataset/phenotypes_toy_dataset.txt
