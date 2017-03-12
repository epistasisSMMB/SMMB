#!/bin/bash

# Compilation if needed
make

# SMMB command line : 
# ./SMMB <path_to_genotypes> <path_to_phenotypes>
./SMMB toy_dataset/GENOTYPES_01p_0005h_05m_2snps_4000ind.txt toy_dataset/PHENOTYPES_4000ind.txt