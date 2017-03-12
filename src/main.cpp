#include <iostream>
#include <stdlib.h>
#include <ctime>

#include "Parameters_file_parsing.hpp"
#include "CSVParser.hpp"
#include "Services.hpp"
#include "Smmb_usecase.hpp"
#include "Smmb.hpp"
#include "Permutations_adapt.hpp"

#include "G2_conditional_test_indep.hpp"
#include "common.h"

using namespace std;

int main(int argc, char *argv[])
{

    if(argc < 3)
    {
        cerr << "Missing parameter :\n"
             << "\t./executable <path_to_genotypes> <path_to_phenotypes>"
             << endl;
        exit(-1);
    }

    // Arguments
    string genos_file = argv[1];
    string phenos_file = argv[2];

//  PARAMETERS
    Parameters_file_parsing params;
    params.list_parameters();

    params.genos_file = genos_file;      // genotype file name (path)
    params.phenos_file = phenos_file;    // phenotype file name (path)
    int header = params.header;
    char separator = params.separator;

//  DATA IMPORTATION
    CSVParser<int> genos_csv(genos_file, separator, header);
    CSVParser<int> phenos_csv(phenos_file, separator, header);
    blas_matrix genos = genos_csv.data();
    blas_matrix phenos_m = phenos_csv.data();
    blas_column phenos(phenos_m, 0);

    params.update_subset_size_large(genos.size2());

    cout << endl << "Data imported : " << genos.size1() << " individuals X " << genos.size2() << " SNPs" << endl;

//  GENERATION OF PERMUTED PHENOTYPES
    cout << "Generation of permuted phenotypes ... ";
    blas_column geno(genos,1);
    Permutations_adapt permut(geno, phenos, genos, params.alpha, params.precision);
    unsigned n_permut = permut.get_n();
    cout << "number of permutations : " << n_permut << endl;
    blas_matrix permuted_phenos(phenos.size(), n_permut);
    cout << "Permuted phenotypes generated." << endl;
    Services::generate_permutations(phenos, permuted_phenos);
	
//  SMMB usecase
    Smmb_usecase sm2b_usecase(genos, phenos, params, permuted_phenos);

    Services::frame_cout("END OF PROGRAM");

    return 0;
}
