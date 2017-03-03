#ifndef SMMB_USECASE_H
#define SMMB_USECASE_H

#include "common.h"
#include "Parameters_file_parsing.hpp"

class Smmb_usecase
{

public:
    Smmb_usecase(blas_matrix & genos, blas_column & phenos, Parameters_file_parsing & params, blas_matrix & permuted_phenos);
    void run(blas_matrix & permuted_phenos);

protected:
    blas_matrix & _genos;
    blas_column & _phenos;
    Parameters_file_parsing & _params;
};

#endif // SMMB_USECASE_H
