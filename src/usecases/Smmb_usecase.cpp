#include "Smmb_usecase.hpp"
#include "Smmb.hpp"

#include <chrono>

using namespace std;
using namespace std::chrono;

Smmb_usecase::Smmb_usecase(blas_matrix & genos, blas_column & phenos, Parameters_file_parsing & params, blas_matrix & permuted_phenos)
    : _genos(genos), _phenos(phenos), _params(params)
{
    run(permuted_phenos);
}


void Smmb_usecase::run(blas_matrix & permuted_phenos)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    Smmb sm2b(_genos, _phenos, _params);
    sm2b.run();
    sm2b.make_consensus(permuted_phenos);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    double duration = duration_cast<milliseconds>(t2-t1).count();
    sm2b.write_result(duration);
}
