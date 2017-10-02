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
//    cout << "constructor of sm2b ok\n";
    sm2b.run();
//    cout << "run of smb2b ok\n";
    sm2b.make_consensus(permuted_phenos);
//    cout << "make_consensus of sm2b ok\n";

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    double duration = duration_cast<milliseconds>(t2-t1).count();
    sm2b.write_result(duration);
}
