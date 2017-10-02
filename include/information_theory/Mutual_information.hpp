#ifndef MUTUAL_INFORMATION_HPP
#define MUTUAL_INFORMATION_HPP

#include "common.h"
#include "Contingency.hpp"
#include <vector>
#include <list>

class Mutual_information
{
public:
    Mutual_information(Contingency const& c);
    Mutual_information(blas_column const& genos, blas_column const& phenos);
    Mutual_information(blas_column const& genos, blas_vector const& phenos);
    Mutual_information(blas_column const& genos,
                       blas_column const& phenos,
                       std::list<unsigned> const& cond_genos_indexes);


    void run(Contingency const& c);
    void run();
    double mi() const;

protected:
    double _mi;
    double _total_obs;
    std::vector<Contingency> _contingencies;
};

#endif // MUTUAL_INFORMATION_HPP
