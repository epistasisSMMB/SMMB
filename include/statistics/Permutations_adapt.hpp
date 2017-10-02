#ifndef PERMUTATIONS_ADAPT_HPP
#define PERMUTATIONS_ADAPT_HPP

#include"Permutations.hpp"
#include "common.h"
#include <string>
#include <list>

class Permutations_adapt : public Permutations
{
public:
    Permutations_adapt(blas_column const& var, blas_column const& phenos, blas_matrix & permuted_phenos, double alpha, double precision, std::string _mode);
    Permutations_adapt(blas_column const& var,
                       blas_column const& phenos,
                       std::list<unsigned> & conditioning_set,
                       blas_matrix & permuted_phenos,
                       double alpha,
                       double precision,
                       std::string _mode);

    Permutations_adapt(blas_column const& var, blas_column const& phenos, blas_matrix & permuted_phenos, double alpha, double precision);

    Permutations_adapt(blas_column const& var, blas_column const& phenos, blas_matrix & permuted_phenos, int r, int n, std::string _mode);


    static int choose_r(double alpha, double precision);
    static int choose_n(double alpha, double precision);
    void run();
    void run(std::list<unsigned> const& l);
    double correction();
    unsigned get_r() const;

protected:
    /* This function choose_r determines the number of test
        statistics that should be sampled in adaptive permutation
        sampling such that p-values at the desired level of
        significance (alpha) will be sampled with a 68% CI
        (about 1 standard error) contained within a specified
        level of precision (c).*/
    void choose_r(double alpha);

    blas_matrix & _permuted_phenos;
    double _precision;    // Desired precision (advised to set in a range from 0.1 to 0.3 (more to less permutations respectively)
    unsigned _Ri; // count number of permuted t_stat exceedind observed t_stat (obs_t_stat)
    unsigned _Bi; // count permutations (needed for final correction)
    std::string _mode;
    unsigned _r;  /* Adaptative number of permutations.
                     Number of permutations with permuted t_stat > observed t_stat to reach for stopping*/
};

#endif // PERMUTATIONS_ADAPT_HPP
