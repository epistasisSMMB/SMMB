#ifndef PERMUTATIONS_HPP
#define PERMUTATIONS_HPP

#include <queue>

#include "common.h"
#include "G2_test_indep.hpp"

class Permutations
{

public:
    /*!
     * \brief Constructor
     * \param var The variable to test for independence. It is variable S in the test "indep(S,P)?" (e.g. a snp)
     * \param phenos This variable will be permuted. It is variable P in the test "indep(S,P)?" (e.g. phenotype)
     * \param n Number of permutations to run
     */
    Permutations(blas_column const& var, blas_column const& phenos, unsigned n);

    /*!
     * \brief Constructor n°2. Number of permutations to run is automatically calculated thanks to (alpha, precision)
     * \param var
     * \param phenos
     * \param alpha Statistical significance level. Needed to compute the number of permutations
     * \param precision Advised to set between 0.1 and 0.3 (more to less permutations respectively)
     */
    Permutations(blas_column const& var, blas_column const& phenos, double alpha, double precision);

    /*!
     * \brief Correct a given p_values by the permutations (previously done)
     * \param pval
     * \return
     */
    double correction(double pval);

    /*!
     * \brief run permutations
     */
    void run();

    unsigned get_n() const;

protected:
    /* */
    /*!
     * \brief choose_b This function choose_b determines the maximal number of
                        permutations for either adaptive or standard permutation.
                        It is a function of alpha (the p-value you would like to
                        estimate) and c (a desired precision level of p-value
                        estimation), where the standard error of the estimated
                        p-value is c*alpha.
     * \param alpha
     * \param precision
     */
    void choose_b(double alpha, double precision);

    blas_column const& _var;                // variable to test for association
    blas_column const& _original_phenos;    // variable to permute (a backup of its original permutation)
    unsigned _n;                            // number of permutations
    std::priority_queue<double, std::vector<double>, std::greater<double> > _pvalues_dist; // Distribution of p-values built from permutations
};

#endif // PERMUTATIONS_HPP
