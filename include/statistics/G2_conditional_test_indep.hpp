#ifndef G2_CONDITIONAL_TEST_INDEP_HPP
#define G2_CONDITIONAL_TEST_INDEP_HPP

#include "G2_test_indep.hpp"
#include "Contingency.hpp"

#include <vector>
#include <list>

class G2_conditional_test_indep : public G2_test_indep
{
public:
    G2_conditional_test_indep(blas_column const& genos,
                              blas_column const& phenos,
                              blas_column const& cond_genos);

    G2_conditional_test_indep(blas_column const& genos,
                              blas_column const& phenos,
                              blas_matrix const& cond_genos_vector);

    G2_conditional_test_indep(blas_column const& genos,
                              blas_column const& phenos,
                              std::list<unsigned> const& cond_genos_indexes, bool print_contingency=false);

    G2_conditional_test_indep(Contingency const& c, unsigned number_of_sub_contigencies);

    void run(bool verbose=false);
    void print_contingencies();

protected:
    std::vector<Contingency> _contingencies;
};

#endif // G2_CONDITIONAL_TEST_INDEP_HPP
