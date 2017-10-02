#ifndef G2_TEST_INDEP_HPP
#define G2_TEST_INDEP_HPP

#include "Contingency.hpp"
#include <boost/numeric/ublas/matrix.hpp>

namespace blas=boost::numeric::ublas;

typedef blas::matrix<int> blas_matrix;
typedef blas::matrix_column<blas::matrix<int> > blas_column;
typedef blas::matrix_row<blas::matrix<int> > blas_row;

class G2_test_indep
{
public:
    G2_test_indep(Contingency const& contingency);
    G2_test_indep(blas_column const& var, blas_column const& phenos);
    G2_test_indep(blas_column const& var, blas_vector const& phenos);
    G2_test_indep();
    void run(Contingency const& c);
    void set_pval(double pval);
    double g2() const;
    double pval() const;
    bool is_reliable() const;



protected:

    void compute_g2(Contingency const& observeds_contingency, Contingency const& expected_contingency);
    bool reliable_test(Contingency const& contingency);
    bool deprecated_reliable_test(Contingency const& expected_contingency);

    double _pval;
    double _g2;
    bool _reliable;
    int _df;

};

#endif // G2_TEST_INDEP_HPP
