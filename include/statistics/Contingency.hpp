#ifndef CONTINGENCY_HPP
#define CONTINGENCY_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "common.h"

namespace blas=boost::numeric::ublas;

//typedef blas::matrix<int> blas_matrix;
//typedef blas::matrix<double> blas_dmatrix;
//typedef blas::matrix_column<blas::matrix<int> > blas_column;
//typedef blas::matrix_row<blas::matrix<int> > blas_row;

class Contingency : public blas_dmatrix
{
public:
    Contingency();
    Contingency(int a, int b);
    Contingency(blas_matrix const& m);
    Contingency(Contingency const& m);
    Contingency(blas_column const& phenos,
                blas_matrix const& genos,
                unsigned index);
    Contingency(blas_column const& phenos, blas_column const& genos);
    Contingency(blas_vector const& phenos, blas_column const& genos);

    double sum() const;
    double sum_col(int index) const;
    double sum_row(int index) const;
    double average() const;
};

#endif // CONTINGENCY_HPP
