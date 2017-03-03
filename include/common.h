#ifndef COMMON_H
#define COMMON_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp> // operations on ublas::matrix (column & rows...)
#include <boost/numeric/ublas/io.hpp>           // Allow to call: cout << ublas::matrix

#include <queue>
#include <vector>
#include "Scoring.hpp"

namespace blas=boost::numeric::ublas;

typedef blas::matrix<int> blas_matrix;
typedef blas::matrix<double> blas_dmatrix;
typedef blas::matrix_column<blas::matrix<int> > blas_column;
typedef blas::matrix_row<blas::matrix<int> > blas_row;

typedef blas::vector<int> blas_vector;
typedef blas::vector<unsigned> blas_uvector;
typedef blas::vector<double> blas_dvector;

typedef std::priority_queue<Scoring, std::vector<Scoring>, std::greater<Scoring> > scoring_queue;

#endif // COMMON_H
