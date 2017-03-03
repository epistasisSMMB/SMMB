#include "Permutations.hpp"

#include <iostream>
#include <algorithm>

using namespace std;

Permutations::Permutations(blas_column const& var, blas_column const& phenos, unsigned n) : _var(var), _original_phenos(phenos), _n(n) {}

Permutations::Permutations(blas_column const& var, blas_column const& phenos, double alpha, double c): _var(var), _original_phenos(phenos)
{
    choose_b(alpha,c);
//    cout << "number of permutations to do: " << _n << endl;
}

void Permutations::run()
{
    cout << "Permutations." << endl;
    for(unsigned i=0; i<_n; ++i)
    {
        blas_vector permuted_phenos(_original_phenos.size());
        // copy column in p -> p can be modified without modifying _original_phenos
        copy(_original_phenos.begin(), _original_phenos.end(), permuted_phenos.begin());
        random_shuffle(permuted_phenos.begin(), permuted_phenos.end());

        G2_test_indep g2(_var, permuted_phenos);
        if(g2.is_reliable())
            _pvalues_dist.push(g2.pval());
    }
}

double Permutations::correction(double pval)
{
    _pvalues_dist.push(pval);
    int i=0;
    while(pval != _pvalues_dist.top())
    {
        i++;
        _pvalues_dist.pop();
    }
    return (double)(i+1)/(_n+1);
}

void Permutations::choose_b(double alpha, double c)
{
    double error = alpha * c;
    int b = alpha*(1 - alpha) / (error*error);
    _n = b;
}

unsigned Permutations::get_n() const
{
    return _n;
}
