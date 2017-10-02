#include "Permutations.hpp"
#include <list>
#include <iostream>
#include <algorithm>

#include "G2_test_indep.hpp"
#include "G2_conditional_test_indep.hpp"

using namespace std;

//-----------------------------------------
// Constructors
//-----------------------------------------
Permutations::Permutations(blas_column const& var, blas_column const& phenos, unsigned n)
    : _var(var), _original_phenos(phenos), _n(n)
{}

Permutations::Permutations(blas_column const& var, blas_column const& phenos, list<unsigned> & conditioning_set, unsigned n)
    : _var(var), _original_phenos(phenos), _n(n), _conditioning_set(conditioning_set)
{}

Permutations::Permutations(blas_column const& var, blas_column const& phenos, double alpha, double c)
    : _var(var), _original_phenos(phenos)
{
    choose_b(alpha,c);
//    cout << "number of permutations to do: " << _n << endl;
}

//-----------------------------------------
// Permutations : run
//-----------------------------------------
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

//-----------------------------------------
// Permutations : run
//-----------------------------------------
void Permutations::run(blas_matrix & permuted_phenos)
{
    // size1 = number of individuals
    // size2 = number of permuted phenos
    for(unsigned i=0; i<permuted_phenos.size2(); ++i)
    {
        blas_column permuted_i(permuted_phenos,i);
        if(! _conditioning_set.empty())
        {
            G2_conditional_test_indep g2_cond(_var, permuted_i,_conditioning_set);
            if(g2_cond.is_reliable())
                _pvalues_dist.push(g2_cond.pval());
        }
    }
}

//-----------------------------------------
// Permutations : correction
//-----------------------------------------
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

//-----------------------------------------
// Permutations : choose_b
//-----------------------------------------
void Permutations::choose_b(double alpha, double c)
{
    double error = alpha * c;
    int b = alpha*(1 - alpha) / (error*error);
    _n = b;
}

//-----------------------------------------
// Permutations : get_n
//-----------------------------------------
unsigned Permutations::get_n() const
{
    return _n;
}
