#include "Mutual_information.hpp"
#include <cmath>

using namespace std;
// 2 var
Mutual_information::Mutual_information(Contingency const& c) : _mi(0.0)
{
    _total_obs = c.sum();
    run(c);
}

Mutual_information::Mutual_information(blas_column const& genos, blas_column const& phenos) : _mi(0.0)
{
    _total_obs = (double)phenos.size();
    Contingency c(phenos, genos);
    run(c);
}

Mutual_information::Mutual_information(blas_column const& genos, blas_vector const& phenos) : _mi(0.0)
{
    _total_obs = (double)phenos.size();
    Contingency c(phenos, genos);
    run(c);
}

Mutual_information::Mutual_information(blas_column const& genos,
                                       blas_column const& phenos,
                                       list<unsigned> const& cond_genos_indexes) : _mi(0.0)
{
    _total_obs = (double)genos.size();
    unsigned n_cond_genos = cond_genos_indexes.size();
    unsigned n_contingencies = pow(3, n_cond_genos);

    _contingencies = vector<Contingency>(n_contingencies);

    if(!cond_genos_indexes.empty())
    {
        blas::matrix_reference<blas_matrix> ref_genos_matrix = genos.data(); // get matrix from a column
        for(unsigned i=0; i<(unsigned)_total_obs; ++i)
        {
            // Put the current observation in the correct contingency table
            unsigned contingency_index = 0;
            unsigned j=0;
            for(list<unsigned>::const_iterator it=cond_genos_indexes.begin(); it!=cond_genos_indexes.end(); ++it, ++j)
                contingency_index += pow(2, j) * ref_genos_matrix(i, *it);
            Contingency& c = _contingencies[contingency_index];
            unsigned cr = phenos(i);
            unsigned cc = genos(i);
            c(cr, cc) += 1;
        }
    }
    else
    {
        for(unsigned i=0; i<(unsigned)_total_obs; ++i)
        {
            Contingency& c = _contingencies[0];
            unsigned cr = phenos(i);
            unsigned cc = genos(i);
            c(cr, cc) += 1;
        }
    }
    run();
}

void Mutual_information::run()
{
    double mi_summed = 0.0;
    for(unsigned i=0; i<_contingencies.size(); ++i)
    {
        _mi = 0;
        Contingency& c = _contingencies[i];
        run(c);
        double p_z = c.sum() / _total_obs;
        _mi *= p_z;
        mi_summed += _mi;
    }
    _mi = mi_summed;
}


void Mutual_information::run(Contingency const& c)
{
    for(unsigned x=0; x<c.size1(); x++)
    {
        double total_x = c.sum_row(x);
        for(unsigned y=0; y<c.size2(); y++)
	{
	    double p_xy = c(x,y) / _total_obs;
	    double p_x = total_x / _total_obs;
	    double p_y = c.sum_col(y) / _total_obs;

            if(p_xy == 0)
                _mi += 0;
            else
                _mi += p_xy * log2(p_xy / (p_x*p_y));
        }
    }
}

double Mutual_information::mi() const
{
    return _mi;
}
