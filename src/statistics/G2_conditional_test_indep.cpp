#include "G2_conditional_test_indep.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

using namespace std;

//----------------------------------------------------
// G2_conditional_test_indep : Constructor 1
//----------------------------------------------------

// Genotypes are coded 0,1,2. Phenotypes are coded 0,1

//--------- Deprecated -----------//
G2_conditional_test_indep::G2_conditional_test_indep(blas_column const& genos,
                                                     blas_column const& phenos,
                                                     blas_column const& cond_genos)
{
    G2_test_indep();
    _df = 6; // 2*1*3 -> N_Class( (geno-1) * (pheno-1) * cond_geno )
    _contingencies = vector<Contingency>(3);
    unsigned n_obs = cond_genos.size();
    for(unsigned i=0; i<n_obs; ++i)
    {
        // Put the current observation in the correct contingency table
        unsigned contingency_index = cond_genos(i);
        Contingency& c = _contingencies[contingency_index];
        unsigned cr = phenos(i); // index of contingency row
        unsigned cc = genos(i);  // index of contingency column
        c(cr, cc) += 1;
    }
    run();
}

//----------------------------------------------------
// G2_conditional_test_indep : Constructor 2
//----------------------------------------------------
//--------- Deprecated -----------//
G2_conditional_test_indep::G2_conditional_test_indep(blas_column const& genos,
                                                     blas_column const& phenos,
                                                     blas_matrix const& v_cond_genos)
{
    G2_test_indep();
    if(v_cond_genos.size2() != 0)
    {
        unsigned n_obs = v_cond_genos.size1();
        unsigned n_cond_genos = v_cond_genos.size2();
        unsigned n_contingencies = pow(3, n_cond_genos);
//        cout << "There are " << n_contingencies << " contingency tables" << endl;
        _df = 2*3*n_cond_genos;
        _contingencies = vector<Contingency>(n_contingencies);

        for(unsigned i=0; i<n_obs; ++i)
        {
            // Put the current observation in the correct contingency table
            unsigned contingency_index = 0;
            for(unsigned j=0; j<n_cond_genos; ++j)
                contingency_index += pow(3, j) * v_cond_genos(i,j);
            Contingency& c = _contingencies[contingency_index];
            unsigned cr = phenos(i);
            unsigned cc = genos(i);
            c(cr, cc) += 1;
        }
        run();
    }
    else
        cerr << "Error::G2_conditional_test_indep: empty conditional set of variables." << endl;
}
//-------- END Deprecated --------//


//----------------------------------------------------
// G2_conditional_test_indep : Constructor 3
//----------------------------------------------------
G2_conditional_test_indep::G2_conditional_test_indep(blas_column const& genos,
                                                     blas_column const& phenos,
                                                     std::list<unsigned> const& cond_genos_indexes,
                                                     bool do_print_contingency)
{
//    cout << "METHOD G2_conditional_test_indep" << endl;
    unsigned n_obs = genos.size();
    unsigned n_cond_genos = cond_genos_indexes.size();
    unsigned n_contingencies = pow(3, n_cond_genos);
//    cout << "There are " << n_contingencies << " contingency tables" << endl;
    _df = 2;
    if(n_cond_genos != 0)
        _df *= 3*n_cond_genos;

    _contingencies = vector<Contingency>(n_contingencies);

    // Fill contingency table (one or multiple)
    if(!cond_genos_indexes.empty())
    {
        blas::matrix_reference<blas_matrix> ref_genos_matrix = genos.data(); // get matrix from a column
        for(unsigned i=0; i<n_obs; ++i)
        {
            // Put the current observation in the correct contingency table
            unsigned contingency_index = 0;
            unsigned j=0;
            for(list<unsigned>::const_iterator it=cond_genos_indexes.begin(); it!=cond_genos_indexes.end(); ++it, ++j)
                contingency_index += pow(3, j) * ref_genos_matrix(i, *it);
            Contingency& c = _contingencies[contingency_index];
            unsigned cr = phenos(i);
            unsigned cc = genos(i);
            c(cr, cc) += 1;
        }
    }
    else
    {
        for(unsigned i=0; i<n_obs; ++i)
        {
            Contingency& c = _contingencies[0];
            unsigned cr = phenos(i);
            unsigned cc = genos(i);
            c(cr, cc) += 1;
        }
    }

    run(do_print_contingency);
//    cout << "METHOD G2_conditional_test_indep finished" << endl;
}

//----------------------------------------------------
// G2_conditional_test_indep : Constructor 4
//----------------------------------------------------
G2_conditional_test_indep::G2_conditional_test_indep(Contingency const& c, unsigned number_of_sub_contigencies)
{
    unsigned nrows = c.size1();
    unsigned n_cols_by_contingency = c.size2() / number_of_sub_contigencies;
    cout << "n_cols_by_contingency = " << n_cols_by_contingency << endl;

    _contingencies = vector<Contingency>(number_of_sub_contigencies);
    for(unsigned i=0; i<number_of_sub_contigencies; i++)
    {
        Contingency & current_c = _contingencies[i];
        current_c.resize(nrows, n_cols_by_contingency);
        for(unsigned j=0; j<nrows; j++)
        {
            for(unsigned k=i*n_cols_by_contingency, l=0; k<(i+1)*n_cols_by_contingency;k++,l++)
            {
                current_c(j,l) = c(j,k);
            }
        }
    }
    _df = 3;
    run(true);
}

//-----------------------------------------
// G2_conditional_test_indep : run
//-----------------------------------------
void G2_conditional_test_indep::run(bool verbose)
{
    _g2 = 0;
    for(unsigned i=0; i<_contingencies.size(); ++i)
    {
        G2_test_indep g2(_contingencies[i]);
        if(!g2.is_reliable())
        {
            _g2 = 0;
            _pval = 1;
            _reliable = false;
            break;
        }
        _reliable = true;
        _g2 += g2.g2();
    }
    if(verbose)
    {
        cout << "g2_stat\t" << _g2 << "\ndf\t" << _df << endl;
        print_contingencies();
    }
    boost::math::chi_squared_distribution<double> chi2_dist(_df);
    _pval = 1 - boost::math::cdf(chi2_dist, _g2);
    if(_pval == 0)
        _pval = 2.0e-16;
}

//-----------------------------------------
// G2_conditional_test_indep : print_contingencies
//-----------------------------------------
void G2_conditional_test_indep::print_contingencies()
{
    cout << "Contingencies:" << endl;
    for(unsigned i=0; i<_contingencies.size(); i++)
    {
        Contingency& c = _contingencies[i];
        cout << c << endl;
    }
    cout << "end" << endl;
}
