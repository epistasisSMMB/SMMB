#include "G2_test_indep.hpp"

#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;
namespace blas = boost::numeric::ublas;

//-----------------------------------------
// G2_test_indep : Constructors
//-----------------------------------------
G2_test_indep::G2_test_indep(Contingency const& c) : _pval(1), _g2(0), _reliable(true)
{
    _df = (c.size1()-1)*(c.size2()-1);
    run(c);
}

G2_test_indep::G2_test_indep(blas_column const& var, blas_column const& phenos) : _pval(1), _g2(0), _reliable(true)
{
    Contingency c(phenos, var);

    _df = (c.size1()-1) * (c.size2()-1);
    run(c);
}

G2_test_indep::G2_test_indep(blas_column const& var, blas_vector const& phenos) : _pval(1), _g2(0), _reliable(true)
{
    Contingency c(phenos, var);

    _df = (c.size1()-1) * (c.size2()-1);
    run(c);
}

G2_test_indep::G2_test_indep(): _pval(1), _g2(0), _reliable(true), _df(0) {}

//-----------------------------------------
// G2_test_indep : run
//-----------------------------------------
void G2_test_indep::run(Contingency const& c)
{
    if(!reliable_test(c))
    {
        _reliable = false;
        return;
    }
    Contingency e(c);
    int n = c.sum();

    // Expected contingency filling
    for(unsigned i=0; i<c.size1(); ++i)
    {
        for(unsigned j=0; j<c.size2(); ++j)
            e(i,j) = (double)(c.sum_row(i) * c.sum_col(j)) / n;
    }

/*    Do not compute g2 and pval if the test is not reliable
//    if(! reliable_test(e))
//    {
//        _reliable = false;
//        return;
//    }*/
    compute_g2(c, e);

//    cout << "g2(non-cond):" << _g2 << endl;
    boost::math::chi_squared_distribution<double> chi2_dist(_df);
    _pval = 1 - boost::math::cdf(chi2_dist, _g2);
}

void G2_test_indep::set_pval(double pval)
{
    _pval = pval;
}

//-----------------------------------------
// G2_test_indep : g2
//-----------------------------------------
double G2_test_indep::g2() const
{
    return _g2;
}

//-----------------------------------------
// G2_test_indep : pval
//-----------------------------------------
double G2_test_indep::pval() const
{
    return _pval;
}

//-----------------------------------------
// G2_test_indep : is_reliable
//-----------------------------------------
bool G2_test_indep::is_reliable() const
{
    return _reliable;
}

//-----------------------------------------
// G2_test_indep : compute_g2
//-----------------------------------------
void G2_test_indep::compute_g2(Contingency const& c, Contingency const& e)
{
//    _g2 = 0;      // initialized in constructor and this function is called in constructor
    for(unsigned i=0; i<c.size1(); ++i)
    {
        for(unsigned j=0; j<c.size2(); ++j)
        {
            if(c(i,j) != 0)
            {
                double div = (double) c(i,j) / e(i,j);
                _g2 += c(i,j) * log(div);
            }
        }
    }
    _g2 *= 2;
}

//-----------------------------------------
// G2_test_indep : reliable_test
//-----------------------------------------

bool G2_test_indep::reliable_test(Contingency const& c)
{
	for(unsigned i=0; i<c.size1(); ++i)
	{
	   for(unsigned j=0; j<c.size2(); ++j)
	   {
		   if(c(i,j) == 0)
			   return false;
	   }
	}
	return true;
}


bool G2_test_indep::deprecated_reliable_test(Contingency const& e)
{
    unsigned ncells = e.size1()*e.size2();
    int count_inf_5 = 0;
    for(unsigned i=0; i<e.size1(); ++i)
    {
        for(unsigned j=0; j<e.size2(); ++j)
        {
            if(e(i,j) < 0 || e(i,j)!=e(i,j)) // test for nan
            {
                return false;
            }
            if(e(i,j) < 5)
            {
                count_inf_5 ++;
                if((double)count_inf_5 / ncells > 0.2)
                {
//                    cout << "Not reliable test, ";
//                    cout << "expected: " << e << endl;
                    return false;
                }
            }
        }
    }
    return true;
}
