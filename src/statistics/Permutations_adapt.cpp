#include "Permutations_adapt.hpp"

#include <cmath>
#include <boost/math/distributions/negative_binomial.hpp>
#include <omp.h>
#include <sstream>

#include "Mutual_information.hpp"
#include "G2_conditional_test_indep.hpp"

using namespace std;

Permutations_adapt::Permutations_adapt(blas_column const& var, blas_column const& phenos, blas_matrix & permuted_phenos, double alpha, double c, string mode)
    : Permutations(var, phenos, alpha, c), _permuted_phenos(permuted_phenos), _precision(c), _Ri(0), _Bi(0), _mode(mode)
{
    this->choose_r(alpha);
}

Permutations_adapt::Permutations_adapt(blas_column const& var,
                   blas_column const& phenos,
                   std::list<unsigned> & conditioning_set,
                   blas_matrix & permuted_phenos,
                   double alpha,
                   double c,
                   std::string mode)
    : Permutations(var, phenos, alpha, c), _permuted_phenos(permuted_phenos), _precision(c), _Ri(0), _Bi(0), _mode(mode)
{
    _conditioning_set = conditioning_set;
    this->choose_r(alpha);
}

Permutations_adapt::Permutations_adapt(blas_column const& var, blas_column const& phenos, blas_matrix & permuted_phenos, double alpha, double c)
    : Permutations(var, phenos, alpha, c), _permuted_phenos(permuted_phenos), _precision(c), _Ri(0), _Bi(0)
{
    this->choose_r(alpha);
}

Permutations_adapt::Permutations_adapt(blas_column const& var, blas_column const& phenos, blas_matrix & permuted_phenos, int r, int n, string mode)
    : Permutations(var, phenos, n), _permuted_phenos(permuted_phenos)
{
    _r = r;
    _Ri = 0;
    _Bi = 0;
    _mode = mode;
}


void Permutations_adapt::choose_r(double alpha)
{
    double error = alpha * _precision;
    int r = 0;
    bool found_r = false;
    while(found_r == false)
    {
        r += 100;
        double p1 = 0.1586553;  // P(X>1) in normal dist
        double p2 = 0.8413447;  // P(X<=1) in normal dist

        boost::math::negative_binomial_distribution<double> nbinom(r, alpha);

        double q1 = boost::math::quantile(nbinom, p1);
        double q2 = boost::math::quantile(nbinom, p2);

        double pval1 = (double)r/(r+q1);
        double pval2 = (double)r/(r+q2);

        double diff1 = abs(pval1 - alpha);
        double diff2 = abs(pval2 - alpha);

        double diff = std::max(diff1, diff2);
//        stringstream s; s << diff << " < " << error << "?\n";
//        cout << s.str();
        if(diff < error)
          found_r = true;
    }
    _r = r;
}

int Permutations_adapt::choose_r(double alpha, double precision)
{
    double error = alpha * precision;
    int r = 0;
    bool found_r = false;
    while(found_r == false)
    {
        r += 100;
        double p1 = 0.1586553;  // P(X>1) in normal dist
        double p2 = 0.8413447;  // P(X<=1) in normal dist

        boost::math::negative_binomial_distribution<double> nbinom(r, alpha);

        double q1 = boost::math::quantile(nbinom, p1);
        double q2 = boost::math::quantile(nbinom, p2);

        double pval1 = (double)r/(r+q1);
        double pval2 = (double)r/(r+q2);

        double diff1 = abs(pval1 - alpha);
        double diff2 = abs(pval2 - alpha);

        double diff = std::max(diff1, diff2);
//        stringstream s; s << diff << " < " << error << "?\n";
//        cout << s.str();
        if(diff < error)
          found_r = true;
    }
    return r;
}

int Permutations_adapt::choose_n(double alpha, double c)
{
    double error = alpha * c;
    int n = alpha*(1 - alpha) / (error*error);
    return n;
}

unsigned Permutations_adapt::get_r() const
{
    return _r;
}

void Permutations_adapt::run()
{
    unsigned permutation_counter = 0;
    if(_mode == "mi")
    {
        //compute observed mutual_information
        Mutual_information mutual_info(_var, _original_phenos);
        double obs_mi = mutual_info.mi();

//        cout << "Permutations (adaptative)." << endl;

        #pragma omp parallel for
        for(unsigned Bi=0 ; Bi<_n; Bi++)
        {
            if(_Ri < _r)
            {
                _Bi++;
                blas_column p(_permuted_phenos, permutation_counter);

                // permuted mutual information calculation
                Mutual_information mi(_var, p);
                if(obs_mi < mi.mi())
                {
                    _Ri++;
                }
                permutation_counter++;
            }
        }

//        while(_Ri < _r && _Bi < _n)
//        {
//            _Bi++;
//            // permut & test
//            blas_column p(_permuted_phenos, permutation_counter);

//            // permuted mutual information calculation
//            Mutual_information mi(_var, p);
//            if(obs_mi < mi.mi())
//                _Ri++;
//            permutation_counter++;
//        }
    }

    else if(_mode == "g2")
    {
        //compute observed t_stat
        G2_test_indep g2_obs(_var, _original_phenos);
        if(!g2_obs.is_reliable())
        {
            cout << "Permutations asked for a non reliable test on observed data. Exit" << endl;
            return;
        }
        double obs_t_stat = g2_obs.g2();

//        cout << "Permutations (adaptative)." << endl;

        /*
        #pragma omp parallel for
        for(unsigned Bi=0 ; Bi<_n; Bi++)
        {
            if(_Ri < _r)
            {
                _Bi++;
                // permuted G2 stat calculations
                G2_test_indep g2(_var, p);
                if(g2.is_reliable())
                {
                    if(obs_t_stat < g2.g2())
                        _Ri++;
                }
                permutation_counter++;
            }
        }*/

        while(_Ri < _r && _Bi < _n && permutation_counter < _permuted_phenos.size2())
        {
            _Bi++;
            // permut & test
            blas_column p(_permuted_phenos, permutation_counter);

            // permuted G2 stat calculations
            G2_test_indep g2(_var, p);
            if(g2.is_reliable())
            {
                if(obs_t_stat < g2.g2())
                    _Ri++;
            }
            permutation_counter++;
        }
    }
    else
        cout << "Permutation mode invalid. No permutations done." << endl;
}

void Permutations_adapt::run(list<unsigned> const& conditional_indexes)
{
//    stringstream s;
//    s << "We are in RUN(), Max number of permut : " << _n << endl;
//    cout << s.str();
    unsigned permutation_counter = 0;
    if(_mode == "mi")
    {
        //compute observed mutual_information
        Mutual_information mutual_info(_var, _original_phenos, conditional_indexes);
        double obs_mi = mutual_info.mi();

//        cout << "Permutations (adaptative)." << endl;
//        #pragma omp parallel for
//        for(unsigned Bi=0 ; Bi<_n; Bi++)
//        {
//            if(_Ri < _r)
//            {
//                _Bi++;
//                blas_column p(_permuted_phenos, permutation_counter);

//                // permuted mutual information calculation
//                Mutual_information mi(_var, p);
//                if(obs_mi < mi.mi())
//                    _Ri++;
//                permutation_counter++;
//            }
//        }

        while(_Ri < _r && _Bi < _n)
        {
            _Bi++;
            // permut & test
            blas_column p(_permuted_phenos, permutation_counter);

            // permuted mutual information calculation
            Mutual_information mi(_var, p, conditional_indexes);
            if(obs_mi < mi.mi())
                _Ri++;
            permutation_counter++;
        }
    }

    else if(_mode == "g2")
    {
        //compute observed t_stat
        G2_conditional_test_indep g2_obs(_var, _original_phenos, conditional_indexes);
        if(!g2_obs.is_reliable())
        {
//            cout << "Permutations asked for a non reliable test on observed data. Exit" << endl;
            return;
        }
        double obs_t_stat = g2_obs.g2();

//        stringstream s; s << permutation_counter << " vs " << _permuted_phenos.size2() << endl;
//        cout << s.str();

        while(_Ri < _r && _Bi < _n && permutation_counter < _permuted_phenos.size2())
        {
            _Bi++;
            // permut & test
            blas_column p(_permuted_phenos, permutation_counter);

            // permuted G2 stat calculations
            G2_conditional_test_indep g2(_var, p, conditional_indexes);
            if(g2.is_reliable())
            {
                if(obs_t_stat < g2.g2())
                    _Ri++;
            }
            permutation_counter++;
        }
    }

    else
        cout << "Permutation mode invalid. No permutations done." << endl;
}

double Permutations_adapt::correction()
{
    double corrected_pv=0;
    if(_Bi < _n)
    {
        corrected_pv = (double)_r/_Bi;
//        stringstream s; s << "correction : _r = " << _r << " & _Bi = " << _Bi << " -> _r/Bi -> pvalue_corr\n";
//        cout << s.str();
    }
    else
    {
        corrected_pv = (double)(_Ri+1)/(_n+1);
//        stringstream s; s << "correction : _Ri+1 = " << _Ri+1 << " & _n+1 = " << _n+1 << " -> (_Ri+1)/(_n+1) -> pvalue_corr\n";
//        cout << s.str();
    }
    return corrected_pv;
}
