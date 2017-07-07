#ifndef SMMB_HPP
#define SMMB_HPP

#include <list>
#include <map>
#include <string>
#include <random>
#include <ctime>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>

#include "common.h"
#include "Parameters_file_parsing.hpp"

class Smmb
{
public:
    Smmb(boost::numeric::ublas::matrix<int> & genos, blas::matrix_column<blas::matrix<int> > & phenos, Parameters_file_parsing const& params);
    ~Smmb();

    void run();
    void run(std::list<unsigned> & snp_indexes);

    void learn_mb(std::list<unsigned> & mb, std::list<unsigned> & geno_indexes);
    void forward(std::list<unsigned> & mb, std::list<unsigned> & geno_indexes);
    void backward(std::list<unsigned> & mb, std::list<unsigned> & geno_indexes);



    void print_mbs() const;
    void make_consensus(blas_matrix & permuted_phenos);

    void list_map_keys(std::map<std::list<unsigned>, double> const& m) const;
    void list_map_keys_signif(std::map<std::list<unsigned>, double> const& m) const;

    void write_result(double duration);

protected:
    boost::numeric::ublas::matrix<int> & _genotypes;
    blas::matrix_column<blas::matrix<int> > & _phenotypes;
    Parameters_file_parsing _params;

    /* Parameters imported in _params attribute :
    const unsigned _var_subset_size;            // sqrt(number of SNPs), by default
    const unsigned _k_to_draw;                  // 3, by default
    const unsigned _n_mbs;                      // number of Markov blankets to learn, before to make a consensus. 100 by default
    const unsigned _n_trials_to_learn_mbs;      // exit condition if not enough MB are learned in n tries
    const unsigned _n_trials_to_learn_1_mb;     // exit condition if a MB is not learned from a subset K after multiple iterations
    std::string _output_id;                     // id dependent from input genotype file name
    double _precision;
    double _alpha;                              // type I error rate
    */
    std::list<std::list<unsigned> > _mbs;       // list of all learned MBs
    std::list<unsigned> _consensus;             // consensus of all learned MBs
    std::mt19937 _rng;                          // randomness generator
    std::ofstream _results_handler;             // output file handler
    std::ofstream _trace_handler;               // trace file handler
    unsigned _n_tests;
    std::map<std::list<unsigned>, double> _consensus_recorded_tests;
    unsigned _smmb_exec_counter;                // count number of executions of smmb

};

#endif // SMMB_HPP
