#include "Smmb.hpp"
#include "Services.hpp"
#include "G2_conditional_test_indep.hpp"
#include "Mutual_information.hpp"
#include "Permutations_adapt.hpp"

#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <array>
#include <map>
#include <fstream>
#include <ctime>
#include <omp.h>
#include <libgen.h>

#include <boost/numeric/ublas/vector.hpp>

using namespace std;

//-----------------------------------------
// Smmb : Constructor
//-----------------------------------------
Smmb::Smmb(boost::numeric::ublas::matrix<int> & genos, blas::matrix_column<blas::matrix<int> > & phenos, Parameters_file_parsing const& params)
     :  _genotypes(genos), _phenotypes(phenos), _params(params)
{
    // seed for randomness
    _rng.seed(std::time(0));

    // Management of output files (traces, results).
    // File handlers will be open until SMMB destructor.
    string file_basename = basename((char*)_params.genos_file.c_str());
    string result_filename = "outputs/RESULT_" + file_basename;
    _results_handler.open(result_filename.c_str(), ios::trunc);

    if(!_results_handler)
    {
        std::cerr << "Error while opening output.txt (by writing access) !\n";
        exit(-1);
    }

    // Counter for the total number of tests of independance during a SMMB execution.
    _n_tests = 0;

    _smmb_exec_counter = params.n_smmb_runs;
}

//-----------------------------------------
// Smmb : Destructor
//-----------------------------------------
Smmb::~Smmb()
{
    // Close file handlers
    _results_handler.close();
}

//-----------------------------------------
// Smmb : run
//-----------------------------------------
void Smmb::run()
{
//    cout << "Running SMMB" << endl;
    unsigned nsnps = _genotypes.size2();

//    map<list<unsigned>, double> recorded_tests;

    // Learning multiple MBs
    // Stop when 100 (i.e. _params.n_mbs) MBs are learned
    // Each iteration this loop is guaranteed to learn a non-empty MB (mandatory for omp parallelization)
    #pragma omp parallel for
    for(unsigned N_MBs=0 ; N_MBs < _params.n_mbs ; N_MBs++)
    {
//        cout << "########### MB " << N_MBs << " #############\n";
        bool we_cant_go_out = true;
        do  // loop to re-sample K snps when no convergence is reached with the previous sampling
        {
            unsigned cpt = 0; // counter of trials to learn a MB (condition to break the sample/draw again in the current iteration)

            // draw a subset of snps (size = params.subset_size_large)
            // indexes of subset of size K will be stored here (default, K = sqrt(Nsnps))
            list<unsigned> snps_subset;
            Services::draw_without_replacement(0, nsnps-1, snps_subset, _params.subset_size_large, _rng);    // snps_subset given by access (reference)

//            cout << "\n\tDraw_without_replacement\n\t";
//            Services::print_list(snps_subset);

            list<unsigned> mb;     // mb to learn (empty for now)
            do  // we stay in the current high-level iteration while we learn empty MBs. We also break this loop if non-convergence is detected.
            {
                learn_mb(mb, snps_subset);
//                cout << "METHOD learn_mb finished\n";
//                cout << "mb learned: "; Services::print_list(mb);
                cpt ++;
            } while (mb.empty() &&  cpt < _params.n_trials_to_learn_mbs);

            if(!mb.empty())
            {
                _mbs.push_back(mb);
                we_cant_go_out = false;
            }
        } while (we_cant_go_out);
    }
}

//-----------------------------------------
// Smmb : run
//-----------------------------------------
void Smmb::run(list<unsigned> & snp_indexes)
{
//    cout << "Running SMMB" << endl;
    unsigned nsnps = snp_indexes.size();
    vector<unsigned> snp_indexesV;
    for(list<unsigned>::const_iterator it = snp_indexes.begin(); it != snp_indexes.end(); ++it)
        snp_indexesV.push_back(*it);

//    map<list<unsigned>, double> recorded_tests;

    // Learning multiple MBs
    // Stop when 100 (i.e. _params.n_mbs) MBs are learned
    // Each iteration this loop is guaranteed to learn a non-empty MB (mandatory for omp parallelization)
    #pragma omp parallel for
    for(unsigned N_MBs=0 ; N_MBs < _params.n_mbs ; N_MBs++)
    {
        bool we_cant_go_out = true;
        do  // loop to re-sample K snps when no convergence is reached with the previous sampling
        {
            unsigned cpt = 0; // counter of trials to learn a MB (condition to break the sample/draw again in the current iteration)

            // draw a subset of snps (size = params.subset_size_large)
            // indexes of subset of size K will be stored here (default, K = sqrt(Nsnps))
            list<unsigned> position_in_snp_indexes;
            Services::draw_without_replacement(0, nsnps-1, position_in_snp_indexes, _params.subset_size_large, _rng);    // snps_subset given by access (reference)
            list<unsigned> snps_subset;

//            cout << "snps_subset before";
            for(list<unsigned>::const_iterator it=position_in_snp_indexes.begin(); it != position_in_snp_indexes.end() ;it++)
                snps_subset.push_back(snp_indexesV.at(*it));

            list<unsigned> mb;     // mb to learn (empty for now)
            do  // we stay in the current high-level iteration while we learn empty MBs. We also break this loop if non-convergence is detected.
            {
                learn_mb(mb, snps_subset);
                cpt ++;
            } while (mb.empty() &&  cpt < _params.n_trials_to_learn_mbs);

            if(!mb.empty())
            {
//                cout << "Resulting non-empty MB :"; Services::print_list(mb);
                _mbs.push_back(mb);
                we_cant_go_out = false;
            }
        } while (we_cant_go_out);
    }
}

//-----------------------------------------
// Smmb : learn_mb
//-----------------------------------------
void Smmb::learn_mb(list<unsigned> & mb, list<unsigned> & geno_indexes)
{
//    cout << "\nMETHOD learn_mb\n";
    unsigned cpt2 = 0;  // iteration counter for 1 MB learning
    list<unsigned> mem_mb; // useful to know if any changes happen during current iteration
    do  // learning of 1 MB
    {
        mem_mb = mb;
        forward(mb, geno_indexes);
        backward(mb, geno_indexes);
        cpt2++;
    }while ((cpt2 < _params.n_trials_to_learn_1_mb && mb.empty()) && (mem_mb != mb));
    backward(mb, geno_indexes);
}

//-----------------------------------------
// Smmb : forward
//-----------------------------------------
void Smmb::forward(list<unsigned> & mb, list<unsigned> & geno_indexes)
{
//    list<unsigned> snps_drawn;                                          // indexes of subset of size k (by default, k=3)
    vector<unsigned> snps_drawnV;
    Services::random_subset(geno_indexes, snps_drawnV, _params.subset_size_small, _rng); // snps_drawn is given by access (reference)
    std::sort(snps_drawnV.begin(), snps_drawnV.end());                                                  // set sorted -> all subsets (powerset) will also be sorted.

    // A vector containing all subsets (including empty set) of snps_drawn
//    vector<list<unsigned> > allsubs;
    vector<vector<unsigned> > allsubsV;
    Services::powerset(snps_drawnV, allsubsV);

//    Services::print_powerset(allsubsV);

    // Analysis of the powerset of snps_drawn
    unsigned best_subset_index = 100;
    double  best_subset_pvalue = 1.1;

    for(unsigned i=0; i<allsubsV.size(); ++i) // Start at index 1 because of [0] --> empty list (due to powerset output)
    {
        vector<unsigned> const& current_combinV = allsubsV[i];   // reference to current subset
        list<unsigned> current_combin;                         // make a list from this vector
        for(unsigned const& i: current_combinV)
            current_combin.push_back(i);

        /*  Explanation of below "for block"
         *  current_combination = {snp1, snp2, snp3}
         *  mb = {snp10}
         *  current_geno = snp1
         *  tmp_mb = {snp10, snp2, snp3}
         *  Indep_test (snp1, pheno | {snp10,snp2,snp3})    ===    Indep_test (current_geno, pheno | tmp_mb)
        */
        unsigned it_index = 0;
        for(list<unsigned>::const_iterator it = current_combin.begin(); it != current_combin.end(); it++)
        {
            unsigned current_geno = *it;
            list<unsigned> current_combin_trunc = current_combin;

            list<unsigned>::iterator iteratr = current_combin_trunc.begin();
            for(unsigned i=0; i<it_index; i++)
                iteratr ++;
            current_combin_trunc.erase(iteratr);

            // temporary merge mb and subset
            list<unsigned> tmp_mb = mb;
            Services::append_list_by_copy(tmp_mb, current_combin_trunc);

            // Conditional independance test
//            cout << "Ind_test(" << current_geno << ", pheno | " ; Services::print_list(tmp_mb,false); cout << ")\n";
            G2_conditional_test_indep cond_g2(blas_column(_genotypes, current_geno), _phenotypes, tmp_mb);
            _n_tests ++;

//            cout << "   reliable: " << cond_g2.is_reliable() << endl;
//            cout << "   stat    : " << cond_g2.g2() << endl;
//            cout << "   p_val   : " << cond_g2.pval() << endl;

            // A G2 independency test is reliable when there are enough observations in each cell of the contingency table.
            // See the method Contingency::is_reliable() for details.
            if(cond_g2.is_reliable())
            {
                if(cond_g2.pval() < best_subset_pvalue)
                {
                    best_subset_pvalue = cond_g2.pval();
                    best_subset_index = i;
                }
            }
            it_index ++;
        }
    }
    // Append the best subset if alpha < threshold
    if(best_subset_pvalue < _params.alpha)
    {
        vector<unsigned> best_subset = allsubsV[best_subset_index];  //copy
        Services::append_vector_to_list(mb, best_subset);
//        cout << "ADD : [ "; for(auto const& elem: best_subset) cout << elem << ", "; cout << "] => " << best_subset_pvalue << endl;
        for(vector<unsigned>::iterator it=best_subset.begin(); it != best_subset.end(); ++it)
            geno_indexes.remove(*it);
    }
}

//-----------------------------------------
// Smmb : backward
//-----------------------------------------
void Smmb::backward(list<unsigned> & mb, list<unsigned> & geno_indexes)
{
    list<unsigned>::iterator it = mb.begin();
    for(unsigned i=0; i<mb.size(); i++)
    {
        unsigned current = *it;
        vector<unsigned> mb_troncated;
        for(unsigned const& i: mb)
        {
            if(i != current)
                mb_troncated.push_back(i);
        }
//        mb_troncated.remove(current); // was useful when mb_troncated was a list<unsigned> (deprecated)

        bool is_erased = false;
        vector<vector<unsigned> > allsubsets;
        Services::powerset(mb_troncated, allsubsets);
        for(unsigned j=0; j<allsubsets.size(); ++j)
        {
            list<unsigned> current_subset;
            for(unsigned const& elem: allsubsets[j])
                current_subset.push_back(elem);
            G2_conditional_test_indep cond_g2(blas_column(_genotypes, current), _phenotypes, current_subset);
            _n_tests++;

            if(cond_g2.pval() > _params.alpha) // non-reliable tests have assigned a p-value equal to 1
            {
//                cout << "Backward removal !\n";
                geno_indexes.push_front(*it);
                it = mb.erase(it);  // erase the current SNP from MB. Implicit iterator incrementation when erasing
                is_erased = true;
                break;              // if SNP is removed, we don't test it against other subsets of the Markov blanket
            }
        }
        if(!is_erased)
            ++it;
    }
}

//-----------------------------------------
// Smmb : make_consensus
//-----------------------------------------
void Smmb::make_consensus(blas_matrix & permuted_phenos)
{
//    cout << "Running SMMB.make_consensus()" << endl;
    _consensus.clear();

    // Count occurences of each Markov blanket and make them unique (mapping MB_list => number_of_occurences)
    map<list<unsigned>, unsigned> mb_occurence;
    for(list<list<unsigned> >::iterator it=_mbs.begin(); it!=_mbs.end(); ++it)
    {
        list<unsigned>& current_mb = *it;
        current_mb.sort();
        auto mb_it = mb_occurence.find(current_mb);
        if(mb_it == mb_occurence.end())
            mb_occurence[current_mb] = 1;
        else
            mb_occurence[current_mb] += 1;
    }

    // Make a consensus (list object) of unique SNP ids
    for(map<list<unsigned>, unsigned>::const_iterator it = mb_occurence.begin(); it != mb_occurence.end(); ++it)
    {
        list<unsigned> const& li = it->first;
        Services::append_list_by_copy(_consensus, li);
    }

    // No duplicate SNP in the consensus
    _consensus.sort();
    _consensus.unique();

//    print_mbs();
    cout << "MBs:\n";
    for(map<list<unsigned>, unsigned>::const_iterator it = mb_occurence.begin(); it != mb_occurence.end(); ++it)
    {
        list<unsigned> const& li = it->first;
        cout << "{ ";
        for (list<unsigned>::const_iterator l_it=li.begin(); l_it != li.end(); ++l_it)
            cout << *l_it << " ";
        cout << "} => " << it->second << "\n" ;
    }

    if(_smmb_exec_counter > 1)
    {
        cout << "Consensus at 1st SMMB exec: "; Services::print_list(_consensus);
        cout << "Size of consensus : " << _consensus.size() << endl;

        _smmb_exec_counter--;
        _mbs.clear();
        run(_consensus);
        make_consensus(permuted_phenos);
    }
    else
    {
//        cout << "Corrections of p-values with permutations : this may take several minutes." << endl;
        int r_permut_adapt = Permutations_adapt::choose_r(_params.alpha, _params.precision);
        int n_permut_adapt = permuted_phenos.size2();

        // BACKWARD PHASE
        list<unsigned> mem_consensus;
        do
        {
            mem_consensus = _consensus;
            list<unsigned> removed_from_consensus;

            #pragma omp parallel for shared(removed_from_consensus)
            for(unsigned i=0; i<_consensus.size(); i++)
            {
//                cout << "Backward analysis of " << i << endl;
                list<unsigned> updated_consensus = _consensus;
                Services::remove_listA_from_listB(removed_from_consensus, updated_consensus);

                list<unsigned>::iterator it = _consensus.begin();
                advance(it, i); // access with iterator to element to test for deletion
                unsigned current = *it;
                updated_consensus.remove(current);

                vector<unsigned> updated_consensusV;
                for(unsigned const& elem: updated_consensus)
                    updated_consensusV.push_back(elem);

                vector<vector<unsigned> > allsubsets;
                Services::powerset(updated_consensusV, allsubsets);

                for(unsigned j=0; j<allsubsets.size(); ++j)
                {
                    list<unsigned> current_subset;
                    for(unsigned const& elem: allsubsets[j])
                        current_subset.push_back(elem);

//                    cout << "Ind_test(" << current << ", pheno | " ; Services::print_list(current_subset,false); cout << ")\n";
                    G2_conditional_test_indep cond_g2(blas_column(_genotypes, current), _phenotypes, current_subset);
//                    cout << "   reliable: " << cond_g2.is_reliable() << endl;
//                    cout << "   stat    : " << cond_g2.g2() << endl;
//                    cout << "   p_val   : " << cond_g2.pval() << endl;

                    if(cond_g2.is_reliable())
                    {
                        Permutations_adapt p_adapt(blas_column(_genotypes,current), _phenotypes, permuted_phenos, r_permut_adapt, n_permut_adapt, "g2");
                        p_adapt.run(current_subset);
                        cond_g2.set_pval(p_adapt.correction());
//                        cout << "Permutations corrected p-value: " << cond_g2.pval() << endl;
                    }

                    list<unsigned> keymap; keymap.push_back(*it);
                    Services::append_list_by_copy(keymap, current_subset);
                    _consensus_recorded_tests[keymap] = cond_g2.pval();

                    #pragma omp atomic
                    _n_tests ++;

                    if(cond_g2.is_reliable() && cond_g2.pval() > _params.alpha)
                    {
                        #pragma omp critical
                        removed_from_consensus.push_front(current);
                        break;
                    }
                }
            }
            Services::remove_listA_from_listB(removed_from_consensus, _consensus);

        } while (_consensus != mem_consensus);

        cout << "Consensus after backward : "; Services::print_list(_consensus);
    }
}

//-----------------------------------------
// Smmb : write_result
//-----------------------------------------
void Smmb::write_result(double duration)
{
    _results_handler << "Duration\t" << duration << " milliseconds" << endl;
    _results_handler << "Input_file\t" << _params.genos_file << endl;
    _results_handler << "N_statistical_tests\t" << _n_tests << endl;
    _results_handler << "Markov_blanket consensus\t";
    list<unsigned>::const_iterator before_end=_consensus.end(); --before_end;
    _results_handler << "{ ";
    for(list<unsigned>::const_iterator it=_consensus.begin(); it!=before_end; ++it)
        _results_handler << *it << " ";
    _results_handler << _consensus.back() << " }" << endl;

    _results_handler << endl << "### Occurrences of weak learner MBs ####" << endl;
    map<list<unsigned>, unsigned> mb_occurence;
    for(list<list<unsigned> >::iterator it=_mbs.begin(); it!=_mbs.end(); ++it)
    {
        list<unsigned>& current_mb = *it;
        current_mb.sort();
        auto mb_it = mb_occurence.find(current_mb);
        if(mb_it == mb_occurence.end())
            mb_occurence[current_mb] = 1;
        else
            mb_occurence[current_mb] += 1;
    }

    unsigned total_number_mbs = 0;
    for(map<list<unsigned>, unsigned>::const_iterator it = mb_occurence.begin(); it != mb_occurence.end(); ++it)
    {
        list<unsigned> const& li = it->first;
        _results_handler << "{ ";
        for (list<unsigned>::const_iterator l_it=li.begin(); l_it != li.end(); ++l_it)
            _results_handler << *l_it << " ";
        _results_handler << "} => " << it->second << "\n" ;
        total_number_mbs += it->second;
    }
    _results_handler << endl;
    _results_handler << "Total number of MBs learned : " << total_number_mbs << endl;

    _results_handler << endl << "### p_values ###" << endl;
    for(map<list<unsigned>, double>::const_iterator it = _consensus_recorded_tests.begin(); it != _consensus_recorded_tests.end(); ++it)
    {
        if(it->second <= _params.alpha)
        {
            list<unsigned> const& current_combin = it->first;
            list<unsigned>::const_iterator l_it=current_combin.begin();
            _results_handler << "Indep( " << *l_it << ", phenotype ";
            if(current_combin.size() > 1)
            {
                _results_handler << "| ";
                advance(l_it, 1);
                for (; l_it != current_combin.end(); ++l_it)
                    _results_handler << *l_it << " ";
                _results_handler << ") => " << it->second << "\n" ;
            }
            else
            {
                _results_handler << ")\n";
            }

        }
    }
}

//-----------------------------------------
// Smmb : print_mbs
//-----------------------------------------
void Smmb::print_mbs() const
{
    cout << _mbs.size() << " MBs:\n";
    cout << "{\n";
    list<list<unsigned> >::const_iterator it;
    for(it = _mbs.begin(); it != _mbs.end(); ++it)
    {
//        list<unsigned>& current = *it;
        cout << "\t"; Services::print_list(*it); cout << endl;
    }
    cout << "}\n";
}

//-----------------------------------------
// Smmb : list_map_keys
//-----------------------------------------
void Smmb::list_map_keys(map<list<unsigned>, double> const& m) const
{
    for(map<list<unsigned>, double>::const_iterator it = m.begin(); it != m.end(); ++it)
    {
        list<unsigned> const& li = it->first;
        Services::print_list(li);
    }
}

//-----------------------------------------
// Smmb : list_map_keys_signif
//-----------------------------------------
void Smmb::list_map_keys_signif(map<list<unsigned>, double> const& m) const
{
    for(map<list<unsigned>, double>::const_iterator it = m.begin(); it != m.end(); ++it)
    {
        if(it->second < _params.alpha)
        {
            list<unsigned> const& li = it->first;
            Services::print_list(li); cout << " => " << it->second << endl;
	}
    }	
}
