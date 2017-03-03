#ifndef SERVICES_HPP
#define SERVICES_HPP

#include <string>
#include <vector>
#include <list>
#include <random>

#include "common.h"

class Services
{
public:
    Services();

    static void frame_cout(std::string);
    static void print_line();
    static void print_list(std::list<unsigned> const& l, bool carriage_return=true);

//    static void draw_without_replacement(unsigned begin, unsigned end, std::list<unsigned> & result, unsigned n_to_draw);     // DEPRECATED
    static void draw_without_replacement(unsigned begin, unsigned end, std::list<unsigned> & result, unsigned n_to_draw, std::mt19937 & rng);

//    static void random_subset(std::list<unsigned> const& snps_subset, std::list<unsigned> & snps_drawn, unsigned n_to_draw);  // DEPRECATED
    static void random_subset(std::list<unsigned> const& snps_subset, std::list<unsigned> & snps_drawn, unsigned n_to_draw, std::mt19937 & rng);

    // Draw one element in the input vector
    static unsigned draw_one(std::vector<unsigned> const& vec, std::mt19937 & rng);
    static int head_or_tail(std::mt19937 & rng);

    static void append_list_by_copy(std::list<unsigned> & l1, std::list<unsigned> const& l2);
    static void remove_listA_from_listB(std::list<unsigned> const& lA, std::list<unsigned> & lB);


    /* POWERSET example:
     * Powerset({1,2,3}, output_vector)
     * returns : [ {}, {1}, {2}, {3}, {1,2}, {1,3}, {2,3}, {1,2,3} ]
    */
    static void powerset(std::list<unsigned> & l, std::vector<std::list<unsigned> > & allsubs);

    static void append_file(std::string to_write, std::string filename,bool endline=false);
    static void append_file(int to_write, std::string filename, bool endline=false);
    static void append_file(double to_write, std::string filename, bool endline=false);

    // TODO expliquer pour le paramètre de nombre de permutations est absent
    static void generate_permutations(blas_column const& vector, blas_matrix & output_permutations);


};

#endif // SERVICES_HPP
