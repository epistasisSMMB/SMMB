#include "Services.hpp"

#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <algorithm>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

using namespace std;

Services::Services(){}

//-----------------------------------------
// Services : frame_cout
//-----------------------------------------
void Services::frame_cout(string s)
{
    cout << endl;
    for (size_t i=0; i<s.length()+4; i++)
        cout << "#";
    cout << endl << "  " << s << endl;
    for (size_t i=0; i<s.length()+4; i++)
        cout << "#";
    cout << endl;
}

//-----------------------------------------
// Services : print_line
//-----------------------------------------
void Services::print_line()
{
    cout << "\n_________________________________________________________\n" << endl;
}

//-----------------------------------------
// Services : print_list
//-----------------------------------------
void Services::print_list(std::list<unsigned> const& l, bool carriage_return)
{
    cout << "{ ";
    for (list<unsigned>::const_iterator it=l.begin(); it != l.end(); ++it)
    {
        cout << *it << " ";
    }
    cout << "}";
    if(carriage_return)
        cout << endl;
}

//-----------------------------------------
// Services : print_powerset
//-----------------------------------------
void Services::print_powerset(vector<vector<unsigned> > const& vectors)
{
    cout << "[" << endl;
    for(unsigned i=0; i<vectors.size(); i++)
    {
        vector<unsigned> const& currentV = vectors[i];
        cout << "\t[ ";
        for(unsigned j=0; j < currentV.size(); j++)
            cout << currentV[j] << " ";
        cout << "]" << endl;
    }
    cout << "]" << endl;
}

//-----------------------------------------
// Services : draw_without_replacement DEPRECATED
//-----------------------------------------
/*void Services::draw_without_replacement(unsigned begin, unsigned end, list<unsigned> & result, unsigned n_to_draw)
{
    boost::random::mt19937 generator;
    generator.seed(static_cast<unsigned int>(std::time(0)));
    list<unsigned> allnumbs;
    for(unsigned i=begin; i<end+1; i++)
        allnumbs.push_back(i);

    while (result.size() < n_to_draw)
    {
        boost::random::uniform_int_distribution<> distribution(begin, allnumbs.size()-1);
        unsigned number = distribution(generator);
        allnumbs.remove(number);
        result.push_back(number);
//        if (find(result.begin(), result.end(), number) == result.end())
//            result.push_back(number);
    }
}*/

//-----------------------------------------
// Services : draw_without_replacement
//-----------------------------------------
void Services::draw_without_replacement(unsigned begin, unsigned end, std::list<unsigned> & result, unsigned n_to_draw, mt19937 & rng)
{
    list<unsigned> allnumbs;
    for(unsigned i=begin; i<=end; i++)
        allnumbs.push_back(i);

    while (result.size() < n_to_draw)
    {
        uniform_int_distribution<int> distribution(begin, allnumbs.size()-1);
        unsigned number = distribution(rng);

        list<unsigned>::const_iterator it = allnumbs.begin();
        advance(it, number);

        allnumbs.erase(it);
        result.push_back(*it);
//        if (find(result.begin(), result.end(), number) == result.end())
//            result.push_back(number);
    }
}

//-----------------------------------------
// Services : random_subset DEPRECATED
//-----------------------------------------
/*void Services::random_subset(list<unsigned> const& set, list<unsigned> & o_subset, unsigned n_to_draw)
{
    boost::random::mt19937 generator;
    generator.seed(static_cast<unsigned int>(std::time(0)));
    boost::random::uniform_int_distribution<> distribution(0, set.size()-1);
    list<unsigned> index_list;
    while (o_subset.size() < n_to_draw)
    {
        unsigned index = distribution(generator);
        if (find(index_list.begin(), index_list.end(), index) == index_list.end())
        {
            // deplace in list until index;
            list<unsigned>::const_iterator it = set.begin();
            advance(it, index);

            unsigned number = *it;
            o_subset.push_back(number);
            index_list.push_back(index);
        }
    }
}*/

//-----------------------------------------
// Services : random_subset
//-----------------------------------------
void Services::random_subset(list<unsigned> const& set, list<unsigned> & o_subset, unsigned n_to_draw, std::mt19937 & rng)
{
    std::uniform_int_distribution<int> distribution(0, set.size()-1);
//    boost::random::uniform_int_distribution<> distribution(0, set.size());
    list<unsigned> index_list;
//    Services::print_list(set);
    while (o_subset.size() < n_to_draw)
    {
        unsigned index = distribution(rng);
//        cout << "index " << index << "\n";
        if (find(index_list.begin(), index_list.end(), index) == index_list.end())
        {
            // deplace in list until index;
            list<unsigned>::const_iterator it = set.begin();
            advance(it, index);

            unsigned number = *it;
            o_subset.push_back(number);
            index_list.push_back(index);
//            cout << "pioche effectuee -> taille = " << o_subset.size() << endl;
        }
    }
}

// with vector as output than list (same implementation)
void Services::random_subset(list<unsigned> const& set, vector<unsigned> & o_subset, unsigned n_to_draw, std::mt19937 & rng)
{
    std::uniform_int_distribution<int> distribution(0, set.size()-1);
//    boost::random::uniform_int_distribution<> distribution(0, set.size());
    list<unsigned> index_list;
//    Services::print_list(set);
    while (o_subset.size() < n_to_draw)
    {
        unsigned index = distribution(rng);
//        cout << "index " << index << "\n";
        if (find(index_list.begin(), index_list.end(), index) == index_list.end())
        {
            // deplace in list until index;
            list<unsigned>::const_iterator it = set.begin();
            advance(it, index);

            unsigned number = *it;
            o_subset.push_back(number);
            index_list.push_back(index);
//            cout << "pioche effectuee -> taille = " << o_subset.size() << endl;
        }
    }
}

//-----------------------------------------
// Services : draw_one
//-----------------------------------------
unsigned Services::draw_one(vector<unsigned> const& vec, std::mt19937 & rng)
{
    std::uniform_int_distribution<int> distribution(0, vec.size()-1);
    unsigned index = distribution(rng);
    return vec[index];
}

//-----------------------------------------
// Services : head_or_tail
//-----------------------------------------
int Services::head_or_tail(std::mt19937 & rng)
{
    std::uniform_int_distribution<int> distribution(0, 1);
    return distribution(rng);
}

//-----------------------------------------
// Services : append_list_by_copy
//-----------------------------------------
void Services::append_list_by_copy(list<unsigned> & l1, list<unsigned> const& l2)
{
    list<unsigned> copy = l2;
    l1.insert(l1.end(), copy.begin(), copy.end());
}

//-----------------------------------------
// Services : append_vector_to_list
//-----------------------------------------
void Services::append_vector_to_list(list<unsigned> & l, vector<unsigned> const& v)
{
    for(unsigned const& i: v)
        l.push_back(i);
}

//-----------------------------------------
// Services : powerset      DEPRECATED
//-----------------------------------------
/*void Services::powerset(list<unsigned> & l, vector<list<unsigned> > & allsubs)
{
    for (int counter = 0; counter < (1 << l.size()); ++counter)
    {
        list<unsigned> combination;
        for (unsigned i=0; i < l.size(); ++i)
        {
            if (counter & (1 << i))
            {
                list<unsigned>::const_iterator it = l.begin();
                advance(it, i);
                combination.push_back(*it);
            }
        }

        allsubs.push_back(combination);
    }
}*/

//-----------------------------------------
// Services : powerset
//-----------------------------------------
void Services::powerset(vector<unsigned> & arr, vector<vector<unsigned> > & output)
{
//    cout << "POWERSET of : [ ";
//    for (unsigned const& elem: arr)
//        cout << elem << " ";
//    cout << "]\n";

    // compute all combinations for all set sizes

    unsigned max_order = 5;
    for(unsigned i=1; i<=max_order; i++)
    {
//        cout << "Powerset :: Cominations of order " << i << " ..." << endl;
        record_combin_order(arr, i, output);
    }
}

//-----------------------------------------
// Services : record_combin_order
//-----------------------------------------
// The main function that prints all combinations of size r
// in arr. This function mainly uses combinationUtil()
void Services::record_combin_order(vector<unsigned> & arr, unsigned r, vector<vector<unsigned> > & output)
{
    // A temporary array to store all combinations one by one
    vector<unsigned> data; data.reserve(r);
    unsigned n = arr.size();

    // Print all combination using temprary array 'data'
    combinations_order_recursive(arr, data, 0, n-1, 0, r, output);
}

//-----------------------------------------
// Services : combinations_order_recursive
//-----------------------------------------
/* arr  ---> Input Array
   data ---> Temporary array to store current combination
   start & end ---> Starting and Ending indexes in arr
   index  ---> Current index in data
   r ---> Size of a combination to be printed */
void Services::combinations_order_recursive(vector<unsigned> & arr, vector<unsigned> & data, unsigned start, unsigned end, unsigned index, unsigned r, vector<vector<unsigned> > & output)
{
    // Current combination is ready to be recorded, record it
    if (index == r)
    {
        vector<unsigned> v;
        for (unsigned j=0; j<r; j++)
            v.push_back(data[j]);
        output.push_back(v);

        return;
    }

    // replace index with all possible elements. The condition
    // "end-i+1 >= r-index" makes sure that including one element
    // at index will make a combination with remaining elements
    // at remaining positions
    for (unsigned i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = arr[i];
        combinations_order_recursive(arr, data, i+1, end, index+1, r, output);
    }
}

//-----------------------------------------
// Services : remove_listA_from_listB
//-----------------------------------------
void Services::remove_listA_from_listB(list<unsigned> const& lA, list<unsigned> & lB)
{
    for(list<unsigned>::const_iterator it=lA.begin(); it != lA.end(); ++it)
    {
        unsigned current = *it;
        lB.remove(current);
    }
}

//-----------------------------------------
// Services : append_file
//-----------------------------------------
void Services::append_file(string to_write, string filename, bool endline)
{
    std::ofstream file(filename.c_str(), ios_base::app);
    if(file)
    {
        file << to_write;
        if(endline)
            file << endl;
        file.close();
    }
    else
    {
        std::cerr << "Impossible d'ouvrir le fichier en écriture !\n";
    }
}

//-----------------------------------------
// Services : append_file
//-----------------------------------------
void Services::append_file(int to_write, string filename, bool endline)
{
    std::ofstream file(filename.c_str(), ios_base::app);
    if(file)
    {
        file << to_write;
        if(endline)
            file << endl;
        file.close();
    }
    else
    {
        std::cerr << "Impossible d'ouvrir le fichier en écriture !\n";
    }
}

//-----------------------------------------
// Services : append_file
//-----------------------------------------
void Services::append_file(double to_write, string filename, bool endline)
{
    std::ofstream file(filename.c_str(), ios_base::app);
    if(file)
    {
        file << to_write;
        if(endline)
            file << endl;
        file.close();
    }
    else
    {
        std::cerr << "Impossible d'ouvrir le fichier en écriture !\n";
    }
}

//-----------------------------------------
// Services : generate_permutations
//-----------------------------------------
void Services::generate_permutations(blas_column const& vector, blas_matrix & output_permutations)
{
    #pragma omp parallel for
    for(unsigned i=0; i<output_permutations.size2(); ++i)
    {
        blas_vector permuted_vec(vector.size());
        copy(vector.begin(), vector.end(), permuted_vec.begin());
        random_shuffle(permuted_vec.begin(), permuted_vec.end());

        blas_column column(output_permutations, i);
        for(unsigned j=0; j<column.size(); j++)
            column(j) = permuted_vec(j);
    }
}
