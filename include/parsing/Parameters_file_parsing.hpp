#ifndef PARAMETERS_FILE_PARSING_H
#define PARAMETERS_FILE_PARSING_H

#include <string>
#include <vector>

class Parameters_file_parsing
{
public:
    Parameters_file_parsing();
    void import_line(std::string const& line);
    void list_parameters() const;
    void update_subset_size_large(unsigned const& n_genos);

// Parameters given in the OPTIONS.txt file
// Reachable from any class that include the current header (Option_file_parsing.hpp)
    int header;
    char separator;
    double alpha;
    double precision;
    unsigned n_smmb_runs;
    unsigned subset_size_large;
    unsigned subset_size_small;
    unsigned n_mbs;
    unsigned n_trials_to_learn_mbs;
    unsigned n_trials_to_learn_1_mb;


    std::string genos_file;
    std::string phenos_file;


private:
    std::vector<std::string> split(std::string const& s, char delim);

};

#endif // PARAMETERS_FILE_PARSING_H
