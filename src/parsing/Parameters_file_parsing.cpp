#include "Parameters_file_parsing.hpp"

#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;


//=================================================
// Constructor
//=================================================
Parameters_file_parsing::Parameters_file_parsing()
{
    ifstream file("./PARAMETERS_SMMB.txt");
    if(file)
    {
        string line;
        while (!file.eof())
        {
            getline(file, line);
            if (line.length() != 0 && line[0] != '#')
            {
                import_line(line);
            }
        }
    }
    else
    {
        std::cerr << "Error while opening PARAMETERS_SMMB.txt !\n";
    }
}

//=================================================
// Parameters_file_parsing : import_line
//=================================================
void Parameters_file_parsing::import_line(string const& line)
{
    vector<string> token = this->split(line, ' ');
    string const& key = token[0];
    string & value = token[1];

    if(key == "header")
        header = atoi(value.c_str());

    else if(key == "separator")
    {
        if(value == "\t")
            separator = '\t';
        else
            separator = value.at(0);
    }
    else if(key == "alpha")
        alpha = atof(value.c_str());

    else if(key == "precision")
        precision = atof(value.c_str());

    else if(key == "n_smmb_runs")
        n_smmb_runs = atoi(value.c_str());

    else if(key == "subset_size_large")
    {
        if(value == "sqrt")
           subset_size_large = 0;   // updated later in the main with method Parameters_file_parsing::update_subset_size_large
        else
            subset_size_large = atoi(value.c_str());
    }

    else if(key == "subset_size_small")
        subset_size_small = atoi (value.c_str());

    else if(key == "n_mbs")
        n_mbs = atoi(value.c_str());

    else if(key == "n_trials_to_learn_mbs")
        n_trials_to_learn_mbs = atoi(value.c_str());

    else if(key == "n_trials_to_learn_1_mb")
        n_trials_to_learn_1_mb = atoi(value.c_str());

    else if(key == "gfile")
        genos_file = value;

    else if(key == "pfile")
        phenos_file = value;

    else {}

}

//=================================================
// Parameters_file_parsing : split
//=================================================
vector<string> Parameters_file_parsing::split(string const& s, char delim)
{
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim))
        tokens.push_back(item);
    return tokens;
}

//=================================================
// Parameters_file_parsing : list_parameters
//=================================================
void Parameters_file_parsing::list_parameters() const
{
    cout << "########### PARAMETERS ###########\n" << "header => " << header << endl
    << "separator => " << separator << endl
    << "alpha => " << alpha << endl
    << "precision => " << precision << endl
    << "subset_size_large => " << subset_size_large << endl
    << "subset_size_small => " << subset_size_small << endl
    << "n_trials_to_learn_mbs => " << n_trials_to_learn_mbs << endl
    << "n_trials_to_learn_1_mb => " << n_trials_to_learn_1_mb << endl
    << "#################################" << endl;
}

//=================================================
// Parameters_file_parsing : update_subset_size_large
//=================================================
void Parameters_file_parsing::update_subset_size_large(unsigned const& n_genos)
{
    if(subset_size_large == 0)
        subset_size_large = sqrt(n_genos);
}
