#ifndef CSVParser_HPP
#define CSVParser_HPP

#include "CSVRow.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>


/*!
 * \class CSVParser
 * \brief Parses a CSV file and import data in a Matrix. Templated.
 */
template <class T>
class CSVParser
{

public:
    /*!
     * \brief Constructor
     * \param[in] i_filename Path to the CSV file to parse
     */
    CSVParser(const std::string & i_filename, char sep=',', const unsigned int header_nrows=0, bool transpose=false);

    /*!
     * \brief Get the imported matrix
     * \return Imported Matrix, templated
     */
    boost::numeric::ublas::matrix<T> &data();

    /*!
     * \brief Count rows in the CSV file
     * \return Row count
     */
    unsigned int count_rows();

    /*!
     * \brief Count cols in the CSV file
     * \return Column count
     */
    unsigned int count_cols();

    /*!
     * \brief Import the CSV data in _matrix (Matrix<T> class)
     */
    void read();

    /*!
     * \brief Transpose and import the CSV data in _matrix (Matrix<T> class)
     */
    void read_transpose();


private:
    std::string _filename;
    char _field_separator;
    unsigned int _header_nrows;
    unsigned int _nrows;
    unsigned int _ncols;
    boost::numeric::ublas::matrix<T> _matrix;
};


/* ###################################################################### IMPLEMENTATION below this point ###################################################################### */

//==============================================================
// Constructors
//==============================================================
template <class T>
CSVParser<T>::CSVParser(const std::string &i_filename, char sep, const unsigned int header_nrows, bool transpose) : _filename(i_filename), _field_separator(sep), _header_nrows(header_nrows)
{
    _nrows = count_rows() - _header_nrows;
    _ncols = count_cols();

    if(!transpose)
        read();
    else
        read_transpose();
}

//==============================================================
// CSVParser : data
//==============================================================
template <class T>
boost::numeric::ublas::matrix<T>& CSVParser<T>::data()
{
    return _matrix;
}

//==============================================================
// CSVParser : count_rows
//==============================================================
template <class T>
unsigned int CSVParser<T>::count_rows()
{
    unsigned int nrows=0;
    std::string buffer;
    std::ifstream file(_filename.c_str(), std::ifstream::in);

    if(file)
    {
        while(getline(file, buffer))
            nrows++;
    }
    else
        std::cout << "CSVParser<T>::count_rows() => Erreur lecture fichier\n";

    return nrows;
}

//==============================================================
// CSVParser : count_cols
//==============================================================
template <class T>
unsigned int CSVParser<T>::count_cols()
{
    unsigned int ncols=0;
    std::ifstream file(_filename.c_str());
    if(file)
    {
        std::string line;
        std::getline(file, line);

        for(std::string::iterator it=line.begin(); it != line.end(); it++)
        {
            if(*it==_field_separator)
                ncols++;
        }
        ncols++;
    }
    else
        std::cout << "CSVParser<T>::count_cols() => Erreur lecture fichier\n";

    return ncols;
}

//==============================================================
// CSVParser : read
//==============================================================
template <class T>
void CSVParser<T>::read()
{
    _matrix = boost::numeric::ublas::matrix<T>(_nrows, _ncols);

    int i=0;
    std::ifstream file(_filename.c_str());
    if(file)
    {
        // Fill Matrix _matrix
        CSVRow<T> row(_field_separator);
        for(unsigned k=0; k<_header_nrows; k++)  // Skip header
            file >> row;
        while(file >> row)
        {
            for(std::size_t j=0; j<row.size(); j++)
                _matrix(i, j) = row[j];
            i++;
        }
//        std::cout << "File imported.\n";
    }
    else
    {
        std::cerr << "error\n";
        exit(-1);
    }
}

//==============================================================
// CSVParser : read_transpose
//==============================================================
template <class T>
void CSVParser<T>::read_transpose()
{
    _matrix.resize(_ncols, _nrows); // Transpose

    int j=0;
    std::ifstream file(_filename.c_str());
    if(file)
    {
        // Fill Matrix _matrix
        CSVRow<T> row(_field_separator);
        for(unsigned k=0; k<_header_nrows; k++)  // Skip header
            file >> row;
        while(file >> row)
        {
            for(std::size_t i=0; i<row.size(); i++)
                _matrix(i, j) = row[i];
            j++;
        }
//        std::cout << "File imported.\n";
    }
    else
    {
        std::cerr << "error\n";
        exit(-1);
    }
}

#endif // CSVParser_HPP
