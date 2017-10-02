#ifndef CSVROW_HPP
#define CSVROW_HPP

/*!
 * \file   CSVRow.hpp
 * \brief  Contains current line in the read CSV file.
 * \author Clément Niel
 * \date   22/04/2015
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>


/*!
 * \class CSVRow
 * \brief Contains current line in the read CSV file. Templated.
 */
template<class T>
class CSVRow
{
    public:
        CSVRow(char sep);

        /*!
         * \brief Retrieves the stored data at the index position of the row.
         * \param[in] index Index position to retrieve
         * \return Information stored in data for a given index
         */
        T const& operator[](std::size_t index) const { return _data[index]; }

        /*!
         * \brief Retrieves the stored data at the index position of the row.
         * \param[in] index Index position to retrieve
         * \return Information stored in data for a given index
         */
        T const at(std::size_t index) const { return _data.at(index); }

        std::size_t size() const { return _data.size(); }

        /*!
         * \brief Call this function to read the next row (if any) ad update the current stored data with this row.
         * \param[in] str Stream from which data is extracted.
         */
        inline void readNextRow(std::istream& str);

    private:
        std::vector<T> _data;
        char _field_separator;
};

/* ###################################################################### IMPLEMENTATION below this point ###################################################################### */
//=========================================
// Constructor
//=========================================
template <class T>
CSVRow<T>::CSVRow(char sep) : _field_separator(sep)
{}

//=========================================
// operator >>
//=========================================
template<class T>
std::istream& operator>>(std::istream& str,CSVRow<T>& data)
{
    data.readNextRow(str);
    return str;
}

//=========================================
// CSVRow : readNextRow
//=========================================
template <class T>
inline void CSVRow<T>::readNextRow(std::istream& str)
{
    std::string line;
    std::getline(str,line);
    line.erase(std::remove(line.begin(), line.end(), '\"'), line.end());

    std::stringstream lineStream(line);
    T cell;

    _data.clear();
    while(lineStream >> cell)
    {
        _data.push_back(cell);
        if(lineStream.peek() == _field_separator)
            lineStream.ignore();
    }
}

//================================================
// CSVRow : readNextRow -> override with strings
//================================================
template <>
inline void CSVRow<std::string>::readNextRow(std::istream& str)
{
    std::string line;
    std::getline(str,line);
    line.erase(std::remove(line.begin(), line.end(), '\"'), line.end());

    std::stringstream lineStream(line);
    std::string cell;

    _data.clear();
    while(std::getline(lineStream, cell, _field_separator))
    {
        _data.push_back(cell);
    }
}


#endif // CSVROW_HPP
