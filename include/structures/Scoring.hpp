#ifndef SCORING_HPP
#define SCORING_HPP

#include <iostream>

class Scoring
{
public:
    Scoring(double p_pvalue, int p_index): p_value(p_pvalue), index(p_index) {}

    bool operator>(const Scoring& other) const
    {
        return p_value > other.p_value;
    }

    bool operator<(const Scoring& other) const
    {
        return p_value < other.p_value;
    }

    void print() const
    {
        std::cout << "score: " << p_value << "\t index : " << index << std::endl;
    }

    double p_value;
    int index;
};

#endif // SCORING_HPP
