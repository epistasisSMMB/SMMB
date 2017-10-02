#ifndef SCORING2_HPP
#define SCORING2_HPP

class Scoring2
{
public:
    Scoring2(double p_pvalue, int p_index1, int p_index2): p_value(p_pvalue)
    {
        if(p_index1 <= p_index2)
        {
            index1 = p_index1;
            index1 = p_index2;
        }
        else
        {
            index1 = p_index2;
            index2 = p_index1;
        }
    }

    bool operator>(const Scoring2& other) const
    {
        return p_value > other.p_value;
    }

    bool operator<(const Scoring2& other) const
    {
        return p_value < other.p_value;
    }

    double p_value;
    int index1;
    int index2;
};

#endif // SCORING2_HPP
