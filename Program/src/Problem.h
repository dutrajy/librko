#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "Data.h"

class Problem
{
public:
    int dimension; // Size of the random keys vector
    virtual ~Problem() = default;
    virtual const double decode(TSol &solution) const = 0;
    virtual void loadFromFile(char *path) = 0;
};

#endif // _PROBLEM_H
