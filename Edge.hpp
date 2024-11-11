#ifndef EDGE_H
#define EDGE_H

#include "Point.hpp"

struct Edge
{
    Point *u, *v;
    double weight;
    bool operator<(Edge const& other)
    {
        return weight < other.weight;
    }
};


#endif // EDGE_H
