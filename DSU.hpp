#ifndef DSU_H
#define DSU_H

#include <map>
#include "Point.hpp"

struct Dsu
{
    std::map<Point*, Point*> dset;

    void makeSet(Point* v)
    {
        if (!dset.count(v))
            dset[v] = v;
    }
    Point* findSet(Point* v)
    {
        if (!dset.count(v))
            return nullptr;
        if (dset[v] == v)
            return v;
        return dset[v] = findSet(dset[v]);
    }
    void unionSet(Point* a, Point* b)
    {
        a = findSet(a);
        b = findSet(b);
        if (a != b && a && b)
            dset[b] = a;
    }
};

#endif // DSU_H