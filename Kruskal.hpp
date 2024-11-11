#ifndef KRUSKAL_H
#define KRUSKAL_H

#include <vector>
#include <algorithm>
#include "DSU.hpp"
#include "Edge.hpp"

int kruskal(int n, std::vector<Edge>& edges, std::vector<Edge>& result)
{
    int cost = 0;

    Dsu dsu;

    for (auto& e : edges)
    {
        dsu.makeSet(e.u);
        dsu.makeSet(e.v);
    }

    std::sort(edges.begin(), edges.end());

    for (Edge e : edges)
    {
        if (dsu.findSet(e.u) != dsu.findSet(e.v))
        {
            cost += e.weight;
            result.push_back(e);
            dsu.unionSet(e.u, e.v);
        }
    }

    return cost;
}

#endif // KRUSKAL_H