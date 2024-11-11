#ifndef UTILS_H
#define UTILS_H

#include "Point.hpp"
#include "KDNode.hpp"

bool compare(const Point& p1, const Point& p2, int depth)
{
    return p1[depth % p1.size()] < p2[depth % p1.size()];
}


double getDiameter(KDNode* node)
{
    if (!node)
        return -1.0;

    double ans = 0.0;
    for (int i = 0; i < node->point->size(); i++)
        ans = std::max(ans, node->mx[i] - node->mn[i]);
    return ans;
}


#endif // UTILS_H