#ifndef GEOMETRYUTILS_H
#define GEOMETRYUTILS_H

#include "Point.hpp"
#include "KDNode.hpp"
#include <cmath>

double nodeDistance(KDNode *n1, KDNode *n2)
{
    for (int d = 0; d < n1->point->size(); ++d)
    {
        if (n1->mn[d] > n2->mx[d] || n2->mn[d] > n1->mx[d])
        {
            double rsqr = 0;
            for (int dd = d; dd < n1->point->size(); ++dd)
            {
                double tmp = std::max(n1->mn[dd] - n2->mx[dd], n2->mn[dd] - n1->mx[dd]);
                tmp = std::max(tmp, (double)0);
                rsqr += tmp * tmp;
            }
            return std::sqrt(rsqr);
        }
    }
    return 0;
}

double distance(Point& a, Point& b)
{
    double xx = 0;
    for (int i = 0; i < a.size(); ++i)
    {
        double yy = (a[i] - b[i]);
        xx += yy * yy;
    }

    return sqrt(xx);
}

bool wellSeparatedEuclidean(KDNode* u, KDNode* v)
{
    const double s = 2;

    double circleDiam_u = 0;
    double circleDiam_v = 0;
    double circleDistance = 0;
    for (int d = 0; d < u->point->size(); d++)
    {
        double uTmpDiff = u->mx[d] - u->mn[d];
        double vTmpDiff = v->mx[d] - v->mn[d];
        double uTmpAvg = (u->mx[d] + u->mn[d])/2;
        double vTmpAvg = (v->mx[d] + v->mn[d])/2;
        circleDistance += (uTmpAvg - vTmpAvg) * (uTmpAvg - vTmpAvg);
        circleDiam_u += uTmpDiff * uTmpDiff;
        circleDiam_v += vTmpDiff * vTmpDiff;
    }
    circleDiam_u = std::sqrt(circleDiam_u);
    circleDiam_v = std::sqrt(circleDiam_v);

    double myRadius = std::max(circleDiam_u, circleDiam_v)/2;
    circleDistance = std::sqrt(circleDistance) - circleDiam_u/2 - circleDiam_v/2;

    return circleDistance >= (s * myRadius);
}

bool wellSeparatedHDBSCAN(KDNode* u, KDNode* v)
{
    double circleDiam_u = 0;
    double circleDiam_v = 0;
    double circleDistance = 0;
    for (int d = 0; d < u->point->size(); d++)
    {
        double uTmpDiff = u->mx[d] - u->mn[d];
        double vTmpDiff = v->mx[d] - v->mn[d];
        double uTmpAvg = (u->mx[d] + u->mn[d])/2;
        double vTmpAvg = (v->mx[d] + v->mn[d])/2;
        circleDistance += (uTmpAvg - vTmpAvg) * (uTmpAvg - vTmpAvg);
        circleDiam_u += uTmpDiff * uTmpDiff;
        circleDiam_v += vTmpDiff * vTmpDiff;
    }
    circleDiam_u = std::sqrt(circleDiam_u);
    circleDiam_v = std::sqrt(circleDiam_v);

    double myRadius = std::max(circleDiam_u, circleDiam_v)/2;
    double myDiam = std::max(2*myRadius, u->cdMx);
    myDiam = std::max(myDiam, v->cdMx);

    circleDistance = sqrt(circleDistance) - circleDiam_u/2 - circleDiam_v/2;
    bool geoSep = circleDistance >= 2 * myRadius;
    circleDistance = std::max(circleDistance, u->cdMn);
    circleDistance = std::max(circleDistance, v->cdMn);

    if (circleDistance >= myDiam)
	    return true || geoSep;
    else
	    return false || geoSep;
}

#endif // GEOMETRYUTILS_H