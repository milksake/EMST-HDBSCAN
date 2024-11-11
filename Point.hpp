#ifndef POINT_H
#define POINT_H

#include <vector>

struct Point
{
    std::vector<double> coords;

    Point(std::initializer_list<double> values) : coords(values) {}

    double operator[](size_t index) const
    {
        return coords[index];
    }

    double& operator[](size_t index)
    {
        return coords[index];
    }

    size_t size() const
    {
        return coords.size();
    }
};

#endif // POINT_H
