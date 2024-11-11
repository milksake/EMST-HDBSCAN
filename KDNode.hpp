#ifndef KDNODE_H
#define KDNODE_H

#include <vector>
#include <algorithm>
#include <thread>
#include <iostream>
#include "Point.hpp"
#include "globals.h"

struct KDNode
{
    Point* point;
    KDNode* left;
    KDNode* right;
    size_t size;

    Point mn, mx;
    double cdMn, cdMx;

    KDNode(Point& p) : point(&p), left(nullptr), right(nullptr), size(1), mn(p), mx(p), cdMn(std::numeric_limits<double>::max()), cdMx(-1) {}
};

void buildKDTreeParallel(std::vector<Point*>& points, KDNode*& node, int depth = 0)
{
    if (points.size() == 1)
    {
        node = new KDNode(*points[0]);
        return;
    }
    if (points.empty())
    {
        node = nullptr;
        return;
    }

    int axis = depth % points[0]->size();
    std::sort(points.begin(), points.end(), [axis](Point* p1, Point* p2) {
        return (*p1)[axis] < (*p2)[axis];
    });

    int median = points.size() / 2;
    node = new KDNode(*points[median]);

    std::vector<Point*> leftPoints(points.begin(), points.begin() + median);
    std::vector<Point*> rightPoints(points.begin() + median, points.end());

    if (node->size > MIN_PARALLEL_RECURSION_SIZE)
    {
        std::thread leftThread(buildKDTreeParallel, std::ref(leftPoints), std::ref(node->left), depth + 1);
        std::thread rightThread(buildKDTreeParallel, std::ref(rightPoints), std::ref(node->right), depth + 1);

        leftThread.join();
        rightThread.join();
    }
    else
    {
        buildKDTreeParallel(std::ref(leftPoints), std::ref(node->left), depth + 1);
        buildKDTreeParallel(std::ref(rightPoints), std::ref(node->right), depth + 1);
    }

    size_t leftSize = node->left ? node->left->size : 0;
    size_t rightSize = node->right ? node->right->size : 0;
    node->size = leftSize + rightSize;

    if (!node->left && !node->right)
        return;
    if (!node->left)
    {
        node->mn = node->right->mn;
        node->mx = node->right->mx;
        return;
    }
    if (!node->right)
    {
        node->mn = node->left->mn;
        node->mx = node->left->mx;
        return;
    }

    for (int i = 0; i < points[0]->size(); i++)
    {
        node->mn[i] = std::min(node->left->mn[i], node->right->mn[i]);
        node->mx[i] = std::max(node->left->mx[i], node->right->mx[i]);
    }
}

void printKDTree(KDNode* root, int depth = 0)
{
    if (root == nullptr) return;

    printKDTree(root->left, depth + 1);

    std::cout << std::string(depth * 2, ' ') << "(";
    for (size_t i = 0; i < root->point->size(); ++i)
    {
        std::cout << (*root->point)[i];
        if (i != root->point->size() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << "), Size: " << root->size << ' ' << root->point << std::endl;

    printKDTree(root->right, depth + 1);
}

#endif // KDNODE_H