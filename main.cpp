// #define _GLIBCXX_DEBUG 1

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <random>
#include <string>
#include <thread>
#include <mutex>
#include <limits>
#include <queue>
#include <fstream>

#include "Point.hpp"
#include "KDNode.hpp"
#include "globals.h"
#include "GeometryUtils.hpp"
#include "DSU.hpp"
#include "utils.hpp"
#include "Edge.hpp"
#include "Kruskal.hpp"

// Variables globales para el procesamiento
std::vector<std::pair<KDNode *, KDNode *>> pairs; // Almacena pares de nodos KDNode que cumplen con las condiciones de separación
std::mutex mutex;                                 // Mutex para asegurar acceso seguro en entornos multihilo
int beta;                                         // Parámetro de control para la cantidad de pares procesados por iteración
double phi;                                       // Distancia mínima entre pares de nodos en ciertos algoritmos
std::vector<Edge> bccps;                          // Conjunto de aristas del bosque de componentes conexas
std::vector<double> coreDist;                     // Distancia del núcleo para cada punto
Point *ref;                                       // Referencia a los puntos en el conjunto de datos
int cost;                                         // Costo acumulado en la construcción de árboles de expansión mínima (MST)
int n;                                            // Número total de puntos

// funcion para encontrar pares bien separados utilizando la recursión paralela
void findPair(KDNode *a, KDNode *b, bool (*wellSeparated)(KDNode *, KDNode *))
{
    if (wellSeparated(a, b))
    {
        mutex.lock();
        pairs.emplace_back(a, b);
        mutex.unlock();
        return;
    }

    if (getDiameter(a) < getDiameter(b))
        std::swap(a, b);

    if (a->size > MIN_PARALLEL_RECURSION_SIZE)
    {
        std::thread leftThread(findPair, a->left, b, wellSeparated);
        std::thread rightThread(findPair, a->right, b, wellSeparated);

        leftThread.join();
        rightThread.join();
    }
    else
    {
        findPair(a->left, b, wellSeparated);
        findPair(a->right, b, wellSeparated);
    }
}

// funcion para dividir el trabajo en múltiples hilos recursivamente usando `wspd` (Well-Separated Pair Decomposition)
void wspd(KDNode *node, bool (*wellSeparated)(KDNode *, KDNode *))
{
    if (node->size > 1)
    {
        if (node->size > MIN_PARALLEL_RECURSION_SIZE)
        {
            std::thread leftThread(wspd, node->left, wellSeparated);
            std::thread rightThread(wspd, node->right, wellSeparated);

            leftThread.join();
            rightThread.join();
        }
        else
        {
            wspd(node->left, wellSeparated);
            wspd(node->right, wellSeparated);
        }

        findPair(node->left, node->right, wellSeparated);
    }
}

// funcion para calcular las distancias entre nodos utilizando diferentes metricas (Euclideana, HDBSCAN)
void bccpEuclidean(KDNode *n1, KDNode *n2, double &currDist, Point *&a, Point *&b)
{
    if (nodeDistance(n1, n2) > currDist)
        return;

    if (n1->size == 1 && n2->size == 1)
    {
        double nwDist = distance(*n1->point, *n2->point);
        if (nwDist < currDist)
        {
            a = n1->point;
            b = n2->point;
            currDist = nwDist;
        }
    }
    else
    {
        if (n1->size == 1)
        {
            if (nodeDistance(n1, n2->left) < nodeDistance(n1, n2->right))
            {
                bccpEuclidean(n1, n2->left, currDist, a, b);
                bccpEuclidean(n1, n2->right, currDist, a, b);
            }
            else
            {
                bccpEuclidean(n1, n2->right, currDist, a, b);
                bccpEuclidean(n1, n2->left, currDist, a, b);
            }
        }
        else if (n2->size == 1)
        {
            if (nodeDistance(n2, n1->left) < nodeDistance(n2, n1->right))
            {
                bccpEuclidean(n1->left, n2, currDist, a, b);
                bccpEuclidean(n1->right, n2, currDist, a, b);
            }
            else
            {
                bccpEuclidean(n1->right, n2, currDist, a, b);
                bccpEuclidean(n1->left, n2, currDist, a, b);
            }
        }
        else
        {
            std::pair<KDNode *, KDNode *> ordering[4];
            ordering[0] = std::make_pair(n1->left, n2->left);
            ordering[1] = std::make_pair(n1->left, n2->right);
            ordering[2] = std::make_pair(n1->right, n2->left);
            ordering[3] = std::make_pair(n1->right, n2->right);

            auto cmp = [&](std::pair<KDNode *, KDNode *> p1, std::pair<KDNode *, KDNode *> p2)
            {
                return nodeDistance(p1.first, p1.second) < nodeDistance(p2.first, p2.second);
            };
            std::sort(ordering, ordering + 4, cmp);

            for (int o = 0; o < 4; ++o)
                bccpEuclidean(ordering[o].first, ordering[o].second, currDist, a, b);
        }
    }
}

void bccpHDBSCAN(KDNode *n1, KDNode *n2, double &currDist, Point *&a, Point *&b)
{
    if (nodeDistance(n1, n2) > currDist)
        return;

    if (n1->size == 1 && n2->size == 1)
    {
        double nwDist = std::max(distance(*n1->point, *n2->point), coreDist[n1->point - ref]);
        nwDist = std::max(nwDist, coreDist[n2->point - ref]);
        if (nwDist < currDist)
        {
            a = n1->point;
            b = n2->point;
            currDist = nwDist;
        }
    }
    else
    {
        if (n1->size == 1)
        {
            if (nodeDistance(n1, n2->left) < nodeDistance(n1, n2->right))
            {
                bccpHDBSCAN(n1, n2->left, currDist, a, b);
                bccpHDBSCAN(n1, n2->right, currDist, a, b);
            }
            else
            {
                bccpHDBSCAN(n1, n2->right, currDist, a, b);
                bccpHDBSCAN(n1, n2->left, currDist, a, b);
            }
        }
        else if (n2->size == 1)
        {
            if (nodeDistance(n2, n1->left) < nodeDistance(n2, n1->right))
            {
                bccpHDBSCAN(n1->left, n2, currDist, a, b);
                bccpHDBSCAN(n1->right, n2, currDist, a, b);
            }
            else
            {
                bccpHDBSCAN(n1->right, n2, currDist, a, b);
                bccpHDBSCAN(n1->left, n2, currDist, a, b);
            }
        }
        else
        {
            std::pair<KDNode *, KDNode *> ordering[4];
            ordering[0] = std::make_pair(n1->left, n2->left);
            ordering[1] = std::make_pair(n1->left, n2->right);
            ordering[2] = std::make_pair(n1->right, n2->left);
            ordering[3] = std::make_pair(n1->right, n2->right);

            auto cmp = [&](std::pair<KDNode *, KDNode *> p1, std::pair<KDNode *, KDNode *> p2)
            { return nodeDistance(p1.first, p1.second) < nodeDistance(p2.first, p2.second); };
            std::sort(ordering, ordering + 4, cmp);

            for (int o = 0; o < 4; ++o)
            {
                bccpHDBSCAN(ordering[o].first, ordering[o].second, currDist, a, b);
            }
        }
    }
}

void splitParallelHelper(const int *first, const int *last, bool (*cond)(int), std::vector<int> &l, std::vector<int> &r)
{
    std::vector<int> ll, rr;
    for (const int *p = first; p < last; p++)
    {
        if (cond(*p))
            ll.push_back(*p);
        else
            rr.push_back(*p);
    }

    mutex.lock();
    l.insert(l.end(), ll.begin(), ll.end());
    r.insert(r.end(), rr.begin(), rr.end());
    mutex.unlock();
}

void splitParallel(const std::vector<int> &arr, bool (*cond)(int), std::vector<int> &l, std::vector<int> &r)
{
    size_t numThreads = std::thread::hardware_concurrency();
    size_t sz = arr.size() / numThreads;
    std::vector<std::thread *> threads(numThreads);

    for (int i = 0; i < numThreads - 1; i++)
    {
        const int *p = arr.data() + i * sz;
        threads[i] = new std::thread(splitParallelHelper, p, p + sz, cond, std::ref(l), std::ref(r));
    }
    threads[numThreads - 1] = new std::thread(splitParallelHelper, arr.data() + (numThreads - 1) * sz, arr.data() + arr.size(), cond, std::ref(l), std::ref(r));

    for (int i = 0; i < numThreads; i++)
    {
        threads[i]->join();
        delete threads[i];
    }
}

void parallelKruskal(std::vector<int> &sl1, std::vector<Edge> &result, Dsu &dsu)
{
    for (auto x : sl1)
    {
        auto &e = bccps[x];
        dsu.makeSet(e.u);
        dsu.makeSet(e.v);
    }

    std::sort(sl1.begin(), sl1.end(), [](int x, int y)
              { return bccps[x] < bccps[y]; });

    for (auto x : sl1)
    {
        auto &e = bccps[x];
        auto fs = dsu.findSet(e.u);
        if (fs == dsu.findSet(e.v) && fs)
            continue;

        cost += e.weight;
        result.push_back(e);
        dsu.unionSet(e.u, e.v);
    }
}

int parallelGeoFilterKruskal(std::vector<Edge> &result, void (*bccp)(KDNode *n1, KDNode *n2, double &currDist, Point *&a, Point *&b))
{
    Dsu unionFind;

    beta = 2;
    cost = 0;

    while (result.size() < (n - 1))
    {
        std::vector<int> ind(pairs.size());
        for (int i = 0; i < ind.size(); i++)
            ind[i] = i;

        std::vector<int> sl, su;
        splitParallel(ind, [](int x)
                      { return pairs[x].first->size + pairs[x].second->size <= beta; }, sl, su);

        phi = std::numeric_limits<double>::max();
        for (int i = 0; i < su.size(); i++)
        {
            int x = su[i];
            phi = std::min(phi, nodeDistance(pairs[x].first, pairs[x].second));
        }

        bccps.resize(sl.size());
        std::vector<int> ind2(sl.size());
        for (int i = 0; i < sl.size(); i++)
        {
            bccps[i].weight = std::numeric_limits<double>::max();
            bccp(pairs[sl[i]].first, pairs[sl[i]].second, bccps[i].weight, bccps[i].u, bccps[i].v);
            ind2[i] = i;
        }

        std::vector<int> sl1, sl2;
        splitParallel(ind2, [](int x)
                      { return bccps[x].weight <= phi; }, sl1, sl2);

        parallelKruskal(sl1, result, unionFind);

        std::vector<std::pair<KDNode *, KDNode *>> newPairs;

        for (int i = 0; i < sl2.size() + su.size(); i++)
        {
            std::pair<KDNode *, KDNode *> *p;
            if (i < sl2.size())
                p = &pairs[sl[sl2[i]]];
            else
                p = &pairs[su[i - sl2.size()]];
            auto fs = unionFind.findSet(p->first->point);
            if (fs == unionFind.findSet(p->second->point) && fs)
                continue;

            newPairs.push_back(*p);
        }

        pairs = newPairs;

        beta *= 2;
    }

    return cost;
}

void knnHelper(KDNode *node, Point &query, int k, int depth, std::priority_queue<std::pair<double, Point *>> &maxHeap)
{
    if (node == nullptr)
        return;

    if (node->size == 1)
    {
        double dist = distance(query, *node->point);
        if (maxHeap.size() < k)
            maxHeap.push({dist, node->point});
        else if (dist < maxHeap.top().first)
        {
            maxHeap.pop();
            maxHeap.push({dist, node->point});
        }
    }

    int axis = depth % query.coords.size();

    KDNode *first = nullptr;
    KDNode *second = nullptr;

    if (query.coords[axis] < node->point->coords[axis])
    {
        first = node->left;
        second = node->right;
    }
    else
    {
        first = node->right;
        second = node->left;
    }

    knnHelper(first, query, k, depth + 1, maxHeap);

    if (maxHeap.size() < k || std::abs(query.coords[axis] - node->point->coords[axis]) < maxHeap.top().first)
        knnHelper(second, query, k, depth + 1, maxHeap);
}

double knnDistance(Point &query, int k, KDNode *root)
{
    std::priority_queue<std::pair<double, Point *>> maxHeap;

    knnHelper(root, query, k, 0, maxHeap);

    while (maxHeap.empty() > 1)
    {
        maxHeap.pop();
    }

    return maxHeap.top().first;
}

void processCoreDistance(KDNode *nd)
{
    if (nd->size == 1)
    {
        if (coreDist[nd->point - ref] > nd->cdMx)
            nd->cdMx = coreDist[nd->point - ref];
        if (coreDist[nd->point - ref] < nd->cdMn)
            nd->cdMn = coreDist[nd->point - ref];
    }
    else
    {
        if (nd->size > MIN_PARALLEL_RECURSION_SIZE)
        {
            std::thread leftThread(processCoreDistance, nd->left);
            std::thread rightThread(processCoreDistance, nd->right);

            leftThread.join();
            rightThread.join();
        }
        else
        {
            processCoreDistance(nd->left);
            processCoreDistance(nd->right);
        }
        nd->cdMx = std::max(nd->left->cdMx, nd->right->cdMx);
        nd->cdMn = std::min(nd->left->cdMn, nd->right->cdMn);
    }
}

struct DendroNode
{
    size_t a, b;
    double c;
    size_t d;
};

void generateDendrogram(std::vector<Edge> &edges, std::vector<DendroNode> &dendro)
{
    std::sort(edges.begin(), edges.end());

    Dsu uf;

    for (auto &e : edges)
    {
        uf.makeSet(e.u);
        uf.makeSet(e.v);
    }

    size_t idx = n;

    std::vector<size_t> idxMap(n);

    std::vector<size_t> sizes(n);

    for (int i = 0; i < n; i++)
    {
        idxMap[i] = i;
        sizes[i] = 1;
    }

    dendro.resize(edges.size());

    for (size_t i = 0; i < n - 1; ++i)
    {
        auto u = uf.findSet(edges[i].u);
        auto v = uf.findSet(edges[i].v);
        dendro[i] = {idxMap[u - ref], idxMap[v - ref], edges[i].weight, sizes[u - ref] + sizes[v - ref]};
        uf.unionSet(u, v);
        auto newIdx = uf.findSet(u);
        idxMap[newIdx - ref] = idx;
        sizes[newIdx - ref] = sizes[u - ref] + sizes[v - ref];
        idx++;
    }
}

int main()
{
    std::vector<Point> points;
    std::ifstream file_("points.txt");
    double x, y;

    while (file_ >> x >> y)
    {
        points.push_back({x, y});
    }

    file_.close();

    // Verificación de los primeros puntos
    for (size_t i = 0; i < 10 && i < points.size(); ++i)
    {
        std::cout << "Point " << i + 1 << ": (" << points[i].coords[0] << ", " << points[i].coords[1] << ")\n";
    }

    std::vector<Point *> pointPtrs(points.size());
    for (int i = 0; i < points.size(); i++)
        pointPtrs[i] = &points[i];

    n = points.size();

    KDNode *root = nullptr;
    buildKDTreeParallel(pointPtrs, root);

    // std::cout << "KD-Tree (In-order traversal):" << std::endl;
    // printKDTree(root);

    wspd(root, wellSeparatedEuclidean);

    // for (auto& p : pairs)
    // {
    //     std::cout << '(';
    //     for (auto c : p.first->point->coords)
    //         std::cout << c << ", ";
    //     std::cout << ") - (";
    //     for (auto c : p.second->point->coords)
    //         std::cout << c << ", ";
    //     std::cout << ")\n";
    // }

    std::vector<Edge> result;
    parallelGeoFilterKruskal(result, bccpEuclidean);

    std::cout << "Edges of the EMST:\n";

    std::ofstream file("emst.txt");

    for (auto &p : points)
    {
        file << "0 ";
        for (auto x : p.coords)
            file << x << ' ';
        file << '\n';
    }

    ref = points.data();
    for (auto &e : result)
    {
        /*
        std::cout << '(';
        for (int i = 0; i < e.u->coords.size(); i++)
            std::cout << e.u->coords[i] << ((i == e.u->coords.size() - 1) ? ") " : ", ");
        std::cout << "- (";
        for (int i = 0; i < e.v->coords.size(); i++)
            std::cout << e.v->coords[i] << ((i == e.v->coords.size() - 1) ? ") " : ", ");
        std::cout << "- " << e.weight << '\n';
        */

        file << "1 " << (e.u - ref) << ' ' << (e.v - ref) << '\n';
    }

    file.close();

    int minPts = 2;
    coreDist.resize(points.size());

    for (int i = 0; i < coreDist.size(); i++)
    {
        coreDist[i] = knnDistance(points[i], minPts, root);
    }

    processCoreDistance(root);

    pairs.clear();
    wspd(root, wellSeparatedHDBSCAN);

    result.clear();
    bccps.clear();
    parallelGeoFilterKruskal(result, bccpHDBSCAN);

    std::cout << "Edges of the HDBSCAN:\n";

    for (auto &e : result)
    {
        /*
        std::cout << '(';
        for (int i = 0; i < e.u->coords.size(); i++)
            std::cout << e.u->coords[i] << ((i == e.u->coords.size() - 1) ? ") " : ", ");
        std::cout << "- (";
        for (int i = 0; i < e.v->coords.size(); i++)
            std::cout << e.v->coords[i] << ((i == e.v->coords.size() - 1) ? ") " : ", ");
        std::cout << "- " << e.weight << '\n';
        */
    }

    std::vector<DendroNode> dendro;
    generateDendrogram(result, dendro);

    std::ofstream file2("dendro.txt");

    for (auto &x : dendro)
        file2 << x.a << ' ' << x.b << ' ' << x.c << ' ' << x.d << '\n';

    file2.close();

    return 0;
}