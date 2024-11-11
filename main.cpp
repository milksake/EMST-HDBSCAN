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

#define MIN_PARALLEL_RECURSION_SIZE 2000

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

bool compare(const Point& p1, const Point& p2, int depth)
{
    return p1[depth % p1.size()] < p2[depth % p1.size()];
}

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

double getDiameter(KDNode* node)
{
    if (!node)
        return -1.0;

    double ans = 0.0;
    for (int i = 0; i < node->point->size(); i++)
        ans = std::max(ans, node->mx[i] - node->mn[i]);
    return ans;
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

struct Edge
{
    Point *u, *v;
    double weight;
    bool operator<(Edge const& other)
    {
        return weight < other.weight;
    }
};

std::vector<std::pair<KDNode*, KDNode*>> pairs;
std::mutex mutex;
int beta;
double phi;
std::vector<Edge> bccps;
std::vector<double> coreDist;
Point* ref;
int cost;
int n;

void findPair(KDNode* a, KDNode* b, bool(*wellSeparated)(KDNode*, KDNode*))
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

void wspd(KDNode* node, bool(*wellSeparated)(KDNode*, KDNode*))
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

void bccpEuclidean(KDNode* n1, KDNode* n2, double& currDist, Point*& a, Point*& b)
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
            std::pair<KDNode*, KDNode*> ordering[4];
            ordering[0] = std::make_pair(n1->left, n2->left);
            ordering[1] = std::make_pair(n1->left, n2->right);
            ordering[2] = std::make_pair(n1->right, n2->left);
            ordering[3] = std::make_pair(n1->right, n2->right);

            auto cmp = [&](std::pair<KDNode*,KDNode*> p1, std::pair<KDNode*,KDNode*> p2) {
                return nodeDistance(p1.first, p1.second) < nodeDistance(p2.first, p2.second);
            };
            std::sort(ordering, ordering + 4, cmp);

            for (int o = 0; o < 4; ++o)
                bccpEuclidean(ordering[o].first, ordering[o].second, currDist, a, b);
        }
    }
}

void bccpHDBSCAN(KDNode* n1, KDNode* n2, double& currDist, Point*& a, Point*& b)
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
            std::pair<KDNode*, KDNode*> ordering[4];
            ordering[0] = std::make_pair(n1->left, n2->left);
            ordering[1] = std::make_pair(n1->left, n2->right);
            ordering[2] = std::make_pair(n1->right, n2->left);
            ordering[3] = std::make_pair(n1->right, n2->right);

            auto cmp = [&](std::pair<KDNode*,KDNode*> p1, std::pair<KDNode*,KDNode*> p2) {
                    return nodeDistance(p1.first, p1.second) < nodeDistance(p2.first, p2.second);};
            std::sort(ordering, ordering + 4, cmp);

            for (int o=0; o<4; ++o) {
                bccpHDBSCAN(ordering[o].first, ordering[o].second, currDist, a, b);}
        }
    }
}

void splitParallelHelper(const int *first, const int *last, bool(*cond)(int), std::vector<int>& l, std::vector<int>& r)
{
    std::vector<int> ll, rr;
    for (const int* p = first; p < last; p++)
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

void splitParallel(const std::vector<int>& arr, bool(*cond)(int), std::vector<int>& l, std::vector<int>& r)
{
    size_t numThreads = std::thread::hardware_concurrency();
    size_t sz = arr.size() / numThreads;
    std::vector<std::thread*> threads(numThreads);

    for (int i = 0; i < numThreads-1; i++)
    {
        const int *p = arr.data() + i * sz;
        threads[i] = new std::thread(splitParallelHelper, p, p + sz, cond, std::ref(l), std::ref(r));
    }
    threads[numThreads-1] = new std::thread(splitParallelHelper, arr.data() + (numThreads-1) * sz, arr.data() + arr.size(), cond, std::ref(l), std::ref(r));
    
    for (int i = 0; i < numThreads; i++)
    {
        threads[i]->join();
        delete threads[i];
    }
}

void parallelKruskal(std::vector<int>& sl1, std::vector<Edge>& result, Dsu& dsu)
{
    for (auto x : sl1)
    {
        auto& e = bccps[x];
        dsu.makeSet(e.u);
        dsu.makeSet(e.v);
    }

    std::sort(sl1.begin(), sl1.end(), [](int x, int y) {
        return bccps[x] < bccps[y];
    });

    for (auto x : sl1)
    {
        auto& e = bccps[x];
        auto fs = dsu.findSet(e.u);
        if (fs == dsu.findSet(e.v) && fs)
            continue;
            
        cost += e.weight;
        result.push_back(e);
        dsu.unionSet(e.u, e.v);
    }
}

int parallelGeoFilterKruskal(std::vector<Edge>& result, void (*bccp)(KDNode* n1, KDNode* n2, double& currDist, Point*& a, Point*& b))
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
        splitParallel(ind, [](int x) {
            return pairs[x].first->size + pairs[x].second->size <= beta;
        }, sl, su);

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
        splitParallel(ind2, [](int x) {
            return bccps[x].weight <= phi;
        }, sl1, sl2);

        parallelKruskal(sl1, result, unionFind);

        std::vector<std::pair<KDNode*, KDNode*>> newPairs;

        for (int i = 0; i < sl2.size() + su.size(); i++)
        {
            std::pair<KDNode*, KDNode*>* p;
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

void knnHelper(KDNode* node, Point& query, int k, int depth, std::priority_queue<std::pair<double, Point*>>& maxHeap)
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

    KDNode* first = nullptr;
    KDNode* second = nullptr;

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

double knnDistance(Point& query, int k, KDNode* root)
{
    std::priority_queue<std::pair<double, Point*>> maxHeap;

    knnHelper(root, query, k, 0, maxHeap);

    while (maxHeap.empty() > 1) {
        maxHeap.pop();
    }

    return maxHeap.top().first;
}

void processCoreDistance(KDNode* nd)
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

void generateDendrogram(std::vector<Edge>& edges, std::vector<DendroNode>& dendro)
{
    std::sort(edges.begin(), edges.end());

    Dsu uf;

    for (auto& e : edges)
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

    for(size_t i = 0; i < n-1; ++i)
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
    std::vector<Point> points = {
        {3, 6}, {17, 15}, {13, 15}, {6, 12}, {9, 1}, {2, 7}, {10, 19},
        {4, 5}, {12, 8}
    };
    std::vector<Point*> pointPtrs(points.size());
    for (int i = 0; i < points.size(); i++)
        pointPtrs[i] = &points[i];

    n = points.size();
    
    KDNode* root = nullptr;
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

    for (auto& p : points)
    {
        file << "0 ";
        for (auto x : p.coords)
            file << x << ' ';
        file << '\n';
    }

    ref = points.data();
    for (auto& e : result)
    {
        std::cout << '(';
        for (int i = 0; i < e.u->coords.size(); i++)
            std::cout << e.u->coords[i] << ((i == e.u->coords.size()-1) ? ") " : ", ");
        std::cout << "- (";
        for (int i = 0; i < e.v->coords.size(); i++)
            std::cout << e.v->coords[i] << ((i == e.v->coords.size()-1) ? ") " : ", ");
        std::cout << "- " << e.weight << '\n';

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

    for (auto& e : result)
    {
        std::cout << '(';
        for (int i = 0; i < e.u->coords.size(); i++)
            std::cout << e.u->coords[i] << ((i == e.u->coords.size()-1) ? ") " : ", ");
        std::cout << "- (";
        for (int i = 0; i < e.v->coords.size(); i++)
            std::cout << e.v->coords[i] << ((i == e.v->coords.size()-1) ? ") " : ", ");
        std::cout << "- " << e.weight << '\n';
    }

    std::vector<DendroNode> dendro;
    generateDendrogram(result, dendro);

    std::ofstream file2("dendro.txt");

    for (auto& x : dendro)
        file2 << x.a << ' ' << x.b << ' ' << x.c << ' ' << x.d << '\n';

    file2.close();

    return 0;
}
