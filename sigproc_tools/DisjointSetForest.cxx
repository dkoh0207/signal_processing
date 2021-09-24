#ifndef __SIGPROC_TOOLS_DISJOINTSETFOREST_CXX__
#define __SIGPROC_TOOLS_DISJOINTSETFOREST_CXX__

#include "DisjointSetForest.h"


int sigproc_tools::DisjointSetForest::Find(const int x) 
{
    if (x == parent[x]) return x;
    else {
        int rep = Find(parent[x]);
        parent[x] = rep; // Path compression
        return rep;
    }
}

void sigproc_tools::DisjointSetForest::MakeSet()
{
    for (int i=0; i< (int) size; ++i) {
        parent[i] = i;
    }
    return;
}

void sigproc_tools::DisjointSetForest::MakeSet(const std::vector<int>& strongEdges)
{
    if (strongEdges.size() < 1) {
        std::string msg = "When constructing disjoint set parent with list of root "
                          "set members, list must contain at least one entry. Returning";
        std::cout << msg << std::endl;
        return;
    }

    for (int i=0; i< (int) size; ++i) {
        parent[i] = i;
    }

    for (const int& x : strongEdges) {
        parent[x] = rootNode;
    }
    return;
}

void sigproc_tools::DisjointSetForest::Union(const int x, const int y)
{
    int repX = Find(x);
    int repY = Find(y);

    if (repX == repY) return;

    else if (repX == rootNode) parent[repY] = repX;

    else if (repY == rootNode) parent[repX] = repY;

    else {
        int rankX = rank[repX];
        int rankY = rank[repY];

        if (rankX < rankY) {
            parent[repX] = repY;
            // std::cout << repX << " -> " << repY << std::endl;
        }
        else if (rankX > rankY) {
            parent[repY] = repX;
            // std::cout << repY << " -> " << repX << std::endl;
        }
        else {
            parent[repX] = repY;
            rank[repY] = rank[repY] + 1;
            // std::cout << repX << " -> " << repY << std::endl;
        }
    }
} 



#endif
