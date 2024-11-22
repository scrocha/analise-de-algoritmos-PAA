#include <iostream>
#include "structs.hpp"

using namespace std;

int main()
{
    WeightedAdjList* adjList = new WeightedAdjList(5);
    adjList->addEdge(0, 1, 1);
    adjList->addEdge(0, 2, 2);
    adjList->addEdge(0, 3, 3);
    adjList->addEdge(0, 4, 4);
    adjList->addEdge(1, 2, 5);
    adjList->addEdge(1, 3, 6);
    adjList->addEdge(1, 4, 7);
    adjList->addEdge(2, 3, 8);
    adjList->addEdge(2, 4, 9);
    adjList->addEdge(3, 4, 10);
    adjList->print();

    delete adjList;

    return 0;
}
