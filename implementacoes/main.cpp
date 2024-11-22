#include <iostream>
#include "structs.hpp"

using namespace std;

int main()
{
    AdjList* adjList = new AdjList(5);
    adjList->addEdge(0, 1);
    adjList->addEdge(0, 2);
    adjList->addEdge(0, 3);
    adjList->addEdge(0, 4);
    adjList->addEdge(1, 2);
    adjList->addEdge(1, 3);
    adjList->addEdge(1, 4);
    adjList->addEdge(2, 3);
    adjList->addEdge(2, 4);
    adjList->addEdge(3, 4);
    adjList->print();

    delete adjList;

    return 0;
}
