#include "structs.hpp"
#include <iostream>
#include <vector>
#include <queue>
#include <climits>

using namespace std;

#define INFPOS INT_MAX

// O(V + E)
void canReach(AdjList* adjList, Vertex start, bool* visited)
{
    int numVertices = adjList->V;
    if (start < 0 || start >= numVertices) { return; }
    // O(V)
    for (int i = 0; i < numVertices; i++) { visited[i] = false; }

    queue<Vertex> queue;
    queue.push(start);
    visited[start] = true;
    // O(V + E)
    while (!queue.empty())
    {
        Vertex currentVertex = queue.front();
        queue.pop();

        Node* adjNode = adjList->adj[currentVertex].head;
        while (adjNode)
        {
            if (!visited[adjNode->vertex])
            {
                queue.push(adjNode->vertex);
                visited[adjNode->vertex] = true;
            }
            adjNode = adjNode->next;
        }
    }
}

// O(V + E)
int descendentCount(AdjList* adjList, Vertex vertex, List* L2)
{
    int numVertices = adjList->V;
    if (vertex < 0 || vertex >= numVertices) { return 0; }

    int count = 0;
    bool* descendents = new bool[numVertices];
    // O(V + E)
    canReach(adjList, vertex, descendents);
    // O(V)
    for (int i = 0; i < numVertices; i++)
    {
        if (descendents[i] && L2->contains(i)) { count++; }
    }
    return count;
}

// O(V * (V + E))
void dfsDescendants(AdjList* adjList, List* L2, bool* inL2)
{
    int numVertices = adjList->V;
    bool* noDescendents = new bool[numVertices];
    // O(V)
    for (int i = 0; i < numVertices; i++) { noDescendents[i] = false; }
    // O(V * (V + E))
    while (true)
    {
        bool addedVertexL2 = false;
        for (int i = 0; i < numVertices; i++)
        {
            if (noDescendents[i] || inL2[i]) { continue; }

            // O(V + E)
            int count = descendentCount(adjList, i, L2);

            if (count == 0) { noDescendents[i] = true; }
            else if (count % 2 == 0 && !L2->contains(i))
            {
                L2->add(i);
                inL2[i] = true;
                addedVertexL2 = true;
            }
        }
        if (!addedVertexL2) { break; }
    }
}

// Questão 1 a
// O(V * (V + E))
List* computeL2_a(AdjList* adjList, List* L1)
{
    if (L1->head == nullptr) { return nullptr; }
    List* L2 = new List();

    bool* inL2 = new bool[adjList->V];
    for (int i = 0; i < adjList->V; i++) { inL2[i] = false; }

    Node* temp = L1->head;
    while (temp)
    {
        L2->add(temp->vertex);
        inL2[temp->vertex] = true;
        temp = temp->next;
    }

    // O(V * (V + E))
    dfsDescendants(adjList, L2, inL2);

    return L2;
}

// O(V * (V + E))
void ancestorCount(AdjList* adjList, List* L2, int* count)
{
    int numVertices = adjList->V;
    for (int i = 0; i < numVertices; i++) { count[i] = 0; }

    Node* currentNode = L2->head;
    // O(V * (V + E))
    while (currentNode)
    {
        int currentVertex = currentNode->vertex;
        bool* visited = new bool[numVertices];
        // O(V + E)
        canReach(adjList, currentVertex, visited);
        // O(V)
        for (int i = 0; i < numVertices; i++)
        {
            if (visited[i]) { count[i]++; }
        }
    }
}

// O(V^2 * (V + E))
void dfsAncestors(AdjList* adjList, List* L2, bool* inL2)
{
    int numVertices = adjList->V;
    bool* noAncestors = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { noAncestors[i] = false; }

    // O(V^2 * (V + E))
    while (true)
    {
        bool addedVertexL2 = false;
        int* count = new int[numVertices];
        // O(V * (V + E))
        ancestorCount(adjList, L2, count);
        // O(V)
        for (int i = 0; i < numVertices; i++)
        {
            if (noAncestors[i] || inL2[i]) { continue; }

            if (count[i] == 0) { noAncestors[i] = true; }
            else if (count[i] % 2 == 0 && !L2->contains(i))
            {
                L2->add(i);
                inL2[i] = true;
                addedVertexL2 = true;
            }
        }
        delete[] count;
        if (!addedVertexL2) { break; }
    }
}

// Questão 1 b
// O(V^2 * (V + E))
List* computeL2_b(AdjList* adjList, List* L1)
{
    if (L1->head == nullptr) { return nullptr; }
    List* L2 = new List();

    bool* inL2 = new bool[adjList->V];
    // O(V)
    for (int i = 0; i < adjList->V; i++) { inL2[i] = false; }

    Node* temp = L1->head;
    while (temp)
    {
        L2->add(temp->vertex);
        inL2[temp->vertex] = true;
        temp = temp->next;
    }
    // O(V^2 * (V + E))
    dfsAncestors(adjList, L2, inL2);

    return L2;
}

int findRoot(AdjList* adjList)
{
    int numVertices = adjList->V;
    bool* hasAncestors = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { hasAncestors[i] = false; }

    int root = -1;
    while (true)
    {
        bool foundRoot = false;
        for (int i = 0; i < numVertices; i++)
        {
            if (hasAncestors[i]) { continue; }

            bool* visited = new bool[numVertices];
            canReach(adjList, i, visited);
            int count = 0;

            for (int j = 0; j < numVertices; j++)
            {
                if (visited[j])
                {
                    hasAncestors[j] = true;
                    count++;
                }
            }

            delete[] visited;

            if (count == numVertices - 2)
            {
                root = i;
                foundRoot = true;
                break;
            }
        }
        if (foundRoot) { break; }
    }
    delete[] hasAncestors;

    return root;
}

int olderCousinCount(AdjList* adjList, int vertex, List* L2)
{
    int count = 0;

    int parent = -1;
    for (int i = 0; i < adjList->V; i++)
    {
        Node* temp = adjList->adj[i].head;
        while (temp) {
            if (temp->vertex == vertex)
            {
                parent = i;
                break;
            }
            temp = temp->next;
        }
        if (parent != -1) break;
    }

    if (parent == -1) return 0;

    Node* sibling = adjList->adj[parent].head;
    while (sibling)
    {
        if (sibling->vertex != vertex && sibling->vertex < vertex && L2->contains(sibling->vertex))
        {
            count++;
        }
        sibling = sibling->next;
    }

    return count;
}


void dfsOlderCousins(AdjList* adjList, List* L2) {
    int numVertices = adjList->V;
    bool* inL2 = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { inL2[i] = false; }

    while (true) {
        bool added = false;
        for (int i = 0; i < numVertices; i++) {
            cout << "i: " << i << "| older cousins: " << olderCousinCount(adjList, i, L2) << endl;
            if (!inL2[i] && olderCousinCount(adjList, i, L2) % 2 != 0) {
                if (!L2->contains(i)) {
                    L2->add(i);
                    inL2[i] = true;
                    added = true;
                }
            }
        }
        if (!added) { break; }
    }

    delete[] inL2;
}

// Questão 1 c
List* computeL2_c(AdjList* adjList, List* L1) {
    List* L2 = new List();

    Node* temp = L1->head;

    while (temp) {
        L2->add(temp->vertex);
        temp = temp->next;
    }
    if (L2->head == nullptr) { return L2; }

    dfsOlderCousins(adjList, L2);

    return L2;
}

// O(V + E)
List* bfs(AdjList* adjList, Vertex start, Vertex end)
{
    if (start == end || start < 0 || end < 0 || start >= adjList->V || end >= adjList->V) { return nullptr; }

    List* path = new List();
    bool* visited = new bool[adjList->V];
    int* parent = new int[adjList->V];
    // O(V)
    for (int i = 0; i < adjList->V; i++)
    {
        visited[i] = false;
        parent[i] = -1;
    }

    queue<Vertex> queue;
    queue.push(start);
    visited[start] = true;

    bool hasEnded = false;
    // O(V + E)
    while (!queue.empty())
    {
        // O(1)
        Vertex currentVertex = queue.front();
        queue.pop();

        Node* adjNode = adjList->adj[currentVertex].head;
        while (adjNode)
        {
            if (!visited[adjNode->vertex])
            {
                queue.push(adjNode->vertex);
                visited[adjNode->vertex] = true;
                parent[adjNode->vertex] = currentVertex;

                if (adjNode->vertex == end)
                {
                    hasEnded = true;
                    break;
                }
            }
            adjNode = adjNode->next;
        }
        if (hasEnded) { break; }
    }

    if (!visited[end] || !hasEnded) { return nullptr; }
    // O(V)
    for (Vertex at = end; at != -1; at = parent[at]) { path->add(at); }

    delete[] visited;
    delete[] parent;

    return path;
}

// Questão 2
// O(V * (V + E))
List* computePath(AdjList* adjList, Vertex start, Vertex end, List* L)
{
    List* path = new List();
    Node* currentNode = L->head;

    if (start == end || start < 0 || end < 0 || start >= adjList->V || end >= adjList->V || currentNode == nullptr) { return path; }

    Node* nextNode = currentNode->next;

    // O((V + E) * V)
    while (nextNode)
    {
        if (currentNode->vertex == nextNode->vertex) { continue; }

        // O(V + E)
        List* temp = bfs(adjList, currentNode->vertex, nextNode->vertex);
        if (temp == nullptr) { return path; }

        // O(V)
        Node* tempNode = temp->head;
        while (tempNode)
        {
            path->add(tempNode->vertex);
            tempNode = tempNode->next;
        }
        currentNode = nextNode;
        nextNode = nextNode->next;
    }
    
    return path;
}

// O(V)
void initializations(Vertex* parent, int* distance, bool* checked, int numVertices)
{
    for (Vertex v = 0; v < numVertices; v++)
    {
        parent[v] = -1;
        distance[v] = INFPOS;
        checked[v] = false;
    }
}

// O(V)
void seekLowest(const int* distance, const bool* checked, int &minDistance, Vertex &v1, int numVertices)
{
    for (Vertex i = 0; i < numVertices; i++)
    {
        if (checked[i]) continue;
        if (distance[i] < minDistance)
        {
            minDistance = distance[i];
            v1 = i;
        }
    }
}

// O(V^2)
void Dijkstra(WeightedAdjList* adjList, Vertex start, int *distance, Vertex *parent)
{
    if (start < 0 || start >= adjList->V) { return; }

    int numVertices = adjList->V;
    bool checked[numVertices];
    initializations(parent, distance, checked, numVertices);
    parent[start] = start;
    distance[start] = 0;

    // O(V * V)
    while (true)
    {
        int minDistance = INT_MAX;
        Vertex currentVertex = -1;
        seekLowest(distance, checked, minDistance, currentVertex, numVertices);
        if (minDistance == INT_MAX) { break; }

        WeightedNode* adjNode = adjList->adj[currentVertex].head;
        // O(V)
        while (adjNode)
        {
            Vertex adjVertex = adjNode->vertex;

            if (!checked[adjVertex])
            {
                int cost = adjNode->weight;
                if (distance[currentVertex] + cost < distance[adjVertex])
                {
                    parent[adjVertex] = currentVertex;
                    distance[adjVertex] = distance[currentVertex] + cost;
                }
            }
            adjNode = adjNode->next;
        }
        checked[currentVertex] = true;
    }
    delete[] checked;
}

// O((V + E) * log(V))
void heapDijsktra(WeightedAdjList* adjList, Vertex start, int* dist, Vertex* parent)
{
    if (start < 0 || start >= adjList->V) { return; }

    int numVertices = adjList->V;
    bool* visited = new bool[numVertices];
    initializations(parent, dist, visited, numVertices);
    dist[start] = 0;

    MinHeap* heap = new MinHeap(numVertices);
    heap->insert(start, 0);

    // O((V + E) * log(V))
    while (!heap->isEmpty())
    {
        MinHeapNode* fromNode = heap->extractMin();
        int fromVertex = fromNode->v;
        if (dist[fromVertex] == INFPOS) { break; }
        visited[fromVertex] = true;

        WeightedNode* toNode = adjList->adj[fromVertex].head;
        while (toNode)
        {
            int toVertex = toNode->vertex;
            if (!visited[toVertex] && dist[fromVertex] != INFPOS && dist[fromVertex] + toNode->weight < dist[toVertex])
            {
                dist[toVertex] = dist[fromVertex] + toNode->weight;
                parent[toVertex] = fromVertex;
                // O(log(V))
                heap->insertOrUpdate(toVertex, dist[toVertex]);
            }
            toNode = toNode->next;
        }
    }
    delete[] heap;
    delete[] visited;
}

// Questão 3
// O(V^3)
Vertex findMinimaxVertex(WeightedAdjList* adjList, List* L)
{
    int numVertices = adjList->V;
    int* maxDist = new int[numVertices];
    for (int i = 0; i < numVertices; i++) { maxDist[i] = INFPOS; }

    Node* currentNode = L->head;
    // O(V)
    while (currentNode)
    {
        int* dist = new int[numVertices];
        int* parent = new int[numVertices];
        // O(V^2)
        Dijkstra(adjList, currentNode->vertex, dist, parent);
        // O(V)
        for (int i = 0; i < numVertices; i++)
        {
            if (dist[i] < INFPOS && dist[i] > maxDist[i])
            {
                maxDist[i] = dist[i];
            }
        }

        delete[] dist;
        delete[] parent;
        currentNode = currentNode->next;
    }

    Vertex minimaxVertex = -1;
    int minimaxValue = INFPOS;
    // O(V)
    for (int i = 0; i < numVertices; i++)
    {
        if (maxDist[i] < minimaxValue)
        {
            minimaxValue = maxDist[i];
            minimaxVertex = i;
        }
    }
    delete[] maxDist;

    return minimaxVertex;
}

// Questão 4
// O(V^3)
List* findCheapestPath(WeightedAdjList* adjList, List* C, int X)
{
    if (C->head == nullptr) { return nullptr; }

    int numVertices = adjList->V;

    bool* inSubgraph = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { inSubgraph[i] = false; }

    Node* currentNode = C->head;
    int start = currentNode->vertex;
    // O(V^3)
    while (currentNode)
    {
        int* dist = new int[numVertices];
        int* parent = new int[numVertices];
        // O(V^2)
        Dijkstra(adjList, currentNode->vertex, dist, parent);
        // O(V)
        for (int i = 0; i < numVertices; i++)
        {
            if (dist[i] <= X)
            {
                inSubgraph[i] = true;
            }
        }
        delete[] dist;
        delete[] parent;

        currentNode = currentNode->next;
    }
    int end = currentNode->vertex;

    WeightedAdjList* subgraph = new WeightedAdjList(numVertices);
    // O(V + E)
    for (int i = 0; i < numVertices; i++)
    {
        if (inSubgraph[i])
        {
            WeightedNode* adjNode = adjList->adj[i].head;
            while (adjNode)
            {
                if (inSubgraph[adjNode->vertex])
                {
                    subgraph->adj[i].add(adjNode->vertex, adjNode->weight);
                }
                adjNode = adjNode->next;
            }
        }
    }
    delete[] inSubgraph;

    int* dist = new int[numVertices];
    int* parent = new int[numVertices];
    // O(V^2)
    Dijkstra(subgraph, start, dist, parent);

    List* path = new List();
    for (int at = end; at != -1 || at != start; at = parent[at]) { path->add(at); }
    path->add(start);

    return path;
}

// O((V + E) * log(V))
void MST(WeightedAdjList* adjList, int* parent, int* key, Vertex start = 0)
{
    int numVertices = adjList->V;
    bool* inMST = new bool[numVertices];
    // O(V)
    initializations(parent, key, inMST, numVertices);

    key[start] = 0;
    MinHeap* heap = new MinHeap(numVertices);
    heap->insert(start, 0);

    // O((V + E) * log(V))
    while (!heap->isEmpty())
    {
        // O(log(V))
        MinHeapNode* currentNode = heap->extractMin();
        int currentVertex = currentNode->v;
        inMST[currentVertex] = true;

        WeightedNode* adjNode = adjList->adj[currentVertex].head;
        while (adjNode)
        {
            int adjVertex = adjNode->vertex;
            int weight = adjNode->weight;

            if (!inMST[adjVertex] && key[adjVertex] > weight)
            {
                key[adjVertex] = weight;
                // O(log(V))
                heap->insertOrUpdate(adjVertex, key[adjVertex]);
                parent[adjVertex] = currentVertex;
            }
            adjNode = adjNode->next;
        }
    }

    delete[] inMST;
    delete[] heap;
}

// O((V + E) * log(V))
int computeMSTCost(WeightedAdjList* adjList)
{
    int numVertices = adjList->V;
    int* parent = new int[numVertices];
    int* key = new int[numVertices];
    // O((V + E) * log(V))
    MST(adjList, parent, key);

    int mstCost = 0;
    // O(V)
    for (int i = 0; i < numVertices; i++)
    {
        if (parent[i] != -1)
        {
            mstCost += key[i];
        }
    }
    delete[] key;
    delete[] parent;

    return mstCost;
}

// Questão 5
// O(E * (V + E) * log(V))
pair<Vertex, Vertex> findCriticalEdge(WeightedAdjList* adjList)
{
    int numVertices = adjList->V;
    int maxIncrease = 0;
    pair<Vertex, Vertex> criticalEdge = {-1, -1};

    // O(E * (V + E) * log(V))
    for (int i = 0; i < numVertices; i++)
    {
        WeightedNode* adjNode = adjList->adj[i].head;
        while (adjNode)
        {
            int currentVertex = i;
            int adjVertex = adjNode->vertex;
            int weight = adjNode->weight;

            adjList->removeEdge(currentVertex, adjVertex);
            adjList->removeEdge(adjVertex, currentVertex);
            
            // O((V + E) * log(V))
            int newMSTCost = computeMSTCost(adjList);

            adjList->addEdge(currentVertex, adjVertex, weight);
            adjList->addEdge(adjVertex, currentVertex, weight);

            // O((V + E) * log(V))
            int increase = newMSTCost - computeMSTCost(adjList);
            if (increase > maxIncrease)
            {
                maxIncrease = increase;
                criticalEdge = {currentVertex, adjVertex};
            }
            adjNode = adjNode->next;
        }
    }

    return criticalEdge;
}

int main()
{

    return 0;
}