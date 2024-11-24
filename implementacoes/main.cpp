#include "structs.hpp"
#include <iostream>
#include <stack>
#include <vector>
#include <queue>
#include <climits>

using namespace std;

#define INFPOS INT_MAX

int descendentCount(AdjList* adjList, int vertex, List* L2)
{
    int count = 0;
    bool* visited = new bool[adjList->V];
    for (int i = 0; i < adjList->V; i++) { visited[i] = false; }
    
    stack<int> stack;
    stack.push(vertex);
    visited[vertex] = true;

    while (!stack.empty())
    {
        int currentVertex = stack.top();
        stack.pop();
        
        Node* temp = adjList->adj[currentVertex].head;
        while (temp)
        {
            if (!visited[temp->vertex])
            {
                stack.push(temp->vertex);
                visited[temp->vertex] = true;

                if (L2->contains(temp->vertex)) { count++; }
            }
            
            temp = temp->next;
        }
    }
    delete[] visited;
    return count;
}

void dfsDescendants(AdjList* adjList, List* L2)
{
    int numVertices = adjList->V;
    bool* inL2 = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { inL2[i] = false; }
    
    while (true)
    {
        bool added = false;
        for (int i = 0; i < numVertices; i++)
        {
            cout << "i: " << i << "| desc: " << descendentCount(adjList, i, L2) << endl;
            if (!inL2[i] && descendentCount(adjList, i, L2) % 2 != 0)
            {
                if (!L2->contains(i))
                {
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

// Questão 1 a
List* computeL2_a(AdjList* adjList, List* L1)
{
    List* L2 = new List();
    
    Node* temp = L1->head;
    while (temp)
    {
        L2->add(temp->vertex);
        temp = temp->next;
    }
    if (L2->head == nullptr) { return L2; }

    dfsDescendants(adjList, L2);

    return L2;
}


int ancestorCount(AdjList* adjList, int vertex, List* L2)
{
    int count = 0;

    Node* temp = L2->head;
    while (temp)
    {
        int current = temp->vertex;
        stack<int> stack;

        bool* visited = new bool[adjList->V];
        for (int i = 0; i < adjList->V; i++) visited[i] = false;

        stack.push(current);
        visited[current] = true;

        while (!stack.empty())
        {
            int currentVertex = stack.top();
            stack.pop();

            Node* temp2 = adjList->adj[currentVertex].head;
            while (temp2)
            {
                if (temp2->vertex == vertex)
                {
                    count++;
                    break;
                }
                if (!visited[temp2->vertex])
                {
                    stack.push(temp2->vertex);
                    visited[temp2->vertex] = true;
                }
                temp2 = temp2->next;
            }
        }

        delete[] visited;
        temp = temp->next;
    }
    return count;
}

void dfsAncestors(AdjList* adjList, List* L2)
{
    int numVertices = adjList->V;
    bool* inL2 = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { inL2[i] = false; }
    
    while (true)
    {
        bool added = false;
        for (int i = 0; i < numVertices; i++)
        {
            cout << "i: " << i << "| anc: " << ancestorCount(adjList, i, L2) << endl;
            if (!inL2[i] && ancestorCount(adjList, i, L2) % 2 != 0)
            {
                if (!L2->contains(i))
                {
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

// Questão 1 b
List* computeL2_b(AdjList* adjList, List* L1)
{
    List* L2 = new List();
    
    Node* temp = L1->head;
    while (temp)
    {
        L2->add(temp->vertex);
        temp = temp->next;
    }
    if (L2->head == nullptr) { return L2; }

    dfsAncestors(adjList, L2);

    return L2;
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