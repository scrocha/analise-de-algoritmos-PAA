#include "structs.hpp"
#include <iostream>
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

    QueueNode* queue = new QueueNode(start);
    visited[start] = true;
    // O(V + E)
    while (!queue->empty())
    {
        Vertex currentVertex = queue->pop();

        Node* adjNode = adjList->adj[currentVertex].head;
        while (adjNode)
        {
            if (!visited[adjNode->vertex])
            {
                queue->push(adjNode->vertex);
                visited[adjNode->vertex] = true;
            }
            adjNode = adjNode->next;
        }
    }
    delete[] queue;
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
    delete[] descendents;
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
    delete[] noDescendents;
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
    delete[] inL2;

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
        delete[] visited;
    }
}

// O(V * (V + E))
void dfsAncestors(AdjList* adjList, List* L2, bool* inL2)
{
    int numVertices = adjList->V;
    bool* noAncestors = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { noAncestors[i] = false; }

    // O(V * (V + E))
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
    delete[] noAncestors;
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
    delete[] inL2;

    return L2;
}

int findSons(AdjList* adjList, Vertex root, Vertex* parent)
{
    int numVertices = adjList->V;

    QueueNode* queue = new QueueNode(root);
    parent[root] = root;
    int count = 0;

    while (!queue->empty())
    {
        Vertex currentVertex = queue->pop();
        Node* adjNode = adjList->adj[currentVertex].head;
        while (adjNode)
        {
            if (parent[adjNode->vertex] == -1)
            {
                parent[adjNode->vertex] = currentVertex;
                queue->push(adjNode->vertex);
                count++;
            }
            adjNode = adjNode->next;
        }
    }
    delete[] queue;
    return count;
}

Vertex findTreeParents(AdjList* adjList, Vertex* parent)
{
    int numVertices = adjList->V;
    bool* hasParent = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { parent[i] = -1; hasParent[i] = false; }

    Vertex root = -1;
    while (true)
    {
        bool foundRoot = false;
        for (int i = 0; i < numVertices; i++)
        {
            if (hasParent[i]) { continue; }

            int count = findSons(adjList, i, parent);
            if (count >= numVertices - 2) { root = i; foundRoot = true; break; }
            else { hasParent[i] = true; parent[i] = -1; }
        }
        if (foundRoot) { break; }
    }
    delete[] hasParent;
    return root;
}

void olderCousinCount(AdjList* adjList, List* L2, int* count)
{
    int numVertices = adjList->V;
    Vertex* parent = new Vertex[numVertices];
    int* layer = new int[numVertices];
    for (int i = 0; i < numVertices; i++) { count[i] = 0; layer[i] = INT_MIN ;}

    Vertex root = findTreeParents(adjList, parent);

    layer[root] = 0;
    QueueNode* queue = new QueueNode(root);

    while (!queue->empty())
    {
        Vertex currentVertex = queue->pop();
        Node* adjNode = adjList->adj[currentVertex].head;
        while (adjNode)
        {
            layer[adjNode->vertex] = layer[currentVertex] + 1;
            queue->push(adjNode->vertex);
            adjNode = adjNode->next;
        }
    }

    for (int i = 0; i < numVertices; i++)
    {
        if (layer[i] == 0) { continue; }

        Node* temp = L2->head;
        while (temp)
        {
            if (layer[temp->vertex] == layer[i] - 1)
            {
                count[i]++;
            }
            temp = temp->next;
        }
    }

}
/*

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
*/
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

    QueueNode* queue = new QueueNode(start);
    visited[start] = true;

    bool hasEnded = false;
    // O(V + E)
    while (!queue->empty())
    {
        // O(1)
        Vertex currentVertex = queue->pop();

        Node* adjNode = adjList->adj[currentVertex].head;
        while (adjNode)
        {
            if (!visited[adjNode->vertex])
            {
                queue->push(adjNode->vertex);
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
    for (Vertex at = end; at != -1 || at != start ; at = parent[at]) { path->add(at); }
    path->add(start);

    delete[] visited;
    delete[] parent;
    delete[] queue;

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
    path->add(start);
    // O((V + E) * V)
    while (nextNode)
    {
        if (currentNode->vertex == nextNode->vertex) { continue; }

        // O(V + E)
        List* temp = bfs(adjList, currentNode->vertex, nextNode->vertex);
        if (temp == nullptr) { return path; }

        // O(V)
        Node* tempNode = temp->head->next;
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
    bool* checked = new bool[numVertices];
    initializations(parent, distance, checked, numVertices);
    parent[start] = start;
    distance[start] = 0;

    // O(V * V)
    while (true)
    {
        int minDistance = INT_MAX;
        Vertex currentVertex = -1;
        // O(V)
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
    bool* inList = new bool[numVertices];
    for (int i = 0; i < numVertices; i++)
    {
        maxDist[i] = INFPOS;
        inList[i] = false;
    }

    Node* tempNode = L->head;
    while (tempNode)
    {
        inList[tempNode->vertex] = true;
        tempNode = tempNode->next;
    }

    Node* currentNode = L->head;
    // O(V * V^2)
    while (currentNode)
    {
        int* dist = new int[numVertices];
        int* parent = new int[numVertices];
        // O(V^2)
        Dijkstra(adjList, currentNode->vertex, dist, parent);
        // O(V)
        for (int i = 0; i < numVertices; i++)
        {
            if (inList[i] || dist[i] >= INFPOS) { continue; }
            else if (dist[i] < maxDist[i]) { maxDist[i] = dist[i]; }
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
        if (!inList[i] && minimaxValue < maxDist[i])
        {
            minimaxValue = maxDist[i];
            minimaxVertex = i;
        }
    }
    delete[] maxDist;
    delete[] inList;

    return minimaxVertex;
}

// Questão 4
// O(V^3)
void findValidVertex(WeightedAdjList* adjList, List* C, int X, bool* valids, bool* inC)
{
    int numVertices = adjList->V;
    // O(V^3)
    for (int i = 0; i < numVertices; i++)
    {
        if (inC[i] || valids[i]) { continue; }

        int* dist = new int[numVertices];
        Vertex* parent = new Vertex[numVertices];
        // O(V^2)
        Dijkstra(adjList, i, dist, parent);

        int minDist = INFPOS;
        Vertex nearest = -1;
        // O(V)
        for (int j = 0; j < numVertices; j++)
        {
            if (inC[j] && dist[j] < minDist) { minDist = dist[j]; Vertex nearest = j; }
        }

        if (minDist <= X)
        {
            valids[i] = true;
            for (Vertex j = nearest; j != -1 || j != i; j = parent[j]) { valids[j] = true; }
        }

        delete[] dist;
        delete[] parent;
    }
    delete[] inC;
}
// O(V^3)
List* findCheapestPath(WeightedAdjList* adjList, List* C, int X)
{
    List* path = new List();
    int numVertices = adjList->V;
    bool* valids = new bool[numVertices];
    bool* inC = new bool[numVertices];
    // O(V)
    for (int i = 0; i < numVertices; i++) { valids[i] = false; inC[i] = false; }
    
    Node* temp = C->head;
    Vertex start = temp->vertex;
    while (temp)
    {
        inC[temp->vertex] = true;
        valids[temp->vertex] = true;
        temp = temp->next;
    }
    Vertex end = temp->vertex;
    // O(V^3)
    findValidVertex(adjList, C, X, valids, inC);

    WeightedAdjList* subgraph = new WeightedAdjList(numVertices);
    // O(V + E)
    for (int i = 0; i < numVertices; i++)
    {
        if (!valids[i]) { continue; }

        WeightedNode* adjNode = adjList->adj[i].head;
        while (adjNode)
        {
            if (valids[adjNode->vertex])
            {
                subgraph->addEdge(i, adjNode->vertex, adjNode->weight);
            }
            adjNode = adjNode->next;
        }
    }

    int* dist = new int[numVertices];
    Vertex* parent = new Vertex[numVertices];
    // O(V^2)
    Dijkstra(subgraph, 0, dist, parent);

    for (Vertex i = end; i != -1 || i != start; i = parent[i]) { path->add(i); }
    path->add(start);

    delete[] valids;
    delete[] inC;
    delete[] dist;
    delete[] parent;
    delete[] subgraph;

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

// O(V + E)
pair<pair<Vertex, Vertex>, int> findMinEdge(WeightedAdjList* adjList, AdjList* mst, Vertex start, Vertex end)
{
    int numVertices = adjList->V;
    int minWeight = INFPOS;
    Vertex minEdge1 = -1;
    Vertex minEdge2 = -1;

    // O(V)
    mst->removeEdgeDual(start, end);
    bool* component1 = new bool[numVertices];
    // O(V + E)
    canReach(mst, start, component1);
    mst->addEdgeDual(start, end);

    // O(V + E)
    for (int i = 0; i < numVertices; i++)
    {
        if (component1[i])
        {
            WeightedNode* adjNode = adjList->adj[i].head;
            while (adjNode)
            {
                if (!component1[adjNode->vertex] && adjNode->weight < minWeight)
                {
                    minWeight = adjNode->weight;
                    minEdge1 = i;
                    minEdge2 = adjNode->vertex;
                }
                adjNode = adjNode->next;
            }
        }
    }
    delete[] component1;

    pair<pair<Vertex, Vertex>, int> minEdgePair = {{minEdge1, minEdge2}, minWeight};

    return minEdgePair;
}

// Questão 5
// O((V + E) * (log(V) + V))
pair<Vertex, Vertex> findCriticalEdge(WeightedAdjList* adjList)
{
    int numVertices = adjList->V;
    int maxIncrease = 0;
    pair<Vertex, Vertex> criticalEdge = {-1, -1};

    int* parent = new int[numVertices];
    int* key = new int[numVertices];
    // O((V + E) * log(V))
    MST(adjList, parent, key);
    
    AdjList* mst = new AdjList(numVertices);
    // O(V)
    for (int i = 0; i < numVertices; i++)
    {
        if (parent[i] != -1) {  mst->addEdgeDual(i, parent[i]); }
    }
    // O(V * (V + E))
    for (int i = 0; i < numVertices; i++)
    {
        if (parent[i] == -1) { continue; }
        int originalWeight = key[i];

        adjList->removeEdgeDual(i, parent[i]);
        // O(V + E)
        pair<pair<Vertex, Vertex>, int> newMinEdge = findMinEdge(adjList, mst, i, parent[i]);

        adjList->addEdgeDual(i, parent[i], originalWeight);

        int increase = newMinEdge.second; - originalWeight;
        if (increase > maxIncrease)
        {
            maxIncrease = increase;
            criticalEdge = {i, parent[i]};
        }
    }
    delete[] mst;
    delete[] parent;
    delete[] key;

    return criticalEdge;
}

int main()
{

    return 0;
}