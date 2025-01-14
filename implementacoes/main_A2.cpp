#include "structs.hpp"
#include <iostream>
#include <climits>

using namespace std;

#define INFPOS INT_MAX

// O(V + E)
Vertex findRoot(AdjList* adjList, Vertex* parent)
{
    int numVertices = adjList->V;
    bool* hasParent = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { parent[i] = -1; hasParent[i] = false; }

    Vertex root = -1;
    // O(V + E)
    for (int i = 0; i < numVertices; i++)
    {
        Node* adjNode = adjList->adj[i].head;
        while (adjNode)
        {
            parent[adjNode->vertex] = i;
            hasParent[adjNode->vertex] = true;
            adjNode = adjNode->next;
        }
    }
    // O(V)
    for (int i = 0; i < numVertices; i++)
    {
        if (!hasParent[i]) { root = i; break; }
    }
    delete[] hasParent;

    return root;
}

// O(V + E)
List* computeL2_a(AdjList* adjList, List* L1)
{
    int numVertices = adjList->V;
    List* L2 = new List();
    bool* inL2 = new bool[numVertices];
    bool* hasChild = new bool[numVertices];
    bool* visited = new bool[numVertices];
    int* descendentsCount = new int[numVertices];
    // O(V)
    for (int i = 0; i < numVertices; i++)
    {
        inL2[i] = false;
        hasChild[i] = false;
        visited[i] = false;
        descendentsCount[i] = 0;
    }

    Node* current = L1->head;
    // O(V)
    while (current)
    {
        L2->add(current->vertex);
        inL2[current->vertex] = true;
        descendentsCount[current->vertex]++;
        current = current->next;
    }

    Vertex* parent = new Vertex[numVertices];
    // O(V + E)
    Vertex root = findRoot(adjList, parent);

    for (int i = 0; i < numVertices; i++)
    {
        if (parent[i] != -1) { hasChild[parent[i]] = true; }
    }
    hasChild[root] = true;

    QueueNode* queue = new QueueNode(-1);
    // O(V)
    for (int i = 0; i < numVertices; i++) { if (!hasChild[i]) { queue->push(i); } }
    queue->pop();
    // O(V + E)
    while (!queue->empty())
    {
        Vertex currentVertex = queue->pop();
        visited[currentVertex] = true;

        if (!inL2[currentVertex] && descendentsCount[currentVertex] % 2 == 1)
        {
            L2->add(currentVertex);
            inL2[currentVertex] = true;
            descendentsCount[currentVertex]++;
        }

        Vertex parentVertex = parent[currentVertex];

        if (parentVertex == -1) { break; }

        descendentsCount[parentVertex] += descendentsCount[currentVertex];
        if (!visited[parentVertex]) { queue->push(parentVertex); }
    }
    delete[] descendentsCount;
    delete[] parent;
    delete[] inL2;
    delete[] hasChild;
    delete[] visited;
    delete[] queue;

    return L2;
}

// O(V + E)
List* computeL2_b(AdjList* adjList, List* L1)
{
    int numVertices = adjList->V;
    List* L2 = new List();
    bool* inL2 = new bool[numVertices];
    int* ancestorsCount = new int[numVertices];
    // O(V)
    for (int i = 0; i < numVertices; i++) { inL2[i] = false; ancestorsCount[i] = 0; }

    Node* current = L1->head;
    while (current)
    {
        L2->add(current->vertex);
        inL2[current->vertex] = true;
        ancestorsCount[current->vertex]++;
        current = current->next;
    }

    Vertex* parent = new Vertex[numVertices];
    // O(V + E)
    Vertex root = findRoot(adjList, parent);
    QueueNode* queue = new QueueNode(root);

    // O(V + E)
    while (!queue->empty())
    {
        Vertex currentVertex = queue->pop();
        Node* adjNode = adjList->adj[currentVertex].head;

        while (adjNode)
        {
            ancestorsCount[adjNode->vertex] += ancestorsCount[currentVertex];
            queue->push(adjNode->vertex);

            if (!inL2[adjNode->vertex] && ancestorsCount[adjNode->vertex] % 2 == 1)
            {
                L2->add(adjNode->vertex);
                inL2[adjNode->vertex] = true;
                ancestorsCount[adjNode->vertex]++;
            }
            adjNode = adjNode->next;
        }
    }
    delete[] ancestorsCount;
    delete[] parent;
    delete[] inL2;
    delete[] queue;

    return L2;
}

// O(V + E)
void dfsPrePostOrder(AdjList* adjList, Vertex vertex, bool* visited, int* preorder, int* postorder, int& preIndex, int& postIndex)
{
    visited[vertex] = true;

    preorder[vertex] = preIndex++;

    Node* child = adjList->adj[vertex].head;
    while (child)
    {
        if (!visited[child->vertex])
        {
            dfsPrePostOrder(adjList, child->vertex, visited, preorder, postorder, preIndex, postIndex);
        }
        child = child->next;
    }

    postorder[vertex] = postIndex++;
}

// O(V + E)
Vertex computePrePostOrder(AdjList* adjList, int* preorder, int* postorder, int* parent)
{
    int numVertices = adjList->V;
    bool* visited = new bool[numVertices];

    for (int i = 0; i < numVertices; i++)
    {
        visited[i] = false;
        preorder[i] = -1;
        postorder[i] = -1;
        parent[i] = -1;
    }
    // O(V + E)
    Vertex root = findRoot(adjList, parent);

    int preIndex = 0;
    int postIndex = 0;
    int layerCount = 0;
    // O(V + E)
    dfsPrePostOrder(adjList, root, visited, preorder, postorder, preIndex, postIndex);

    delete[] visited;
    return root;
}

// O(V + E)
List* computeL2_c(AdjList* adjList, List* L1)
{
    int numVertices = adjList->V;
    List* L2 = new List();
    bool* inL2 = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { inL2[i] = false; }

    Node* current = L1->head;
    while (current)
    {
        L2->add(current->vertex);
        inL2[current->vertex] = true;
        current = current->next;
    }

    int* preorder = new int[numVertices];
    int* postorder = new int[numVertices];
    int* parent = new int[numVertices];
    // O(V + E)
    computePrePostOrder(adjList, preorder, postorder, parent);
    // O(V + E)
    while (true)
    {
        bool added = false;
        for (int i = 0; i < numVertices; i++)
        {
            if (inL2[i]) { continue; }
            int count = 0;

            Node* current = L2->head;
            while (current)
            {
                if (preorder[current->vertex] < preorder[i]) { count++; }
                current = current->next;
            }

            if (count % 2 == 1)
            {
                L2->add(i);
                inL2[i] = true;
                added = true;
            }
        }
        if (!added) { break; }
    }
    delete[] inL2;
    delete[] parent;
    delete[] preorder;
    delete[] postorder;

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

// O(V^2)
void DijkstraMultiSource(WeightedAdjList* adjList, List* L, int *distance, Vertex *parent)
{
    int numVertices = adjList->V;
    bool* checked = new bool[numVertices];
    initializations(parent, distance, checked, numVertices);
    
    Node* current = L->head;
    while (current)
    {
        distance[current->vertex] = 0;
        parent[current->vertex] = current->vertex;

        current = current->next;
    }

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
void heapDijsktraMultiSource(WeightedAdjList* adjList, List* L, int* dist, Vertex* parent, Vertex& endVertex)
{
    int numVertices = adjList->V;
    bool* visited = new bool[numVertices];
    initializations(parent, dist, visited, numVertices);
    
    MinHeap* heap = new MinHeap(numVertices);
    Node* current = L->head;
    // O(V)
    while (current)
    {
        dist[current->vertex] = 0;
        parent[current->vertex] = current->vertex;
        heap->insert(current->vertex, 0);
        current = current->next;
    }
    endVertex = current->vertex;

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

    Vertex minimaxVertex = -1;
    int minimaxValue = INFPOS;
    // O(V^3)
    for (int i = 0; i < numVertices; i++)
    {
        if (!L->contains(i)) { continue; }
        int* dist = new int[numVertices];
        Vertex* parent = new Vertex[numVertices];
        // O(V^2)
        Dijkstra(adjList, i, dist, parent);

        int maxDist = 0;
        Node* current = L->head;
        while (current)
        {
            if (dist[current->vertex] > maxDist) { maxDist = dist[current->vertex]; }
            current = current->next;
        }

        if (maxDist < minimaxValue)
        {
            minimaxValue = maxDist;
            minimaxVertex = i;
        }
        delete[] dist;
        delete[] parent;
    }

    return minimaxVertex;
}

// Questão 4
// O((V + E) * log(V))
List* findCheapestPath(WeightedAdjList* adjList, List* C, int X)
{
    int numVertices = adjList->V;
    int* dist = new int[numVertices];
    Vertex* parent = new Vertex[numVertices];
    Vertex startVertex = C->head->vertex;
    Vertex endVertex = -1;
    // O((V + E) * log(V))
    heapDijsktraMultiSource(adjList, C, dist, parent, endVertex);

    WeightedAdjList* adjListCopy = new WeightedAdjList(numVertices);

    // O(V + E)
    for (int i = 0; i < numVertices; i++)
    {
        if (dist[i] > X) { continue; }
        WeightedNode* adjNode = adjList->adj[i].head;
        while (adjNode)
        {
            if (dist[adjNode->vertex] > X) { adjNode = adjNode->next; continue; }
            
            adjListCopy->addEdge(i, adjNode->vertex, adjNode->weight);
            adjNode = adjNode->next;
        }
    }

    List* path = new List();
    // O((V + E) * log(V))
    heapDijsktra(adjListCopy, startVertex, dist, parent);

    for (Vertex at = endVertex; at != -1 || at != startVertex; at = parent[at]) { path->add(at); }
    path->add(startVertex);

    return path;
}

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
