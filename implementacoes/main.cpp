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

List* bfs(AdjList* adjList, Vertex start, Vertex end)
{
    if (start == end || start < 0 || end < 0 || start >= adjList->V || end >= adjList->V) { return nullptr; }

    List* path = new List();
    bool* visited = new bool[adjList->V];
    int* parent = new int[adjList->V];

    for (int i = 0; i < adjList->V; i++)
    {
        visited[i] = false;
        parent[i] = -1;
    }

    queue<Vertex> queue;
    queue.push(start);
    visited[start] = true;

    bool hasEnded = false;
    while (!queue.empty())
    {
        Vertex currentVertex = queue.front();
        queue.pop();

        Node* temp = adjList->adj[currentVertex].head;
        while (temp)
        {
            if (!visited[temp->vertex])
            {
                queue.push(temp->vertex);
                visited[temp->vertex] = true;
                parent[temp->vertex] = currentVertex;

                if (temp->vertex == end)
                {
                    hasEnded = true;
                    break;
                }
            }
            temp = temp->next;
        }
        if (hasEnded) { break; }
    }

    if (!visited[end] || !hasEnded) { return nullptr; }

    for (Vertex at = end; at != -1; at = parent[at]) { path->add(at); }

    delete[] visited;
    delete[] parent;

    return path;
}

List* computePath(AdjList* adjList, Vertex start, Vertex end, List* L)
{
    List* path = new List();
    Node* current = L->head;

    if (start == end || start < 0 || end < 0 || start >= adjList->V || end >= adjList->V || current == nullptr) { return path; }

    Node* next = current->next;

    while (next)
    {
        if (current->vertex == next->vertex) { continue; }

        List* temp = bfs(adjList, current->vertex, next->vertex);
        if (temp == nullptr) { return path; }

        Node* tempNode = temp->head;
        while (tempNode)
        {
            path->add(tempNode->vertex);
            tempNode = tempNode->next;
        }
        current = next;
        next = next->next;
    }
    
    return path;
}

void Djisktra(WeightedAdjList* adjList, Vertex start, int* dist, int* parent)
{
    if (start < 0 || start >= adjList->V) { return; }
    
    bool* visited = new bool[adjList->V];
    for (int i = 0; i < adjList->V; i++)
    {
        dist[i] = INFPOS;
        parent[i] = -1;
        visited[i] = false;
    }
    dist[start] = 0;

    MinHeap* heap = new MinHeap(adjList->V);
    heap->insert(start, 0);

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
                heap->insertOrUpdate(toVertex, dist[toVertex]);
            }
            toNode = toNode->next;
        }
    }
    delete[] heap;
    delete[] visited;
}

Vertex findMinimaxVertex(WeightedAdjList* adjList, List* L)
{
    int numVertices = adjList->V;
    int* maxDist = new int[numVertices];
    for (int i = 0; i < numVertices; i++) { maxDist[i] = INFPOS; }

    Node* temp = L->head;
    while (temp)
    {
        int* dist = new int[numVertices];
        int* parent = new int[numVertices];
        Djisktra(adjList, temp->vertex, dist, parent);

        for (int i = 0; i < numVertices; i++)
        {
            if (dist[i] < INFPOS && dist[i] > maxDist[i])
            {
                maxDist[i] = dist[i];
            }
        }

        delete[] dist;
        delete[] parent;
        temp = temp->next;
    }

    Vertex minimaxVertex = -1;
    int minimaxValue = INFPOS;
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

List* findCheapestPath(WeightedAdjList* adjList, List* C, int X)
{
    if (C->head == nullptr) { return nullptr; }

    int numVertices = adjList->V;

    bool* inSubgraph = new bool[numVertices];
    for (int i = 0; i < numVertices; i++) { inSubgraph[i] = false; }

    Node* temp = C->head;
    int start = temp->vertex;
    while (temp)
    {
        int* dist = new int[numVertices];
        int* parent = new int[numVertices];
        Djisktra(adjList, temp->vertex, dist, parent);

        for (int i = 0; i < numVertices; i++)
        {
            if (dist[i] <= X)
            {
                inSubgraph[i] = true;
            }
        }
        delete[] dist;
        delete[] parent;

        temp = temp->next;
    }
    int end = temp->vertex;

    WeightedAdjList* subgraph = new WeightedAdjList(numVertices);

    for (int i = 0; i < numVertices; i++)
    {
        if (inSubgraph[i])
        {
            WeightedNode* temp = adjList->adj[i].head;
            while (temp)
            {
                if (inSubgraph[temp->vertex])
                {
                    subgraph->adj[i].add(temp->vertex, temp->weight);
                }
                temp = temp->next;
            }
        }
    }
    delete[] inSubgraph;

    int* dist = new int[numVertices];
    int* parent = new int[numVertices];
    Djisktra(subgraph, start, dist, parent);

    List* path = new List();
    for (int at = end; at != -1 || at != start; at = parent[at]) { path->add(at); }
    path->add(start);

    return path;
}

int main()
{

    return 0;
}