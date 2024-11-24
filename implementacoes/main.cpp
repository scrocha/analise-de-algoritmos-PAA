#include "structs.hpp"
#include <iostream>
#include <stack>
#include <vector>
#include <queue>
#include <limits>

using namespace std;

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
    List* path = new List();

    if (start == end || start < 0 || end < 0 || start >= adjList->V || end >= adjList->V) { return nullptr; }

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




int main()
{

    return 0;
}