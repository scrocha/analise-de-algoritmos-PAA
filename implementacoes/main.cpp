#include <iostream>
#include "structs.hpp"
#include <stack>

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


int main()
{

    return 0;
}