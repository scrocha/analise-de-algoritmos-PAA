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


int main()
{

    return 0;
}