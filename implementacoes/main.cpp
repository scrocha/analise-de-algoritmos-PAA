#include <iostream>
using namespace std;

struct Node
{
    int vertex;
    Node* next;

    Node(int val) : vertex(val), next(nullptr) {}
};

void addEdge(Node* adjacencyList[], int v, int w)
{
    Node* newNodeV = new Node(w);
    if (adjacencyList[v] == nullptr) { adjacencyList[v] = newNodeV; }
    else
    {
        Node* temp = adjacencyList[v];
        while (temp->next != nullptr) { temp = temp->next; }
        temp->next = newNodeV;
    }

    Node* newNodeW = new Node(v);
    if (adjacencyList[w] == nullptr) { adjacencyList[w] = newNodeW; }
    else
    {
        Node* temp = adjacencyList[w];
        while (temp->next != nullptr) { temp = temp->next; }
        temp->next = newNodeW;
    }
}

void displayGraph(Node* adjacencyList[], int numVertices)
{
    for (int i = 0; i < numVertices; i++)
    {
        cout << "Vértice " << i << ": ";
        Node* temp = adjacencyList[i];
        while (temp != nullptr)
        {
            cout << temp->vertex << " ";
            temp = temp->next;
        }
        cout << endl;
    }
}

void freeGraph(Node* adjList[], int numVertices)
{
    for (int i = 0; i < numVertices; i++)
    {
        Node* temp = adjList[i];
        while (temp != nullptr)
        {
            Node* nextNode = temp->next;
            delete temp;
            temp = nextNode;
        }
    }
}

int main()
{
    int numVertices = 4;
    Node* adjList[numVertices];

    for (int i = 0; i < numVertices; i++)
    {
        adjList[i] = nullptr;
    }

    addEdge(adjList, 0, 1);
    addEdge(adjList, 0, 2);
    addEdge(adjList, 1, 3);
    addEdge(adjList, 2, 3);

    cout << "Grafo (listas de adjacência):" << endl;
    displayGraph(adjList, numVertices);

    freeGraph(adjList, numVertices);

    return 0;
}
