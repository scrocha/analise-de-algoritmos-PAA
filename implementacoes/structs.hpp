#include <iostream>

#ifndef STRUCTS_H
#define STRUCTS_H

using namespace std;

#define Vertex int
struct Node
{
    Vertex vertex;
    Node* next;

    Node(Vertex val) : vertex(val), next(nullptr) {}
};

struct List
{
    Node* head;

    List() : head(nullptr) {}

    void add(Vertex vertex)
    {
        Node* newNode = new Node(vertex);
        newNode->next = head;
        head = newNode;
    }

    ~List()
    {
        Node* current = head;
        while (current != nullptr)
        {
            Node* next = current->next;
            delete current;
            current = next;
        }
    }

    void print()
    {
        Node* current = head;
        cout << "[";
        while (current != nullptr)
        {
            cout << " " << current->vertex << ",";
            current = current->next;
        }
        cout << "]" << endl;
    }

    bool contains(Vertex vertex)
    {
        Node* current = head;
        while (current != nullptr)
        {
            if (current->vertex == vertex) { return true; }
            current = current->next;
        }
        return false;
    }

    void remove(Vertex vertex)
    {
        Node* current = head;
        Node* previous = nullptr;
        while (current != nullptr)
        {
            if (current->vertex == vertex)
            {
                if (previous == nullptr) { head = current->next; }
                else { previous->next = current->next; }
                delete current;
                return;
            }
            previous = current;
            current = current->next;
        }
    }
};

struct AdjList
{
    List* adj;
    int V;
    int E;

    AdjList(int V) : V(V), E(0) { adj = new List[V]; }

    ~AdjList()
    {
        delete[] adj;
    }

    void addEdge(Vertex v, Vertex w)
    {
        adj[v].add(w);
        E++;
    }

    void addEdgeDual(Vertex v, Vertex w)
    {
        adj[v].add(w);
        adj[w].add(v);
        E += 2;
    }

    bool containsEdge(Vertex v, Vertex w)
    {
        return adj[v].contains(w);
    }

    void print()
    {
        for (int v = 0; v < V; v++)
        {
            cout << v << ": ";
            adj[v].print();
        }
    }

    void removeEdge(Vertex v, Vertex w)
    {
        adj[v].remove(w);
        E--;
    }

    void removeEdgeDual(Vertex v, Vertex w)
    {
        adj[v].remove(w);
        adj[w].remove(v);
        E -= 2;
    }
};

struct Edge
{
    Vertex v;
    Vertex w;

    Edge(Vertex v, Vertex w) : v(v), w(w) {}
};

struct WeightedEdge
{
    Vertex v;
    Vertex w;
    double weight;

    WeightedEdge(Vertex v, Vertex w, double weight) : v(v), w(w), weight(weight) {}
};

#endif