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

struct WeightedNode
{
    Vertex vertex;
    int weight;
    WeightedNode* next;

    WeightedNode(Vertex val, int w) : vertex(val), weight(w), next(nullptr) {}
};

struct WeightedList
{
    WeightedNode* head;

    WeightedList() : head(nullptr) {}

    void add(Vertex vertex, int weight)
    {
        WeightedNode* newNode = new WeightedNode(vertex, weight);
        newNode->next = head;
        head = newNode;
    }

    ~WeightedList()
    {
        WeightedNode* current = head;
        while (current != nullptr)
        {
            WeightedNode* next = current->next;
            delete current;
            current = next;
        }
    }

    void print()
    {
        WeightedNode* current = head;
        cout << "[";
        while (current != nullptr)
        {
            cout << " " << current->vertex << " (" << current->weight << "),";
            current = current->next;
        }
        cout << "]" << endl;
    }

    bool contains(Vertex vertex)
    {
        WeightedNode* current = head;
        while (current != nullptr)
        {
            if (current->vertex == vertex) { return true; }
            current = current->next;
        }
        return false;
    }

    void remove(Vertex vertex)
    {
        WeightedNode* current = head;
        WeightedNode* previous = nullptr;
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

struct WeightedAdjList
{
    WeightedList* adj;
    int V;
    int E;

    WeightedAdjList(int V) : V(V), E(0) { adj = new WeightedList[V]; }

    ~WeightedAdjList()
    {
        delete[] adj;
    }

    void addEdge(Vertex v, Vertex w, int weight)
    {
        adj[v].add(w, weight);
        E++;
    }

    void addEdgeDual(Vertex v, Vertex w, int weight)
    {
        adj[v].add(w, weight);
        adj[w].add(v, weight);
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

struct QueueNode
{
    Vertex vertex;
    QueueNode* next;

    QueueNode(Vertex vertex) : vertex(vertex), next(nullptr) {}

    ~QueueNode() { delete next; }
    // O(n)
    void push(Vertex vertex)
    {
        if (next == nullptr) { next = new QueueNode(vertex); }
        else { next->push(vertex); }
    }
    // O(1)
    bool empty() const { return next == nullptr; }
    // O(1)
    Vertex pop()
    {
        Vertex vertex = next->vertex;

        QueueNode* temp = next;
        next = next->next;
        temp->next = nullptr;
        delete temp;

        return vertex;
    }
    // O(1)
    Vertex front() const { return next->vertex; }
    // O(n)
    Vertex back() const
    {
        QueueNode* current = next;
        while (current->next != nullptr) { current = current->next; }
        return current->vertex;
    }
};

struct MinHeapNode
{
    Vertex v;
    int dist;

    MinHeapNode(Vertex v, int dist) : v(v), dist(dist) {}
};

struct MinHeap
{
    MinHeapNode** array;
    int* pos;
    int size;
    int capacity;

    MinHeap(int capacity) : size(0), capacity(capacity)
    {
        array = new MinHeapNode*[capacity];
        pos = new int[capacity];
        for (int i = 0; i < capacity; i++) { pos[i] = -1; }
    }

    ~MinHeap()
    {
        delete[] array;
        delete[] pos;
    }
    // O(1)
    void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b)
    {
        MinHeapNode* t = *a;
        *a = *b;
        *b = t;
    }
    // O(log(n))
    void minHeapify(int idx)
    {
        int smallest = idx;
        int left = 2 * idx + 1;
        int right = 2 * idx + 2;

        if (left < size && array[left]->dist < array[smallest]->dist)
            smallest = left;

        if (right < size && array[right]->dist < array[smallest]->dist)
            smallest = right;

        if (smallest != idx)
        {
            MinHeapNode* smallestNode = array[smallest];
            MinHeapNode* idxNode = array[idx];

            pos[smallestNode->v] = idx;
            pos[idxNode->v] = smallest;

            swapMinHeapNode(&array[smallest], &array[idx]);

            minHeapify(smallest);
        }
    }
    // O(1)
    bool isEmpty() const
    {
        return size == 0;
    }
    // O(log(n))
    MinHeapNode* extractMin()
    {
        if (isEmpty()) { return nullptr; }

        MinHeapNode* root = array[0];

        MinHeapNode* lastNode = array[size - 1];
        array[0] = lastNode;

        pos[root->v] = size - 1;
        pos[lastNode->v] = 0;

        --size;
        minHeapify(0);

        return root;
    }
    // O(log(n))
    void decreaseKey(Vertex v, int dist)
    {
        int i = pos[v];
        array[i]->dist = dist;

        while (i && array[i]->dist < array[(i - 1) / 2]->dist)
        {
            pos[array[i]->v] = (i - 1) / 2;
            pos[array[(i - 1) / 2]->v] = i;
            swapMinHeapNode(&array[i], &array[(i - 1) / 2]);

            i = (i - 1) / 2;
        }
    }
    // O(log(n))
    void insert(Vertex v, int dist)
    {
        if (size == capacity) { return; }

        MinHeapNode* node = new MinHeapNode(v, dist);

        array[size] = node;
        pos[v] = ++size;

        decreaseKey(v, dist);
    }
    // O(log(n))
    void insertOrUpdate(Vertex v, int dist)
    {
        if (pos[v] == -1) { insert(v, dist); }
        else { decreaseKey(v, dist); }
    }
};

#endif