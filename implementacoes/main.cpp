#include <iostream>
#include <unordered_map>
#include <vector>

using namespace std;

#define swap(v, i, j) { int temp = v[i]; v[i] = v[j]; v[j] = temp; }

// class Node
// {
//     public:
//         Node(int key, char data):
//             m_key(key),
//             m_data(data),
//             m_leftNode(nullptr),
//             m_rightNode(nullptr),
//             m_parentNode(nullptr) {}
        
//         Node & leftNode() const { return * m_leftNode; }
//         void setLeftNode(Node * node) { m_leftNode = node;}

//         Node & rightNode() const { return * m_rightNode; }
//         void setRightNode(Node * node) { m_rightNode = node; }

//         Node & parentNode() const { return * m_parentNode; }
//         void setParentNode(Node * node) { m_parentNode =
// node; }
// private:
// int m_key;
// char m_data;
// Node * m_leftNode;
// Node * m_rightNode;
// Node * m_parentNode;
// };

int partition(int arr[], int inicio, int fim)
{
    int pivo = arr[fim];
    int i = inicio - 1;

    for (int j = inicio; j < fim; j++)
    {
        if (arr[j] < pivo)
        {
            i++;
            swap(arr, i, j);
        }
    }

    swap(arr, i + 1, fim);
    return i + 1;
}

int quickSelect(int arr[], int inicio, int fim, int k)
{
    if (inicio <= fim)
    {
        int pivo = partition(arr, inicio, fim);

        if (pivo == k)
        {
            return arr[pivo];
        }

        if (pivo > k)
        {
            return quickSelect(arr, inicio, pivo - 1, k);
        }

        return quickSelect(arr, pivo + 1, fim, k);
    }

    return arr[inicio];
}

int quest_7(int arr[], int n, int k)
{
    unordered_map<int, int> hash;
    
    for (int i = 0; i < n; i++)
    {
        hash[arr[i]]++;
    }

    vector<int> freq;
    
    for (const auto& pair : hash)
    {
        freq.push_back(pair.second);
    }

    int k_freq = quickSelect(freq.data(), 0, freq.size() - 1, freq.size() - k);

    for (const auto& pair : hash)
    {
        if (pair.second == k_freq)
        {
            return pair.first;
        }
    }

    return -1;
}



int main() {
    int arr[] = {10, 15, 3, 7, 8};
    int n = sizeof(arr) / sizeof(arr[0]);
    int x = 13;

    pair<int, int> result = quest_8(arr, n, x);

    if (result.first != -1 && result.second != -1) {
        cout << "Par encontrado nos índices: " << result.first << " e " << result.second << endl;
    } else {
        cout << "Nenhum par encontrado." << endl;
    }

    return 0;
}