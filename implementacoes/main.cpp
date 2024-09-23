#include <iostream>
#include <unordered_map>
#include <vector>
#include <climits>

using namespace std;

#define swap(v, i, j) { int temp = v[i]; v[i] = v[j]; v[j] = temp; }

typedef struct Node
{
    int data;
    Node* left;
    Node* right;

    Node(int data) : data(data), left(nullptr), right(nullptr) {}
} Node;

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

int quickSelectMOM(int arr[], int inicio, int fim, int k)
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
            return quickSelectMOM(arr, inicio, pivo - 1, k);
        }

        return quickSelectMOM(arr, pivo + 1, fim, k);
    }

    return arr[inicio];
}

void merge(int arr[], int inicio1, int inicio2, int fim2)
{
    int* r = (int*) malloc((fim2 - inicio1) * sizeof(int));
    if (!r) return;

    int a = inicio1, b = inicio2, i = 0;
    while (a < inicio2 && b < fim2)
    {
        r[i++] = arr[a] < arr[b] ? arr[a++] : arr[b++];
    }
    while (a < inicio2) r[i++] = arr[a++];
    while (b < fim2) r[i++] = arr[b++];

    for (a = inicio1; a < fim2; a++)
    {
        arr[a] = r[a - inicio1];
    }

    free(r);
}

void mergeSort(int arr[], int inicio, int fim)
{
    if (inicio < fim - 1)
    {
        int meio = (inicio + fim) / 2;
        mergeSort(arr, inicio, meio);
        mergeSort(arr, meio, fim);
        merge(arr, inicio, meio, fim);
    }
}

void inOrder(Node* root)
{
    if (root)
    {
        inOrder(root->left);
        cout << root->data << " ";
        inOrder(root->right);
    }
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

    int k_freq = quickSelectMOM(freq.data(), 0, freq.size() - 1, freq.size() - k);

    for (const auto& pair : hash)
    {
        if (pair.second == k_freq)
        {
            return pair.first;
        }
    }

    return -1;
}

pair<int, int> quest_8(int arr[], int n, int x)
{
    unordered_map<int, int> hash;

    for (int i = 0; i < n; i++)
    {
        int complemento = x - arr[i];

        if (hash.find(complemento) != hash.end())
        {
            return {hash[complemento], i};
        }

        hash[arr[i]] = i;
    }

    return {-1, -1};
}

int quest_9(int** const arr, int n, int m)
{
    unordered_map<int, int> hash;
    int count = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            hash[arr[i][j]]++;
        }

    }

    for (int i = 0; i < n; i++)
    {
        bool incomum = true;
        for (int j = 0; j < m; j++)
        {
            if (hash[arr[i][j]] > 1)
            {
                incomum = false;
                break;
            }
        }

        if (incomum)
        {
            count++;
        }
    }

    return count;
}

int quest_10(int arr[], int n)
{
    int* counting = (int*) malloc((n + 1) * sizeof(int));
    if (!counting) return -1;

    for (int i = 0; i < n; i++)
    {
        arr[i] < n ? counting[arr[i]]++ : counting[n]++;
    }

    int total = 0;

    for (int x = n; x >= 0; x--)
    {
        total += counting[x];
        cout << total << ' ' << counting[x] << endl;
        if (total == x)
        {
            free(counting);
            return x;
        }
    }

    return -1;
}

int BSTsearch(Node* current, int prev, int minDiff)
{
    if (current)
    {
        minDiff = BSTsearch(current->left, prev, minDiff);

        if (prev != -1)
        {
            minDiff = min(minDiff, abs(current->data - prev));
            cout << current->data << ' ' << prev << ' ' << minDiff << endl;
        }
        prev = current->data;

        minDiff = BSTsearch(current->right, prev, minDiff);
    }

    return minDiff;
}

int quest_11(Node* root)
{
    int prev = -1;
    int minDiff = INT_MAX;
    if (!root) return -1;
    else if (!root->left && !root->right) return -1;
    else if (!root->left && root->right) minDiff = root->right->data - root->data;
    else if (root->left && !root->right) minDiff = root->data - root->left->data;

    minDiff = BSTsearch(root, prev, minDiff);
    return minDiff;
}


int* quest_12(int arr[], int n, int x, int k)
{
    int* result = (int*) malloc(k * sizeof(int));
    if (!result) return nullptr;

    int* diff = (int*) malloc(n * sizeof(int));
    if (!diff)
    {
        free(result);
        return nullptr;
    }

    unordered_map<int, int> hash;

    for (int i = 0; i < n; i++)
    {
        diff[i] = abs(arr[i] - x);
        hash[arr[i]] = diff[i];
    }

    int value = quickSelectMOM(diff, 0, n - 1, k - 1);
    int count = 0;

    for (int i = 0; i < n && count < k; i++)
    {
        if (hash[arr[i]] <= value)
        {
            result[count++] = arr[i];
        }
    }

    free(diff);

    return result;
}

vector<int> quest_13(int arr[], int n, int x)
{
    vector<int> result;

    mergeSort(arr, 0, n);

    int maxSize = 1;
    int size = 1;

    int bestStart = 0;
    int start = 0;

    for (int end = 1; end < n; end++)
    {
        if (arr[end] - arr[end - 1] <= x)
        {
            size++;

            if (size > maxSize)
            {
                maxSize = size;
                bestStart = start;
            }
        }
        else
        {
            start = end;
            size = 1;
        }
    }

    for (int i = bestStart; i < bestStart + maxSize; i++)
    {
        result.push_back(arr[i]);
    }

    return result;
}

pair<int, int> quest_14(int arr[], int n)
{
    return {0, 0};
}

Node* arrayToBST(int arr[], int inicio, int fim)
{
    if (inicio > fim) return nullptr;

    int meio = (inicio + fim) / 2;
    Node* root = new Node(arr[meio]);

    root->left = arrayToBST(arr, inicio, meio - 1);
    root->right = arrayToBST(arr, meio + 1, fim);

    return root;
}

Node* quest_15(int arr[], int n)
{
    mergeSort(arr, 0, n);
    
    return arrayToBST(arr, 0, n - 1);
}

int main()
{
    int arr[] = {1, 5, 10, 15, 20, 25, 27};
    int n = sizeof(arr) / sizeof(arr[0]);
    int x = 25;
    int k = 3;

    int* arr2 = quest_12(arr, n, x, k);
    for (int i = 0; i < k; i++)
    {
        cout << arr2[i] << ' ';
    }
    free(arr2);
    return 0;
}