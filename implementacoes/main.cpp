#include <iostream>
#include <unordered_map>
#include <queue>
#include <vector>
#include <climits>

using namespace std;

#define swap(v, i, j) { int temp = v[i]; v[i] = v[j]; v[j] = temp; }

#define swap2(v, i, j) { pair<int, int> temp = v[i]; v[i] = v[j]; v[j] = temp; }

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

void heapify(int arr[], int n, int i)
{
    int index = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    
    if (right < n && arr[right] > arr[index])
    {
        index = right;
    }

    else if (left < n && arr[left] > arr[index])
    {
        index = left;
    }

    if (index != i)
    {
        swap(arr, i, index);
        heapify(arr, n, index);
    }

}

void buildHeap(int arr[], int n)
{
    for (int i = n / 2 - 1; i >= 0; i--)
    {
        heapify(arr, n, i);
    }
}

void heapSort(int arr[], int n)
{
    buildHeap(arr, n);

    for (int i = n - 1; i > 0; i--)
    {
        swap(arr, 0, i);
        heapify(arr, i, 0);
    }
}

void merge(int arr[], int inicio, int meio, int fim)
{
    int* r = (int*) malloc((fim - inicio) * sizeof(int));
    if (!r) return;

    int a = inicio, b = meio, i = 0;
    while (a < meio && b < fim)
    {
        r[i++] = arr[a] < arr[b] ? arr[a++] : arr[b++];
    }
    while (a < meio) r[i++] = arr[a++];
    while (b < fim) r[i++] = arr[b++];

    for (a = inicio; a < fim; a++)
    {
        arr[a] = r[a - inicio];
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

int partition2(int arr[], int inicio, int fim, int pivo)
{
    int i = inicio - 1;

    for (int j = inicio; j < fim; j++)
    {
        if (arr[j] < pivo)
        {
            i++;
            swap(arr, i, j);
        }
    }
    return i + 1;
}

int medianOf(int arr[], int n)
{
    mergeSort(arr, 0, n);
    return arr[n / 2];
}

int selectMOM(int arr[], int inicio, int fim, int k)
{
    int n = fim - inicio + 1;
    if (k <= 0 || k > n) return -1;

    int grupos = (n + 4) / 5;
    int* median = new int[grupos];
    
    int i = 0;
    int pos = inicio;

    while (pos < fim)
    {
        int size = min(5, fim - pos + 1);
        median[i++] = medianOf(&arr[pos], size);
        pos += 5;
    }

    int mom = (i == 1) ? median[i - 1] : selectMOM(median, 0, i - 1, i / 2);

    int j = partition2(arr, inicio, fim, mom);

    if (j - inicio == k - 1) return arr[j];
    if (j - inicio > k - 1) return selectMOM(arr, inicio, j - 1, k);
    return selectMOM(arr, j + 1, fim, k - j + inicio - 1);
}

int quickSelectMOM(int arr[], int inicio, int fim, int k)
{
    if (inicio <= fim)
    {
        int pivo = selectMOM(arr, inicio, fim, k);
        return pivo;
    }

    return -1;
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

    int k_freq = quickSelectMOM(freq.data(), 0, freq.size() - 1, freq.size() + 1 - k);

    for (const auto& pair : hash)
    {
        if (pair.second == k_freq) return pair.first;
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
    unordered_map<int, int> freq;

    for (int i = 0; i < n; i++)
    {
        arr[i] < n ? freq[arr[i]]++: freq[n]++;
    }

    int count = 0;
    for (int i = n; i >= 0; i--)
    {
        if (freq[i] == 0} continue;
        count += freq[i];
        if (count >= i) return i;
    }

    return -1;
}

void BSTtoVector(Node* root, vector<int>& result)
{
    if (root == nullptr) return;
    
    BSTtoVector(root->left, result);

    result.push_back(root->data);

    BSTtoVector(root->right, result);
}

int quest_11(Node* root)
{
    vector<int> result;
    BSTtoVector(root, result);

    int n = result.size();
    int best = INT_MAX;

    for (int i = 1; i < n; i++)
    {
        int diff = result[i] - result[i - 1];
        if (diff < best)
        {
            best = diff;
        }
    }

    return best;
}

void heapify(pair<int, int>* arr, int n, int i)
{
    int index = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    
    if (left < n && arr[left].first > arr[index].first)
    {
        index = left;
    }
    
    if (right < n && arr[right].first > arr[index].first)
    {
        index = right;
    }

    if (index != i)
    {
        swap2(arr, i, index);
        heapify(arr, n, index);
    }

}

void buildHeap(pair<int, int>* arr, int n)
{
    for (int i = n / 2 - 1; i >= 0; i--)
    {
        heapify(arr, n, i);
    }
}

int* quest_12(int arr[], int n, int x, int k)
{
    if (k > n || k <= 0) return nullptr;

    int* result = (int*) malloc(k * sizeof(int));
    if (!result) return nullptr;

    pair<int, int>* maxHeap = (pair<int, int>*) malloc(k * sizeof(pair<int, int>));
    if (!maxHeap) { free(result); return nullptr; }
    
    for (int i = 0; i < k; ++i)
    {
        int distance = abs(arr[i] - x);
        maxHeap[i] = {distance, arr[i]};
    }

    buildHeap(maxHeap, k);

    for (int i = k; i < n; ++i)
    {
        int distance = abs(arr[i] - x);

        if (distance < maxHeap[0].first)
        {
            maxHeap[0] = {distance, arr[i]};
            heapify(maxHeap, k, 0);
        }
    }

    for (int i = 0; i < k; ++i)
    {
        result[k-1 - i] = maxHeap[i].second;
    }

    free(maxHeap);
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
    int i = 0, j = n - 1, maixmo = 0;
    pair<int, int> result = {0, 0};

    while (i < j)
    {
        int produto = (j - i) * min(arr[i], arr[j]);

        if (produto > maixmo)
        {
            maixmo = produto;
            result = {i, j};
        }
        
        if (arr[i] < arr[j]) { i++; }
        else { j--; }
    }

    return result;
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
    return 0;
}
