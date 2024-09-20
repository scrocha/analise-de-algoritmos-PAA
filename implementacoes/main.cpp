#include <iostream>
#include <unordered_map>
#include <vector>

using namespace std;

#define swap(v, i, j) { int temp = v[i]; v[i] = v[j]; v[j] = temp; }

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
    if (inicio < fim -1)
    {
        int meio = (inicio + fim) / 2;
        mergeSort(arr, inicio, meio);
        mergeSort(arr, meio, fim);
        merge(arr, inicio, meio, fim);
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
        if (total >= x)
        {
            return x;
        }
    }

    return -1;
}

int quest_11(int arvore[], int n)
{
    return 0;
}

int* quest_12(int arr[], int n, int x, int k)
{
   return nullptr; 
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


int main() {
    int arr[] = {1, 7, 4, 9, 10, 15, 32};
    int n = sizeof(arr) / sizeof(arr[0]);
    int x = 2;

    vector<int> result = quest_13(arr, n, x);

    cout << "O maior subconjunto onde a diferença é <= " << x << " é: ";
    for (int i = 0; i < result.size(); i++) {
        cout << result[i] << " ";
    }
    cout << endl;

    return 0;
}