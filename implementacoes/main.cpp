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
    int* proximos = (int*) malloc(k*sizeof(int));
    if (!proximos) return NULL;

    int* dists = (int*) malloc(n*sizeof(int));
    if (!dists) { free(proximos); return NULL; }
    
    unordered_map<int, int> hash;

    for (int i=0; i < n; i++)
    {
        int dist = abs(arr[i] - x);
        dists[i] = dist;
        hash[dist] = i;
    }

    for (int i=0; i < k; i++)
    {
        int distLoc = quickSelect(dists, i, n-1, i);
        proximos[i] = arr[hash[distLoc]];
    }

    return proximos;
}


int main()
{
    int arr[] = {1, 3, 7, 10, 15, 20, 25};
    int n = sizeof(arr) / sizeof(arr[0]);
    int x = 12;  // Valor de referência
    int k = 3;   // Número de elementos mais próximos a buscar

    int* result = quest_12(arr, n, x, k);

    if (result)
    {
        printf("Os %d elementos mais próximos de %d são:\n", k, x);
        for (int i = 0; i < k; i++)
        {
            printf("%d ", result[i]);
        }
        printf("\n");

        // Libera a memória alocada
        free(result);
    }
    else
    {
        printf("Falha na alocação de memória.\n");
    }


    return 0;
}