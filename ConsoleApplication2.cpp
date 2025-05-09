#include <iostream>
#include <vector>
#include <omp.h>
#include <climits>
#include <fstream>
using namespace std;
const int INF = 1e9;
vector<vector<int>> generateConnectedGraph(int N, int minWeight, int maxWeight)
{
    const int INF = 1e9;
    vector<vector<int>> adjMatrix(N, vector<int>(N, INF));
    srand(time(0));
    for (int i = 0; i < N; ++i) adjMatrix[i][i] = 0;
    for (int i = 1; i < N; ++i)
    {
        int weight = minWeight + rand() % (maxWeight - minWeight + 1);
        adjMatrix[i][i - 1] = weight;
        adjMatrix[i - 1][i] = weight;
    }
    int extraEdges = rand() % (N * (N - 1) / 2 - (N - 1));
    for (int i = 1; i <= extraEdges; ++i)
    {
        int u = rand() % N;
        int v = rand() % N;
        if ((u != v) && (abs(u - v) != 1) && (adjMatrix[u][v] == INF))
        {
            int weight = minWeight + rand() % (maxWeight - minWeight + 1);
            adjMatrix[u][v] = adjMatrix[v][u] = weight;
        }
    }
    return adjMatrix;
}
void adjMatrixToFile(const vector<vector<int>>& matrix, const string& filename)
{
    int N = matrix.size();
    ofstream fout(filename);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j) fout << matrix[i][j] << " ";
        fout << "\n";
    }
    fout.close();
}
void answersToFile(const vector<vector<int>>& matrix, const string& filename)
{
    int N = matrix.size();
    ofstream fout(filename);
    for (int i = 0; i < N; ++i)
    {
        for (int j = i + 1; j < N; ++j) fout << "The shortest distance between " << i << " and " << j << " is: " << matrix[i][j] << "\n";
    }
    fout.close();
}
void floydWarshall(std::vector<std::vector<int>>& dist, int N, int num_threads) 
{
    omp_set_num_threads(num_threads);
    for (int k = 0; k < N; ++k) 
    {
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) 
        {
            for (int j = 0; j < N; ++j) 
            {
                if (dist[i][k] < 1e9 && dist[k][j] < 1e9 && dist[i][k] + dist[k][j] < dist[i][j]) dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }
}
int main() 
{
    int N;
    int num_threads;
    cout << "Enter the number of vertexes: ";
    cin >> N;
    cout << "Enter the number of threads: ";
    cin >> num_threads;
    vector<vector<int>> Matrix;
    Matrix = generateConnectedGraph(N, 1, 10);
    adjMatrixToFile(Matrix, "adjacency_matrix.txt");
    double start_time = omp_get_wtime();
    floydWarshall(Matrix, N, num_threads);
    double end_time = omp_get_wtime();
    answersToFile(Matrix, "shortest_paths.txt");
    std::cout << "Execution time: " << (end_time - start_time) << " seconds" << "\n";
    return 0;
}





