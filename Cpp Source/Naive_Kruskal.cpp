// Language: C++
//
// File: Naive_Kruskal.cpp
// Description: Naive implementation of Kruskal's algorithm for computing the Minimum Spanning Tree (MST)
//              of a randomly generated undirected weighted graph. This implementation uses insertion sort
//              for edge sorting and a basic disjoint-set (union-find) structure without path optimization.
//              Designed for comparing algorithm performance and educational purposes.
// Author: [Gova]

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>  // For high resolution time measurement

using namespace std;

// Structure representing an edge in the graph which connects two vertices with a specific weight.
struct Edge {
    int src, dest, weight;
};

// Structure representing a subset for union-find, used to detect cycles.
// 'parent' holds the representative of the set, and 'rank' helps optimize union operations.
struct Subset {
    int parent, rank;
};

class Graph {
private:
    vector<Edge> edges;   // Collection of all edges in the graph.
    int numVertices;      // Total number of vertices in the graph.

public:
    /**
     * Constructor for Graph.
     *
     * @param V Number of vertices in the graph.
     */
    Graph(int V) : numVertices(V) {}

    /**
     * Adds an undirected edge between two vertices.
     *
     * @param src Source vertex index.
     * @param dest Destination vertex index.
     * @param weight Weight of the edge.
     */
    void addEdge(int src, int dest, int weight) {
        // Save the new edge to be considered in the MST construction.
        edges.push_back({src, dest, weight});
    }

    /**
     * Finds the set for an element using path compression.
     * This efficiently detects cycles and ensures the union-find tree remains flat.
     *
     * @param subsets Array of subsets representing the disjoint sets.
     * @param i The element to find the set for.
     * @return The representative of the set which contains element i.
     */
    int find(Subset subsets[], int i) {
        if (subsets[i].parent != i)
            subsets[i].parent = find(subsets, subsets[i].parent);
        return subsets[i].parent;
    }

    /**
     * Unites two sets using union by rank.
     * This minimizes the depth of the resulting trees during union operations.
     *
     * @param subsets Array of subsets.
     * @param x First set representative.
     * @param y Second set representative.
     */
    void Union(Subset subsets[], int x, int y) {
        int xroot = find(subsets, x);
        int yroot = find(subsets, y);

        if (subsets[xroot].rank < subsets[yroot].rank)
            subsets[xroot].parent = yroot;
        else if (subsets[xroot].rank > subsets[yroot].rank)
            subsets[yroot].parent = xroot;
        else {
            subsets[yroot].parent = xroot;
            subsets[xroot].rank++;
        }
    }

    /**
     * Sorts the edges in non-decreasing order based on their weights using insertion sort.
     * Although suboptimal for large datasets, insertion sort is used here to illustrate the naive approach.
     */
    void insertionSort() {
        int n = edges.size();
        for (int i = 1; i < n; i++) {
            Edge key = edges[i];
            int j = i - 1;
            // Shift elements that are larger than key one position to the right.
            while (j >= 0 && edges[j].weight > key.weight) {
                edges[j + 1] = edges[j];
                j = j - 1;
            }
            edges[j + 1] = key;
        }
    }

    /**
     * Implements Kruskal's algorithm to compute the MST.
     * This function sorts the edges and then iteratively adds the smallest edge that does not create a cycle.
     *
     * @return A vector containing the edges included in the MST.
     */
    vector<Edge> kruskalMST() {
        vector<Edge> result;  // This will store the resultant MST.

        // Sort edges based on non-decreasing edge weight.
        insertionSort();

        // Allocate memory for V subsets, initially making each vertex its own subset.
        Subset* subsets = new Subset[numVertices];
        for (int i = 0; i < numVertices; i++) {
            subsets[i].parent = i;
            subsets[i].rank = 0;
        }

        int i = 0; // Index to iterate over sorted edges.
        // Continue until MST contains (numVertices - 1) edges.
        while (result.size() < numVertices - 1) {
            // Consider the next edge with the smallest weight.
            Edge next_edge = edges[i++];

            int x = find(subsets, next_edge.src);
            int y = find(subsets, next_edge.dest);

            // Check if including this edge would create a cycle.
            if (x != y) {
                result.push_back(next_edge);
                Union(subsets, x, y);  // Combine the sets.
            }
        }

        // Free allocated memory for subsets.
        delete[] subsets;

        return result;
    }

    /**
     * Generates a complete graph with random edge weights.
     * This simulates a random undirected graph for testing the MST algorithm.
     *
     * @param numPoints Number of vertices.
     * @param maxWeight Maximum possible weight for any edge.
     */
    void generateRandomGraph(int numPoints, int maxWeight) {
        // Seed the random number generator to ensure varied graph structures.
        srand(static_cast<unsigned int>(time(nullptr)));
        for (int i = 0; i < numPoints; ++i) {
            for (int j = i + 1; j < numPoints; ++j) {
                int weight = rand() % maxWeight + 1; // Generate a weight between 1 and maxWeight.
                addEdge(i, j, weight);
            }
        }
    }

    /**
     * Computes the total cost of the MST.
     *
     * @param MST Vector of edges that constitute the MST.
     * @return Total weight (cost) of the MST.
     */
    int calculateMSTCost(const vector<Edge>& MST) {
        int cost = 0;
        for (const Edge& edge : MST) {
            cost += edge.weight;  // Summing edge weights for the overall MST cost.
        }
        return cost;
    }
};

int main() {
    // Define the graph size and weight range. Adjust these for scaling tests.
    const int numPoints = 500;   // Number of vertices in the graph.
    const int maxWeight = 1000;  // Maximum weight for an edge.

    // Create the graph and generate a random graph instance.
    Graph g(numPoints);
    g.generateRandomGraph(numPoints, maxWeight);

    // Measure the execution time of the MST computation.
    struct timeval start, stop;
    gettimeofday(&start, NULL);
    vector<Edge> MST = g.kruskalMST();
    gettimeofday(&stop, NULL);

    // Calculate elapsed time in microseconds.
    long long elapsedTime = (stop.tv_sec - start.tv_sec) * 1000000LL +
                            (stop.tv_usec - start.tv_usec);

    // Output the MST edges with clear formatting for easy debugging.
    cout << "Edges in MST:" << "\n";
    for (const Edge& edge : MST) {
        cout << edge.src << " - " << edge.dest << " : " << edge.weight << "\n";
    }

    // Calculate and display the MST cost.
    int cost = g.calculateMSTCost(MST);
    cout << "Cost of MST: " << cost << "\n";
    cout << "Running time: " << elapsedTime << " microseconds\n";

    return 0;
}
// End of Naive_Kruskal.cpp implementation.