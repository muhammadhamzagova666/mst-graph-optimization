// Language: C++
//
// File: Efficient_Kruskal.cpp
// Description: Efficient implementation of Kruskal's algorithm to compute the Minimum Spanning Tree (MST)
//              for a graph. This module generates a random undirected weighted graph and measures the runtime
//              performance of the MST computation. Designed for researchers and developers comparing algorithmic
//              efficiencies in network design and optimization.
//
// Author: [Gova]

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>

using namespace std;

// Structure representing an edge in the graph, connecting two vertices with a specified weight.
struct Edge {
    int src, dest, weight;
};

// Structure for union-find subsets used in Kruskal's algorithm, supporting path compression and union by rank.
struct Subset {
    int parent, rank;
};

class Graph {
private:
    vector<Edge> edges;  // List of all potential edges in the graph.
    int numVertices;     // Total number of vertices in the graph.

public:
    /**
     * Constructor for Graph.
     * @param V: Number of vertices in the graph.
     */
    Graph(int V) : numVertices(V) {}

    /**
     * Adds an undirected edge between two vertices with a given weight.
     * @param src: Source vertex index.
     * @param dest: Destination vertex index.
     * @param weight: Weight of the edge.
     */
    void addEdge(int src, int dest, int weight) {
        // Store the edge for later use in MST computation.
        edges.push_back({src, dest, weight});
    }

    /**
     * Find the subset of an element i using path compression.
     * This helps flatten the structure of the tree, improving future queries.
     *
     * @param subsets: Array of subsets representing disjoint sets.
     * @param i: The element whose root subset is to be found.
     * @return The root of the subset containing i.
     */
    int find(Subset subsets[], int i) {
        if (subsets[i].parent != i)
            subsets[i].parent = find(subsets, subsets[i].parent);
        return subsets[i].parent;
    }

    /**
     * Merges two subsets x and y using union by rank.
     * This minimizes the height of the resulting tree.
     *
     * @param subsets: Array of subsets.
     * @param x: First element index.
     * @param y: Second element index.
     */
    void Union(Subset subsets[], int x, int y) {
        int xroot = find(subsets, x);
        int yroot = find(subsets, y);

        // Attach smaller rank tree under root of the higher rank tree.
        if (subsets[xroot].rank < subsets[yroot].rank)
            subsets[xroot].parent = yroot;
        else if (subsets[xroot].rank > subsets[yroot].rank)
            subsets[yroot].parent = xroot;
        else {
            // If ranks are same, choose one as root and increment its rank.
            subsets[yroot].parent = xroot;
            subsets[xroot].rank++;
        }
    }

    /**
     * Implements Kruskal's algorithm to compute the MST.
     * The algorithm sorts edges by weight and adds them one by one while avoiding cycles.
     *
     * @return A vector of edges representing the MST.
     */
    vector<Edge> kruskalMST() {
        vector<Edge> result;    // Holds the resulting MST edges.

        // Sort edges in non-decreasing order to use the greedy approach.
        sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
            return a.weight < b.weight;
        });

        // Initialize subsets: each vertex is initially in its own subset.
        Subset* subsets = new Subset[numVertices];
        for (int i = 0; i < numVertices; i++) {
            subsets[i].parent = i;
            subsets[i].rank = 0;
        }

        int i = 0;  // Index to iterate through sorted edges.

        // Continue until the MST has (numVertices - 1) edges.
        while (result.size() < numVertices - 1) {
            // Select the next minimum weight edge.
            Edge next_edge = edges[i++];

            // Determine the subsets of the vertices connected by this edge.
            int x = find(subsets, next_edge.src);
            int y = find(subsets, next_edge.dest);

            // If including this edge doesn't form a cycle, add it to the MST.
            if (x != y) {
                result.push_back(next_edge);
                Union(subsets, x, y);
            }
        }

        // Clean up dynamic memory.
        delete[] subsets;
        return result;
    }

    /**
     * Generates a complete graph with random weights.
     * This simulates a random undirected graph with all possible edges between points.
     *
     * @param numPoints: Number of vertices (points) to generate.
     * @param maxWeight: The maximum weight for any edge.
     */
    void generateRandomGraph(int numPoints, int maxWeight) {
        // Seed the random number generator to ensure varied results.
        srand(static_cast<unsigned int>(time(nullptr)));

        // For each unique pair, add an edge with a random weight.
        for (int i = 0; i < numPoints; ++i) {
            for (int j = i + 1; j < numPoints; ++j) {
                int weight = rand() % maxWeight + 1; // Ensure weight is between 1 and maxWeight.
                addEdge(i, j, weight);
            }
        }
    }

    /**
     * Computes the total cost (sum of edge weights) of the MST.
     *
     * @param MST: Vector containing the edges that form the MST.
     * @return Total weight/cost of the MST.
     */
    int calculateMSTCost(const vector<Edge>& MST) {
        int cost = 0;
        for (const Edge& edge : MST) {
            cost += edge.weight;  // Accumulate each edge's weight.
        }
        return cost;
    }
};

int main() {
    // Define the number of vertices and the maximum possible edge weight
    const int numPoints = 500;  // Adjustable parameter for graph size.
    const int maxWeight = 1000; // Maximum weight for any edge in the graph.

    // Instantiate the graph and generate its random weighted edges.
    Graph g(numPoints);
    g.generateRandomGraph(numPoints, maxWeight);

    // Record the starting time before MST computation.
    struct timeval start;
    gettimeofday(&start, NULL);

    // Compute the MST using the efficient Kruskal algorithm.
    vector<Edge> MST = g.kruskalMST();

    // Record the ending time and calculate the duration in microseconds.
    struct timeval stop;
    gettimeofday(&stop, NULL);
    long long duration = (stop.tv_sec - start.tv_sec) * 1000000LL + stop.tv_usec - start.tv_usec;

    // Informative output detailing the MST edges, total cost, and execution time.
    cout << "Minimum Spanning Tree Edges:" << "\n";
    for (const Edge& edge : MST) {
        // Using formatted output to help with debugging and clarity.
        cout << "Vertex " << edge.src << " -- Vertex " << edge.dest 
             << " | Edge Weight: " << edge.weight << "\n";
    }

    int cost = g.calculateMSTCost(MST);
    cout << "Total MST Cost: " << cost << "\n";
    cout << "Execution Time: " << duration << " microseconds\n";

    return 0;
}
// End of Efficient Kruskal algorithm implementation.