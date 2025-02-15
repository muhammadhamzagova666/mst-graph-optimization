// Language: C++
//
// File: Naive_Prim.cpp
// Description: Naive implementation of Prim's algorithm for computing the Minimum Spanning Tree (MST)
//              from a randomly generated undirected weighted graph. The algorithm uses a simple linear
//              search to select the next vertex, making it suitable for educational purposes and small graphs.
// Author: [Gova]

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h> // For measuring execution time
#include <math.h>
#include <bits/stdc++.h>

using namespace std;

/// Structure representing a 2D point.
/// Points are used to generate vertices with coordinates, where the edge weight is based on Euclidean distance.
struct Point {
    int x, y;
};

/// Node in the adjacency list representing an edge.
/// Contains the destination vertex, edge weight, and a pointer to the next node.
struct AdjListNode {
    int dest;
    int weight;
    struct AdjListNode* next;
};

/// Adjacency list for a vertex.
/// Holds the head pointer for a linked list of AdjListNode.
struct AdjList {
    struct AdjListNode* head;
};

/// Graph structure represented as an array of adjacency lists.
/// V: Total number of vertices.
struct Graph {
    int V;                 // Number of vertices
    struct AdjList* array; // Array of adjacency lists (one per vertex)
};

/// Creates a new adjacency list node for a given destination and weight.
///
/// @param dest Destination vertex index.
/// @param weight Weight assigned to the edge.
/// @return A pointer to the newly created AdjListNode.
struct AdjListNode* newAdjListNode(int dest, int weight) {
    struct AdjListNode* newNode = (struct AdjListNode*)malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->weight = weight;
    newNode->next = NULL;
    return newNode;
}

/// Allocates and initializes a new graph with V vertices.
///
/// @param V Number of vertices in the graph.
/// @return A pointer to the newly created Graph.
struct Graph* createGraph(int V) {
    struct Graph* graph = (struct Graph*)malloc(sizeof(struct Graph));
    graph->V = V;
    graph->array = (struct AdjList*)malloc(V * sizeof(struct AdjList));

    // Initialize each adjacency list as empty.
    for (int i = 0; i < V; ++i)
        graph->array[i].head = NULL;

    return graph;
}

/// Adds an undirected edge between two vertices with the given weight.
/// The function updates the adjacency lists of both vertices.
///
/// @param graph Pointer to the Graph.
/// @param src Source vertex index.
/// @param dest Destination vertex index.
/// @param weight Weight of the edge.
void addEdge(struct Graph* graph, int src, int dest, int weight) {
    // Insert node for destination in the source vertex's list.
    struct AdjListNode* newNode = newAdjListNode(dest, weight);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;

    // Since the graph is undirected, add an edge from destination to source.
    newNode = newAdjListNode(src, weight);
    newNode->next = graph->array[dest].head;
    graph->array[dest].head = newNode;
}

/// Calculates the Euclidean distance between two points.
/// The distance is used as the weight for edges between vertices.
///
/// @param p1 First point.
/// @param p2 Second point.
/// @return The computed weight between p1 and p2.
int calculateWeight(struct Point p1, struct Point p2) {
    return (int)sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

/// Prints the entire graph's adjacency list representation.
/// For each vertex, outputs all connecting edges with weights.
///
/// @param graph Pointer to the Graph.
void printGraph(struct Graph* graph) {
    int V = graph->V;
    printf("Randomly Generated Graph:\n");
    for (int i = 0; i < V; ++i) {
        struct AdjListNode* pCrawl = graph->array[i].head;
        printf("Vertex %d: ", i);
        // The following print also includes a recalculated weight as a sanity check.
        while (pCrawl != NULL) {
            printf("(%d, %d, %d) ", pCrawl->dest, pCrawl->weight,
                   calculateWeight({i, 0}, {pCrawl->dest, 0}));
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}

/// Prints the MST (represented via a parent array) in a formatted manner.
///
/// @param arr Parent array where each index represents a vertex and its value the parent vertex.
/// @param n Number of vertices.
void printArr(int arr[], int n) {
    printf("Edges in Minimum Spanning Tree:\n");
    for (int i = 1; i < n; ++i)
        printf("%d - %d\n", arr[i], i);
}

/// Generates random points for graph vertices.
/// Each point is assigned random x and y coordinates in the range [0, 100].
///
/// @param points Array to be populated with generated points.
/// @param V Number of points (vertices) to generate.
void generateRandomPoints(struct Point points[], int V) {
    srand(time(NULL)); // Seed random number generator for varied output
    for (int i = 0; i < V; ++i) {
        points[i].x = rand() % 100;
        points[i].y = rand() % 100;
    }
}

/// Creates a complete graph from an array of points.
/// For every pair of points, computes the edge weight and adds an undirected edge.
///
/// @param graph Pointer to the Graph to be built.
/// @param points Array of points (vertices).
/// @param V Total number of points (vertices).
void createGraphFromPoints(struct Graph* graph, struct Point points[], int V) {
    for (int i = 0; i < V; ++i) {
        for (int j = i + 1; j < V; ++j) {
            int weight = calculateWeight(points[i], points[j]);
            addEdge(graph, i, j, weight);
        }
    }
}

/// Computes the total cost of the MST based on the parent relationships.
/// Iterates through each vertex (except the root) and sums the selected edge weights.
///
/// @param parent Array representing the MST; parent[i] is the parent of vertex i.
/// @param graph Pointer to the Graph.
/// @return The total cost (sum of edge weights) of the MST.
int calculateMSTCost(int parent[], struct Graph* graph) {
    int cost = 0;
    for (int i = 1; i < graph->V; ++i) {
        struct AdjListNode* pCrawl = graph->array[i].head;
        while (pCrawl != NULL) {
            if (pCrawl->dest == parent[i])
                cost += pCrawl->weight;
            pCrawl = pCrawl->next;
        }
    }
    return cost;
}

/// Implements the naive version of Prim's algorithm to compute the MST.
/// Uses linear search to pick the next minimal edge, updating keys and parent information.
///
/// @param graph Pointer to the Graph.
void PrimMST(struct Graph* graph) {
    int V = graph->V;
    int parent[V]; // Array to store the constructed MST
    int key[V];    // Minimum weight to reach each vertex
    int inMST[V];  // Boolean array to track vertices included in MST

    // Initialize all keys as infinite and marks all vertices as not yet included.
    for (int v = 0; v < V; ++v) {
        key[v] = INT_MAX;
        inMST[v] = 0;
    }

    key[0] = 0;    // Start from the first vertex, make its key 0 to pick it first.
    parent[0] = -1; // Root node has no parent.

    // The MST will have V-1 edges.
    for (int count = 0; count < V - 1; ++count) {
        int u = -1;

        // Choose the vertex with the minimum key value not yet included in MST.
        for (int v = 0; v < V; ++v) {
            if (!inMST[v] && (u == -1 || key[v] < key[u]))
                u = v;
        }

        inMST[u] = 1; // Include vertex u in MST.

        // Update key and parent for adjacent vertices of the picked vertex.
        struct AdjListNode* pCrawl = graph->array[u].head;
        while (pCrawl != NULL) {
            int v = pCrawl->dest;
            // Only update if v is not yet in MST and the new weight is lower than current key.
            if (!inMST[v] && pCrawl->weight < key[v]) {
                key[v] = pCrawl->weight;
                parent[v] = u;
            }
            pCrawl = pCrawl->next;
        }
    }

    // Output the edges of the MST.
    printArr(parent, V);

    // Calculate and print the total cost of the MST.
    int mstCost = calculateMSTCost(parent, graph);
    printf("Cost of Minimum Spanning Tree: %d\n", mstCost);
}

/// Main function drives the MST computation.
/// Reads the number of nodes from the user, generates a random graph from points, applies Prim's algorithm,
/// and reports the resulting MST along with the execution time.
int main() {
    int V;  // Number of vertices (points) in the graph.
    cout << "Enter number of nodes: ";
    cin >> V;
    
    // Create a graph with V vertices.
    struct Graph* graph = createGraph(V);

    // Generate random points for vertices and build a complete graph.
    struct Point points[V];
    generateRandomPoints(points, V);
    createGraphFromPoints(graph, points, V);

    // Optionally, print the generated graph structure
    // printGraph(graph);

    // Measure execution time of Prim's algorithm using microsecond resolution.
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // Compute the MST using Prim's algorithm.
    PrimMST(graph);

    gettimeofday(&end, NULL);

    // Compute elapsed time in microseconds.
    long seconds = end.tv_sec - start.tv_sec;
    long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);

    printf("Execution Time: %ld microseconds\n", micros);

    return 0;
}
// End of Naive Prim's Algorithm implementation.