// Language: C++
//
// File: Efficient_Prim.cpp
// Description: Efficient implementation of Prim's algorithm to compute the Minimum Spanning Tree (MST)
//              for a randomly generated weighted graph. The code calculates the MST cost and tracks execution time,
//              providing insight into the performance of the algorithm.
// Author: [Gova]

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <bits/stdc++.h>

using namespace std;

// Structure representing a 2D point used to compute edge weights from Euclidean distance.
struct Point {
    int x, y;
};

// Node in the adjacency list for a graph vertex.
// Each node represents an edge to a neighboring vertex with a weight.
struct AdjListNode {
    int dest;                       // Destination vertex of the edge.
    int weight;                     // Weight of the edge (e.g., Euclidean distance).
    struct AdjListNode* next;       // Pointer to the next node in the adjacency list.
};

// Linked list to hold all adjacency list nodes for a given graph vertex.
struct AdjList {
    struct AdjListNode* head;       // Head pointer for the linked list.
};

// Graph structure represented as an array of adjacency lists.
// V indicates the number of vertices.
struct Graph {
    int V;                          // Total number of vertices in the graph.
    struct AdjList* array;          // Array of adjacency lists.
};

// Structure representing a node in the minimum heap used for Prim's algorithm.
// Each node contains a vertex (v) and its key value for MST computation.
struct MinHeapNode {
    int v;                          // Vertex index.
    int key;                        // Current minimum weight to connect the vertex.
};

// Structure representing the minimum heap data structure.
// Used to quickly extract the vertex with the minimum key value.
struct MinHeap {
    int size;                       // Current number of elements in the heap.
    int capacity;                   // Maximum capacity of the heap.
    int* pos;                       // Position map to track indices of vertices in the heap array.
    struct MinHeapNode** array;     // Array of pointers to heap nodes.
};

/**
 * Creates a new 2D point.
 *
 * @param x Horizontal coordinate.
 * @param y Vertical coordinate.
 * @return Pointer to the newly created Point.
 */
struct Point* createPoint(int x, int y) {
    struct Point* point = (struct Point*)malloc(sizeof(struct Point));
    point->x = x;
    point->y = y;
    return point;
}

/**
 * Creates a new adjacency list node for a given destination and weight.
 *
 * @param dest Destination vertex index.
 * @param weight Weight associated with the edge.
 * @return Pointer to the new adjacency list node.
 */
struct AdjListNode* newAdjListNode(int dest, int weight) {
    struct AdjListNode* newNode = (struct AdjListNode*)malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->weight = weight;
    newNode->next = NULL;
    return newNode;
}

/**
 * Creates a graph with V vertices.
 *
 * @param V Number of vertices.
 * @return Pointer to the newly created Graph.
 */
struct Graph* createGraph(int V) {
    struct Graph* graph = (struct Graph*)malloc(sizeof(struct Graph));
    graph->V = V;
    graph->array = (struct AdjList*)malloc(V * sizeof(struct AdjList));

    // Initialize each adjacency list as empty.
    for (int i = 0; i < V; ++i)
        graph->array[i].head = NULL;

    return graph;
}

/**
 * Adds an undirected edge to the graph by updating the adjacency lists of both source and destination.
 *
 * @param graph Pointer to the Graph.
 * @param src Source vertex.
 * @param dest Destination vertex.
 * @param weight Weight of the edge.
 */
void addEdge(struct Graph* graph, int src, int dest, int weight) {
    // Insert node for destination in source list.
    struct AdjListNode* newNode = newAdjListNode(dest, weight);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;

    // Since the graph is undirected, add an edge from destination to source as well.
    newNode = newAdjListNode(src, weight);
    newNode->next = graph->array[dest].head;
    graph->array[dest].head = newNode;
}

/**
 * Creates a new minimum heap node.
 *
 * @param v Vertex index.
 * @param key Initial key value.
 * @return Pointer to the new MinHeapNode.
 */
struct MinHeapNode* newMinHeapNode(int v, int key) {
    struct MinHeapNode* minHeapNode = (struct MinHeapNode*)malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->key = key;
    return minHeapNode;
}

/**
 * Creates a minimum heap with a specified capacity.
 *
 * @param capacity Maximum capacity for the heap.
 * @return Pointer to the new MinHeap.
 */
struct MinHeap* createMinHeap(int capacity) {
    struct MinHeap* minHeap = (struct MinHeap*)malloc(sizeof(struct MinHeap));
    minHeap->pos = (int*)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (struct MinHeapNode**)malloc(capacity * sizeof(struct MinHeapNode*));
    return minHeap;
}

/**
 * Swaps two nodes of the min heap.
 *
 * @param a Pointer to the first node.
 * @param b Pointer to the second node.
 */
void swapMinHeapNode(struct MinHeapNode** a, struct MinHeapNode** b) {
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

/**
 * Maintains the heap property by heapifying at a given index.
 *
 * @param minHeap Pointer to the min heap.
 * @param idx Index at which to maintain the heap property.
 */
void minHeapify(struct MinHeap* minHeap, int idx) {
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    // Compare left child with current smallest.
    if (left < minHeap->size && minHeap->array[left]->key < minHeap->array[smallest]->key)
        smallest = left;

    // Compare right child with current smallest.
    if (right < minHeap->size && minHeap->array[right]->key < minHeap->array[smallest]->key)
        smallest = right;

    // If the smallest is not the current node, swap and continue heapifying.
    if (smallest != idx) {
        struct MinHeapNode* smallestNode = minHeap->array[smallest];
        struct MinHeapNode* idxNode = minHeap->array[idx];

        // Update positions in the heap.
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);
        minHeapify(minHeap, smallest);
    }
}

/**
 * Checks if the min heap is empty.
 *
 * @param minHeap Pointer to the min heap.
 * @return true if empty, false otherwise.
 */
bool isEmpty(struct MinHeap* minHeap) {
    return minHeap->size == 0;
}

/**
 * Extracts the minimum key node from the heap.
 *
 * @param minHeap Pointer to the min heap.
 * @return The minimum key heap node.
 */
struct MinHeapNode* extractMin(struct MinHeap* minHeap) {
    if (isEmpty(minHeap))
        return NULL;

    // Remove the root of the heap (min element).
    struct MinHeapNode* root = minHeap->array[0];
    struct MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update positions after extraction.
    minHeap->pos[root->v] = minHeap->size - 1;
    minHeap->pos[lastNode->v] = 0;

    --minHeap->size;
    minHeapify(minHeap, 0);
    return root;
}

/**
 * Decreases the key value of a given vertex in the min heap.
 * This function is crucial to update the cost of adding a vertex to the MST.
 *
 * @param minHeap Pointer to the min heap.
 * @param v Vertex for which the key is to be decreased.
 * @param key New key value.
 */
void decreaseKey(struct MinHeap* minHeap, int v, int key) {
    int i = minHeap->pos[v];
    minHeap->array[i]->key = key;

    // Traverse up the heap tree until the heap property is satisfied.
    while (i && minHeap->array[i]->key < minHeap->array[(i - 1) / 2]->key) {
        minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
        minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

/**
 * Determines whether a given vertex is in the min heap.
 *
 * @param minHeap Pointer to the min heap.
 * @param v Vertex index.
 * @return true if the vertex is in the heap, false otherwise.
 */
bool isInMinHeap(struct MinHeap* minHeap, int v) {
    return minHeap->pos[v] < minHeap->size;
}

/**
 * Prints the MST stored in the parent array.
 * This function outputs the edge between the parent and vertex.
 *
 * @param arr Parent array where index represents a vertex and value its parent.
 * @param n Number of vertices.
 */
void printArr(int arr[], int n) {
    for (int i = 1; i < n; ++i)
        printf("%d - %d\n", arr[i], i);
}

/**
 * Calculates the total cost (sum of edge weights) of the MST.
 *
 * @param parent Array representing the MST structure; parent[i] gives the parent vertex of i.
 * @param graph Pointer to the graph.
 * @return The total cost of the MST.
 */
int calculateMSTCost(int parent[], struct Graph* graph) {
    int cost = 0;
    // Start from vertex 1 since vertex 0 is the root of the MST.
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

/**
 * Generates a random graph by creating random points and adding edges with weights
 * computed from the Euclidean distance between point pairs.
 *
 * @param graph Pointer to the Graph.
 * @param numPoints Number of vertices (points) to generate.
 */
void generateRandomGraph(struct Graph* graph, int numPoints) {
    struct Point** points = (struct Point**)malloc(numPoints * sizeof(struct Point*));

    // Generate random points within a bounded coordinate plane.
    for (int i = 0; i < numPoints; ++i) {
        points[i] = createPoint(rand() % 100, rand() % 100);
    }

    // For each unique pair of points, compute the Euclidean distance and add an undirected edge.
    for (int i = 0; i < numPoints; ++i) {
        for (int j = i + 1; j < numPoints; ++j) {
            int dist = (int)sqrt(pow(points[i]->x - points[j]->x, 2) + 
                                  pow(points[i]->y - points[j]->y, 2));
            addEdge(graph, i, j, dist);
        }
    }

    // Free memory allocated for the temporary points.
    for (int i = 0; i < numPoints; ++i) {
        free(points[i]);
    }
    free(points);
}

/**
 * Prints the entire graph as an adjacency list.
 *
 * @param graph Pointer to the Graph.
 */
void printGraph(struct Graph* graph) {
    for (int i = 0; i < graph->V; ++i) {
        struct AdjListNode* pCrawl = graph->array[i].head;
        printf("Adjacency list of vertex %d:\n", i);
        while (pCrawl) {
            printf("-> %d(%d) ", pCrawl->dest, pCrawl->weight);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}

/**
 * Executes Prim's algorithm to compute the MST of the graph.
 * It uses a min heap to efficiently select the next vertex with the minimum connecting edge.
 *
 * @param graph Pointer to the Graph.
 */
void PrimMST(struct Graph* graph) {
    int V = graph->V;
    int parent[V];      // Array to store the constructed MST.
    int key[V];         // Array to store key values used to pick minimum weight edge.

    struct MinHeap* minHeap = createMinHeap(V);

    // Initialize min heap with all vertices. Set key value of vertex 0 to 0 to pick it first.
    for (int v = 1; v < V; ++v) {
        parent[v] = -1;             // No parent initially assigned.
        key[v] = INT_MAX;           // Set infinite initial key value.
        minHeap->array[v] = newMinHeapNode(v, key[v]);
        minHeap->pos[v] = v;
    }
    key[0] = 0;
    minHeap->array[0] = newMinHeapNode(0, key[0]);
    minHeap->pos[0] = 0;
    minHeap->size = V;

    // Continue until all vertices are included in the MST.
    while (!isEmpty(minHeap)) {
        struct MinHeapNode* minHeapNode = extractMin(minHeap);
        int u = minHeapNode->v;

        // Traverse all adjacent vertices of u.
        struct AdjListNode* pCrawl = graph->array[u].head;
        while (pCrawl != NULL) {
            int v = pCrawl->dest;

            /* If v is still in heap and the weight of edge u-v is less than the current key value of v,
             update key value and parent accordingly */
            if (isInMinHeap(minHeap, v) && pCrawl->weight < key[v]) {
                key[v] = pCrawl->weight;
                parent[v] = u;
                decreaseKey(minHeap, v, key[v]);
            }
            pCrawl = pCrawl->next;
        }
    }

    printf("Edges of Minimum Spanning Tree:\n");
    printArr(parent, V);

    int mstCost = calculateMSTCost(parent, graph);
    printf("Cost of Minimum Spanning Tree: %d\n", mstCost);
}

/**
 * Main function to drive the MST computation.
 * It reads the number of nodes from the user, generates a random graph, prints the graph,
 * and displays the MST along with the execution time.
 */
int main() {
    int numPoints;  // Desired number of nodes in the graph.
    cout << "Enter number of nodes: ";
    cin >> numPoints;
    
    // Create a graph with the specified number of vertices.
    struct Graph* graph = createGraph(numPoints);

    // Generate random weighted edges based on Euclidean distances.
    generateRandomGraph(graph, numPoints);

    // Debug: Output the generated graph structure.
    printf("Generated Graph:\n");
    printGraph(graph);
    printf("\n");

    // Track the execution time of Prim's algorithm using microsecond resolution.
    struct timeval start, end;
    gettimeofday(&start, NULL);
    PrimMST(graph);
    gettimeofday(&end, NULL);

    long seconds = end.tv_sec - start.tv_sec;
    long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);

    printf("\nExecution Time: %ld microseconds\n", micros);
    return 0;
}
// End of Efficient Prim's algorithm implementation.