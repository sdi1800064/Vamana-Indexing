#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <math.h>

#define INF 1e9         // A large number to represent "infinity"

// Struct for a single edge
typedef struct {
    int to;  // Index of the point this edge connects to
} Edge;

// Struct for a single point
typedef struct {
    int index;              // Index of the point
    float *coordinates;     // Vector coordinates
    int edge_count;         // Current number of Outgoing edges
    Edge *edges;            // Array of Outgoing edges (connections to other points)
    int N_in_count;         // Current number of incoming edges
    Edge *N_in;             // Array of incoming edges
} Point;

// Struct for the whole graph
typedef struct {
    Point *points; // Array of points
    int num_dimensions; // Number of dimensions of each point
    int num_vectors; // Number of vectors in the graph
} Graph;


// Function prototypes
Graph* create_random_graph(float **base_vectors, int base_num_dimensions, int max_edges, int base_num_vectors);
void add_random_edges(Graph* graph, int max_edges);
void fprint_graph_coordinates(Graph* graph, FILE *outputfd);
void fprint_graph(Graph* graph, FILE *outputfd);
void robust_prune(Graph* graph, int point_index, int *V, int *V_size, float a, int R);

void GreedySearch(Graph *graph, int **V, int *V_size, int *l, int L, float *Xq, int k, int Xs);
void insert_closest(int *l, float *distances, int new_node, float new_distance, int L);
int arrayContains(int *V, int V_size, int node);
double squared_euclidean_distance(float *p, float *q, int n);

int calculate_medoid(Graph *graph, int *sample_point_indexes, int num_sample_points);
int* sample_points(int num_sample_points);



#endif