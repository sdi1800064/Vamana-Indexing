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
    int index;      // Index of the point
    float *coordinates;  // Vector coordinates
    int edge_count; // Current number of edges
    Edge *edges; // Array of edges (connections to other points)
} Point;

// Struct for the whole graph
typedef struct {
    Point *points; // Array of points
} Graph;


// Function prototypes
Graph* create_random_graph(float **base_vectors, int base_num_dimensions, int max_edges, int base_num_vectors);
void add_random_edges(Graph* graph, int max_edges, int base_num_vectors);
void fprint_graph_coordinates(Graph* graph, int base_num_vectors, int base_num_coordinates, FILE *outputfd);
void fprint_graph(Graph* graph, int base_num_vectors, FILE *outputfd);
void robust_prune(Graph* graph, int point_index, int *V, int *V_size, float a, int R, int base_num_dimensions);

void GreedySearch(Graph *graph, int dimensions, int **V, int *V_size, int *l, int L, float *Xq, int k, int Xs);
void insert_closest(int *l, float *distances, int new_node, float new_distance, int L);
int V_contains(int **V, int V_size, int node);
void greedy_search(Graph* graph, int base_num_vectors, int base_num_dimensions, int start_index, int query_index, int k, int L, int **V, int *V_size, int *l);
void add_candidate(int* candidate_list, int* candidate_count, int candidate_index);
void prune_V(Graph* graph, int** V, int* V_size, int base_num_dimensions, int L, int query_index);


#endif