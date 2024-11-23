#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#define INF 1e9         // A large number to represent "infinity"

// ==================== Structures ==================== //

// Define a struct for the result
typedef struct {
    uint32_t num_vectors;
    float category;  // Array for discretized categorical attribute C
    float timestamp;  // Array for normalized timestamp attribute T
    float vectors[100];     // Array for the remaining 100-dimensional vectors
} DatasetInfo;

// Define the QueryInfo structure
typedef struct {
    uint32_t num_queries;    // Number of queries
    float query_type;        // Query type (0, 1, 2, or 3)
    float v;                 // Specific query value for categorical attribute (or -1)
    float l;                 // Lower bound for timestamp attribute (or -1)
    float r;                 // Upper bound for timestamp attribute (or -1)
    float query_vector[100]; // 100-dimensional query vector
} QueryInfo;

typedef struct{
    int *filters[2];
    int num_indexes;
    int filters_size;
} filterInfo;

// Struct for a single point
typedef struct {
    int index;              // Index of the point
    float *coordinates;     // Vector coordinates
    int edge_count;         // Current number of Outgoing edges
    int *edges;            // Array of Outgoing edges (connections to other points)
} Point;

// Struct for the whole graph
typedef struct {
    Point *points; // Array of points
    int num_dimensions; // Number of dimensions of each point
    int num_points; // Number of vectors in the graph
} Graph;

// ==================== Structures ==================== //

// Function prototypes
double squared_euclidean_distance(float *p, float *q, int n);

Graph* create_random_graph(float **base_vectors, int base_num_dimensions, int max_edges, int base_num_vectors);

void add_random_edges(Graph* graph, int max_edges);
void fprint_graph_coordinates(Graph* graph, FILE *outputfd);
void fprint_graph(Graph* graph, FILE *outputfd);

void addEdge(Point *point, int toIndex);
int edgeExists(Point *point, int toIndex);
void robustPrune(Graph *graph, int p_index, int *V, int V_size, float a, int R);

int *get_the_difference(int *Lamda, int Lamda_size, int *V, int V_size, int *Lamda_minus_V_size);

void swap(int *a, int *b);
void swap_float(float *a, float *b);
void sort_array(Graph *graph, int *l_temp, int l_temp_size, float *Xq);
void sort_filter_array(int *array[2], int size);
void printArray(int *array, int array_size);
int arrayContains(int *V, int V_size, int node);
void add_to_dynamic_array(int **array, int *size, int element);


void greedy_search(Graph *graph, float *Xq, int start_index, int **V, int *V_size, int **Lamda, int *Lamda_size, int L);

int calculate_medoid(Graph *graph, int *sample_point_indexes, int num_sample_points);
int* sample_points(int max, int num_sample_points);

void check_for_duplicates(int *array, int size);
void vamana_indexing(Graph *graph, int L, float a, int R);

#endif