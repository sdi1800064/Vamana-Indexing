#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "../headers/structs.h"

#define INF 1e9         // A large number to represent "infinity"

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