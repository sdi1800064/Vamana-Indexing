#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "../headers/structs.h"

#define INF 1e9         // A large number to represent "infinity"

// Function prototypes
double squared_euclidean_distance(float *p, float *q, int n);

Graph initialise_graph(DatasetInfo* dataset, int max_edges);
void free_graph(Graph graph);
void fprint_graph_coordinates(Graph* graph, FILE *outputfd);
void fprint_graph(Graph* graph, FILE *outputfd);

void addEdge(Point *point, int toIndex);
int edgeExists(Point *point, int toIndex);

int *get_the_difference(int *Lamda, int Lamda_size, int *V, int V_size, int *Lamda_minus_V_size);

void swap(int *a, int *b);
void swap_float(float *a, float *b);
void sort_array(Graph *graph, int *l_temp, int l_temp_size, float *Xq);
void printArray(int *array, int array_size);
int arrayContains(int *V, int V_size, int node);
void add_to_dynamic_array(int **array, int *size, int element);
void printArray(int *array, int array_size);
void check_for_duplicates(int *array, int size);
int sort_array_based_on_dataset(DatasetInfo *dataset, int *array, int array_size, float *Xq);
int arrayContainsForRobustRrune(int *V, int V_size, int node);


// Functions for the medoid calculation
int calculate_medoid(Graph *graph, int *sample_point_indexes, int num_sample_points);
int* sample_points(Graph graph, int num_sample_points);

int filtered_Robust_prune(Graph *graph, int p_index, int *V, int V_size, float a, int R);
void filtered_greedy_search(Graph *graph, float *Xq, int* start_index, int start_index_size, int **V, int *V_size, int **Lamda, int *Lamda_size, int L, int query_filter);

Graph create_random_graph(DatasetInfo dataset, int base_num_dimensions, int max_edges);
void robustPrune(Graph *graph, int p_index, int *V, int V_size, float a, int R);
void greedy_search(Graph *graph, float *Xq, int start_index, int **V, int *V_size, int **Lamda, int *Lamda_size, int L);
Graph vamana_indexing(DatasetInfo dataset, int L, float a, int R);
Graph* stitched_vamana_indexing(DatasetInfo* dataset, int L_small, float a, int R_small);
int calculate_medoid(Graph *graph, int *sample_point_indexes, int num_sample_points);
FilteredMethoidList * get_filtered_medoids(Graph *graph, int *t, filterInfo *filterInfo);
int find_medoid_for_point(FilteredMethoidList* filteredMedoids, Point* point, int medoid_index);

FilteredMedoid * findClosestDataPoints(Point **groupedData, filterInfo *groupSizes, int numCategories, int t);



// functions that can be used if we want to use multiple filters
// int compute_intersection(int *set1, int set1_size, int *set2, int set2_size, int *intersection);
// bool is_not_subset(int *intersection, int intersection_size, int *set3, int set3_size);
void check_for_duplicates(int *array, int size);
void generate_random_permutation(int *perm, int n);

Graph filtered_vamana_indexing(DatasetInfo* dataset, int L, float a, int R,filterInfo *filterinfo);

#endif