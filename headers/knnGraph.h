// knnGraph.h

#ifndef KNNGRAPH_H
#define KNNGRAPH_H

// Structure to store the distance and the index of a neighbor
typedef struct {
    double distance;
    int index;
} Neighbor;

double squared_euclidean_distance(float *v1, float *v2, int dim);
int compare_neighbors(const void *a, const void *b);
void bf_find_k_nearest_neighbors(float **vectors, int num_vectors, int dim, int vector_idx, Neighbor neighbors[]);
void bf_build_knn_graph(int **knn_graph, float **vectors, int num_vectors, int dim, int k, FILE *outputfd);
void bf_find_k_nearest_neighbors_for_point(float **vectors, int num_vectors, int dim, float *query_point, int k, int **knn_graph);

#endif