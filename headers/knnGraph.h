// knnGraph.h

#ifndef KNNGRAPH_H
#define KNNGRAPH_H

// Structure to store the distance and the index of a neighbor
typedef struct {
    double distance;
    int index;
} Neighbor;

int compare_neighbors(const void *a, const void *b);
void find_k_nearest_neighbors(float **vectors, int num_vectors, int dim, int vector_idx, Neighbor neighbors[]);
void build_knn_graph(int **knn_graph, float **vectors, int num_vectors, int dim, int k, FILE *outputfd);
void find_k_nearest_neighbors_for_point(float **vectors, int num_vectors, int dim, float *query_point, int k, int **knn_graph);

#endif