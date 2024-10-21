#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../headers/knnGraph.h"


// Function to calculate squared Euclidean distance between two vectors
double squared_euclidean_distance(float *p, float *q, int n) {
    // printf("Coordinates located at %d and point %d\n", &p, &q);
    float sum = 0.0f;
    for (int i = 0; i < n; i++) {
        float diff = p[i] - q[i];  // Calculate the difference
        sum += diff * diff;        // Square the difference and add to sum
    }
    return sqrtf(sum);  // Use sqrtf for floats
}

// Function to compare two neighbors (used for sorting)
int compare_neighbors(const void *a, const void *b) {
    Neighbor *neighbor1 = (Neighbor *)a;
    Neighbor *neighbor2 = (Neighbor *)b;
    if (neighbor1->distance < neighbor2->distance) return -1;
    else if (neighbor1->distance > neighbor2->distance) return 1;
    else return 0;
}

// Function to find the K nearest neighbors of a given vector
void bf_find_k_nearest_neighbors(float **vectors, int num_vectors, int dim, int vector_idx, Neighbor neighbors[]) {
    // Calculate distances from vector_idx to all other vectors
    for (int i = 0; i < num_vectors; i++) {
        if (i == vector_idx) {
            neighbors[i].distance = INFINITY;  // Ignore self
        } else {
            neighbors[i].distance = squared_euclidean_distance(vectors[vector_idx], vectors[i], dim);   // Calculate the euclidean distance
        }
        neighbors[i].index = i;  // Store the index of the neighbor
    }

    // Sort the neighbors by distance
    qsort(neighbors, num_vectors, sizeof(Neighbor), compare_neighbors);
}

// Function to build the KNN graph
void bf_build_knn_graph(int **knn_graph, float **vectors, int num_vectors, int dim, int k, FILE *outputfd) {
    Neighbor *neighbors = (Neighbor *)malloc(num_vectors * sizeof(Neighbor));


    fprintf(outputfd, "%d-Nearest Neighbors for the each vector : \n", k);
    // For each vector, find the K nearest neighbors
    for (int i = 0; i < num_vectors; i++) {
        bf_find_k_nearest_neighbors(vectors, num_vectors, dim, i, neighbors);

        // Store the indices of the K nearest neighbors in the graph
        for (int j = 0; j < k; j++) {
            knn_graph[i][j] = neighbors[j].index;
        }

    }

    // Free the memory for neighbors
    free(neighbors);

}

// Function to find K nearest neighbors of a given coordinate
void bf_find_k_nearest_neighbors_for_point(float **vectors, int num_vectors, int dim, float *query_point, int k, int **knn_graph) {
    // Array to hold distances
    double *distances = (double *)malloc(num_vectors * sizeof(double));
    int *indices = (int *)malloc(num_vectors * sizeof(int));

    // Calculate distances from the query point to all vectors
    for (int i = 0; i < num_vectors; i++) {
        distances[i] = squared_euclidean_distance(vectors[i], query_point, dim);
        indices[i] = i;  // Store original index
    }

    // Sort distances and corresponding indices
    for (int i = 0; i < num_vectors - 1; i++) {
        for (int j = 0; j < num_vectors - i - 1; j++) {
            if (distances[j] > distances[j + 1]) {
                // Swap distances
                double temp_dist = distances[j];
                distances[j] = distances[j + 1];
                distances[j + 1] = temp_dist;

                // Swap corresponding indices
                int temp_index = indices[j];
                indices[j] = indices[j + 1];
                indices[j + 1] = temp_index;
            }
        }
    }

    // Store the K nearest neighbors in the knn_graph
    for (int i = 0; i < k; i++) {
        knn_graph[0][i] = indices[i];  // Store indices of the nearest neighbors
    }

    // Free the allocated memory
    free(distances);
    free(indices);
}
