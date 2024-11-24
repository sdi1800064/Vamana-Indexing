#include <float.h>
#include "../headers/fvecs.h"
#include "stdio.h"
#include "string.h"
#include "structs.h"
#include "stdlib.h"
#include "math.h"
#define K 100

typedef struct {
    int index;
    float distance;
} Neighbor;


int compare_neighbors(const void *a, const void *b) {
    Neighbor *neighborA = (Neighbor *)a;
    Neighbor *neighborB = (Neighbor *)b;
    return (neighborA->distance > neighborB->distance) - (neighborA->distance < neighborB->distance);
}

Neighbor* find_closest_neighbors(DatasetInfo *dataset_info, QueryInfo *query_info, int *actual_neighbors_count) {
    QueryPoint query;
    for (int i = 0; i < query_info->num_queries; i++) {
        if (query_info->queries[i].v != -1) {
            query = query_info->queries[i];
            break;
        }
    }
    Neighbor *neighbors = (Neighbor *)malloc(K * sizeof(Neighbor));
    for (int i = 0; i < K; i++) {
        neighbors[i].index = -1;
        neighbors[i].distance = FLT_MAX;
    }

    int count = 0;
    for (int i = 0; i < dataset_info->num_vectors; i++) {
        if (dataset_info->datapoints[i].category == query.v) {
            float distance = squared_euclidean_distance(query.query_vector, dataset_info->datapoints[i].vectors, 100);
            if (distance < neighbors[K-1].distance) {
                neighbors[K-1].index = i;
                neighbors[K-1].distance = distance;
                qsort(neighbors, K, sizeof(Neighbor), compare_neighbors);
                if (count < K) {
                    count++;
                }
            }
        }
    }

    *actual_neighbors_count = count;
    return neighbors;
}

void print_neighbors(Neighbor *neighbors) {
    for (int i = 0; i < K; i++) {
        printf("Index: %d, Distance: %f\n", neighbors[i].index, neighbors[i].distance);
    }
}

int main(int argc, char *argv[])
{

    float** base_vectors;
    uint32_t base_num_vectors;
    int base_num_dimensions;

    // Variables for Query file
    float** query_vectors;
    int query_num_vectors;
    int query_num_dimensions;

    // Variables for GroundTruth
    int** groundtruth_vectors;
    int groundtruth_num_vectors;
    int groundtruth_num_dimensions;
    filterInfo *filters = NULL;


    // Initialize variables
    char *base_file_name = NULL;
    char *query_file_name = NULL;
    char *groundtruth_file_name = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            base_file_name = argv[i + 1];
            i++;  // Move past the flag and value
        } else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
            query_file_name = argv[i + 1];
            i++;
        } else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
            groundtruth_file_name = argv[i + 1];
            i++;
        }
    }

    DatasetInfo *dataset_info = read_dataset(base_file_name, &base_num_vectors, filters);
    QueryInfo *query_info = read_query_dataset(query_file_name, &query_num_vectors);\
    //print_query_dataset(query_info);
    print_dataset(dataset_info);
    int actual_neighbors_count;
    Neighbor *neighbors = find_closest_neighbors(dataset_info, query_info, &actual_neighbors_count);
    print_neighbors(neighbors);
    free_dataset(dataset_info);
    return 0;
}