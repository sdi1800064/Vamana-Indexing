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


void save_neighbors_to_file(const char *filename, Neighbor **all_neighbors, int num_queries, int num_neighbors,
                            int num_dimensions) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        perror("Error opening file");
        return;
    }

    // Write the number of dimensions
    if (fwrite(&num_dimensions, sizeof(int), 1, file) != 1) {
        perror("Error writing number of dimensions");
        fclose(file);
        return;
    }

    // Write each set of neighbors to the file
    for (int q = 0; q < num_queries; q++) {
        for (int i = 0; i < num_neighbors; i++) {
            // Write the index
            if (fwrite(&all_neighbors[q][i].index, sizeof(int), 1, file) != 1) {
                perror("Error writing neighbor index");
                fclose(file);
                return;
            }

            // Write the distance
            if (fwrite(&all_neighbors[q][i].distance, sizeof(float), 1, file) != 1) {
                perror("Error writing neighbor distance");
                fclose(file);
                return;
            }
        }
    }

    fclose(file);
}

int compare_neighbors(const void *a, const void *b) {
    Neighbor *neighborA = (Neighbor *) a;
    Neighbor *neighborB = (Neighbor *) b;
    return (neighborA->distance > neighborB->distance) - (neighborA->distance < neighborB->distance);
}

Neighbor **find_closest_neighbors(DatasetInfo *dataset_info, QueryInfo *query_info, int *actual_neighbors_count) {
    Neighbor **all_neighbors = (Neighbor **) malloc(query_info->num_queries * sizeof(Neighbor *));
    for (int q = 0; q < query_info->num_queries; q++) {
        all_neighbors[q] = (Neighbor *) malloc(K * sizeof(Neighbor));
    }
    int count = 0;

    for (int q = 0; q < query_info->num_queries; q++) {
        QueryPoint query = query_info->queries[q];

        // Initialize neighbors for each query
        for (int i = 0; i < K; i++) {
            all_neighbors[q][i].index = -1;
            all_neighbors[q][i].distance = FLT_MAX;
        }

        for (int i = 0; i < dataset_info->num_vectors; i++) {
            if (query.query_type == 1) {
                if (dataset_info->datapoints[i].category == query.v) {
                    float distance = squared_euclidean_distance(query.query_vector, dataset_info->datapoints[i].vectors,
                                                                100);
                    if (distance < all_neighbors[q][K - 1].distance) {
                        all_neighbors[q][K - 1].index = i;
                        all_neighbors[q][K - 1].distance = distance;
                        qsort(all_neighbors[q], K, sizeof(Neighbor), compare_neighbors);
                        if (count < K) {
                            count++;
                        }
                    }
                }
            } else {
                float distance = squared_euclidean_distance(query.query_vector, dataset_info->datapoints[i].vectors,
                                                            100);
                if (distance < all_neighbors[q][K - 1].distance) {
                    all_neighbors[q][K - 1].index = i;
                    all_neighbors[q][K - 1].distance = distance;
                    qsort(all_neighbors[q], K, sizeof(Neighbor), compare_neighbors);
                    if (count < K) {
                        count++;
                    }
                }
            }
        }
    }

    *actual_neighbors_count = count;
    return all_neighbors;
}

void print_neighbors(Neighbor **all_neighbors, int num_queries) {
    for (int q = 0; q < num_queries; q++) {
        printf("Query %d:\n", q);
        for (int i = 0; i < K; i++) {
            printf("Index: %d, Distance: %f\n", all_neighbors[q][i].index, all_neighbors[q][i].distance);
        }
    }
}

int main(int argc, char *argv[]) {

    uint32_t base_num_vectors;
    int query_num_vectors;
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
    QueryInfo *query_info = read_query_dataset(query_file_name, &query_num_vectors);
    print_query_dataset(query_info);
    print_dataset(dataset_info);
    int actual_neighbors_count;
    Neighbor **all_neighbors = find_closest_neighbors(dataset_info, query_info, &actual_neighbors_count);
    print_neighbors(all_neighbors, query_info->num_queries);
    save_neighbors_to_file("neighbors.ivecs", all_neighbors, query_info->num_queries, actual_neighbors_count, 1);
    free_dataset(dataset_info);
    free_query_dataset(query_info);

    return 0;
}