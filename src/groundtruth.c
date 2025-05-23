#include <float.h>
#include "../headers/dataset.h"
#include "stdio.h"
#include "string.h"
#include "structs.h"
#include "stdlib.h"
#include "time.h"
#define K 100

typedef struct {
    int index;
    float distance;
} Neighbor;

void read_ivecs(const char* filename, int*** vectors, int* num_vectors, int* dimension) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        printf("Error opening file.\n");
        exit(EXIT_FAILURE);
    }

    // Get the size of the file
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);

    // Read the first 4 bytes to get the dimensionality (d)
    if(fread(dimension, sizeof(int), 1, file) != 1){
        printf("Error reading dimension.\n");
        exit(EXIT_FAILURE);
    }

    // Calculate the number of vectors
    *num_vectors = file_size / ((*dimension + 1) * sizeof(int));

    // Allocate memory for the vectors
    *vectors = (int**)malloc((*num_vectors) * sizeof(int*));
    for (int i = 0; i < *num_vectors; i++) {
        (*vectors)[i] = (int*)malloc((*dimension) * sizeof(int));
    }

    // Set file pointer back to the beginning
    fseek(file, 0, SEEK_SET);

    // Read all vectors from the file
    for (int i = 0; i < *num_vectors; i++) {
        int dim;
        if(fread(&dim, sizeof(int), 1, file) != 1){  // Read the dimension (should always match *dimension)
            printf("Error reading dimension.\n");
            exit(EXIT_FAILURE);
        }

        if (dim != *dimension) {
            printf("Error: Dimensionality mismatch.\n");
            exit(EXIT_FAILURE);
        }

        // Read the vector components
        if(fread((*vectors)[i], sizeof(int), *dimension, file) != 1){
            printf("Error reading vector.\n");
            exit(EXIT_FAILURE);
        }
        printf("\n");
    }

    fclose(file);
}


void save_neighbors_to_file(const char *filename, Neighbor **neighbours, int num_vectors, int dimension) {
    printf("Saving neighbors to file %s: ", filename);
    fflush(stdout);
    FILE *file = fopen(filename, "wb");
    if (!file) {
        perror("Error opening file");
        return;
    }
    for (int i = 0; i < num_vectors; i++) {
        // Write the dimension of the vector
        // printf("%d: ", i);

        // Write the indices of the neighbors for this query as the vector
        for (int j = 0; j < dimension; j++) {
            if(&(neighbours[i][j].index) != -1){
                if (fwrite(&(neighbours[i][j].index), sizeof(int), 1, file) != 1) {
                    perror("Error writing neighbor index");
                    fclose(file);
                    return;
                }
            }
        }
    }

    fclose(file);
    printf("Done.\n");
    fflush(stdout);
}

int compare_neighbors(const void *a, const void *b) {
    Neighbor *neighborA = (Neighbor *) a;
    Neighbor *neighborB = (Neighbor *) b;
    return (neighborA->distance > neighborB->distance) - (neighborA->distance < neighborB->distance);
}

Neighbor **find_closest_neighbors(DatasetInfo *dataset_info, QueryInfo *query_info, int *actual_neighbors_count) {
    printf("Calculating groundtruth: ");
    fflush(stdout);
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
                    float distance = squared_euclidean_distance(query.query_vector, dataset_info->datapoints[i].vectors,100);

                    if (distance < all_neighbors[q][K - 1].distance) {
                        all_neighbors[q][K - 1].index = i;
                        all_neighbors[q][K - 1].distance = distance;
                        qsort(all_neighbors[q], K, sizeof(Neighbor), compare_neighbors);
                        if (count < K) {
                            count++;
                        }
                    }
                }
            } else if(query.query_type == 0) {
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
    printf("Done.\n");
    fflush(stdout);
    return all_neighbors;
}

void print_neighbors(Neighbor **all_neighbors, int num_queries) {
    for (int q = 0; q < num_queries; q++) {
        printf("Query %d:\n", q);
        for (int i = 0; i < K; i++) {
            if(all_neighbors[q][i].index == -1){
                break;
            }
            printf("%d ", all_neighbors[q][i].index);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {

    time_t t = time(NULL);
    // Initialize variables
    char *base_file_name = NULL;
    char *query_file_name = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            base_file_name = argv[i + 1];
            i++;  // Move past the flag and value
        } else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
            query_file_name = argv[i + 1];
            i++;
        }
    }



    DatasetInfo *dataset_info = read_dataset(base_file_name);
    QueryInfo *query_info = read_query_dataset(query_file_name);
    // print_query_dataset(query_info);
    // print_dataset(dataset_info);
    int actual_neighbors_count;
    time_t start_groundtruth_time = time(NULL);
    Neighbor **all_neighbors = find_closest_neighbors(dataset_info, query_info, &actual_neighbors_count);
    time_t end_groundtruth_time = time(NULL);
    // print_neighbors(all_neighbors, query_info->num_queries);
    save_neighbors_to_file("groundtruth.ivecs", all_neighbors, query_info->num_queries,  K);
    

    free_dataset(dataset_info);
    free_query_dataset(query_info);

    printf("Time to calculate groundtruth: %ld seconds || Time taken for the whole program: %ld\n", end_groundtruth_time - start_groundtruth_time, time(NULL) - t);
    printf("Exiting program\n");
    return 0;
}

