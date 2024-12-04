#include <float.h>
#include "../headers/dataset.h"
#include "stdio.h"
#include "string.h"
#include "structs.h"
#include "stdlib.h"
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
    fread(dimension, sizeof(int), 1, file);

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
        fread(&dim, sizeof(int), 1, file);  // Read the dimension (should always match *dimension)

        if (dim != *dimension) {
            printf("Error: Dimensionality mismatch.\n");
            exit(EXIT_FAILURE);
        }

        // Read the vector components
        fread((*vectors)[i], sizeof(int), *dimension, file);
        printf("\n");
    }

    fclose(file);
}


void save_neighbors_to_file(const char *filename, Neighbor **neighbours, int num_vectors, int dimension) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        perror("Error opening file");
        return;
    }
    for (int i = 0; i < num_vectors; i++) {

        // Write the indices of the neighbors for this query as the vector
        for (int j = 0; j < dimension; j++) {
            if (fwrite(&(neighbours[i][j].index), sizeof(int), 1, file) != 1) {
                perror("Error writing neighbor index");
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

    uint32_t base_num_vectors;
    int query_num_vectors;


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
    Neighbor **all_neighbors = find_closest_neighbors(dataset_info, query_info, &actual_neighbors_count);
    
    char* GROUNDTRUTH_FILE_NAME = "neighbors.ivecs";
    save_neighbors_to_file(GROUNDTRUTH_FILE_NAME, all_neighbors, query_info->num_queries,  K);

    int** ground_truth = readGroundTruth(GROUNDTRUTH_FILE_NAME, query_info->num_queries);

    for(int i = 0; i < query_info->num_queries; i++) {
        if(ground_truth[i][0] != -1){
            printf("Query %d of type %d has closest neighbors: ", i, query_info->queries[i].query_type);
            for(int j = 0; j < 100; j++) {
                if(ground_truth[i][j] == -1){
                    break;
                }
                printf("%d ", ground_truth[i][j]);
            }
            printf("\n\n");
        }
    }


//    int** groundtruth_vectors;
//    int groundtruth_num_vectors;
//    int groundtruth_num_dimensions;
//    read_ivecs("neighbors.ivecs", &groundtruth_vectors, &groundtruth_num_vectors, &groundtruth_num_dimensions);
    free_dataset(dataset_info);
    free_query_dataset(query_info);

    return 0;
}

