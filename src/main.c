#include <stdio.h>
#include <stdlib.h>

#include "../headers/fvecs.h"
#include "../headers/knnGraph.h"

#define BASE_FILENAME "testSets/siftsmall/siftsmall_base.fvecs"
#define QUERY_FILENAME "testSets/siftsmall/siftsmall_query.fvecs"
#define GROUNDTRUTH_FILENAME "testSets/siftsmall/siftsmall_groundtruth.ivecs"

int main() {
    float** base_vectors;
    int base_num_vectors;
    int base_num_dimensions;

    float** query_vectors;
    int query_num_vectors;
    int query_num_dimensions;

    float** groundtruth_vectors;
    int groundtruth_num_vectors;
    int groundtruth_num_dimensions;

    // Get Folder name

    // Get user's Base-file

    // Get user's Query-file

    // Get user's GroundTruth


    // Create an output file
    FILE *file_check = fopen("output.txt", "r");
    if (file_check != NULL) {
        // File exists, so close it and delete it
        fclose(file_check);
        if (remove("output.txt") == 0) {
            printf("Existing output.txt file deleted successfully.\n");
        } else {
            printf("Error: Could not delete existing output.txt file!\n");
            return 1;
        }
    }

    // Create a new file (overwrite or create anew)
    FILE *outputfd = fopen("output.txt", "w");
    if (outputfd == NULL) {
        printf("Error: Could not create output.txt file!\n");
        return 1;
    }


    // Read the base_vectors from the file
    read_fvecs(BASE_FILENAME, &base_vectors, &base_num_vectors, &base_num_dimensions);
    fprintf(outputfd, "Base-Vector dimensionality: %d\n", base_num_dimensions);
    fprintf(outputfd, "Number of Base-Vectorts: %d\n", base_num_vectors);

    read_fvecs(QUERY_FILENAME, &query_vectors, &query_num_vectors, &query_num_dimensions);
    fprintf(outputfd, "Query-Vector dimensionality: %d\n", query_num_dimensions);
    fprintf(outputfd, "Number of query-Vectorts: %d\n", query_num_vectors);

    read_fvecs(GROUNDTRUTH_FILENAME, &groundtruth_vectors, &groundtruth_num_vectors, &groundtruth_num_dimensions);
    fprintf(outputfd, "Groundtruth-Vector dimensionality: %d\n", groundtruth_num_dimensions);
    fprintf(outputfd, "Number of Groundtruth-Vectorts: %d\n", groundtruth_num_vectors);


    // Example: Print the first vector's components
    // fprintf(outputfd, "Printing vector coordinates:\n");
    // for (int i = 0; i < base_num_vectors; i++) {
    //     fprintf(outputfd, "%d : ", i+1);
    //     for(int j=0; j < base_num_dimensions; j++){
    //             fprintf(outputfd, "%.1f, ", base_vectors[i][j]);
    //     }    
    //     fprintf(outputfd, "\n");   
    // }
    // fprintf(outputfd, "\n");

    int k = 4;
    
    // Allocate memory for the graph
    int **knn_graph = (int **)malloc(base_num_vectors * sizeof(int *));
    for (int i = 0; i < base_num_vectors; i++) {
        knn_graph[i] = (int *)malloc(k * sizeof(int));
    }

    build_knn_graph(knn_graph, base_vectors, base_num_vectors, base_num_dimensions, k, outputfd);

    // Printing each vector and its k-Nearest Neighbors
    fprintf(outputfd, "Printing KNN Graph:\n");
    for (int i = 0; i < base_num_vectors; i++) {
        fprintf(outputfd, "%d is neighbors with : ", i);
        for(int j=0; j < k; j++){
                fprintf(outputfd, "%d, ", knn_graph[i][j]);
        }    
        fprintf(outputfd, "\n");   
    }
    fprintf(outputfd, "\n");

    printf("Exiting the program..\n");

    // Free the allocated memory of vector matrix
    for (int i = 0; i < base_num_vectors; i++) {
        free(base_vectors[i]);
    }
    free(base_vectors);

    // Close the file
    fclose(outputfd);

    return 0;
}
