#include <stdio.h>
#include <stdlib.h>

#include "../headers/fvecs.h"
#include "../headers/knnGraph.h"
#include "../headers/graph.h"

#define BASE_FILENAME "testSets/siftsmall/siftsmall_base.fvecs"
#define QUERY_FILENAME "testSets/siftsmall/siftsmall_query.fvecs"
#define GROUNDTRUTH_FILENAME "testSets/siftsmall/siftsmall_groundtruth.ivecs"

int main() {
    // Variables for Base file
    float** base_vectors;
    int base_num_vectors;
    int base_num_dimensions;

    // Variables for Query file
    float** query_vectors;
    int query_num_vectors;
    int query_num_dimensions;

    // Variables for GroundTruth
    float** groundtruth_vectors;
    int groundtruth_num_vectors;
    int groundtruth_num_dimensions;

    // Get Folder name

    // Get user's Base-file

    // Get user's Query-file

    // Get user's GroundTruth
    printf("Starting...\n");

    // Create an output file
    FILE *file_check = fopen("output.txt", "r");
    if (file_check != NULL) {
        // File exists, so close it and delete it
        printf("Existing output.txt file found. Deleting...\n");
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

    // Read the query_vectors from the file
    read_fvecs(QUERY_FILENAME, &query_vectors, &query_num_vectors, &query_num_dimensions);
    fprintf(outputfd, "Query-Vector dimensionality: %d\n", query_num_dimensions);
    fprintf(outputfd, "Number of query-Vectorts: %d\n", query_num_vectors);

    // printf("Printing query vector...\n");
    // fprintFloatVectors(query_vectors, query_num_vectors, query_num_dimensions, outputfd);

    // Read the groundtruth_vectors from the file
    read_fvecs(GROUNDTRUTH_FILENAME, &groundtruth_vectors, &groundtruth_num_vectors, &groundtruth_num_dimensions);
    fprintf(outputfd, "Groundtruth-Vector dimensionality: %d\n", groundtruth_num_dimensions);
    fprintf(outputfd, "Number of Groundtruth-Vectorts: %d\n", groundtruth_num_vectors);


    printf("Creating the random graph...\n");
    Graph *base_graph = create_random_graph(base_vectors, base_num_dimensions, 20, base_num_vectors);

    //============== UNCOMMENT THIS TO PRINT THE RANDOM GRAPH =================//
    // fprint_graph(base_graph, base_num_vectors, outputfd);
    
    int *V;                     // V - Visited List
    int L = 10;                 // L - max number of candidates
    int R = 5;                  // R - max radius for RobustPrune
    int V_size = 0;             // V_size - Size of V
    int start_index = 0;        // Start point
    int query_index = 40;       // Query point < query_vectors
    float a = 1;                // a - RobustPrune
    int *l;                     // L~ - k nearest neighbors list
    int k = 4;                  // k - number of nearest neighbor


    printf("Allocating memory for V and l...\n");
    // Allocate memory for V
    V = (int*)malloc((V_size + 1) * sizeof(int));
    if (!V) {
        printf("Memory allocation of V failed!\n");
        exit(1);
    }

    // L~ - L nearest neighbors to Xq list
    l = (int*)calloc(L, sizeof(int));
    if (!l) {
        printf("Memory allocation of l failed!\n");
        exit(1);
    }
    printf("Memory allocated\n");

    printf("V_size before entering GreedySearch: %d\n", V_size);
    GreedySearch(base_graph, base_num_dimensions, &V, &V_size, l, L, query_vectors[query_index], k, start_index);

    // printf("Printing V and l...\n");
    // for( int i = 0; i < V_size; i++){
    //     printf("%d ", V[i]);
    // }
    // printf("\n");
    fprintf(outputfd, "V: \n");
    for (int i = 0; i < V_size; i++) {
        fprintf(outputfd, "%d,", V[i]);
    }
    fprintf(outputfd, "\n");

    printf("Printing L closest to Xq...\n");
    fprintf(outputfd,"Printing L closest to Xq...\n");
    fprintf(outputfd, "Xq: ");
    for (int i = 0; i < base_num_dimensions; i++) {
        fprintf(outputfd, "%.1f ", query_vectors[query_index][i]);
    }
    fprintf(outputfd, "\n");

    for (int i = 0; i < k; i++) {
        fprintf(outputfd, "%d: ", l[i]);
        for (int j = 0; j < base_num_dimensions; j++){
            fprintf(outputfd, "%.1f ", base_graph->points[l[i]].coordinates[j]);
        }
        fprintf(outputfd, "\n");
    }

    printf("Exiting the program..\n");

    // Close the file
    fclose(outputfd);

    return 0;
}
