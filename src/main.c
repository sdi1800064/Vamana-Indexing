#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../headers/fvecs.h"
#include "../headers/graph.h"

int main(int argc, char *argv[]) {
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
    
    printf("Starting...\n");
    // Initialize variables
    char *base_file_name = NULL;
    char *query_file_name = NULL;
    char *groundtruth_file_name = NULL;
    int k = -1;
    int R = -1;
    float a = -1.0;
    int L = -1;

    // Iterate through command line arguments
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
        } else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) {
            k = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "-R") == 0 && i + 1 < argc) {
            R = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "-a") == 0 && i + 1 < argc) {
            a = atof(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "-L") == 0 && i + 1 < argc) {
            L = atof(argv[i + 1]);
            i++;
        }
    }

    // Check if all required arguments are provided
    if (!base_file_name || !query_file_name || !groundtruth_file_name || k == -1 || R == -1 || a == -1.0) {
        fprintf(stderr, "Error: Missing required arguments.\n");
        fprintf(stderr, "Usage: %s -b base_file_name -q query_file_name -g groundtruth_file_name -k (int k) -R (int R) -a (float a) -L (float L)\n", argv[0]);
        return 1;
    }


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
    read_fvecs(base_file_name, &base_vectors, &base_num_vectors, &base_num_dimensions);
    fprintf(outputfd, "Base-Vector dimensionality: %d\n", base_num_dimensions);
    fprintf(outputfd, "Number of Base-Vectorts: %d\n", base_num_vectors);

    // Read the query_vectors from the file
    read_fvecs(query_file_name, &query_vectors, &query_num_vectors, &query_num_dimensions);
    fprintf(outputfd, "Query-Vector dimensionality: %d\n", query_num_dimensions);
    fprintf(outputfd, "Number of query-Vectorts: %d\n", query_num_vectors);

    //============== UNCOMMENT THIS TO PRINT THE QUERY POINTS =================//
    // printf("Printing query vector...\n");
    // fprintFloatVectors(query_vectors, query_num_vectors, query_num_dimensions, outputfd);

    // Read the groundtruth_vectors from the file
    read_fvecs(groundtruth_file_name, &groundtruth_vectors, &groundtruth_num_vectors, &groundtruth_num_dimensions);
    fprintf(outputfd, "Groundtruth-Vector dimensionality: %d\n", groundtruth_num_dimensions);
    fprintf(outputfd, "Number of Groundtruth-Vectorts: %d\n", groundtruth_num_vectors);


    Graph *base_graph = create_random_graph(base_vectors, base_num_dimensions, R, base_num_vectors);

    //============== UNCOMMENT THIS TO PRINT THE RANDOM GRAPH =================//
    // fprint_graph(base_graph, base_num_vectors, outputfd);
    
    int *V;                     // V - Visited List
    int V_size = 0;             // V_size - Size of V
    int start_index = 0;        // Start point
    int query_index = 1;       // Query point < query_vectors
    int *l;                     // L~ - k nearest neighbors list


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

    GreedySearch(base_graph, &V, &V_size, l, L, query_vectors[query_index], k, start_index);

    printf("Query coordinates : \n");
    for (int i = 0; i < base_num_dimensions; i++) {
        printf("%.1f ", query_vectors[query_index][i]);
    }
    printf("\n");

    fprintf(outputfd, "Visited %d nodes: \n", V_size);
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
