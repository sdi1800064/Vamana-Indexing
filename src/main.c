#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../headers/fvecs.h"
#include "../headers/graph.h"

int main(int argc, char *argv[]) {

    // Get the current time
    time_t t = time(NULL);

    // Seed the random number generator with the current time
    srand((unsigned int)t);
    
    // Variables for Base file
    float** base_vectors;
    int base_num_vectors;
    int base_num_dimensions;

    // Variables for Query file
    float** query_vectors;
    int query_num_vectors;
    int query_num_dimensions;

    // Variables for GroundTruth
    int** groundtruth_vectors;
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
            L = atoi(argv[i + 1]);
            i++;
        }
    }

    // Check if all required arguments are provided
    if (!base_file_name || !query_file_name || !groundtruth_file_name || k == -1 || R == -1 || a == -1.0) {
        fprintf(stderr, "Error: Missing required arguments.\n");
        fprintf(stderr, "Usage: %s -b base_file_name -q query_file_name -g groundtruth_file_name -k (int k) -R (int R) -a (float a) -L (int L)\n", argv[0]);
        return 1;
    }

    if(L < k)  {
        fprintf(stderr, "Error: L must be greater than k.\n");
    }


    // Create an output file
    FILE *file_check = fopen("output.txt", "r");
    if (file_check != NULL) {
        // File exists, so close it and delete it
        printf("Existing output.txt file found. Deleting...\n");
        fclose(file_check);
        if (remove("output.txt") == 0) {
            printf("File deleted successfully.\n");
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
    read_ivecs(groundtruth_file_name, &groundtruth_vectors, &groundtruth_num_vectors, &groundtruth_num_dimensions);
    fprintf(outputfd, "Groundtruth-Vector dimensionality: %d\n", groundtruth_num_dimensions);
    fprintf(outputfd, "Number of Groundtruth-Vectorts: %d\n", groundtruth_num_vectors);
    // fprintIntVectors(groundtruth_vectors, groundtruth_num_vectors, groundtruth_num_dimensions, outputfd);


    // ============== CREATING RANDOM GRAPH ================= //
    Graph *base_graph = create_random_graph(base_vectors, 5, R, 20);


    // ------ UNCOMMENT THIS TO PRINT THE RANDOM GRAPH ------ //
    fprintf(outputfd, "Printing random graph in output file...\n");
    fprint_graph(base_graph, outputfd);

    printf("printed graph to the output file\n");

    // ============== VAMANA INDEXING ================= //
    vamana_indexing(base_graph, L, a, R, outputfd);

    float recall = 0.0;
    int *K_CLOSEST[query_num_vectors];
    int *V = NULL;
    int V_size = 0;

    printf("= Calculating recall..\n");

    // How many of the k-returned are in the groundtruth
    for(int i = 0; i < query_num_vectors; i++) {
        int random_index = rand() % base_num_vectors;
        K_CLOSEST[i] = (int*)malloc(L * sizeof(int));
        K_CLOSEST[i] = greedy_search(base_graph, random_index, query_vectors[i], L, &V_size, &V);
    }    
    free(V);

    for(int i = 0; i < query_num_vectors; i++) {
        for(int j = 0; j < k; j++) {
            int predicted_index = K_CLOSEST[i][j];
            for(int l = 0; l < groundtruth_num_dimensions; l++) {
                int groundtruth_index = groundtruth_vectors[i][l];
                if(predicted_index == groundtruth_index) {
                    recall += 1.0;
                    break;
                }
            }
        }
    }

    printf("Recall: %f\n", (recall / query_num_vectors));

    
    printf("Exiting the program..\n");

    // Close the file
    fclose(outputfd);

    return 0;
}
