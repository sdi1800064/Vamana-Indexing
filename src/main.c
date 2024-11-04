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
        exit(1);
    }

    printf("k = %d | L = %d | a = %f | R = %d\n", k, L, a, R);


    // ===================== CREATE OUTPUT FILE IF NEEDED ================== //
    // Create an output file
    
    // FILE *file_check = fopen("output.txt", "r");
    // if (file_check != NULL) {
    //     // File exists, so close it and delete it
    //     printf("Existing output.txt file found. Deleting...\n");
    //     fclose(file_check);
    //     if (remove("output.txt") == 0) {
    //         printf("File deleted successfully.\n");
    //     } else {
    //         printf("Error: Could not delete existing output.txt file!\n");
    //         return 1;
    //     }
    // }

    // // Create a new file (overwrite or create anew)
    // FILE *outputfd = fopen("output.txt", "w");
    // if (outputfd == NULL) {
    //     printf("Error: Could not create output.txt file!\n");
    //     return 1;
    // }

    // Read the base_vectors from the file
    read_fvecs(base_file_name, &base_vectors, &base_num_vectors, &base_num_dimensions);
    printf("Base-Vector dimensionality: %d\n", base_num_dimensions);
    printf("Number of Base-Vectorts: %d\n", base_num_vectors);

    // Read the query_vectors from the file
    read_fvecs(query_file_name, &query_vectors, &query_num_vectors, &query_num_dimensions);
    printf("Query-Vector dimensionality: %d\n", query_num_dimensions);
    printf("Number of Query-Vectorts: %d\n", query_num_vectors);

    //============== UNCOMMENT THIS TO PRINT THE QUERY POINTS =================//
    // printf("Printing query vector...\n");
    // fprintFloatVectors(query_vectors, query_num_vectors, query_num_dimensions, outputfd);

    // Read the groundtruth_vectors from the file
    read_ivecs(groundtruth_file_name, &groundtruth_vectors, &groundtruth_num_vectors, &groundtruth_num_dimensions);
    printf("Number of Groundtruth indexes: %d\n", groundtruth_num_dimensions);
    printf("Groundtruths for queries: %d\n", groundtruth_num_vectors);
    // fprintIntVectors(groundtruth_vectors, groundtruth_num_vectors, groundtruth_num_dimensions, outputfd);


    // ============== CREATING RANDOM GRAPH ================= //
    Graph *base_graph = create_random_graph(base_vectors, base_num_dimensions, R, base_num_vectors);


    // ------ UNCOMMENT THIS TO PRINT THE RANDOM GRAPH ------ //
    // fprintf(outputfd, "Printing random graph in output file...\n");
    // fprint_graph(base_graph, outputfd);

    // printf("printed graph to the output file\n");


    // ============== RECALL BEFORE VAMANA INDEXING ================= //
    
    int *K_CLOSEST[query_num_vectors];
    int *V = NULL;
    int V_size = 0;
    int lamda_size = 0;
    float recall = 0.0;


    int total_count = 0;
    int count = 0;
    int prediction_count = 0;

    // =============== UUNCOMMENT TO PRINT THE RECALL BEFORE VAMANA INDEXING ================== //
    // for(int i = 0; i < query_num_vectors; i++) {
    //     K_CLOSEST[i] = (int*)malloc(sizeof(int));
    //     int random_graph_index = rand() % base_num_vectors;
    //     greedy_search(base_graph, query_vectors[i], random_graph_index, &V, &V_size, &K_CLOSEST[i], &lamda_size, L, k);
    //     for(int j = 0; j < k; j++) {
    //         int predicted_index = K_CLOSEST[i][j];
    //         for(int l = 0; l < groundtruth_num_vectors; l++) {
    //             if(predicted_index == groundtruth_vectors[i][l]) {
    //                 count++;
    //                 break;
    //             }
    //         }
    //     }
    //     // printf("Query %d: Found %d / %d\n", i, count, k);
    //     total_count += count;
    //     recall += count / k;
    //     count = 0;
    // }
    // int prediction_count = k * query_num_vectors;
    // printf("Found %d / %d before\n", total_count, prediction_count);
    // recall = (float)total_count / prediction_count;
    // printf("Recall: %.3f%%\n", recall*100);

    // ============== VAMANA INDEXING ================= //

    vamana_indexing(base_graph, L, a, R);

    // ============ RECALL AFTER VAMANA INDEXING ============= //

    printf("= Calculating recall..\n");

    V = NULL;
    V_size = 0;
    lamda_size = 0;
    recall = 0.0;

    total_count = 0;
    count = 0;
    for(int i = 0; i < query_num_vectors; i++) {
        K_CLOSEST[i] = (int*)malloc(sizeof(int));
        int random_graph_index = rand() % base_num_vectors;
        greedy_search(base_graph, query_vectors[i], random_graph_index, &V, &V_size, &K_CLOSEST[i], &lamda_size, L);
        for(int j = 0; j < k; j++) {
            int predicted_index = K_CLOSEST[i][j];
            for(int l = 0; l < groundtruth_num_vectors; l++) {
                if(predicted_index == groundtruth_vectors[i][l]) {
                    count++;
                    break;
                }
            }
        }
        // printf("Query %d: Found %d / %d\n", i, count, k);
        recall += count / k;
        total_count += count;
        count = 0;
        free(K_CLOSEST[i]);
    }

    free(V);
    prediction_count = k * query_num_vectors;
    recall = (float) total_count / prediction_count;
    printf("Found %d / %d\n", total_count, prediction_count);
    printf("Recall: %.3f%%\n", recall*100);


    printf("Exiting the program..\n");

    // Close the file
    // fclose(outputfd);

    return 0;
}
