#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../headers/dataset.h"
#include "../headers/graph.h"

int main(int argc, char *argv[]) {

    // Get the current time
    time_t t = time(NULL);

    // Seed the random number generator with the current time
    srand((unsigned int)t);
    
    printf("Starting...\n");
    // Initialize variables
    char *base_file_name = NULL;
    // char *query_file_name = NULL;
    // char *groundtruth_file_name = NULL;
    int k = -1;
    int R = -1;
    float a = -1.0;
    int L = -1;

    // Iterate through command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            base_file_name = argv[i + 1];
            i++;  // Move past the flag and value
        // } else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
        //     query_file_name = argv[i + 1];
        //     i++;
        // } else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
        //     groundtruth_file_name = argv[i + 1];
        //     i++;
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
    // if (!base_file_name || !query_file_name || !groundtruth_file_name || k == -1 || R == -1 || a == -1.0) {
    //     fprintf(stderr, "Error: Missing required arguments.\n");
    //     fprintf(stderr, "Usage: %s -b base_file_name -q query_file_name -g groundtruth_file_name -k (int k) -R (int R) -a (float a) -L (int L)\n", argv[0]);
    //     return 1;
    // }

    if (!base_file_name || k == -1 || R == -1 || a == -1.0) {
        fprintf(stderr, "Error: Missing required arguments.\n");
        fprintf(stderr, "Usage: %s -b base_file_name -k (int k) -R (int R) -a (float a) -L (int L)\n", argv[0]);
        return 1;
    }

    if(L < k)  {
        fprintf(stderr, "Error: L must be greater than k.\n");
        exit(1);
    }

    printf("k = %d | L = %d | a = %f | R = %d\n", k, L, a, R);

    // ===================== CREATE OUTPUT FILE IF NEEDED ================== //
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

    // ======================================= //


    DatasetInfo* dataSet;
    dataSet = read_dataset(base_file_name);

    // Initialise the base Dataset
    // QueryInfo* querySet;
    // querySet = read_query_dataset(query_file_name);


    // ============== VAMANA INDEXING ================= //
    int R_small = R/2;
    int graph_size = dataSet->filterInfo.num_filters;
    Graph* graph = (Graph*)malloc(graph_size * sizeof(Graph));
    if (graph == NULL) {
        printf("Error: Could not allocate memory for graph!\n");
        return 1;
    }
    // print_dataset(dataSet);
    graph = stitched_vamana_indexing(dataSet, L, a, R_small, R);
    printf("Printing Graphs..\n");
    for(int i = 0; i < graph_size; i++) {
        fprintf(outputfd, "Graph %d:\n", i);
        fprint_graph(&graph[i], outputfd);
        fprintf(outputfd, "====================================\n\n");        
    }

    // ============ RECALL AFTER VAMANA INDEXING ============= //

    // printf("Calculating recall..\n");


    // int *K_CLOSEST[query_num_vectors];
    // int *V = NULL;
    // int V_size = 0;
    // int lamda_size = 0;
    // float recall = 0.0;


    // int total_count = 0;
    // int count = 0;
    // int prediction_count = 0;
    // count = 0;
    // for(int i = 0; i < query_num_vectors; i++) {
    //     K_CLOSEST[i] = (int*)malloc(sizeof(int));
    //     int random_graph_index = rand() % base_num_vectors;
    //     greedy_search(base_graph, query_vectors[i], random_graph_index, &V, &V_size, &K_CLOSEST[i], &lamda_size, L);
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
    //     recall += count / k;
    //     total_count += count;
    //     count = 0;
    //     free(K_CLOSEST[i]);
    // }

    // free(V);
    // prediction_count = k * query_num_vectors;
    // recall = (float) total_count / prediction_count;
    // printf("Found %d / %d\n", total_count, prediction_count);
    // printf("Recall: %.3f%%\n", recall*100);

    for(int i = 0; i < dataSet->filterInfo.num_filters; i++) {
        free_graph(graph[i]);
    }   
    free(graph);

    printf("Exiting the program..\n");

    // Close the file
    // fclose(outputfd);

    return 0;
}
