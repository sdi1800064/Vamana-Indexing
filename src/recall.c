#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
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
    char *graph_file_name = NULL;
    char *query_file_name = NULL;
    char *groundtruth_file_name = NULL;
    int k = -1;
    int L = -1;

    // Iterate through command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            base_file_name = argv[i + 1];
            i++;
        } else if (strcmp(argv[i], "--graph") == 0 && i + 1 < argc) {
            graph_file_name = argv[i + 1];
            i++;
        } else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
            query_file_name = argv[i + 1];
            i++;
        } else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
            groundtruth_file_name = argv[i + 1];
            i++;
        } else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) {
            k = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "-L") == 0 && i + 1 < argc) {
            L = atoi(argv[i + 1]);
            i++;
        }
    }

    // Check if all required arguments are provided
    if (!base_file_name || !graph_file_name || !query_file_name || !groundtruth_file_name || k == -1 || L == -1) {
        fprintf(stderr, "Error: Missing required arguments.\n");
        fprintf(stderr, "Usage: %s -b base_file_name --graph graph_file_name -q query_file_name -g groundtruth_file_name -k (int k) -L (int L)\n", argv[0]);
        return 1;
    }

    if(L < k){
        fprintf(stderr, "Error: L must be greater than k.\n");
        fprintf(stderr, "Usage: %s -b base_file_name --graph graph_file_name -q query_file_name -g groundtruth_file_name -k (int k) -L (int L)\n", argv[0]);
        return 1;
    }

    printf("k = %d || L = %d\n", k, L);

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

    // ======================================= //

    // MUST INCLUDE THE DATASET. MANDATORY FOR SORTING
    DatasetInfo* dataSet;
    dataSet = read_dataset(base_file_name);

    // Initialise the base Dataset
    QueryInfo* querySet;
    querySet = read_query_dataset(query_file_name);

    // Initialise the groundtruth Dataset
    // groundTruthSet[num_queries][100]
    int** groundTruthSet;
    groundTruthSet = readGroundTruth(groundtruth_file_name, querySet->num_queries);
    if(groundTruthSet == NULL) {
        perror("Error reading ground truth");
        exit(EXIT_FAILURE);
    }

    // Initialise the graphs
    Graph* stitchedGraphs;
    int stitchedGraphs_size;
    stitchedGraphs = readGraphs(graph_file_name, &stitchedGraphs_size);


    // ============ RECALL AFTER VAMANA INDEXING ============= //

    printf("Calculating recall: ");
    fflush(stdout);

    float recall = 0.0;

    int total_count = 0;
    int count = 0;
    int prediction_count = 0;
    int ground_truth_count[querySet->num_queries];
    
    // Calculate the number of groundtruths for each query
    for(int i = 0; i < querySet->num_queries; i++) {
        ground_truth_count[i] = 0;
        if(querySet->queries[i].v != 144){
            for(int j = 0; j < 100; j++) {
                if(groundTruthSet[i][j] != -1) {
                    ground_truth_count[i]++;
                }       
                
            }
        }
            
    }

    int *V = NULL;
    int V_size = 0;
    int *K_CLOSEST = (int*)malloc(sizeof(int));
    int k_closest_size = 0;
    time_t start_recall_time = time(NULL);

    for( int i = 0; i < querySet->num_queries; i++) {
        int query_type = querySet->queries[i].query_type;
        if(query_type == 0) {
            // printf("Query %d of type %d\n", i, query_type);        
            int sumOfAllClosest = 0;
            int* arrayOfPredictedIndexes = (int*)malloc(sizeof(int));
            for(int j = 0; j < stitchedGraphs_size; j++) {
                int random_graph_index = rand() % stitchedGraphs[j].num_points;

                greedy_search(&stitchedGraphs[j], querySet->queries[i].query_vector, random_graph_index, &V, &V_size, &K_CLOSEST, &k_closest_size, L);
                
                for( int k = 0; k < k_closest_size; k++) {
                    add_to_dynamic_array(&arrayOfPredictedIndexes, &sumOfAllClosest, stitchedGraphs[j].points[K_CLOSEST[k]].index);
                }
                
                V = NULL;
                V_size = 0;
                k_closest_size = 0;

            }

            // MUST INCLUDE THE DATASET
            sort_array_based_on_dataset(dataSet, arrayOfPredictedIndexes, sumOfAllClosest, querySet->queries[i].query_vector);
            for(int j = 0; j < k; j++) {
                int predicted_index = arrayOfPredictedIndexes[j];
                for(int l = 0; l < 100; l++) {
                    if(groundTruthSet[i][l] == -1){
                        break;
                    }
                    if(predicted_index == groundTruthSet[i][l]) {
                        count++;
                        break;
                    }
                }
            }
            // printf("Query %d: Found %d / %d\n", i, count, k);
            prediction_count += k;
            recall += count / k;
            total_count += count;
            count = 0;
        
        } else if(query_type == 1 && querySet->queries[i].v != 144) {
            // printf("Query %d of type %d and filter %d\n", i, query_type, querySet->queries[i].v);

            int query_filter = querySet->queries[i].v;
            // printf("Filter_graph: %d with first index -> %d | filter %d\n", query_filter, stitchedGraphs[query_filter].points[0].index, stitchedGraphs[query_filter].points[0].category);
            int random_graph_index = rand() % stitchedGraphs[query_filter].num_points;

            greedy_search(&stitchedGraphs[query_filter], querySet->queries[i].query_vector, random_graph_index, &V, &V_size, &K_CLOSEST, &k_closest_size, L);
            
            int k_minimum = k;
            if(k_closest_size < k) {
                k_minimum = k_closest_size;
            }
            for(int k = 0; k < k_minimum; k++) {
                int predicted_index = stitchedGraphs[query_filter].points[K_CLOSEST[k]].index;
                for(int l = 0; l < 100; l++) {
                    if(groundTruthSet[i][l] == -1){
                        break;
                    }
                    if(predicted_index == groundTruthSet[i][l]) {
                        count++;
                        break;
                    }
                }                   
            }
            // printf("Query %d: Found %d / %d\n", i, count, k_minimum);
            prediction_count += k_minimum;
            recall += count / k_minimum;
            total_count += count;
            count = 0;
            V = NULL;
            V_size = 0;
            k_closest_size = 0;
        } 
    }
    printf("Done.\n");
    fflush(stdout);
    time_t end_recall_time = time(NULL);

    free(K_CLOSEST);
    free(V);
    recall = (float) total_count / prediction_count;
    printf("Found %d / %d\n", total_count, prediction_count);
    printf("Recall: %.3f%%\n", recall*100);

    // ============== FREE MEMORY ================= //
    
    
    free(dataSet->datapoints);
    free(dataSet);
    free(querySet->queries);
    free(querySet);
    free(groundTruthSet);
    printf("Time taken to compute recall: %ld seconds\n", end_recall_time - start_recall_time);
    printf("Exiting the program..\n");

    // Close the file
    // fclose(outputfd);

    return 0;
}
