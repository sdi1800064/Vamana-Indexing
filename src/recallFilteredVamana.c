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

    // Code to create the graph

    char *query_file_name = NULL;
    char *groundtruth_file_name = NULL;
    int k = -1;
    int L = -1;
    int R = -1;
    float a = -1.0;

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
        }
        else if (strcmp(argv[i], "-R") == 0 && i + 1 < argc) {
            R = atoi(argv[i + 1]);
            i++;
        }
        else if (strcmp(argv[i], "-a") == 0 && i + 1 < argc) {
            a = atof(argv[i + 1]);
            i++;
        }
        else if (strcmp(argv[i], "-L") == 0 && i + 1 < argc) {
            L = atoi(argv[i + 1]);
            i++;
        }
    }

    // Check if all required arguments are provided
    if (!base_file_name || !graph_file_name || !query_file_name || !groundtruth_file_name || k == -1 || L == -1 || a == -1.0 || R == -1) {
        fprintf(stderr, "Error: Missing required arguments.\n");
        fprintf(stderr, "Usage: %s -b base_file_name --graph graph_file_name -q query_file_name -g groundtruth_file_name -k (int k) -L (int L)\n", argv[0]);
        return 1;
    }
    
    
    if(L < k){
        fprintf(stderr, "Error: L must be greater than k.\n");
        fprintf(stderr, "Usage: %s -b base_file_name --graph graph_file_name -q query_file_name -g groundtruth_file_name -k (int k) -L (int L)\n", argv[0]);
        return 1;
    }

    printf("k = %d || L = %d || a = %f || R = %d\n", k, L, a, R);

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

    //Check if the file exists,if it does not create it and write the data else read the data
    char new_groundtruth_file_name[100];
    sprintf(new_groundtruth_file_name, "%s%s%d%s", graph_file_name, "_R", R, ".bin");


    Graph filteredVamanaGraph;
    // Check if the file exists
    if (access(graph_file_name, F_OK) != -1) {
        // File exists, read the data
        filteredVamanaGraph = *readGraph(graph_file_name);
    } else {
        time_t start_vamana = time(NULL);
        filteredVamanaGraph = filtered_vamana_indexing(dataSet, L,a,R,&(dataSet->filterInfo));
        time_t end_vamana = time(NULL);
        double time_vamana = difftime(end_vamana, start_vamana);
        printf("Time to create Filtered graph: %.3f seconds\n", time_vamana);
        writeVamanaGraph(&filteredVamanaGraph, graph_file_name);
    }



    // ============ RECALL AFTER VAMANA INDEXING ============= //

    printf("Calculating recall: ");
    fflush(stdout);

    float recall = 0.0;

    int total_count = 0;
    int filtered_total_count = 0;
    int filtered_prediction_count = 0;
    int count = 0;
    int filtered_total_relevant = 0;
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

            int* V = NULL;
            int *lamda = NULL;
            int lamda_size = 0;
            int V_size = 0;
            int startIndex = filteredVamanaGraph.medoid;
            filtered_greedy_search(&filteredVamanaGraph, querySet->queries[i].query_vector, &startIndex, 1, &V,&V_size,&lamda,&lamda_size,k, querySet->queries[i].v);
//             MUST INCLUDE THE DATASET
            for( int a = 0; a < k; a++) {
                    add_to_dynamic_array(&arrayOfPredictedIndexes, &sumOfAllClosest,lamda[a]);
            }
//            sort_array_based_on_dataset(dataSet, arrayOfPredictedIndexes, sumOfAllClosest, querySet->queries[i].query_vector);
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
            int sumOfAllClosest = 0;
            int* arrayOfPredictedIndexes = (int*)malloc(sizeof(int));

            int* V = NULL;
            int *lamda = NULL;
            int lamda_size = 0;
            int V_size = 0;
            int filteredmedoid;
            for(int a =0; a<=filteredVamanaGraph.filteredMedoids.size; a++){
                if(filteredVamanaGraph.filteredMedoids.metoids[a].category == querySet->queries[i].v){
                    filteredmedoid = filteredVamanaGraph.filteredMedoids.metoids[a].index;
                    break;
                }
            }
            filtered_greedy_search(&filteredVamanaGraph, querySet->queries[i].query_vector, &filteredmedoid, 1, &V,&V_size,&lamda,&lamda_size,k, querySet->queries[i].v);

            for( int a = 0; a < lamda_size; a++) {
                add_to_dynamic_array(&arrayOfPredictedIndexes, &sumOfAllClosest,lamda[a]);
            }
//            sort_array_based_on_dataset(dataSet, arrayOfPredictedIndexes, sumOfAllClosest, querySet->queries[i].query_vector);


            int k_a = lamda_size;
            if(lamda_size < k) {
                k_a = lamda_size;
            }

            for(int j = 0; j < k_a; j++) {
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
            // printf("Query %d: Found %d / %d\n", i, count, k_minimum);
            filtered_total_relevant += ground_truth_count[i]<k ? ground_truth_count[i] : k;
            prediction_count += k_a;
            recall += count / 100;
            total_count += count;
            filtered_total_count += count;
            filtered_prediction_count += k_a;
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
    printf("Found only for filtered %d / %d\n", filtered_total_count, filtered_prediction_count);
    float filtered_recall = 0.0f;
    if (filtered_total_relevant > 0) {
        filtered_recall = (float) filtered_total_count / (float) filtered_total_relevant;
    }
    printf("Recall: %.3f%%\n", recall*100);
    printf("Filtered Recall: %.3f%%\n", filtered_recall*100);

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
