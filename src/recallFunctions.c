#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <float.h>

#include "../headers/dataset.h"
#include "../headers/graph.h"
#include "../headers/structs.h"
#include "../headers/threadFunctions.h"
#include "../headers/timeFunctions.h"

void calculateRecallStitched(DatasetInfo* dataSet, QueryInfo* querySet, int** groundTruthSet, Graph* stitchedGraphs, int stitchedGraphs_count, int L, int k, int numOfThreads, FILE* resultFile) {
    
    printf("Calculating recall: ");
    fprintf(resultFile,"Calculating recall: ");
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

    double start_recall_time = get_current_time();
    double total_query_time = 0.0;
    double max_query_time = 0.0;
    double min_query_time = DBL_MAX;

    for( int i = 0; i < querySet->num_queries; i++) {
        double start_query_time = get_current_time();
        int query_type = querySet->queries[i].query_type;
        if(query_type == 0) {
            int sumOfAllClosest = 0;
            int* arrayOfPredictedIndexes = (int*)malloc(sizeof(int));
            for(int j = 0; j < stitchedGraphs_count; j++) {
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
            int query_filter = querySet->queries[i].v;
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
        double end_query_time = get_current_time();
        double query_time = end_query_time - start_query_time;
        total_query_time += query_time;
        if (query_time > max_query_time) {
            max_query_time = query_time;
        }
        if (query_time < min_query_time) {
            min_query_time = query_time;
        }
    }
    double elapsed = get_elapsed_time(start_recall_time);
    printf("Done.\n");
    fflush(stdout);
    printf("Time taken: %.3f seconds\n", elapsed);

    free(K_CLOSEST);
    free(V);
    recall = (float) total_count / prediction_count;
    printf("Found %d / %d\n", total_count, prediction_count);
    printf("Recall: %.3f%%\n", recall*100);
    printf("Average query time: %.3f seconds\n", total_query_time / querySet->num_queries);
    printf("Max query time: %.3f seconds\n", max_query_time);
    printf("Min query time: %.3f seconds\n", min_query_time);

    fprintf(resultFile, "Done.\n");
    fprintf(resultFile, "Time taken: %.3f seconds\n", elapsed);
    fprintf(resultFile, "Found %d / %d\n", total_count, prediction_count);
    fprintf(resultFile, "Recall: %.3f%%\n", recall*100);
    fprintf(resultFile, "Average query time: %.3f seconds\n", total_query_time / querySet->num_queries);
    fprintf(resultFile, "Max query time: %.3f seconds\n", max_query_time);
    fprintf(resultFile, "Min query time: %.3f seconds\n", min_query_time);

}