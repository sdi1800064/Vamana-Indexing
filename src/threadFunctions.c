#include "graph.h"
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "dataset.h"

#include "timeFunctions.h"

// Function to compare two integers for sorting
int compare_desc(const void* a, const void* b) {
    return (*(int*)b - *(int*)a);
}

// Function to execute by each thread
void* thread_function_vamana_indexing(void* args) {

    ThreadArgs* thread_args = (ThreadArgs*) args;
    DatasetInfo* dataset = thread_args->dataset;
    int* filterList = thread_args->filterList;
    int numOfFilters = thread_args->numOfFilters;
    int L_small = thread_args->L_small;
    float a = thread_args->a;
    int R_small = thread_args->R_small;

    double start_time = get_current_time();
    double total_time = 0.0;
    double max_time = 0.0;
    double min_time = DBL_MAX;
    double avg_time = 0.0;

    // printf("Thread %d | numofFilters %d:\n", thread_args->id, numOfFilters);
    // for (int i = 0; i < numOfFilters; i++) {
    //     printf("%d ", filterList[i]);
    // }
    // printf("\n");

    // Process the assigned filters
    for (int i = 0; i < numOfFilters; i++) {
        // Perform the tasks for the specific filter
        int filterIndex = filterList[i];

        // Validate filterIndex
        if (filterIndex < 0 || filterIndex >= dataset->filterInfo.num_filters) {
            printf("Thread %d: Invalid filterIndex %d\n", thread_args->id, filterIndex);
            continue; // Skip this iteration if the index is invalid
        }

        double start_filter_time = get_current_time();
        // printf("\nthread %d | i = %d | filterIndex = %d | filterCount = %d\n", thread_args->id, i, filterIndex, dataset->filterInfo.filtersPoints[filterIndex].count);
        // for(int j = 0; j < dataset->filterInfo.filtersPoints[filterIndex].count; j++){
        //     printf("%d  ", dataset->filterInfo.filtersPoints[filterIndex].point_indexes[j]);
        // }
        // printf("\n");
        // printf("Thread %d traversing...\n", thread_args->id);
        // fflush(stdout);
        // printf("Thread %d Indexing filter graph %d with %d points\n",thread_args->id, filterIndex, dataset->filterInfo.filtersPoints[filterIndex].count);
        DatasetInfo filter_dataset;
        filter_dataset.num_vectors = dataset->filterInfo.filtersPoints[filterIndex].count;
        filter_dataset.datapoints = (DataPoint *)calloc(filter_dataset.num_vectors, sizeof(DataPoint));
        if(filter_dataset.datapoints == NULL){
            printf("Memory allocation failed filter_dataset.datapoint!\n");
            exit(1);
        }

        // printf(" Thread %d | number of point in filter %d\n", thread_args->id, filter_dataset.num_vectors);
        // if(filter_dataset.num_vectors ==15 && thread_args->id == 1){
        //     cprint_dataset(dataset, filterIndex);
        // }

        for(int j = 0; j < filter_dataset.num_vectors; j++){
            int filter_point_index = dataset->filterInfo.filtersPoints[filterIndex].point_indexes[j];

            // Validate filter_point_index
            if (filter_point_index < 0 || filter_point_index >= dataset->num_vectors) {
                printf("Thread %d: Invalid point index %d\n", thread_args->id, filter_point_index);
                continue; // Skip this iteration if the index is invalid
            }
            filter_dataset.datapoints[j].category = dataset->datapoints[filter_point_index].category;
            filter_dataset.datapoints[j].point_index = dataset->datapoints[filter_point_index].point_index;
            filter_dataset.datapoints[j].timestamp = dataset->datapoints[filter_point_index].timestamp;
            for(int k = 0; k < 100; k++){
                filter_dataset.datapoints[j].vectors[k] = dataset->datapoints[filter_point_index].vectors[k];
            }
        }
        // For each filter, run the vamana algorithm
        Graph temp_graph = vamana_indexing(filter_dataset, L_small, a, R_small);
        thread_args->filter_graph[i] = temp_graph;
        free(filter_dataset.datapoints); // Free the allocated memory after use

        double end_filter_time = get_current_time();
        double filter_time = end_filter_time - start_filter_time;
        total_time += filter_time;
        if(filter_time > max_time){
            max_time = filter_time;
        }
        if(filter_time < min_time){
            min_time = filter_time;
        }
    }
    double elapsedTime = get_elapsed_time(start_time);
    avg_time = total_time / numOfFilters;
    thread_args->totalTimeTaken = elapsedTime;
    thread_args->min_time = min_time;
    thread_args->max_time = max_time;
    thread_args->avg_time = avg_time;
    // printf("Thread %d finished\n", thread_args->id);
    return NULL;
}


// Modify the threadStitchedVamanaIndexing function to divide filters among threads
Graph* threadStitchedVamanaIndexing(DatasetInfo* dataset, int L_small, float a, int R_small, int numOfThreads) {
    
    if(dataset == NULL){
        printf("Dataset is NULL\n");
        return NULL;
    }
    if(dataset->filterInfo.filtersPoints == NULL || dataset->filterInfo.num_filters == 0){
        printf("No filters present\n");
        return NULL;
    }

    DatasetInfo *filter_dataset = (DatasetInfo *)malloc(dataset->filterInfo.num_filters * sizeof(DatasetInfo));
    if(filter_dataset == NULL){
        printf("Memory allocation failed filter_set!\n");
        exit(1);
    }

    // Create a thread pool with a fixed number of threads
    pthread_t* threads = (pthread_t*) malloc(numOfThreads * sizeof(pthread_t));
    ThreadArgs* thread_args = (ThreadArgs*) malloc(numOfThreads * sizeof(ThreadArgs));

    // Number of Points each thread will have to index
    int threadLoad[numOfThreads];
    memset(threadLoad, 0, sizeof(threadLoad));

    // List of Filters 
    int threadFilterList[numOfThreads][dataset->filterInfo.num_filters];
    memset(threadFilterList, -1, sizeof(threadFilterList));

    // Number of Filters
    int threadFilterCount[numOfThreads];
    memset(threadFilterCount, 0, sizeof(threadFilterCount));

    // List of order
    int threadOrder[dataset->filterInfo.num_filters];
    memset(threadOrder, -1, sizeof(threadOrder));

    if(numOfThreads == 1){
        threadFilterCount[0] = dataset->filterInfo.num_filters;
        for(int i = 0; i < dataset->filterInfo.num_filters; i++){
            threadFilterList[0][i] = i;
        }
    }else{
        // Distibute the filters equaly
        for(int i = 0; i < dataset->filterInfo.num_filters; i++){
            // The thread we are loading
            int min_thread = 0;
            for(int j = 0; j < numOfThreads; j++){
                // Find the thread with the least load
                if(threadLoad[j] < threadLoad[min_thread]){
                    min_thread = j;
                }
            }

            // Add the thread to the order
            threadOrder[i] = min_thread;
            
            // Add the filter
            // printf("Adding filter %d to thread %d\n", i, min_thread);
            threadFilterList[min_thread][threadFilterCount[min_thread]] = i;
            ++threadFilterCount[min_thread];

            // Up the load
            if(dataset->filterInfo.filtersPoints[i].count > 500){
                int loadCount = 2 * dataset->filterInfo.filtersPoints[i].count;
                threadLoad[min_thread] += loadCount;    
            }else{
                threadLoad[min_thread] += dataset->filterInfo.filtersPoints[i].count;
            }
        }
    }

    
    int totalFilters = 0;
    // Assign filters to each thread
    for (int i = 0; i < numOfThreads; i++) {
        // printf("Thread %d had %d load\n", i, threadLoad[i]);
        thread_args[i].id = i;
        thread_args[i].dataset = (DatasetInfo*)malloc(sizeof(DatasetInfo));
        thread_args[i].dataset = dataset;
        // Add the list of Filters to the arguments of the thread "i"
        thread_args[i].filterList = (int*)malloc(threadFilterCount[i] * sizeof(int));

        int j = 0;
        while(threadFilterList[i][j] != -1 && j < dataset->filterInfo.num_filters){
            // printf("Loading thread_args[%d].filterList[%d] <- %d\n", i, j, threadFilterList[i][j]);
            thread_args[i].filterList[j] = threadFilterList[i][j];
            ++j;
        }

        thread_args[i].numOfFilters = threadFilterCount[i];
        totalFilters += thread_args[i].numOfFilters;
        thread_args[i].L_small = L_small;
        thread_args[i].a = a;
        thread_args[i].R_small = R_small;
        thread_args[i].filter_graph = (Graph*) malloc(thread_args[i].numOfFilters * sizeof(Graph));
        thread_args[i].totalTimeTaken = 0.0;
        thread_args[i].avg_time = 0.0;
        thread_args[i].min_time = 0.0;
        thread_args[i].max_time = 0.0;

        pthread_create(&threads[i], NULL, thread_function_vamana_indexing, &thread_args[i]);
    }
    

    //  threadFilterCount  threadFilterList  threadOrder
    // filter_graph[num_filters of the dataset]
    Graph* filter_graph = (Graph*) malloc(dataset->filterInfo.num_filters * sizeof(Graph));

    for (int i = 0; i < numOfThreads; i++) {
        // thread_filter_graph[threadFilterCount[i]]
        pthread_join(threads[i], NULL);
        Graph* thread_filter_graphs = thread_args[i].filter_graph;
        int thread_numOfGraphs = threadFilterCount[i];
        // printf("Number of graphs assigned to thread %d : %d\n", i, thread_numOfGraphs);

        // Print the times of each thread
        printf("Thread %d / %d -> #Filters = %d | min_time=%.2f | max_time=%.2f | avg_time=%.2f | total_time_taken=%.2f |\n", i+1, numOfThreads, thread_numOfGraphs, thread_args[i].min_time, thread_args[i].max_time, thread_args[i].avg_time, thread_args[i].totalTimeTaken);

        // For every filter in the list of filters of the thread
        for(int j = 0; j < thread_numOfGraphs; j++){
            // printf("Inserting Graph %d of category %d | Thread %d | to Position %d \n", j, thread_filter_graphs[j].points[0].category, i, threadFilterList[i][j]);
            filter_graph[threadFilterList[i][j]] = thread_filter_graphs[j];
            // fprint_graph(&filter_graph[threadFilterList[i][j]], stdout);
        }

    }
    // printf("Joined!\n");
    free(threads);
    free(thread_args);

    return filter_graph;
}
