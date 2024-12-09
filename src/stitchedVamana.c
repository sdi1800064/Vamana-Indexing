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
    int R = -1;
    float a = -1.0;
    int L = -1;

    // Iterate through command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            base_file_name = argv[i + 1];
            i++;  // Move past the flag and value
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
    if (!base_file_name || R == -1 || a == -1.0 || L == -1) {
        fprintf(stderr, "Error: Missing required arguments.\n");
        fprintf(stderr, "Usage: %s -b base_file_name -R (int R) -a (float a) -L (int L)\n", argv[0]);
        return 1;
    }

    printf("L = %d | a = %f | R = %d\n", L, a, R);


    DatasetInfo* dataSet;
    dataSet = read_dataset(base_file_name);
    int num_of_graphs = dataSet->filterInfo.num_filters;
    printf("Number of graphs that will be created %d\n", num_of_graphs);
    int Rsmall = R/2;
    Graph* graph = (Graph *)malloc(num_of_graphs * sizeof(Graph));

    time_t start_vamana = time(NULL);
    graph = stitched_vamana_indexing(dataSet, L, a, Rsmall);
    time_t end_vamana = time(NULL);

    double time_vamana = difftime(end_vamana, start_vamana);
    
    char* graph_file_name = "stitchedGraph";
    char new_groundtruth_file_name[100];
    sprintf(new_groundtruth_file_name, "%s%s%d%s", graph_file_name, "_R", R, ".bin");

    writeGraphs(graph, num_of_graphs, graph_file_name);
    for(int i = 0; i < num_of_graphs; i++){
        free_graph(graph[i]);
    }
    printf("Time taken to index all graphs %.2f seconds || Time taken for the whole program %.2f seconds\n", time_vamana, difftime(time(NULL), t));
    free(graph);
    free_dataset(dataSet);
    printf("Exiting program\n");
    return 0;
}