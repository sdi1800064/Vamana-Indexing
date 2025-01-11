#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "../headers/dataset.h"
#include "../headers/graph.h"
#include "../headers/threadFunctions.h"
#include "../headers/timeFunctions.h"

int main(int argc, char *argv[]) {

    // Get the current time
    time_t t = time(NULL);

    // Seed the random number generator with the current time
    srand((unsigned int)t);
        // Initialize variables
    char *base_file_name = NULL;
    int R = -1;
    float a = -1.0;
    int L = -1;
    int numOfThreads = 0;

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
        } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            numOfThreads = atoi(argv[i + 1]);
            i++;
        }
    }

    // Check if all required arguments are provided
    if (!base_file_name || R == -1 || a == -1.0 || L == -1 || numOfThreads <= 0) {
        fprintf(stderr, "Error: Missing required arguments.\n");
        fprintf(stderr, "Usage: %s -b (base_file_name) -R (int R) -a (float a) -L (int L) -t (int numOfThreads)\n", argv[0]);
        return 1;
    }

    printf("\nbase_file = %s | L = %d | a = %f | R = %d | threads = %d\n", base_file_name, L, a, R, numOfThreads);
    printf("---------------------------------------------\n");


    DatasetInfo* dataSet;
    dataSet = read_dataset(base_file_name);
    printf("---------------------------------------------\n");

    int num_of_graphs = dataSet->filterInfo.num_filters;
    int Rsmall = R/2;
    Graph* graph = (Graph *)malloc(num_of_graphs * sizeof(Graph));

    double startVamana = get_current_time();

    graph = threadStitchedVamanaIndexing(dataSet, L, a, Rsmall, numOfThreads);
    printf("---------------------------------------------\n");
    double time_vamana = get_elapsed_time(startVamana);

    char graph_file_name[100]; // Ensure this is large enough
    snprintf(graph_file_name, sizeof(graph_file_name), "%s_%s%d%s", "stitchedGraph", "R", R, "_");
    snprintf(graph_file_name, sizeof(graph_file_name), "%s#%d", dataSet->num_vectors, ".bin");

    writeGraphs(graph, num_of_graphs, graph_file_name);

    for(int i = 0; i < num_of_graphs; i++){
        free_graph(graph[i]);
    }
    printf("Time taken to index all graphs %.2f seconds\n", time_vamana);
    free(graph);
    free_dataset(dataSet);
    printf("Exiting program\n");
    return 0;
}