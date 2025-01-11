#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "../headers/dataset.h"
#include "../headers/graph.h"

int main(int argc, char *argv[]) {

   time_t t = time(NULL);

    // Seed the random number generator with the current time
    srand((unsigned int)t);

    printf("Starting...\n");
    // Initialize variables
    char *base_file_name = NULL;
    char *graph_file_name = NULL;

    // Code to create the graph

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
    if (!base_file_name || !graph_file_name || L == -1 || a == -1.0 || R == -1) {
        fprintf(stderr, "Error: Missing required arguments.\n");
        fprintf(stderr, "Usage: %s -b base_file_name --graph graph_file_name -q query_file_name -g groundtruth_file_name -k (int k) -L (int L)\n", argv[0]);
        return 1;
    }

    printf("L = %d || a = %f || R = %d\n",  L, a, R);

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

    // Initialise the base Datase


    Graph filteredVamanaGraph;
    // Check if the file exists

        time_t start_vamana = time(NULL);
        filteredVamanaGraph = filtered_vamana_indexing(dataSet, L,a,R,&(dataSet->filterInfo));
        time_t end_vamana = time(NULL);
        double time_vamana = difftime(end_vamana, start_vamana);
        printf("Time to create Filtered graph: %.3f seconds\n", time_vamana);
        writeVamanaGraph(&filteredVamanaGraph, graph_file_name);

        return 0;
}