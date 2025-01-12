#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "../headers/dataset.h"
#include "../headers/graph.h"
#include "../headers/threadFunctions.h"
#include "../headers/recallFunctions.h"

#define NUM_THREADS 4

int main(int argc, char *argv[]) {

    // Get the current time
    time_t timet = time(NULL);
    char *filename = "results.txt";

    // Open/create the file
    FILE *resultFile = fopen(filename, "a");
    if (!resultFile) {
        perror("Error opening file");
        exit(1);
    }
    fseek(resultFile, 0, SEEK_END);
    fprintf(resultFile,"---------------------------------------------\n");

    fprintf(resultFile, "Starting recall\n");


    // Seed the random number generator with the current time
    srand((unsigned int)timet);
    
    printf("Starting...\n");
    // Initialize variables
    char *base_file_name = NULL;
    char *graph_file_name = NULL;
    char *query_file_name = NULL;
    char *groundtruth_file_name = NULL;
    int k = -1;
    int L = -1;
    float a = -1.0;
    int R = -1;
    int t = -1;

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
        } else if (strcmp(argv[i], "-a") == 0 && i + 1 < argc) {
            a = atof(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "-R") == 0 && i + 1 < argc) {
            R = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            t = atoi(argv[i + 1]);
            i++;
        }
    }

    // Check if all required arguments are provided
    if (!base_file_name || !graph_file_name || !query_file_name || !groundtruth_file_name || k == -1 || L == -1 || a == -1 || R == -1 || t == -1) {
        fprintf(stderr, "Error: Missing required arguments.\n");
        fprintf(stderr, "Usage: %s -b base_file_name --graph graph_file_name -q query_file_name -g groundtruth_file_name -k (int k) -L (int L) -R (int R) -t (int t)\n", argv[0]);
        return 1;
    }

    if(L < k){
        fprintf(stderr, "Error: L must be greater than k.\n");
        fprintf(stderr, "Usage: %s -b base_file_name --graph graph_file_name -q query_file_name -g groundtruth_file_name -k (int k) -L (int L) -R (int R) -t (int t)\n", argv[0]);
        return 1;
    }

    printf("baseFile = %s || queryFile = %s || groundTruthFile = %s\nk = %d || L = %d || a = %f || R = %d || t = %d\n", base_file_name, query_file_name, groundtruth_file_name, k, L, a, R, t);
    fprintf(resultFile, "baseFile = %s || queryFile = %s || groundTruthFile = %s\nk = %d || L = %d || a = %f || R = %d || t = %d\n", base_file_name, query_file_name, groundtruth_file_name, k, L, a, R, t);

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
    
    int stitchedGraphs_count = 0;
    char new_graph_file_name[100];
    snprintf(new_graph_file_name, sizeof(new_graph_file_name), "%s_R%d_L%d_a%.2f_#%d.bin", 
            "graphs/stitchedGraph", R, L, a, dataSet->num_vectors);
    int Rsmall = R/2;
    if(access(graph_file_name, F_OK) == -1){
        graph_file_name = new_graph_file_name;
    }

    if (access(new_graph_file_name, F_OK) == -1) {
        stitchedGraphs_count = dataSet->filterInfo.num_filters;
        stitchedGraphs = (Graph *)malloc(stitchedGraphs_count * sizeof(Graph));
        time_t start_vamana = time(NULL);
        stitchedGraphs = threadStitchedVamanaIndexing(dataSet, L, a, Rsmall, NUM_THREADS, resultFile);
        time_t end_vamana = time(NULL);
        writeGraphs(stitchedGraphs, stitchedGraphs_count, new_graph_file_name);
        double time_vamana = difftime(end_vamana, start_vamana);
        printf("Time to create graphs: %.3f seconds\n", time_vamana);
    }
    stitchedGraphs = readGraphs(new_graph_file_name, &stitchedGraphs_count);
    
    // ============ RECALL AFTER VAMANA INDEXING ============= //

    calculateRecallStitched(dataSet, querySet, groundTruthSet, stitchedGraphs, stitchedGraphs_count, L, k, NUM_THREADS, resultFile);

    // ============== FREE MEMORY ================= //
    
    free(dataSet->datapoints);
    free(dataSet);
    free(querySet->queries);
    free(querySet);
    free(groundTruthSet);
    printf("Exiting the program..\n");
    fprintf(resultFile,"---------------------------------------------\n");
    fclose(resultFile);

    // Close the file
    // fclose(outputfd);

    return 0;
}
