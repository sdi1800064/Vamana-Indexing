#include <stdio.h>
#include <stdlib.h>

#include "../headers/fvecs.h"
#include "../headers/graph.h"

// Function to sort an array using Bubble Sort
void sort_filter_array(int *array[2], int size) {
    if (size <= 1) return; // No need to sort arrays of size 0 or 1

    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (array[0][j] > array[0][j + 1]) {
                // Swap elements
                int temp0 = array[0][j];
                array[0][j] = array[0][j + 1];
                array[0][j + 1] = temp0;

                int temp1 = array[1][j];
                array[1][j] = array[1][j + 1];
                array[1][j + 1] = temp1;
            }
        }
    }
}

int main() {
    const char *filename = "testSets/dummy-data.bin";
    const char *queryfile = "testSets/dummy-queries.bin";
    uint32_t total_vectors;
    
    uint32_t total_query_vectors;
    filterInfo filters;
    filters.filters[0] = (int*)malloc(sizeof(int));
    filters.filters[1] = (int*)malloc(sizeof(int));
    filters.filters_size = 0;
    DatasetInfo *dataset = read_dataset(filename, &total_vectors, &filters);
    QueryInfo *queryDataSet = read_query_dataset(queryfile, &total_query_vectors);

    if (!dataset) {
        fprintf(stderr, "Failed to read the dataset.\n");
        return EXIT_FAILURE;
    }

    sort_filter_array(filters.filters, filters.filters_size);

    printf("Number of filters: %u\n", filters.filters_size);
    for (int i = 0; i < filters.filters_size; i++) {
        printf("Filter %d: num %d\n", filters.filters[0][i], filters.filters[1][i] );
    }

    cprint_dataset(dataset, 2);
    cprint_query_dataset(queryDataSet, 4);


    // Example: Print the first vector's details
    // print_dataset(dataset);

    free_dataset(dataset);
    return EXIT_SUCCESS;
}
