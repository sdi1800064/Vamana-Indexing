#include <stdio.h>
#include <stdlib.h>

#include "../headers/fvecs.h"
#include "../headers/graph.h"


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
        printf("Filter %d: number of vectors with that filter: %d\n", filters.filters[0][i], filters.filters[1][i] );
    }

    // cprint_dataset(dataset, 4);
    // cprint_query_dataset(queryDataSet, 4);

    // print_dataset(dataset);

    free_dataset(dataset);
    return EXIT_SUCCESS;
}
