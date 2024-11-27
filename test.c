#include <stdio.h>
#include <stdlib.h>

#include "../headers/dataset.h"
#include "../headers/graph.h"


int main() {
    const char *filename = "testSets/dummy-data.bin";
    const char *queryfile = "testSets/dummy-queries.bin";
    uint32_t total_vectors;
    
    DatasetInfo *dataset = read_dataset(filename);
    // QueryInfo *queryDataSet = read_query_dataset(queryfile, &total_query_vectors);

    if (!dataset) {
        fprintf(stderr, "Failed to read the dataset.\n");
        return EXIT_FAILURE;
    }

    // for( int i = 0; i < dataset->filterInfo.num_filters; i++){
    //     printf("Filter %d: number of vectors: %d\n", dataset->filterInfo.filtersPoints[i].filter_index, dataset->filterInfo.filtersPoints[i].count );
    //     for (int j = 0; j < dataset->filterInfo.filtersPoints[i].count; j++) {
    //         printf("%d ", dataset->filterInfo.filtersPoints[i].point_indexes[j]);
    //     }
    //     printf("\n");
    // }

    sort_filter_array(&dataset->filterInfo);
    

    // cprint_dataset(dataset, 4);
    // cprint_query_dataset(queryDataSet, 4);

    // print_dataset(dataset);

    free_dataset(dataset);
    return EXIT_SUCCESS;
}
