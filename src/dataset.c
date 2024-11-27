#include <stdio.h>
#include <stdlib.h>
#include "../headers/dataset.h"
#include "../headers/graph.h"

// Function to dynamically add a filter or update its count
void add_or_increment_filter(filterInfo *filtersInfo, int point_category, int point_index) {

    // The filters array is empty
    if( filtersInfo->num_filters == 0) {
        filtersInfo->filtersPoints[0].filter_index = point_category;
        filtersInfo->filtersPoints[0].count = 1;
        filtersInfo->filtersPoints[0].point_indexes = (int *)malloc(sizeof(int));
        filtersInfo->filtersPoints[0].point_indexes[0] = point_index;
        filtersInfo->num_filters++;
        return;
    }

    // Check if the filter already exists
    for ( int i = 0; i < filtersInfo->num_filters; i++) {
        // If the filter already exists, Add it to the point indexes
        if ( filtersInfo->filtersPoints[i].filter_index == point_category) {
            filtersInfo->filtersPoints[i].count++;
            filtersInfo->filtersPoints[i].point_indexes = (int *)realloc(filtersInfo->filtersPoints[i].point_indexes, (filtersInfo->filtersPoints[i].count) * sizeof(int));
            filtersInfo->filtersPoints[i].point_indexes[filtersInfo->filtersPoints[i].count - 1] = point_index;
            return;
        }
    }

    // If the filter doesn't exist, create it
    filtersInfo->filtersPoints = (filterPoint *)realloc(filtersInfo->filtersPoints, (filtersInfo->num_filters + 1) * sizeof(filterPoint));
    filtersInfo->filtersPoints[filtersInfo->num_filters].filter_index = point_category;
    filtersInfo->filtersPoints[filtersInfo->num_filters].count = 1;
    filtersInfo->filtersPoints[filtersInfo->num_filters].point_indexes = (int *)malloc(sizeof(int));
    filtersInfo->filtersPoints[filtersInfo->num_filters].point_indexes[0] = point_index;
    filtersInfo->num_filters++;

}


filterInfo initialise_filters() {

    filterInfo filtersInfo;
    filtersInfo.num_filters = 0;

    filtersInfo.filtersPoints = (filterPoint *)malloc(sizeof(filterPoint));
    if (!filtersInfo.filtersPoints) {
        perror("Memory allocation error for filters");
        exit(EXIT_FAILURE);
    }

    filtersInfo.filtersPoints[0].filter_index = -1;
    filtersInfo.filtersPoints[0].count = 0;
    filtersInfo.filtersPoints[0].point_indexes = (int *)malloc(sizeof(int));

    return filtersInfo;
}

// Function to read the dataset
DatasetInfo* read_dataset(const char *filename) {
    printf("Reading dataset from %s.\n", filename);
    
    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }

    // Read the number of vectors
    uint32_t num_vectors;
    if (fread(&num_vectors, sizeof(uint32_t), 1, file) != 1) {
        perror("Error reading number of vectors");
        fclose(file);
        return NULL;
    }

    // Allocate memory for the DatasetInfo structure
    DatasetInfo *dataset = (DatasetInfo *)malloc(sizeof(DatasetInfo));
    if (!dataset) {
        perror("Memory allocation error");
        fclose(file);
        return NULL;
    }
    dataset->num_vectors = (int)num_vectors;
    dataset->filterInfo = initialise_filters();

    // Allocate memory for the DataPoint array
    dataset->datapoints = (DataPoint *)malloc(num_vectors * sizeof(DataPoint));
    if (!dataset->datapoints) {
        perror("Memory allocation error");
        free(dataset);
        fclose(file);
        return NULL;
    }

    // Define the vector size
    const int vector_size = (2 + 100) * sizeof(float); // 2 attributes + 100 dimensions

    // Read each vector from the file
    int count = 0;
    for (int i = 0; i < dataset->num_vectors; i++) {
        float buffer[102];

        // Read the entire vector into the buffer
        if (fread(buffer, vector_size, 1, file) != 1) {
            perror("Error reading vector data");
            free(dataset->datapoints);
            free(dataset);
            fclose(file);
            return NULL;
        }

        // Populate the DataPoint structure
        dataset->datapoints[i].category = (int)buffer[0];
        dataset->datapoints[i].timestamp = buffer[1];
        for (int j = 0; j < 100; j++) {
            dataset->datapoints[i].vectors[j] = buffer[j + 2];
        }

        // Update the filters with the category
        // printf("Updating index %d to filter index %d.\n", i, dataset->datapoints[i].category);
        count++;
        add_or_increment_filter(&dataset->filterInfo, dataset->datapoints[i].category, i);
        
    }
    printf("counted inside the dataset %d\n", count);
    fclose(file);
    printf("Finished reading dataset.\n");
    return dataset;
}

// Function to free the allocated dataset
void free_dataset(DatasetInfo *dataset) {
    free(dataset->datapoints);
    free(dataset);
}

void print_dataset(DatasetInfo *dataset) {
    for (int i = 0; i < dataset->num_vectors; i++) {
        printf("Vector %d: Category = %d, Timestamp = %f, Vectors = ", i, dataset->datapoints[i].category, dataset->datapoints[i].timestamp);
        for (int j = 0; j < 100; j++) {
            printf("%f ", dataset->datapoints[i].vectors[j]);
        }
        printf("\n");
    }
}

void cprint_dataset(DatasetInfo *dataset, int target_category) {
    for (int i = 0; i < dataset->num_vectors; i++) {
        if(dataset->datapoints[i].category == target_category){
            printf("Vector %d: Category = %d\n", i, dataset->datapoints[i].category);
        }
    }
}


// Function to read the query set
QueryInfo* read_query_dataset(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }

    // Read the number of queries
    uint32_t num_queries;
    if (fread(&num_queries, sizeof(uint32_t), 1, file) != 1) {
        perror("Error reading number of queries");
        fclose(file);
        return NULL;
    }

    // Allocate memory for the QueryInfo structure
    QueryInfo *query_info = (QueryInfo *)malloc(sizeof(QueryInfo));
    if (!query_info) {
        perror("Memory allocation error");
        fclose(file);
        return NULL;
    }
    query_info->num_queries = (int)num_queries;

    // Allocate memory for the QueryPoint array
    query_info->queries = (QueryPoint *)malloc(num_queries * sizeof(QueryPoint));
    if (!query_info->queries) {
        perror("Memory allocation error");
        free(query_info);
        fclose(file);
        return NULL;
    }

    // Define the query size
    const size_t query_size = (4 + 100) * sizeof(float); // 4 attributes + 100 dimensions

    // Read each query from the file
    for (uint32_t i = 0; i < num_queries; i++) {
        float buffer[104];

        // Read the entire query into the buffer
        if (fread(buffer, query_size, 1, file) != 1) {
            perror("Error reading query data");
            fclose(file);
            free(query_info->queries);
            free(query_info);
            return NULL;
        }

        // Populate the QueryPoint structure
        query_info->queries[i].query_type = (int)buffer[0];
        query_info->queries[i].v = (int)buffer[1];
        query_info->queries[i].l = buffer[2];
        query_info->queries[i].r = buffer[3];
        for (size_t j = 0; j < 100; j++) {
            query_info->queries[i].query_vector[j] = buffer[j + 4];
        }
    }

    fclose(file);
    return query_info;
}

// Function to free the allocated query set
void free_query_dataset(QueryInfo *queries) {
    free(queries->queries);
    free(queries);
}


void print_query_dataset(QueryInfo *dataset) {
    for (int i = 0; i < dataset->num_queries; i++) {
        QueryPoint query = dataset->queries[i];
        printf("Vector %d: Query-Category = %d, Filter = %d, Timestamp = %f, Timestamp = %f, Vectors = ", i, query.query_type, query.v, query.l, query.r);
        for (int j = 0; j < 100; j++) {
            printf("%f ", query.query_vector[j]);
        }
        printf("\n");
    }
}

void cprint_query_dataset(QueryInfo *dataset, int category) {
    for (int i = 0; i < dataset->num_queries; i++) {
        QueryPoint query = dataset->queries[i];
        if(query.v == category){
            printf("Vector %d: Query-Category = %d, Filter = %d, Timestamp = %f, Timestamp = %f, Vectors = ", i, query.query_type, query.v, query.l, query.r);
            printf("\n");
        }
    }
}
