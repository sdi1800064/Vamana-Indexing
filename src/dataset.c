#include <stdio.h>
#include <stdlib.h>
#include "../headers/dataset.h"
#include "../headers/graph.h"

// Function to dynamically add a filter or update its count
void add_or_increment_filter(filterInfo *filters, int category) {
    // Check if the filter already exists
    for (int i = 0; i < filters->filters_size; i++) {
        if (filters->filters[0][i] == category) {
            // Increment the count for this filter
            filters->filters[1][i]++;
            return;
        }
    }

    // If not found, add a new filter
    filters->filters_size++;
    filters->filters[0] = (int *)realloc(filters->filters[0], filters->filters_size * sizeof(int));
    filters->filters[1] = (int *)realloc(filters->filters[1], filters->filters_size * sizeof(int));

    if (!filters->filters[0] || !filters->filters[1]) {
        perror("Memory allocation error for filters");
        exit(EXIT_FAILURE);
    }

    // Add the new filter and initialize its count to 1
    filters->filters[0][filters->filters_size - 1] = category;
    filters->filters[1][filters->filters_size - 1] = 1;
}

// Function to read the dataset
DatasetInfo* read_dataset(const char *filename, filterInfo *filters) {
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
        add_or_increment_filter(filters, dataset->datapoints[i].category);
    }

    fclose(file);
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
