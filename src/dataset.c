#include <stdio.h>
#include <stdlib.h>
#include "../headers/dataset.h"
#include "../headers/graph.h"

void sort_filter_info(DatasetInfo* dataset) {
    filterInfo* filters = &dataset->filterInfo;
    int num_filters = filters->num_filters;

    // Use bubble sort to sort the filters based on filter_index
    for (int i = 0; i < num_filters - 1; i++) {
        for (int j = 0; j < num_filters - i - 1; j++) {
            if (filters->filtersPoints[j].filter_index > filters->filtersPoints[j + 1].filter_index) {
                // Swap the filterPoints
                filterPoint temp = filters->filtersPoints[j];
                filters->filtersPoints[j] = filters->filtersPoints[j + 1];
                filters->filtersPoints[j + 1] = temp;
            }
        }
    }
}

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
    printf("Reading dataset from %s: ", filename);
    fflush(stdout);
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
        dataset->datapoints[i].point_index = i;
        dataset->datapoints[i].category = (int)buffer[0];
        dataset->datapoints[i].timestamp = buffer[1];
        for (int j = 0; j < 100; j++) {
            dataset->datapoints[i].vectors[j] = buffer[j + 2];
        }

        // Update the filters with the category
        // printf("Updating index %d to filter index %d.\n", i, dataset->datapoints[i].category);
        count++;
        add_or_increment_filter(&dataset->filterInfo, dataset->datapoints[i].category, dataset->datapoints[i].point_index);
        
    }

    sort_filter_info(dataset);
    fclose(file);
    printf("Done.\n");
    fflush(stdout);
    return dataset;
}

// Function to free the allocated dataset
void free_dataset(DatasetInfo *dataset) {
    free(dataset->datapoints);
    free(dataset);
}

void print_dataset(DatasetInfo *dataset) {
    for (int i = 0; i < dataset->num_vectors; i++) {
        printf("index %d: Category = %d, Timestamp = %f, Vectors :\n ", dataset->datapoints[i].point_index, dataset->datapoints[i].category, dataset->datapoints[i].timestamp);
        for (int j = 0; j < 100; j++) {
            printf("%f ", dataset->datapoints[i].vectors[j]);
        }
        printf("\n\n");
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
    printf("Reading queryset from %s: ", filename);
    fflush(stdout);
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
    // printf("Number of queries: %d\n", num_queries);

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
    printf("Done.\n");
    fflush(stdout);
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
        // for (int j = 0; j < 100; j++) {
        //     printf("%f ", query.query_vector[j]);
        // }
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


int** readGroundTruth(char* filename, int num_queries) {
    printf("Reading ground truth from file %s: ", filename);
    fflush(stdout);
    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }
    int** ground_truth = (int**)malloc(num_queries * sizeof(int*));
    for(int i = 0; i < num_queries; i++) {
        ground_truth[i] = (int*)malloc(100 * sizeof(int));
        for(int j = 0; j < 100; j++) {
            if(fread(&ground_truth[i][j], sizeof(int), 1, file) != 1) {
                perror("Error reading number of vectors");
                fclose(file);
                return NULL;
            }   
        }
    }
    fclose(file);
    printf("Done.\n");
    fflush(stdout);
    return ground_truth;
}


void writeGraphs(Graph* graph, int num_of_graphs, char* filename) {
    printf("Writing graphs to file %s: ", filename);
    fflush(stdout);
    FILE *file = fopen(filename, "wb");
    if (!file) {
        perror("Error opening file");
        return;
    }

    fwrite(&num_of_graphs, sizeof(int), 1, file);
    for(int i = 0; i < num_of_graphs; i++){
        fwrite(&graph[i].num_points, sizeof(int), 1, file);
        fwrite(&graph[i].num_dimensions, sizeof(int), 1, file);
        for(int j = 0; j < graph[i].num_points; j++) {
            fwrite(&graph[i].points[j].index, sizeof(int), 1, file);
            fwrite(&graph[i].points[j].category, sizeof(int), 1, file);
            fwrite(&graph[i].points[j].edge_count, sizeof(int), 1, file);
            for(int k = 0; k < graph[i].points[j].edge_count; k++) {
                fwrite(&graph[i].points[j].edges[k], sizeof(int), 1, file);
            }
            for(int k = 0; k < graph[i].num_dimensions; k++) {
                fwrite(&graph[i].points[j].coordinates[k], sizeof(float), 1, file);
            }
        }
    }
    fclose(file);
    printf("Done.\n");
    fflush(stdout);
}

Graph* readGraphs(char* filename, int* num_of_graphs) {
    printf("Reading graphs from file %s: ", filename);
    fflush(stdout);

    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Read the number of graphs
    if (fread(num_of_graphs, sizeof(int), 1, file) != 1) {
        perror("Error reading number of graphs");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // printf("Number of graphs to be read: %d\n", *num_of_graphs);
    Graph* graph = (Graph*)malloc(*num_of_graphs * sizeof(Graph));
    if (!graph) {
        perror("Memory allocation error for graphs");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < *num_of_graphs; i++) {
     
        // Read the number of points
        if (fread(&graph[i].num_points, sizeof(int), 1, file) != 1) {
            perror("Error reading number of points");
            fclose(file);
            exit(EXIT_FAILURE);
        }

        // Read the number of dimensions
        if (fread(&graph[i].num_dimensions, sizeof(int), 1, file) != 1) {
            perror("Error reading number of dimensions");
            fclose(file);
            exit(EXIT_FAILURE);
        }

        // Allocate memory for the points
        graph[i].points = (Point *)malloc(graph[i].num_points * sizeof(Point));
        if (!graph[i].points) {
            perror("Memory allocation error for points");
            fclose(file);
            exit(EXIT_FAILURE);
        }

        // Read each point from the file
        for (int j = 0; j < graph[i].num_points; j++) {
            Point *point = &graph[i].points[j];

            // Read the index
            if (fread(&point->index, sizeof(int), 1, file) != 1) {
                perror("Error reading index");
                fclose(file);
                exit(EXIT_FAILURE);
            }

            // Read the category
            if (fread(&point->category, sizeof(int), 1, file) != 1) {
                perror("Error reading category");
                fclose(file);
                exit(EXIT_FAILURE);
            }

            // Read the edge count
            if (fread(&point->edge_count, sizeof(int), 1, file) != 1) {
                perror("Error reading edge count");
                fclose(file);
                exit(EXIT_FAILURE);
            }

            // Allocate memory for the edges
            point->edges = (int *)malloc(point->edge_count * sizeof(int));
            if (!point->edges) {
                perror("Memory allocation error for edges");
                fclose(file);
                exit(EXIT_FAILURE);
            }

            // Read the edges
            for (int k = 0; k < point->edge_count; k++) {
                if (fread(&point->edges[k], sizeof(int), 1, file) != 1) {
                    perror("Error reading edge");
                    fclose(file);
                    exit(EXIT_FAILURE);
                }
            }

            // Allocate memory for the coordinates
            point->coordinates = (float *)malloc(graph[i].num_dimensions * sizeof(float));
            if (!point->coordinates) {
                perror("Memory allocation error for coordinates");
                fclose(file);
                exit(EXIT_FAILURE);
            }

            // Read the coordinates
            for (int k = 0; k < graph[i].num_dimensions; k++) {
                if (fread(&point->coordinates[k], sizeof(float), 1, file) != 1) {
                    perror("Error reading coordinate");
                    fclose(file);
                    exit(EXIT_FAILURE);
                }
            }
        }
   
    }
    fclose(file);
    printf("Done.\n");
    fflush(stdout);
    return graph;
}
