#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdint.h>
#include <pthread.h>

// ==================== Structures ==================== //

typedef struct {
    int filter_index;
    int count;
    int *point_indexes;
} filterPoint;

typedef struct{
    filterPoint *filtersPoints;
    int num_filters;
} filterInfo;

typedef struct {
    int point_index;
    int category;  // Categorical attribute C
    float timestamp;  // Timestamp attribute T
    float *vectors;     // Array for the remaining 100-dimensional vectors
} DataPoint;

// Define a struct for the result
typedef struct {
    int num_vectors;
    int num_dimensions;
    DataPoint *datapoints;
    filterInfo filterInfo;
} DatasetInfo;


typedef struct {
    int query_type; // 0: no attribute, 1: categorical, 2: timestamp, 3: categorical and timestamp
    int v; // Value for categorical attribute (or -1)
    float l; // Lower bound for timestamp attribute (or -1)
    float r; // Upper bound for timestamp attribute (or -1)
    float *query_vector; // 100-dimensional query vector
} QueryPoint;

// Define the QueryInfo structure
typedef struct {
    int num_queries;    // Number of queries
    int num_dimensions;
    QueryPoint *queries;
} QueryInfo;


// Struct for a single point
typedef struct {
    int index;              // Index of the point
    float *coordinates;     // Vector coordinates
    int category;           // Category of the point
    int edge_count;         // Current number of Outgoing edges
    int *edges;            // Array of Outgoing edges (connections to other points)
} Point;

typedef struct {
    int index;
    int category;
} FilteredMedoid;

typedef  struct {
    FilteredMedoid*metoids;
    int size;
} FilteredMethoidList;

// Struct for the whole graph
typedef struct {
    Point *points; // Array of points
    int num_dimensions; // Number of dimensions of each point
    int num_points; // Number of vectors in the graph
    FilteredMethoidList filteredMedoids;
    int  medoid;
} Graph;

typedef struct  {
    int * points;
    int num_points;
    int filter_Category;
}PointsPerFilter;


typedef struct {
    // The parameters you need for do_work:
    int L;
    float a;
    int R;
    Graph *graph;
    int medoid_index;
    FilteredMethoidList *filteredMedoids;

    // The array of per-filter subsets
    PointsPerFilter *pointsPerFilterArray;

    // How many filters total?
    int num_filters;

    // Current filter index (which filter to process next)
    int current_idx;

    // A mutex to protect current_idx
    pthread_mutex_t lock;
} ParallelData;

typedef struct {
    int id;
    DatasetInfo* dataset;
    int *filterList;
    int numOfFilters;
    int L_small;
    float a;
    int R_small;
    Graph* filter_graph;
    double totalTimeTaken;
    double min_time;
    double max_time;
    double avg_time;
} ThreadArgs;


// ==================== Structures ==================== //



#endif