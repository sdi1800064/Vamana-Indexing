#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>

#include "../headers/graph.h"


/**
 * Creates a random graph with the given number of vectors and dimensions
 * and random edges within the given maximum number of edges per point.
 * 
 * @param base_vectors The 2D array of vectors to use as points in the graph
 * @param base_num_dimensions The number of dimensions of each vector
 * @param max_edges The maximum number of edges each point can have
 * @param base_num_points The number of vectors in the graph
 * @return A pointer to the created graph
 */
Graph* create_random_graph(float **base_vectors, int base_num_dimensions, int max_edges, int base_num_points) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    
    if (!graph) {
        printf("Memory allocation failed!\n");
        exit(1);
    }

    graph->points = (Point*)malloc(base_num_points * sizeof(Point));

    if (!graph->points) {
        printf("Memory allocation failed!\n");
        exit(1);
    }

    graph->num_points = base_num_points;
    graph->num_dimensions = base_num_dimensions;

    // Initialize points with vector values
    for (int i = 0; i < base_num_points; i++) {
        graph->points[i].index = i;

        // Allocate memory for coordinates
        graph->points[i].coordinates = (float*)malloc(base_num_dimensions * sizeof(float));
        if (!graph->points[i].coordinates) {
            printf("Memory allocation failed!\n");
            exit(1);
        }

        // Allocate memory for edges
        graph->points[i].edges = (int*)malloc(max_edges * sizeof(int));
        if (!graph->points[i].edges) {
            printf("Memory allocation failed!\n");
            exit(1);
        }

        // Add random edges to every point untill they have max_edges, no duplicates and no self-loops
        graph->points[i].edge_count = 0; // No edges initially
        while(graph->points[i].edge_count < max_edges) {
            int edge = rand() % base_num_points;
            if (edge != i && !arrayContains(graph->points[i].edges, graph->points[i].edge_count, edge)) {
                graph->points[i].edges[graph->points[i].edge_count] = edge;
                graph->points[i].edge_count++;
            }
        }

        // Initialize coordinates
        for(int j = 0; j < base_num_dimensions; j++) {
            graph->points[i].coordinates[j] = base_vectors[i][j];
        }
        
        // if(i % 100 == 0){
        //     printf("Point %d Initialized\n ", i);
        // }
    }

    for (int i = 0; i < base_num_points; i++) {
        free(base_vectors[i]);
    }
    free(base_vectors);

    return graph;
}


/**
 * Function to print the graph coordinates for debugging
 * 
 * @param graph The graph to print
 * @param base_num_points The number of vectors in the graph
 * @param base_num_coordinates The number of coordinates per vector
 * @param outputfd The file descriptor to write to
 */
void fprint_graph_coordinates(Graph* graph, FILE *outputfd) {
    printf("Printing graph coordinates..\n");
    
    // Iterate over each point in the graph
    for (int i = 0; i < graph->num_points; i++) {
        Point p = graph->points[i];
        
        // Print the coordinates of the current point
        fprintf(outputfd, "Point %d: ( ", p.index);
        for(int j = 0; j < graph->num_dimensions; j++) {
            if (j == graph->num_dimensions - 1) {
                fprintf(outputfd, "%f )\n", p.coordinates[j]);
                break;
            }
            fprintf(outputfd, "%.2f, ", p.coordinates[j]);
        }
        
        // Print the edges of the current point
        fprintf(outputfd, " Point %d edges: ", p.index);
        for (int j = 0; j < p.edge_count; j++) {
            fprintf(outputfd, "%d ", p.edges[j]);
        }
        fprintf(outputfd, "\n");
    }
}

/**
 * Prints the edges of the graph for debugging
 * @param graph The graph to print
 * @param base_num_points The number of vectors in the graph
 * @param outputfd The file descriptor to write to
 */
void fprint_graph(Graph* graph, FILE *outputfd) {
    printf("Printing graph..\n");
    for (int i = 0; i < graph->num_points; i++) {
        Point p = graph->points[i];
        fprintf(outputfd, " Point %d edges: ", p.index);
        
        for (int j = 0; j < p.edge_count; j++) {
            fprintf(outputfd, "%d ", p.edges[j]);
        }
        fprintf(outputfd, "\n");
    }
}

// Add an edge to a point
void addEdge(Point *point, int toIndex) {
    for(int i = 0; i < point->edge_count; i++) {
        if (point->edges[i] == toIndex) return;
    }
    // printf("Adding edge from %d to %d | Edge count before adding: %d\n", point->index, toIndex, point->edge_count);
    point->edges[point->edge_count] = toIndex;
    point->edge_count += 1;
    // printf("Edge count after adding: %d\n", point->edge_count);

}

// Function to add an element to a dynamically resized array
void add_to_dynamic_array(int **array, int *size, int element) {
    if(arrayContains(*array, *size, element)) {
        return;
    }
    printf("adding element %d to array of size %d\n", element, *size);
    int old_size = *size;
    (*size)++;

    int *temp = *array;
    *array = realloc(*array, (*size + 1) * sizeof(int));
    if(!(*array)) {
        printf("Memory re-Allocation failed!\n");
        *array = temp;
    }
    (*array)[old_size] = element;
    printf("Added %d to position %d/%d \n", element, old_size, *size);

}


// Check if point `toIndex` is already in point's edges
int edgeExists(Point *point, int toIndex) {
    for (int i = 0; i < point->edge_count; i++) {
        if (point->edges[i] == toIndex) return 1;
    }
    return 0;
}

// The robustPrune function
void robustPrune(Graph *graph, int p_index, int *V, int V_size, float a, int R) {

    // Add to V all the neighbors of p and remove p from V if it exists
    // Set reset edges of p
    // While V is not empty
    // set as P* the closest point in V to p
    // Add P* to the neighbors of p
    // if the number of neighbors of p is equal to R then break
    // for every point p' in V
    // if a * d(p*, p') <= d(p, p') then remove p' from V

    // printf("-- Starting PRUNING --\n");

    // add to V all the neighbors of p and remove p from v if it exists
    printf("Pruning..\n");
    for (int i = 0; i < graph->points[p_index].edge_count; i++) {
        int toIndex = graph->points[p_index].edges[i];
        // Add to V all the neighbors of p
        if (!arrayContains(V, V_size, toIndex)) {
            add_to_dynamic_array(&V, &V_size, toIndex);
        }
    }

    // remove from V the point p if it exists
    int position;
    position = arrayContains(V, V_size, p_index);
    if(position != -1) {
        V[position] = V[V_size - 1];
        V_size--;
    }

    // Set reset edges of p
    for (int i = 0; i < graph->points[p_index].edge_count; i++) {
        graph->points[p_index].edges[i] = -1;
    }
    graph->points[p_index].edge_count = 0;


    Point *p_star = NULL;

    // While V is not empty
    while(V_size > 0) {

        // set as P* the closest point in V to p
        int min_distance = INF;
        int min_index = -1;
        for (int i = 0; i < V_size; i++) {
            p_star = &graph->points[V[i]];
            float distance = squared_euclidean_distance(graph->points[p_index].coordinates, p_star->coordinates, graph->num_dimensions);
            if (distance < min_distance) {
                min_distance = distance;
                min_index = i;
            }
        }
        p_star = &graph->points[V[min_index]];

        addEdge(&graph->points[p_index], p_star->index);

        // if the number of neighbors of p is equal to R then break
        if (graph->points[p_index].edge_count == R) {
            break;
        }

        // for every point p' in V
        for (int i = 0; i < V_size; i++) {
            // if a * d(p*, p') <= d(p, p') then remove p' from V
            // d(p*, p') = squared_euclidean_distance(graph->points[V[i]].coordinates, p_star->coordinates, graph->num_dimensions
            float DISTANCE_PSTAR_PPRIME = squared_euclidean_distance(p_star->coordinates, graph->points[V[i]].coordinates, graph->num_dimensions);
            float DISTANCE_P_PPRIME = squared_euclidean_distance(graph->points[p_index].coordinates, graph->points[V[i]].coordinates, graph->num_dimensions);
            if (a * DISTANCE_PSTAR_PPRIME <= DISTANCE_P_PPRIME) {
                V[i] = V[V_size - 1];
                V_size--;
            }
        }

    }

}


// Function to swap two elements in an array
void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Function to swap two elements in a float array
void swap_float(float *a, float *b) {
    float temp = *a;
    *a = *b;
    *b = temp;
}

// Sort function to sort array based on distance from Xq
void sort_array(Graph *graph, int *array, int array_size, float *Xq) {
    // Create an array to hold distances
    float *distances = (float *)malloc(array_size * sizeof(float));
    if (distances == NULL) {
        perror("Failed to allocate memory for distances");
        exit(EXIT_FAILURE);
    }

    // Calculate distances for each index in array and store in distances array
    for (int i = 0; i < array_size; i++) {
        int point_index = array[i];
        Point *point = &graph->points[point_index];
        distances[i] = squared_euclidean_distance(point->coordinates, Xq, graph->num_dimensions);
    }

    // Sort both the distances and indexes arrays using a basic selection sort
    for (int i = 0; i < array_size - 1; i++) {
        int min_idx = i;
        for (int j = i + 1; j < array_size; j++) {
            if (distances[j] < distances[min_idx]) {
                min_idx = j;
            }
        }
        
        // Swap distances
        swap_float(&distances[i], &distances[min_idx]);
        
        // Swap corresponding indexes in the array
        swap(&array[i], &array[min_idx]);
    }

    // Free the temporary distances array
    free(distances);

}


int *get_the_difference(int *Lamda, int Lamda_size, int *V, int V_size, int *Lamda_minus_V_size) {

    *Lamda_minus_V_size = 0;
    int *Lamda_minus_V = (int *)malloc(sizeof(int));
    if(Lamda_minus_V == NULL) {
        perror("Failed to allocate memory for Lamda_minus_V");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < Lamda_size; i++) {
        if(!arrayContains(V, V_size, Lamda[i]) && !arrayContains(Lamda_minus_V, *Lamda_minus_V_size, Lamda[i])) {
            add_to_dynamic_array(&Lamda_minus_V, Lamda_minus_V_size, Lamda[i]);
        }
    }
    return Lamda_minus_V;
}

void printArray(int *array, int array_size) {
    for (int i = 0; i < array_size; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");
}

// Function to update Lamda with closest L points to Xq
void update_lamda(int *lamda, int *lamda_size, Graph *graph, float *Xq_coords, int L) {
    // Sort Lamda based on distance to Xq and retain only closest L points
    for (int i = 0; i < *lamda_size - 1; i++) {
        for (int j = i + 1; j < *lamda_size; j++) {
            float dist_i = squared_euclidean_distance(graph->points[lamda[i]].coordinates, Xq_coords, graph->num_dimensions);
            float dist_j = squared_euclidean_distance(graph->points[lamda[j]].coordinates, Xq_coords, graph->num_dimensions);
            if (dist_i > dist_j) {
                int temp = lamda[i];
                lamda[i] = lamda[j];
                lamda[j] = temp;
            }
        }
    }
    // Truncate Lamda to retain only the L closest points
    if (*lamda_size > L) {
        *lamda_size = L;
    }
}

int* greedy_search(Graph *graph, int start_index, float *Xq, int L, int *V_size, int **V) {
    // Allocate space for Lamda and V (both are arrays of indexes)
    printf("===========GREEDY SEARCH===========\n");
    int *lamda = malloc(sizeof(int));
    if(lamda == NULL) {
        perror("Failed to allocate memory for lamda");
        exit(EXIT_FAILURE);
    }
    int lamda_size = 0;
    
    if(*V != NULL) {
        free(*V);
    }
    *V = NULL;

    *V = malloc(sizeof(int));
    if(V == NULL) {
        perror("Failed to allocate memory for V");
        exit(EXIT_FAILURE);
    }
    V_size = 0;
    
    // Initialize Lamda with the start point index
    add_to_dynamic_array(&lamda, &lamda_size, start_index);
    
    while (lamda_size > 0) {
        // Find the closest point in Lamda (by index) that is not in V
        printf("Looping..\n");
        int closest_index = -1;
        float min_distance = FLT_MAX;
        
        for (int i = 0; i < lamda_size; i++) {
            int point_index = lamda[i];
            printf("Here 1\n");
            if(arrayContains(*V, *V_size, point_index)){
                printf("Here 2\n");
                continue;
            }
            printf("Here 10\n");
            // Calculate distance to Xq
            float distance = squared_euclidean_distance(graph->points[point_index].coordinates, Xq, graph->num_dimensions);
            printf("Here 11\n");
            if (distance < min_distance) {
                min_distance = distance;
                closest_index = i;
            }
        }
        printf("Here 3\n");
        // If no unvisited points are found, exit the loop
        if (closest_index == -1) break;
        
        // Set P* as the closest point and add it to V
        int P_star_index = lamda[closest_index];
        add_to_dynamic_array(V, V_size, P_star_index);
        printf("Here 4\n");
        // Add neighbors of P* to Lamda
        Point *P_star = &graph->points[P_star_index];
        for (int i = 0; i < P_star->edge_count; i++) {
            int neighbor_index = P_star->edges[i];
            // Check if neighbor is already in Lamda or V
            if (!arrayContains(lamda, lamda_size, neighbor_index) && !arrayContains(*V, *V_size, neighbor_index)) {
                add_to_dynamic_array(&lamda, &lamda_size, neighbor_index);
                printf("Here 5\n");
            }
        }
        
        // If Lamda size exceeds L, keep only the closest L points to Xq
        if (lamda_size > L) {
            sort_array(graph, lamda, lamda_size, Xq);
            lamda_size = L;
            printf("Here 6\n");
        }
    }
    printf("Here 7\n");
    // Set the output result size
    printf("printing lamda before exiting greedy search\n");
    for(int i = 0; i < lamda_size; i++) {
        printf("%d ", lamda[i]);
    }
    printf("\n");
    // Return V as the list of closest points
    printf("==========FINISHED GREEDY SEARCH==================\n");

    return lamda;
}

// void greedy_search(Graph *graph, float *Xq, int start_index, int **V, int *V_size, int **Lamda, int L, int k) {

//     // Initialiase Lamda with the start point
//     // Initialise V to NULL
//     // while there are points in Lamda not yet visited do
//     //  set P* as the closest point from Lamda but not from V to Xq
//     //  add P* to V
//     //  add the neighbors of P* to Lamda
//     //  if the number of points in Lamda > L
//     //      then update Lampda to retain closest L points to Xw
//     // return closest k points from Lamda and V

//     printf("==========GREEDY SEARCH==========\n");


    // if (!graph || !Xq || !V || !V_size || !Lamda ) {
    //     fprintf(stderr, "Invalid input: null pointer\n");
    // }

    // printf("Starting from point %d | V_size %d ", start_index, *V_size);
    // // Allocate initial space for V and Lamda if needed
    // if (*V == NULL || *Lamda == NULL) {
    //     *V = (int *)malloc(sizeof(int));
    //     *Lamda = (int *)malloc(sizeof(int));
    //     if (!*V || !*Lamda) {
    //         perror("Failed to allocate memory for V or Lamda");
    //         free(*V);
    //         free(*Lamda);
    //     }
    // }

    // *V_size = 0;
    // int *Lamda_size = (int*)malloc(sizeof(int));
    // if(Lamda_size == NULL) {
    //     perror("Failed to allocate memory for Lamda_size");
    // }
    // *Lamda_size = 0;

    // bool *visited = (bool *)calloc(graph->num_points, sizeof(bool));
    // if (!visited) {
    //     perror("Failed to allocate memory for visited");
    // }

    // // Initialize Lamda with the start point
    // add_to_dynamic_array(Lamda, Lamda_size, start_index);

    // int *Lamda_minus_V = NULL;
    // int Lamda_minus_V_size = 0;
    // if(!visited[(*Lamda)[0]]){
    //     add_to_dynamic_array(&Lamda_minus_V, &Lamda_minus_V_size, (*Lamda)[0]);
    // }

    // printArray(*Lamda, *Lamda_size);

    // Point *p_star = NULL;

    // // Main loop: Continue until Lamda_minus_V is empty
    // while (Lamda_minus_V_size > 0) {
    //     printf("__________________Looping________________\n");
    //     // Select the point p* from Lamda_minus_V that is closest to the query point Xq
    //     int closest_index = -1;
    //     float min_distance = FLT_MAX;
    //     for (int i = 0; i < Lamda_minus_V_size; i++) {
    //         int current_index = Lamda_minus_V[i];
    //         printf("Current index: %d\n", current_index);
    //         float distance = squared_euclidean_distance(graph->points[current_index].coordinates, Xq, graph->num_dimensions);
    //         if (distance < min_distance) {
    //             min_distance = distance;
    //             closest_index = current_index;
    //         }
    //     }

        
    //     // Set p* as the closest point
    //     printf("Moved to another point: %d\n", closest_index);
    //     p_star = &graph->points[closest_index];
    //     printf("Neighbors of p*: ");
    //     for (int i = 0; i < p_star->edge_count; i++) {
    //         printf("%d ", p_star->edges[i]);
    //     }
    //     printf("\n");

    //     // Add p* to V
    //     visited[closest_index] = true;

    //     // Add the neighbors of p* to Lamda
    //     for (int i = 0; i < p_star->edge_count; i++) {
    //         int neighbor_index = p_star->edges[i];
    //         if(visited[neighbor_index] == false){
    //             add_to_dynamic_array(Lamda, Lamda_size, neighbor_index);
    //             printf("Printing Lamda for every new neighbor added\n");
    //             printArray(*Lamda, *Lamda_size);
    //         }
    //     }

    //     // Check if |Lamda| > L
    //     if (*Lamda_size > L) {
    //         check_for_duplicates(*Lamda, *Lamda_size);
    //         sort_array(graph, *Lamda, *Lamda_size, Xq);
    //         printArray(*Lamda, *Lamda_size);
    //         *Lamda = (int *)realloc(*Lamda, L * sizeof(int));
    //         if (!*Lamda) {
    //             perror("Failed to reallocate memory for Lamda");
    //             free(*V);
    //             free(Lamda_minus_V);
    //         }
    //         *Lamda_size = L;
    //         printf("Printing Lamda after sorting\n");
    //         printArray(*Lamda, *Lamda_size);
    //     }
        
    //     // Update Lamda_minus_V
    //     free(Lamda_minus_V);
    //     Lamda_minus_V = NULL;
    //     Lamda_minus_V = (int*)malloc(sizeof(int));
    //     Lamda_minus_V_size = 0;
    //     for(int i = 0; i < *Lamda_size; i++) {
    //         if(!visited[(*Lamda)[i]]) {
    //             printf("Adding to Lamda_minus_V: %d\n", (*Lamda)[i]);
    //             add_to_dynamic_array(&Lamda_minus_V, &Lamda_minus_V_size, (*Lamda)[i]);
    //         }
    //     }
    //     printf("Printing Lamda_minus_V after finding new points\n");
    //     printArray(Lamda_minus_V, Lamda_minus_V_size);

    //     if (!Lamda_minus_V) {
    //         perror("Failed to allocate memory for Lamda_minus_V");
    //         free(*V);
    //         free(*Lamda);
    //     } 
    // }

    // free(Lamda_minus_V);

    // for(int i = 0; i < graph->num_points; i++) {
    //     if(visited[i]) {
    //         add_to_dynamic_array(V, V_size, i);
    //     }
    // }

    // printArray(*V, *V_size);

    // // Sort V and Lamda
    // sort_array(graph, *V, *V_size, Xq);
    // sort_array(graph, *Lamda, *Lamda_size, Xq);

    // (*Lamda) = (int*)realloc((*Lamda), k * sizeof(int));

    // if(!(*Lamda)) {
    //     perror("Failed to reallocate memory for Lamda");
    //     free(*V);
    // }

    // *Lamda_size = k;

    // printArray(*Lamda, *Lamda_size);
    // exit(1);
// }

/**
 * Checks if a node is in the visited array V
 * 
 * @param V The V array
 * @param V_size The size of the V array
 * @param node The index to search for
 * @return 1 if the node is in the array, 0 otherwise
 */
int arrayContains(int *V, int V_size, int node) {
    // Iterate through the array and check if the node is present
    printf("ArrayContains CALLED\n");
    if(V_size == 0) return 0;
    for (int i = 0; i < V_size; i++) {
        if (V[i] == node) {
            return 1;  // Node found
        }
    }
    return 0;  // Node not found
}


/**
 * Calculates the squared Euclidean distance between two vectors of length n
 * 
 * @param p The first vector
 * @param q The second vector
 * @param n The length of the vectors
 * @return The squared Euclidean distance between the two vectors
 */
double squared_euclidean_distance(float *p, float *q, int n) {
    // Calculate the difference between the elements of the two vectors
    float sum = 0.0f;
    for (int i = 0; i < n; i++) {
        float diff = p[i] - q[i];  // Calculate the difference
        sum += diff * diff;        // Square the difference and add to sum
    }
    return sum;  // Return the squared Euclidean distance
}

/**
 * Calculates the medoid from a set of sampled points
 *
 * The medoid is the point from the set of sampled points which has the smallest
 * sum of distances to all other points in the set.
 *
 * @param graph The graph that contains the sampled points
 * @param sample_point_indexes The indexes of the sampled points in the graph
 * @param num_sample_points The number of sampled points
 * @return The index of the medoid in the original graph
 */
int calculate_medoid(Graph *graph, int *sample_point_indexes, int num_sample_points) {
    
    float *distance_sums = (float *)malloc(num_sample_points * sizeof(float));

    // Initialize the distance sums to zero
    for (int i = 0; i < num_sample_points; i++) {
        distance_sums[i] = 0.0;
    }


    // Loop through all pairs of sampled points and calculate the pairwise distances
    for (int i = 0; i < num_sample_points; i++) {
        Point *point_i = &graph->points[sample_point_indexes[i]];
        for (int j = 0; j < num_sample_points; j++) {
            if (i != j) {
                Point *point_j = &graph->points[sample_point_indexes[j]];
                float dist = squared_euclidean_distance(point_i->coordinates, point_j->coordinates, graph->num_dimensions);
                distance_sums[i] += dist;
            }
        }
    }

    // Find the point with the smallest sum of distances (i.e., the medoid)
    int medoid_index = 0;
    float min_distance_sum = distance_sums[0];
    for (int i = 1; i < num_sample_points; i++) {
        if (distance_sums[i] < min_distance_sum) {
            min_distance_sum = distance_sums[i];
            medoid_index = i;
        }
    }

    // Free the allocated memory
    free(distance_sums);

    // Return the index of the medoid in the original graph
    return sample_point_indexes[medoid_index];
}


/**
 * Function to randomly sample points from the graph
 *
 * @param graph The graph to sample points from
 * @param num_sample_points The number of sample points to select
 * @return An array of indices representing the sampled points
 */
int* sample_points(int max, int num_sample_points) {
    int *sample_point_indexes = (int *)calloc(num_sample_points, sizeof(int));

    // Check if memory allocation was successful
    if (sample_point_indexes == NULL) {
        printf("Memory allocation failed!\n");
        exit(1);
    }

    // Temporary index for random selection
    int s_index = 0;
    int i = 0;
    while (i < num_sample_points) {
        // Randomly select an index from the graph
        s_index = rand() % max;
        // Ensure no duplicate indices are selected
        if (!arrayContains(sample_point_indexes, i, s_index)) {
            sample_point_indexes[i] = s_index;
            i++;
        }
    }
    return sample_point_indexes;
}

// Function to check for duplicates in an array
void check_for_duplicates(int *array, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            if (array[i] == array[j]) {
                printf("Duplicate found: %d\n", array[i]);
                return;
            }
        }
    }
}

void vamana_indexing(Graph *graph, int L, float a, int R, FILE *outputfd) {

    printf("===============Starting Vamana Indexing===============\n");

    printf("L = %d\na = %f\nR = %d\n", L, a, R);

    // get the medoid of the graph
    // traverse the graph in a random way without repetitions
    // for every random point in the graph
    // run greedysearch(medoid, random point, k=1, L) to get the visited list V
    // run robustPrune(random Point, Visited list, a, R)
    // for every neighbor of random point p'
    // if | neighbor(p') U random point | > R then
    // run robustPrune(p', Neighbors of p' U random point, a, R)
    // else add to neighbors of p' the random point
    // end for

    // Calculating the medoid
    int *sample_point_indexes = sample_points(graph->num_points, graph->num_points / 10);
    int medoid_index = calculate_medoid(graph, sample_point_indexes, graph->num_points / 10);
    printf("Medoid index: %d\n", medoid_index);

    fprintf(outputfd,  "\n\nGraph before indexing:\n");
    fprint_graph(graph, outputfd);

    // traverse the graph in a random way without repetitions
    
    bool *shuffled_point_indexes = (bool *)calloc(graph->num_points, sizeof(bool));
    int s_index;
    int i=0;
    
    while(i < graph->num_points) {
        s_index = rand() % graph->num_points;
        if (!shuffled_point_indexes[s_index]) {
            shuffled_point_indexes[s_index] = true;

            // run greedysearch(medoid, random point, k=1, L) to get the visited list V
            int *V = NULL;
            int V_size = 0;
            // int *Lamda = NULL;

            greedy_search(graph, medoid_index, graph->points[s_index].coordinates, L, &V_size, &V);
            // greedy_search(graph, graph->points[s_index].coordinates, medoid_index, &V, &V_size, &Lamda, L, 1);

            // run robustPrune(random Point, Visited list, a, R)
            robustPrune(graph, s_index, V, V_size, a, R);

            int *new_V = NULL;
            int new_V_size = 0;

            // for every neighbor of random point p'
            for(int j = 0; j < graph->points[s_index].edge_count; j++) {

                // if | neighbor(p') U random point | > R then
                Point *P_PRIME = &graph->points[graph->points[s_index].edges[j]];
                
                if(P_PRIME->edge_count + 1 > R) {
                    // run robustPrune(p', Neighbors of p' U random point, a, R)
                    // create a new visited list that contains only the neighbors of p' and the random point
                    new_V = NULL;
                    new_V = (int *)malloc((P_PRIME->edge_count + 1) * sizeof(int));
                    new_V_size = 0;
                    for (int k = 0; k < P_PRIME->edge_count; k++) {
                        if (!arrayContains(new_V, new_V_size, P_PRIME->edges[k])) {
                            add_to_dynamic_array(&new_V, &new_V_size, P_PRIME->edges[k]);
                        }
                    }
                    add_to_dynamic_array(&new_V, &new_V_size, s_index);

                    // run robustPrune(p', new_V, a, R)
                    robustPrune(graph, P_PRIME->index, new_V, new_V_size, a, R);
                    
                } else{
                    // else add to neighbors of p' the random point
                    addEdge(P_PRIME, s_index);
                }
            }
            // Traversed to another point of the graph
            i++;
            printf("Points remaining: %d\n", graph->num_points - i);
            // free(new_V);
        
        }

    }


    fprintf(outputfd, "\n\nGraph after Vamana Indexing\n");
    fprint_graph(graph, outputfd);

    free(shuffled_point_indexes);
}


// Function to return a pointer to a new array containing the row at index `num`
int* get_kclosest(int** kclosest, int num, int k) {
    // Allocate memory for the result array of size k
    int* result = (int*)malloc(sizeof(int) * k);
    
    // Copy the elements from the row `num` into `result`
    for (int i = 0; i < k; i++) {
        result[i] = kclosest[num][i];
    }
    
    return result;
}