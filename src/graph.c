#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#include "../headers/graph.h"
#include "../headers/dataset.h"
#include "../headers/structs.h"



/**
 * Initializes a graph with the given dataset and maximum number of edges
 *
 * @param dataset The dataset to use for initializing the graph
 * @param max_edges The maximum number of edges that each point in the graph can have
 * @return A pointer to the newly created graph
 */
Graph initialise_graph(DatasetInfo* dataset, int max_edges) {
    Graph graph;

    int base_num_points = dataset->num_vectors;
    int base_num_dimensions = 100;

    graph.points = (Point*)malloc(base_num_points * sizeof(Point));

    if (!graph.points) {
        printf("Memory allocation failed!\n");
        exit(1);
    }

    graph.num_points = base_num_points;
    graph.num_dimensions = base_num_dimensions;

    // Initialize points with dataset values
    for (int i = 0; i < base_num_points; i++) {
        graph.points[i].index = i;
        graph.points[i].category = dataset->datapoints[i].category;

        // Allocate memory for coordinates
        graph.points[i].coordinates = (float*)malloc(base_num_dimensions * sizeof(float));
        if (!graph.points[i].coordinates) {
            printf("Memory allocation failed!\n");
            exit(1);
        }

        // Initialize coordinates
        for (int j = 0; j < base_num_dimensions; j++) {
            graph.points[i].coordinates[j] = dataset->datapoints[i].vectors[j];
        }

        // Initialize edges with no connections
        graph.points[i].edges = malloc(max_edges * sizeof(int));
        graph.points[i].edge_count = 0;
    }
    return graph;
}

// Free the memory allocated for the graph
void free_graph(Graph graph) {
    for (int i = 0; i < graph.num_points; i++) {
        free(graph.points[i].coordinates);
        free(graph.points[i].edges);
    }
    free(graph.points);
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
        
        // Print the edges of the current point
        fprintf(outputfd, " Point %d | category: %d | edges: ", p.index, p.category);
        for (int j = 0; j < p.edge_count; j++) {
            fprintf(outputfd, "%d ", p.edges[j]);
        }
        printf("\n");
        // Print the coordinates of the current point
        fprintf(outputfd, "Coordinates: ( ");
        for(int j = 0; j < graph->num_dimensions; j++) {
            if (j == graph->num_dimensions - 1) {
                fprintf(outputfd, "%f )\n\n", p.coordinates[j]);
                break;
            }
            fprintf(outputfd, "%.2f, ", p.coordinates[j]);
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
    for (int i = 0; i < graph->num_points; i++) {
        Point p = graph->points[i];
        fprintf(outputfd, " Point %d | category: %d | edges: ", p.index, p.category);
        
        for (int j = 0; j < p.edge_count; j++) {
            fprintf(outputfd, "%d ", p.edges[j]);
        }
        fprintf(outputfd, "\n");
    }
}


/**
 * Adds an edge from the given point to the point with index `toIndex` if that edge does not already exist.
 * 
 * @param point The point to add the edge from
 * @param toIndex The index of the point to add the edge to
 */
void addEdge(Point *point, int toIndex) {
    for(int i = 0; i < point->edge_count; i++) {
        if (point->edges[i] == toIndex) return;
    }
    // printf("Adding edge from %d to %d | Edge count before adding: %d\n", point->index, toIndex, point->edge_count);
    point->edges[point->edge_count] = toIndex;
    point->edge_count += 1;
    // printf("Edge count after adding: %d\n", point->edge_count);

}

/**
 * Adds an element to a dynamic array if it is not already present.
 *
 * This function checks if the given element is already in the array.
 * If not, it reallocates memory to increase the array size by one and
 * adds the element to the end of the array. If memory reallocation fails,
 * the original array is preserved.
 *
 * @param array A pointer to the dynamic array to be modified.
 * @param size A pointer to the current size of the array, which is updated.
 * @param element The element to be added to the array.
 */

void add_to_dynamic_array(int **array, int *size, int element) {
    // printf("Adding %d to array of size %d\n", element, *size);

    if(arrayContains(*array, *size, element)) {
        // printf("Element %d already in array\n", element);
        return;
    }
    int old_size = *size;
    (*size)++;

    int *temp = *array;
    *array = realloc(*array, (*size + 1) * sizeof(int));
    if(!(*array)) {
        printf("Memory re-Allocation failed!\n");
        *array = temp;
        printf("Reverted to old array\n");
    }
    (*array)[old_size] = element;
    // printf("Array now has size %d\n", *size);

}


// Check if point `toIndex` is already in point's edges
int edgeExists(Point *point, int toIndex) {
    for (int i = 0; i < point->edge_count; i++) {
        if (point->edges[i] == toIndex) return 1;
    }
    return 0;
}


/**
 * Robust Pruning Algorithm
 * 
 * Add to V all the neighbors of p and remove p from V if it exists
 * Set reset edges of p
 * While V is not empty
 * set as P* the closest point in V to p
 * Add P* to the neighbors of p
 * if the number of neighbors of p is equal to R then break
 * for every point p' in V
 * if a * d(p*, p') <= d(p, p') then remove p' from V
 *
 * @param graph The graph containing the point p
 * @param p_index The index of the point p
 * @param V The set of points that are candidates to be neighbors of p
 * @param V_size The size of V
 * @param a The pruning parameter
 * @param R The maximum number of neighbors of p
 */
int filtered_Robust_prune(Graph *graph, int p_index, int *V, int V_size, float a, int R) {

    if(V_size == 0) {
        return 0;
    }

    // add to V all the neighbors of p and remove p from v if it exists
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

    // Reset edges of p
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
            // if Fp' ^ Fp !C Fp* then continue ( since we only have one filter we check
            //                                    if the filter of p' is the same as the filter of p
            //                                    and the same as the filter of p* )
            int F_p_prime = graph->points[V[i]].category;
            int F_p = graph->points[p_index].category;
            int F_p_star = p_star->category;

            if ( F_p_prime == F_p && F_p_prime == F_p_star ){
                continue;
            }

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
    return 1;
}


/**
 * Swaps the values of two int pointers.
 *
 * @param a The first int pointer.
 * @param b The second int pointer.
 */
void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}


/**
 * Swaps the values of two float pointers.
 *
 * @param a The first float pointer.
 * @param b The second float pointer.
 */
void swap_float(float *a, float *b) {
    float temp = *a;
    *a = *b;
    *b = temp;
}


/**
 * Sorts the given array of point indexes in ascending order of their squared Euclidean distances from the given query point Xq.
 *
 * The function first calculates the squared Euclidean distances for each point in the array and stores them in a temporary array.
 * It then uses a basic selection sort to sort both the distances array and the input array of indexes.
 * After sorting, the temporary distances array is freed.
 *
 * @param graph The graph containing the points.
 * @param array The array of point indexes to sort.
 * @param array_size The size of the array.
 * @param Xq The query point.
 */
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

int sort_array_based_on_dataset(DatasetInfo *dataset, int *array, int array_size, float *Xq) {
    if(dataset->datapoints == 0 || dataset == NULL) {
        return 0;
    }
    // Create an array to hold distances
    float *distances = (float *)malloc(array_size * sizeof(float));
    if (distances == NULL) {
        perror("Failed to allocate memory for distances");
        exit(EXIT_FAILURE);
    }

    // Calculate distances for each index in array and store in distances array
    for (int i = 0; i < array_size; i++) {
        int point_index = array[i];
        DataPoint *point = &dataset->datapoints[point_index];
        distances[i] = squared_euclidean_distance(point->vectors, Xq, 100);
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
    return 1;
}



/**
 * Calculates the difference between two arrays of integers, Lamda and V.
 *
 * The returned array, Lamda_minus_V, contains the elements of Lamda that are not in V.
 * The size of Lamda_minus_V is stored in *Lamda_minus_V_size.
 *
 * @param Lamda The first array of integers.
 * @param Lamda_size The size of Lamda.
 * @param V The second array of integers.
 * @param V_size The size of V.
 * @param Lamda_minus_V_size The size of the returned array.
 * @return The returned array.
 */
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


/**
 * Prints the elements of an integer array.
 *
 * This function iterates through the given array and prints each element
 * followed by a space. After printing all elements, it prints a newline.
 *
 * @param array Pointer to the integer array to be printed.
 * @param array_size The number of elements in the array.
 */
void printArray(int *array, int array_size) {
    for (int i = 0; i < array_size; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");
}


/**
 * Performs a greedy search on the graph to find the L closest points to a given query point Xq.
 *
 * The algorithm works by selecting the point in Lamda that is closest to Xq, adding it to V, and then adding the neighbors of that point to Lamda.
 * It continues until all points in Lamda are visited or the limit L is reached.
 *
 * @param graph The graph to search.
 * @param Xq The query point.
 * @param start_index The index of the starting point of the search.
 * @param V The list of visited points.
 * @param V_size The size of V.
 * @param Lamda The list of points to consider.
 * @param Lamda_size The size of Lamda.
 * @param L The limit on the number of points to consider.
 */
void filtered_greedy_search(Graph *graph, float *Xq, int* start_index, int start_index_size, int **V, int *V_size, int **Lamda, int *Lamda_size, int L, int query_filter) {

    // Allocate initial space for V and Lamda
    *V = (int *)malloc(sizeof(int) * graph->num_points);
    *V_size = 0;
    *Lamda = (int *)malloc(sizeof(int) * graph->num_points);
    *Lamda_size = 0;


    // for s ε S : if Fs ^ Fx != 0 then add s to Lamda
    for(int i = 0 ; i < start_index_size; i++) {
        if(graph->points[start_index[i]].category == query_filter || query_filter == -1) {
            // Initialize Lamda with the start point
            add_to_dynamic_array(Lamda, Lamda_size, start_index[i]);
        }
    }


    // Create an int array that contains points that are in Lamda but not in V
    int *Lamda_minus_V;
    int Lamda_minus_V_size = 0;
    Lamda_minus_V = get_the_difference(*Lamda, *Lamda_size, *V, *V_size, &Lamda_minus_V_size);


    // for every point in Lamda_minus_V
    // calculate the distance from xq to that point
    // the one that is closest is p*
    // Take P* 's   Neighbors and add them to Lamda
    // Take P* and add it to V

    // Main loop: Continue until Lamda_minus_V is empty
    while (Lamda_minus_V_size > 0) {
        // Select the point p* from Lamda_minus_V that is closest to the query point Xq
        int closest_index = -1;
        float min_distance = FLT_MAX;
        for (int i = 0; i < Lamda_minus_V_size; i++) {
            int current_index = Lamda_minus_V[i];
            float distance = squared_euclidean_distance(graph->points[current_index].coordinates, Xq, graph->num_dimensions);
            if (distance < min_distance) {
                min_distance = distance;
                closest_index = current_index;
            }
        }
        Point *p_star = &graph->points[closest_index];

        // Add p* to V
        add_to_dynamic_array(V, V_size, p_star->index);

        // Copy the current neighbors of p* to a temp array
        int *temp_neighbors = (int *)malloc(sizeof(int) * p_star->edge_count);
        int temp_neighbors_size = p_star->edge_count;
        memcpy(temp_neighbors, p_star->edges, sizeof(int) * p_star->edge_count); 

        p_star->edge_count = 0;       

        // Keep the neighbors of p* that have at least one common filter with Xq and are not visited
        // Nout'(p*) <- ( p' ε Nout(p*) : Fs ^ Fx != 0, p' !ε V )
        for (int i = 0; i < temp_neighbors_size; i++) {
            Point *neighbor = &graph->points[temp_neighbors[i]];
            if (neighbor->category == query_filter || query_filter == -1) {
                addEdge(p_star, neighbor->index);
            }
        }

        free(temp_neighbors);

        // Add the neighbors of p* to Lamda
        // Lamda <- Lamda U Nout(p*)
        for (int i = 0; i < p_star->edge_count; i++) {
            int neighbor_index = p_star->edges[i];
            add_to_dynamic_array(Lamda, Lamda_size, neighbor_index);
        }

        // Check if |Lamda| > L
        if (*Lamda_size > L) {
            // Update Lamda to contain the L closest points to Xq
            sort_array(graph, *Lamda, *Lamda_size, Xq);
            int *temp = *Lamda;
            *Lamda_size = L;
            *Lamda = (int *)realloc(*Lamda, L * sizeof(int));
            if(*Lamda == NULL) {
                printf("Memory reallocation failed!\n");
                *Lamda = temp;
            }
        }
        
        // Update Lamda_minus_V
        Lamda_minus_V = get_the_difference(*Lamda, *Lamda_size, *V, *V_size, &Lamda_minus_V_size); 
    }
    free(Lamda_minus_V);
}


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
 * @return An array of indices representing the sampled point indexes or NULL if an error occurs
 */
int* sample_points(Graph graph, int num_sample_points) {
    int *sample_point_indexes = (int *)calloc(num_sample_points, sizeof(int));

    if( num_sample_points > graph.num_points) {
        printf("Error: num_sample_points > graph.num_points\n");
        return NULL;
    }

    int indices[graph.num_points];
    for (int i = 0; i < graph.num_points; i++) {
        indices[i] = i;
    }

    srand(time(NULL));
    
    for (int i = graph.num_points - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        // Swap indices[i] and indices[j]
        int temp = indices[i];
        indices[i] = indices[j];
        indices[j] = temp;
    }

    for (int i = 0; i < num_sample_points; i++) {
        sample_point_indexes[i] = indices[i];
    }

    return sample_point_indexes;
}


/**
 * Checks an array for duplicate elements
 *
 * @param array The array to check
 * @param size The size of the array
 * @return void
 *
 * This function checks the given array for any duplicate elements. If a
 * duplicate is found, it prints a message indicating the duplicate value.
 */
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




/**
 * Performs Vamana indexing on a given graph to optimize its structure.
 *
 * get the medoid of the graph
 * traverse the graph in a random way without repetitions
 * for every random point in the graph
 * run greedysearch(medoid, random point, k=1, L) to get the visited list V
 * run robustPrune(random Point, Visited list, a, R)
 * for every neighbor of random point p'
 * if | neighbor(p') U random point | > R then
 * run robustPrune(p', Neighbors of p' U random point, a, R)
 * else add to neighbors of p' the random point
 * end for
 *
 * @param graph The graph to be indexed.
 * @param L The parameter controlling the size of the visited list in the greedy search.
 * @param a The pruning parameter that influences the robustness of pruning.
 * @param R The maximum number of neighbors allowed for a point after pruning.
 */
Graph* filtered_vamana_indexing(DatasetInfo* dataset, int L, float a, int R,filterInfo *filterinfo) {

    printf("Starting Vamana Indexing\n");

    Graph graph = initialise_graph(dataset, R);
    printf("Graph initialised\n");

    // Calculating the medoid
    // ======= CHANGE THE PERCENTAGE OF THE SAMPLES HERE ======== //
    int percentage = 20;
    // ========================================================== //

    int num_sample_points = graph.num_points / percentage;
    int *sample_point_indexes = sample_points(graph, num_sample_points);
    int medoid_index = calculate_medoid(&graph, sample_point_indexes, graph.num_points);
    FilteredMethoidList* filteredMedoids = get_filtered_medoids(&graph, &num_sample_points, filterinfo);
    printf("Medoid index: %d\n", medoid_index);


    // traverse the graph in a random way without repetitions
    bool *shuffled_point_indexes = (bool *)calloc(graph.num_points, sizeof(bool));
    int s_index;
    int i=0;
    int* V = NULL;
    int *lamda = NULL;
    int lamda_size = 0;
    int V_size = 0;

    int starting_points[128];
    for(i = 0; i < 128; i++) {
        starting_points[i] = rand() % graph.num_points;
    }
    
    while(i < graph.num_points) {
        s_index = rand() % graph.num_points;
        if (!shuffled_point_indexes[s_index]) {
            /**
             * Finds the medoid from the filtered medoids that is in the same category as the given datapoint.
             * If the datapoint does not have a category, it returns the global medoid.
             *
             * @param filteredMedoids The list of filtered medoids.
             * @param datapoint The datapoint to find the medoid for.
             * @param medoid_index The global medoid index.
             * @return The index of the medoid in the same category as the datapoint, or the global medoid index if no category.
             */
            shuffled_point_indexes[s_index] = true;
            // printf("Indexing point %d\n", s_index);

            // Set as the starting index the node that is at the starting_points[category]
           int F_Fx_si = find_medoid_for_point(filteredMedoids, &graph.points[s_index], medoid_index);
//            int *S_Fx_si = malloc(sizeof(int));
//            S_Fx_si[0] = graph->points[F_Fx_si];
            int sp_size = 1;

            // =============== GREEDY SEARCH ================ //
            filtered_greedy_search(&graph, graph.points[s_index].coordinates, &F_Fx_si, sp_size, &V, &V_size, &lamda, &lamda_size, L, graph.points[s_index].category);

            // =============== ROBUST PRUNE ================ //
            filtered_Robust_prune(&graph, s_index, V, V_size, a, R);

            int *new_V = NULL;
            int new_V_size = 0;

            // for every neighbor of random point p'
            for(int j = 0; j < graph.points[s_index].edge_count; j++) {

                // if | neighbor(p') U random point | > R then
                Point *P_PRIME = &graph.points[graph.points[s_index].edges[j]];
                
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
                    filtered_Robust_prune(&graph, P_PRIME->index, new_V, new_V_size, a, R);
                    
                } else{
                    // else add to neighbors of p' the random point
                    addEdge(P_PRIME, s_index);
                }
            }
            // Traversed to another point of the graph
            i++;

            V = NULL;
            V_size = 0;
            lamda = NULL;
            lamda_size = 0;
        }

    }
    free(V);
    free(lamda);            
    free(shuffled_point_indexes);
    return &graph;
}


// ==================== Stitched Vamana Indexing ==================== //


/**
 * Creates a graph with the given dataset and maximum number of edges per point.
 * Every point will have a random number of edges between 0 and max_edges, with no duplicates and no self-loops.
 * The graph will be initialized with random coordinates between 0 and 1.
 * 
 * @param dataset The dataset to use for initializing the graph
 * @param base_num_dimensions The number of dimensions of the graph
 * @param max_edges The maximum number of edges for each point in the graph
 * @return A pointer to the newly created graph
 */
Graph create_random_graph(DatasetInfo dataset, int base_num_dimensions, int max_edges) {
    int base_num_points = dataset.num_vectors;
    Graph graph;
    graph.points = (Point*)malloc(base_num_points * sizeof(Point));

    if (!graph.points) {
        printf("Memory allocation failed!\n");
        exit(1);
    }

    graph.num_points = base_num_points;
    graph.num_dimensions = base_num_dimensions;

    // Initialize points with vector values
    for (int i = 0; i < base_num_points; i++) {
        graph.points[i].index = dataset.datapoints[i].point_index;
        graph.points[i].category = dataset.datapoints[i].category;

        // Allocate memory for coordinates
        graph.points[i].coordinates = (float*)malloc(base_num_dimensions * sizeof(float));
        if (!graph.points[i].coordinates) {
            printf("Memory allocation failed!\n");
            exit(1);
        }

        // Allocate memory for edges
        graph.points[i].edges = (int*)malloc(max_edges * sizeof(int));
        if (!graph.points[i].edges) {
            printf("Memory allocation failed!\n");
            exit(1);
        }

        // Add random edges to every point until they have max_edges, no duplicates and no self-loops
        graph.points[i].edge_count = 0; // No edges initially

        int temp_max_edges = max_edges;
        if(graph.num_points <= max_edges){
            temp_max_edges = graph.num_points - 1;
        }

        while(graph.points[i].edge_count < temp_max_edges) {
            if(graph.num_points == 1) break;
            
            int edge = rand() % base_num_points;
            if (edge != i && !arrayContains(graph.points[i].edges, graph.points[i].edge_count, edge)) {
                graph.points[i].edges[graph.points[i].edge_count] = edge;
                graph.points[i].edge_count++;
            }
        }

        // Initialize coordinates
        for(int j = 0; j < base_num_dimensions; j++) {
            graph.points[i].coordinates[j] = dataset.datapoints[i].vectors[j];
        }

    }

    return graph;
}


/**
 * Robust Pruning Algorithm
 * 
 * Add to V all the neighbors of p and remove p from V if it exists
 * Set reset edges of p
 * While V is not empty
 * set as P* the closest point in V to p
 * Add P* to the neighbors of p
 * if the number of neighbors of p is equal to R then break
 * for every point p' in V
 * if a * d(p*, p') <= d(p, p') then remove p' from V
 *
 * @param graph The graph containing the point p
 * @param p_index The index of the point p
 * @param V The set of points that are candidates to be neighbors of p
 * @param V_size The size of V
 * @param a The pruning parameter
 * @param R The maximum number of neighbors of p
 */
void robustPrune(Graph *graph, int p_index, int *V, int V_size, float a, int R) {

    // Add to V all the neighbors of p and remove p from v if it exists
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
    
        if(V[min_index] != p_index){
            p_star = &graph->points[V[min_index]];
            // printf("Adding edge from %d -> %d\n", p_index, V[min_index]);
            addEdge(&graph->points[p_index], V[min_index]);
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
        }else{
            V[min_index] = V[V_size - 1];
            V_size--;
        }

    }
    
}

/**
 * Performs a greedy search on the graph to find the L closest points to a given query point Xq.
 *
 * The algorithm works by selecting the point in Lamda that is closest to Xq, adding it to V, and then adding the neighbors of that point to Lamda.
 * It continues until all points in Lamda are visited or the limit L is reached.
 *
 * @param graph The graph to search.
 * @param Xq The query point.
 * @param start_index The index of the starting point of the search.
 * @param V The list of visited points.
 * @param V_size The size of V.
 * @param Lamda The list of points to consider.
 * @param Lamda_size The size of Lamda.
 * @param L The limit on the number of points to consider.
 */
void greedy_search(Graph *graph, float *Xq, int start_index, int **V, int *V_size, int **Lamda, int *Lamda_size, int L) {

    // Allocate initial space for V and Lamda
    *V = (int *)malloc(sizeof(int));
    *V_size = 0;
    *Lamda = (int *)malloc(sizeof(int));
    *Lamda_size = 0;

    // Initialize Lamda with the start point
    add_to_dynamic_array(Lamda, Lamda_size, start_index);

    // Create an int array that contains points that are in Lamda but not in V
    int *Lamda_minus_V;
    int Lamda_minus_V_size = 0;
    Lamda_minus_V = get_the_difference(*Lamda, *Lamda_size, *V, *V_size, &Lamda_minus_V_size);


    // for every point in Lamda_minus_V
    // calculate the distance from xq to that point
    // the one that is closest is p*
    // Take P* 's   Neighbors and add them to Lamda
    // Take P* and add it to V

    // Main loop: Continue until Lamda_minus_V is empty
    while (Lamda_minus_V_size > 0) {
        // Select the point p* from Lamda_minus_V that is closest to the query point Xq
        int closest_index = -1;
        float min_distance = FLT_MAX;
        for (int i = 0; i < Lamda_minus_V_size; i++) {
            int current_index = Lamda_minus_V[i];
            // printf("Current index: %d\n", current_index);
            float distance = squared_euclidean_distance(graph->points[current_index].coordinates, Xq, graph->num_dimensions);
            if (distance < min_distance) {
                min_distance = distance;
                closest_index = current_index;
            }
        }
        Point *p_star = &graph->points[closest_index];

        // Add p* position to V
        // add_to_dynamic_array(V, V_size, p_star->index);
        add_to_dynamic_array(V, V_size, closest_index);


        // Add the neighbors of p* to Lamda
        for (int i = 0; i < p_star->edge_count; i++) {
            int neighbor_index = p_star->edges[i];
            add_to_dynamic_array(Lamda, Lamda_size, neighbor_index);
        }

        // Check if |Lamda| > L
        if (*Lamda_size > L) {
            // Update Lamda to contain the L closest points to Xq
            sort_array(graph, *Lamda, *Lamda_size, Xq);
            int *temp = *Lamda;
            *Lamda_size = L;
            *Lamda = (int *)realloc(*Lamda, L * sizeof(int));
            if(*Lamda == NULL) {
                printf("Memory reallocation failed!\n");
                *Lamda = temp;
            }
        }
        
        // Update Lamda_minus_V
        Lamda_minus_V = get_the_difference(*Lamda, *Lamda_size, *V, *V_size, &Lamda_minus_V_size); 

    }
    free(Lamda_minus_V);
}

/**
 * Performs Vamana indexing on a given graph to optimize its structure.
 *
 * get the medoid of the graph
 * traverse the graph in a random way without repetitions
 * for every random point in the graph
 * run greedysearch(medoid, random point, k=1, L) to get the visited list V
 * run robustPrune(random Point, Visited list, a, R)
 * for every neighbor of random point p'
 * if | neighbor(p') U random point | > R then
 * run robustPrune(p', Neighbors of p' U random point, a, R)
 * else add to neighbors of p' the random point
 * end for
 *
 * @param dataset The dataset to be indexed.
 * @param L The parameter controlling the size of the visited list in the greedy search.
 * @param a The pruning parameter that influences the robustness of pruning.
 * @param R The maximum number of neighbors allowed for a point after pruning.
 */
Graph vamana_indexing(DatasetInfo dataset, int L, float a, int R) {

    // Initialize a random graph from the dataset
    Graph graph = create_random_graph(dataset, 100, R);
    if(graph.num_points < R){
        return graph;
    }


    // Calculating the medoid
    // ======= CHANGE THE PERCENTAGE OF THE SAMPLES HERE ======== //
    int percentage = 20;
    // ========================================================== //

    int num_sample_points = graph.num_points / percentage;
    int *sample_point_indexes = sample_points(graph, num_sample_points);
    if(sample_point_indexes == NULL) {
        printf("Error: sample_points returned NULL\n");
        exit(1);
    }

    int medoid_index = calculate_medoid(&graph, sample_point_indexes, num_sample_points);       // Position in the graph
    // int temp_medoid_index = graph.points[medoid_index].index;                                   // Actual index

    // traverse the graph in a random way without repetitions
    bool *shuffled_point_indexes = (bool *)calloc(graph.num_points, sizeof(bool));
    int s_index;
    int i=0;
    int* V = NULL;
    int *lamda = NULL;
    int lamda_size = 0;
    int V_size = 0;
    
    while(i < graph.num_points) {
        s_index = rand() % graph.num_points;
        if (!shuffled_point_indexes[s_index]) {
            shuffled_point_indexes[s_index] = true;
            // printf("Indexing graph %d point %d | %d / %d\n", graph.points[s_index].category, s_index, i, graph.num_points);

            // =============== GREEDY SEARCH ================ //
            greedy_search(&graph, graph.points[s_index].coordinates, medoid_index, &V, &V_size, &lamda, &lamda_size, L);
            // =============== ROBUST PRUNE ================ //

            robustPrune(&graph, s_index, V, V_size, a, R);

            int *new_V = NULL;
            int new_V_size = 0;
            
            // for every neighbor of random point p'
            for(int j = 0; j < graph.points[s_index].edge_count; j++) {

                // if | neighbor(p') U random point | > R then
                Point *P_PRIME = &graph.points[graph.points[s_index].edges[j]];
                int P_PRIME_graph_position = graph.points[s_index].edges[j];
                
                if(P_PRIME->edge_count + 1 > R) {
                    // run robustPrune(p', Neighbors of p' U random point, a, R)
                    // create a new visited list that contains only the neighbors of p' and the random point
                    new_V = NULL;
                    new_V = (int *) malloc((P_PRIME->edge_count + 1) * sizeof(int));
                    new_V_size = 0;
                    for (int k = 0; k < P_PRIME->edge_count; k++) {
                        if (!arrayContains(new_V, new_V_size, P_PRIME->edges[k])) {
                            add_to_dynamic_array(&new_V, &new_V_size, P_PRIME->edges[k]);
                        }
                    }
                    add_to_dynamic_array(&new_V, &new_V_size, s_index);

                    // run robustPrune(p', new_V, a, R)
                    robustPrune(&graph, P_PRIME_graph_position, new_V, new_V_size, a, R);
                    
                } else{
                    // else add to neighbors of p' the random point
                    addEdge(P_PRIME, s_index);
                }
            }
            // Traversed to another point of the graph
            i++;

            V = NULL;
            V_size = 0;
            lamda = NULL;
            lamda_size = 0;
        }

    }
    free(V);
    free(lamda);            
    free(shuffled_point_indexes);
    return graph;
}

FilteredMethoidList * get_filtered_medoids(Graph *graph, int *t, filterInfo *filterInfo)
{
    Point **groupedData = malloc(filterInfo->num_filters * sizeof(Point*));
    for (int i = 0; i < filterInfo->num_filters; i++) {
        groupedData[i] = malloc(filterInfo->filtersPoints[i].count * sizeof(Point));
    }

    int *current_position = (int *)calloc(filterInfo->num_filters, sizeof(int));

    for (int i = 0; i < graph->num_points; i++) {
        int category = graph->points[i].category;
        for (int j = 0; j < filterInfo->num_filters; j++) {
            if (filterInfo->filtersPoints[j].filter_index == category) {
                groupedData[j][current_position[j]++] = graph->points[i];
                break;
            }
        }
    }

    free(current_position);


    FilteredMethoidList* filteredMedoids = malloc(sizeof(FilteredMethoidList));

    filteredMedoids->metoids = findClosestDataPoints(groupedData,filterInfo,filterInfo->num_filters,*t);
    filteredMedoids->size=filterInfo->num_filters;


    for (int i = 0; i < filterInfo->num_filters; i++) {
        free(groupedData[i]);
    }
    free(groupedData);

    return filteredMedoids;
}

void generate_random_permutation(int *perm, int n) {
    for (int i = 0; i < n; i++) {
        perm[i] = i;
    }
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = perm[i];
        perm[i] = perm[j];
        perm[j] = temp;
    }
}



FilteredMedoid* findClosestDataPoints(Point **groupedData, filterInfo *filterInfo, int numCategories, int t) {
    FilteredMedoid*closestPoints = malloc(numCategories * sizeof(FilteredMedoid));
    srand(time(NULL)); // Seed for random permutation

    for (int cat = 0; cat < numCategories; cat++) {
        int groupSize = filterInfo->filtersPoints[cat].count;
        if (groupSize == 0) continue;

        // If the category has only one element, directly assign its index as the medoid
        if (groupSize == 1) {
            closestPoints[cat].index = groupedData[cat][0].index;
            closestPoints[cat].category = cat;
            continue;
        }

        // Adjust the permutation size to be the minimum of `groupSize` and `t`
        int permSize = groupSize < t ? groupSize : t;

        int *randomPermutation = malloc(permSize * sizeof(int));
        generate_random_permutation(randomPermutation, groupSize);

        double minDistance = DBL_MAX;
        int closestIndex = -1;

        for (int i = 0; i < permSize; i++) {
            int permIndex = randomPermutation[i];
            float *permVector = groupedData[cat][permIndex].coordinates;

            double totalDistance = 0.0;
            for (int j = 0; j < groupSize; j++) {
                if (j == permIndex) continue;
                float *currentVector = groupedData[cat][j].coordinates;
                totalDistance += squared_euclidean_distance(permVector, currentVector, 100);
            }

            if (totalDistance < minDistance) {
                minDistance = totalDistance;
                closestIndex = permIndex;
            }
        }

        // Store the index and category of the closest DataPoint for this category
        closestPoints[cat].index = groupedData[cat][closestIndex].index;
        closestPoints[cat].category = cat;
        free(randomPermutation);
    }

    return closestPoints;
}

int find_medoid_for_point(FilteredMethoidList* filteredMedoids, Point* datapoint, int medoid_index)
{
    if (datapoint->category == -1) {
        return medoid_index;
    }

    for (int i = 0; i < filteredMedoids->size; i++) {
        if (filteredMedoids->metoids[i].category == datapoint->category) {
            return filteredMedoids->metoids[i].index;
        }
    }

    return medoid_index;
}

/**
 * Stitched Vamana Indexing
 *
 * This algorithm is a combination of multiple Vamana indexing.
 * It takes a dataset and divides it into multiple sub datasets, each one
 * corresponding to a filter.
 * For each sub dataset, it runs the Vamana algorithm, and then it takes the
 * visited list of each filter and runs the filteredRobustPrune algorithm.
 *
 * @param dataset The dataset to be indexed.
 * @param L_small The parameter controlling the size of the visited list in the greedy search.
 * @param a The pruning parameter that influences the robustness of pruning.
 * @param R_small The maximum number of neighbors allowed for a point after pruning.
 * @param R_stitched The maximum number of neighbors allowed for a point after pruning.
 * @return The graph with the indexed points.
 */
Graph* stitched_vamana_indexing(DatasetInfo* dataset, int L_small, float a, int R_small, int R_stitched) {
    
    printf("Stitched Vamana Indexing: ");
    fflush(stdout);
    // Graph stitched_graph = initialise_graph(dataset, R_small);


    Graph *filter_graph = (Graph *)malloc(dataset->filterInfo.num_filters * sizeof(Graph));
    if(filter_graph == NULL){
        printf("Memory allocation failed filter_graph!\n");
        exit(1);
    }

    DatasetInfo *filter_dataset = (DatasetInfo *)malloc(dataset->filterInfo.num_filters * sizeof(DatasetInfo));
    if(filter_dataset == NULL){
        printf("Memory allocation failed filter_set!\n");
        exit(1);
    }

    // Find all the indexed with that filter
    for( int i = 0; i < dataset->filterInfo.num_filters; i++){

        // Find all the points fo each filter
        filter_dataset[i].num_vectors = dataset->filterInfo.filtersPoints[i].count;
        filter_dataset[i].datapoints = (DataPoint *)calloc(filter_dataset[i].num_vectors, sizeof(DataPoint));
        if(filter_dataset[i].datapoints == NULL){
            printf("Memory allocation failed filter_dataset.datapoint!\n");
            exit(1);
        }

        for(int j = 0; j < filter_dataset[i].num_vectors; j++){
            int filter_point_index = dataset->filterInfo.filtersPoints[i].point_indexes[j];
            filter_dataset[i].datapoints[j].category = dataset->datapoints[filter_point_index].category;
            filter_dataset[i].datapoints[j].point_index = dataset->datapoints[dataset->filterInfo.filtersPoints[i].point_indexes[j]].point_index;
            filter_dataset[i].datapoints[j].timestamp = dataset->datapoints[dataset->filterInfo.filtersPoints[i].point_indexes[j]].timestamp;
            for(int k = 0; k < 100; k++){
                filter_dataset[i].datapoints[j].vectors[k] = dataset->datapoints[dataset->filterInfo.filtersPoints[i].point_indexes[j]].vectors[k];
            }

        }

        // For each filter, run the vamana algorithm
        filter_graph[i] = vamana_indexing(filter_dataset[i], L_small, a, R_small);

    }

    // Stitch the graphs

    // for(int i = 0; i < stitched_graph.num_points; i++){
    //     printf("Running filteredRobustPrune for point: %d\n", i);
    //     // Run the filteredRobustPrune algorithm
    //     // void filtered_Robust_prune(Graph *graph, int p_index, int *V, int V_size, float a, int R)
    //     int u_index = stitched_graph.points[i].index;
    //     int *V = stitched_graph.points[u_index].edges;
    //     int V_size = stitched_graph.points[u_index].edge_count;
    //     filtered_Robust_prune(filter_graph, u_index, V, V_size, a, R_stitched);
    //     printf("Finished running filteredRobustPrune for point: %d\n", i);
    // }

    for(int i = 0; i < dataset->filterInfo.num_filters; i++){
        free(filter_dataset[i].datapoints);
    }
    free(filter_dataset);
    printf("Done.\n");
    fflush(stdout);
    // Return array of graphs
    return filter_graph;
}


// ============= Functions that can be used if we have multiple filters ============= //

// /**
//  * Computes the intersection of two sets.
//  *
//  * @param set1 The first set
//  * @param set1_size The size of the first set
//  * @param set2 The second set
//  * @param set2_size The size of the second set
//  * @param intersection The computed intersection
//  * @return The size of the intersection
//  */
// int compute_intersection(int *set1, int set1_size, int *set2, int set2_size, int *intersection) {
//     int k = 0; // Counter for the intersection
//     for (int i = 0; i < set1_size; i++) {
//         if (arrayContains(set2, set2_size, set1[i])) {
//             intersection[k++] = set1[i];
//         }
//     }
//     return k; // Return the size of the intersection
// }


// /**
//  * Checks if any element of the intersection of two sets is not present in the third set
//  * @param intersection The intersection of two sets
//  * @param intersection_size The size of the intersection
//  * @param set3 The third set
//  * @param set3_size The size of set3
//  * @return true if any element of intersection is not in set3, false if all elements of intersection are in set3
//  */
// bool is_not_subset(int *intersection, int intersection_size, int *set3, int set3_size) {
//     for (int i = 0; i < intersection_size; i++) {
//         if (!arrayContains( set3, set3_size, intersection[i])) {
//             return true; // Found an element in intersection that is not in set3
//         }
//     }
//     return false; // All elements of intersection are in set3
// }

// ============= END Functions that can be used if we have multiple filters ============= //

