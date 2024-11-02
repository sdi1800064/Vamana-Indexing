#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>

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

    // Create the graph
    printf("Creating random graph with %d vectors and %d dimensions\n", base_num_points, base_num_dimensions);

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
        graph->points[i].edges = (Edge*)malloc(max_edges * sizeof(Edge));
        if (!graph->points[i].edges) {
            printf("Memory allocation failed!\n");
            exit(1);
        }

        // Initialize coordinates
        for(int j = 0; j < base_num_dimensions; j++) {
            graph->points[i].coordinates[j] = base_vectors[i][j];
        }
        graph->points[i].edge_count = 0; // No edges initially
        // if(i % 100 == 0){
        //     printf("Point %d Initialized\n ", i);
        // }
    }
    printf("All graph-points initiallized\n");

    // Add random edges
    add_random_edges(graph, max_edges);
    

    printf("Freeing base vectors..\n");
    for (int i = 0; i < base_num_points; i++) {
        free(base_vectors[i]);
    }
    free(base_vectors);
    printf("Graph created\n");

    return graph;
}

/**
 * Function to randomly connect points with edges
 * 
 * @param graph The graph to add edges to
 * @param max_edges The maximum number of edges a point can have
 * @param base_num_points The number of vectors in the graph
 */
void add_random_edges(Graph* graph, int max_edges) {
    printf("Adding random edges..\n");
    for (int i = 0; i < graph->num_points; i++) {
        for (int j = 0; j < max_edges; j++) {
            int random_index;
            do {
                // Generate a random index of a point to connect to
                random_index = rand() % graph->num_points;
            } while (random_index == i); // Avoid self-loop
            // printf("Connecting point %d to point %d\n", i, random_index);

            // Add the edge if there is space
            if ((graph->points[i].edge_count) < max_edges) {
                graph->points[i].edges[graph->points[i].edge_count].to = random_index;
                graph->points[i].edge_count++;
            }
            // printf("Point %d now has %d edges\n", i, graph->points[i].edge_count);
        }
        // printf("Point %d has edges\n", i);
    }
    printf("Edges added\n");
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
            fprintf(outputfd, "%d ", p.edges[j].to);
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
            fprintf(outputfd, "%d ", p.edges[j].to);
        }
        fprintf(outputfd, "\n");
    }
}

// Add an edge to a point
void addEdge(Point *point, int toIndex) {
    point->edge_count++;
    point->edges[point->edge_count].to = toIndex;

}

// Check if point `toIndex` is already in point's edges
int edgeExists(Point *point, int toIndex) {
    for (int i = 0; i < point->edge_count; i++) {
        if (point->edges[i].to == toIndex) return 1;
    }
    return 0;
}

// The robustPrune function
void robustPrune(Graph *graph, int p_index, int *V, int V_size, float a, int R) {
    printf("Starting Robust Prune\n");
    Point *p = &graph->points[p_index];

    // Doing Robust Pruning for p_index
    printf("Doing Robust Pruning for point %d with V_size = %d\n", p_index, V_size);
    printf("V_size before adding %d neighbors = %d\n", p->edge_count, V_size);
    
    // Add edges of p to V if not already included
    // For every edge if p 
    for (int i = 0; i < p->edge_count; i++) {
        // Check if neighbor is already in V
        int neighborIndex = p->edges[i].to;
        int found = 0;
        for (int j = 0; j < V_size; j++) {
            if (V[j] == neighborIndex) {
                // printf("Neighbor %d already in V\n", neighborIndex);
                found = 1;
                break;
            }
        }
        // if not found, add to V
        if (!found) {
            printf("resize V from %d to %d\n", V_size, V_size + 1);
            V = realloc(V, (V_size + 1) * sizeof(int));
            if (!V) {
                printf("Memory allocation of V failed!\n");
                exit(1);
            }
            V[V_size] = neighborIndex;
            V_size++;
        }
    }

    // printf("V_size after adding neighbors = %d\n", V_size);

    // Clear edges of p
    p->edges = calloc(R, sizeof(Edge));
    p->edge_count = 0;

    // printf("Resetted the edges of point %d\n", p_index);    

    while (V_size > 0) {

        // Find p* - the point in V closest to p
        float minDistance = FLT_MAX;        // Minimum distance
        int p_star_index = -1;              // Index of p*
        int minIndex = -1;                  // Index of p* in V

        // Iterate over V to add p* to the edges of p
        for (int i = 0; i < V_size; i++) {
            // Skip if p is already in V
            if (V[i] == p_index) {
                continue;
            }
            int candidateIndex = V[i];
            // Find the point in V closest to p
            float dist = squared_euclidean_distance(p->coordinates, graph->points[candidateIndex].coordinates, graph->num_dimensions);
            if (dist < minDistance) {
                minDistance = dist;
                p_star_index = candidateIndex;
                minIndex = i;
            }
        }

        // Add p* to the edges of p
        addEdge(p, p_star_index);
        // printf("Added edge from point %d to point %d\n", p_index, p_star_index);

        // If the number of edges reaches R, break
        if (p->edge_count == R) break;

        // Remove p* from V by swapping with the last element
        V[minIndex] = V[--V_size];
        // printf("Removed point %d from V\n", p_star_index);

        // Prune V based on the distance criterion
        for (int i = 0; i < V_size; ) {
            // p'
            int p_prime_index = V[i];
            // d(p, p')
            float d_pp = squared_euclidean_distance(p->coordinates, graph->points[p_prime_index].coordinates, graph->num_dimensions);
            // d(p*, p')
            float d_pstar_pprime = squared_euclidean_distance(graph->points[p_star_index].coordinates, graph->points[p_prime_index].coordinates, graph->num_dimensions);
            if (a * d_pstar_pprime <= d_pp) {
                // Remove p' from V by swapping with the last element
                V[i] = V[--V_size];
            } else {
                i++;
            }
        }
    }
    printf("Done Robust Pruning\n");
}


// Comparison function for qsort to sort indices in l_temp based on distance to Xq
int compare(const void *a, const void *b, void *param) {
    // Extract graph, query point, and dimensions from parameters
    Graph *graph = ((void **)param)[0];
    float *Xq = ((void **)param)[1];
    int num_dimensions = *((int *)((void **)param)[2]);

    int index_a = *(int *)a;
    int index_b = *(int *)b;

    // Calculate distances from points to the query point Xq
    float distance_a = squared_euclidean_distance(graph->points[index_a].coordinates, Xq, num_dimensions);
    float distance_b = squared_euclidean_distance(graph->points[index_b].coordinates, Xq, num_dimensions);

    // Sort in ascending order based on distance
    if (distance_a < distance_b) return -1;
    if (distance_a > distance_b) return 1;
    return 0;
}

// Function to sort l_temp based on distances to query point Xq
void sort_array(Graph *graph, int *l_temp, int l_temp_size, float *Xq) {
    // Create parameter array to pass multiple values to qsort_r
    void *params[3] = {graph, Xq, &graph->num_dimensions};

    // Sort the array in place using qsort with custom comparator
    qsort_r(l_temp, l_temp_size, sizeof(int), compare, params);
}


// Helper function to check if Lamda \ V is not empty
int exists_in_difference(int *Lamda, int Lamda_size, int *V, int V_size) {
    for (int i = 0; i < Lamda_size; i++) {
        int in_V = 0;
        for (int j = 0; j < V_size; j++) {
            if (Lamda[i] == V[j]) {
                in_V = 1;
                break;
            }
        }
        if (!in_V) return 1; // Found an element in Lamda that is not in V
    }
    return 0;
}

// Function to add an element to a dynamically resized array
void add_to_dynamic_array(int **array, int *size, int element) {
    *array = realloc(*array, (*size + 1) * sizeof(int));
    if (*array == NULL) {
        perror("Failed to allocate memory");
        exit(EXIT_FAILURE);
    }
    (*array)[*size] = element;
    (*size)++;
}

// Helper function to check if an element is in an array
int is_in_array(int *array, int size, int element) {
    for (int i = 0; i < size; i++) {
        if (array[i] == element) {
            return 1;
        }
    }
    return 0;
}

// Greedy search function
void greedy_search(Graph *graph, float *Xq, int start_index, int **V, int *V_size, int **Lamda, int *Lamda_size, int L, int k) {
    printf("Starting Greedy Search\n");
    bool visited[graph->num_points]; // Array to mark visited points
    for (int i = 0; i < graph->num_points; i++) visited[i] = 0; // Initialize visited array

    int current = start_index; // Start from the given point
    visited[current] = 1; // Mark starting point as visited
    add_to_dynamic_array(V, V_size, current); // Add starting point to V

    while (exists_in_difference(*Lamda, *Lamda_size, *V, *V_size)) {
        Point *current_point = &graph->points[current];
        float min_distance = FLT_MAX;
        int closest_neighbor = -1;

        // Iterate over neighbors of the current point
        for (int i = 0; i < current_point->edge_count; i++) {
            int neighbor_index = current_point->edges[i].to;
            Point *neighbor = &graph->points[neighbor_index];

            // Add this neighbor to Lamda if not already visited
            if (!visited[neighbor_index] && !is_in_array(*Lamda, *Lamda_size, neighbor_index)) {
                add_to_dynamic_array(Lamda, Lamda_size, neighbor_index);
            }

            // Calculate distance to Xq
            float distance = squared_euclidean_distance(neighbor->coordinates, Xq, graph->num_dimensions);

            // Update the closest neighbor if this one is closer
            if (distance < min_distance && !visited[neighbor_index]) {
                min_distance = distance;
                closest_neighbor = neighbor_index;
            }
        }

        // Check if the closest neighbor brings us closer to Xq
        float current_distance = squared_euclidean_distance(current_point->coordinates, Xq, graph->num_dimensions);
        if (closest_neighbor == -1 || min_distance >= current_distance) {
            // No unvisited neighbors closer than the current point, stop search
            break;
        }

        // Move to the closest neighbor
        current = closest_neighbor;
        visited[current] = 1; // Mark as visited
        printf("Visited point: %d\n", current);
        add_to_dynamic_array(V, V_size, current); // Add to V

        if (*Lamda_size > L) {
            sort_array(graph, *Lamda, *Lamda_size, Xq);
            *Lamda_size = L;
            *(Lamda) = realloc(*Lamda, *Lamda_size * sizeof(int));
        }
    }
    printf("Finished Greedy Search\n");
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
    int current_index = 0;
    int i = 0;
    while (i < num_sample_points) {
        // Randomly select an index from the graph
        current_index = rand() % max;
        // Ensure no duplicate indices are selected
        if (!arrayContains(sample_point_indexes, i, current_index)) {
            sample_point_indexes[i] = current_index;
            i++;
        }
    }
    return sample_point_indexes;
}

void vamana_indexing(Graph *graph, int k, int L, float a, int R, FILE *outputfd) {

    printf("Starting Vamana Indexing\n");

    printf("k = %d\n L = %d\n a = %f\n R = %d\n", k, L, a, R);
    
    int max = graph->num_points;            // max - Number of points in graph
    int num_sample_points = graph->num_points / 10;        // num_sample_points - Number of sample points     

    int *sample_point_indexes = sample_points(max, num_sample_points);          // sample_point_indexes - Indexes of sampled points
    int medoid_index = calculate_medoid(graph, sample_point_indexes, num_sample_points);            // medoid_index - Index of medoid
    printf("Medoid index: %d\n", medoid_index);

    // Initiallizations
    bool *shuffled_point_indexes = (bool*)calloc(max, sizeof(bool));            // shuffled_point_indexes - Indexes of points in a random order
    int current_index = -1;        // Temporary index for random selection
    int *V;                     // V - Visited List
    int V_size = 0;             // V_size - Size of V
    int *l;                     // L~ - k nearest neighbors list

    // Allocate memory for V
    V = (int*)malloc((V_size + 1) * sizeof(int));
    if (!V) {
        printf("Memory allocation of V failed!\n");
        exit(1);
    }

    // L~ - L nearest neighbors to Xq list
    int lamda_size = L;
    l = (int*)calloc(lamda_size, sizeof(int));
    if (!l) {
        printf("Memory allocation of l failed!\n");
        exit(1);
    }

    printf("Memory allocation of shuffled_point_indexes and V and l completed successfully\n");

    int i = 0;

    // Going through all points of the graph in a random order
    while (i < max) {
        fprintf(outputfd, "\nProcessing %d\n", i);
        printf("Points remaining %d\n", max - i);
        // Randomly select an index from the graph
        current_index = rand() % max;
        // Ensure no duplicate indices are selected
        if (shuffled_point_indexes[current_index] == false) {

            // Mark the current index as visited
            shuffled_point_indexes[current_index] = true;

            Point *q_point = &graph->points[current_index];
            
            // Reset the visited list
            V_size = 0;
            V = NULL;
            V = (int*)malloc(sizeof(int));
            if (!V) {
                printf("Vamana_Indexing: Memory allocation of V failed!\n");
                exit(1);
            }

            printf("V got resetted\n");
            
            // Perform Greedy Search GreedySearch(s, Xsi, 1, L) for the current point (current_index)
            printf("Performing Greedy Search for point %d\n", current_index);
            fprintf(outputfd, "Performing Greedy Search for point %d\n", current_index);
            greedy_search(graph, q_point->coordinates, medoid_index, &V, &V_size, &l, &lamda_size, L, 1);
            // GreedySearch(graph, &V, &V_size, l, L, graph->points[current_index].coordinates, 1, medoid_index);
        
            printf("Performing Robust Prune for point %d\n", current_index);
            fprintf(outputfd, "Performing Robust Prune for point %d\n", current_index);
            robustPrune(graph, current_index, V, V_size, a, R);
            printf("RobustPrune completed\n");


            // Nout(j) <- Nout(j)U{s(i)}
            printf("Updating Nout(j) for point %d\n", current_index);

            for(int j = 0; j < q_point->edge_count; j++) {
                Point *point_j = &graph->points[q_point->edges[j].to];
                
                // Check if |Nout(j)U{s(i)}| > R
                //          RobustPrune(j,Nout(j)U{s(i)},a,R)
                //      else
                //          update Nout(j) <- Nout(j)U{s(i)}
                if(point_j->edge_count + 1 > R) {

                    // Create Temporary array to hold all the points in Nout(j) + s
                    int *Nout_plus_s_j = (int*)malloc((point_j->edge_count + 1) * sizeof(int));
                    for (int k = 0; k < point_j->edge_count; k++) {
                        Nout_plus_s_j[k] = point_j->edges[k].to;
                    }
                    Nout_plus_s_j[point_j->edge_count] = current_index;

                    // RobustPrune(j,Nout(j)U{s(i)},a,R)
                    robustPrune(graph, point_j->index, Nout_plus_s_j, point_j->edge_count + 1, a, R);
                    free(Nout_plus_s_j);
                }
                else {
                    if (point_j->edge_count < R){
                        // Update Nout(j) <- Nout(j)U{s(i)}
                        printf("P* %d had %d edges\n", point_j->index, point_j->edge_count);
                        printf("Adding edge from %d to %d\n", current_index, point_j->index);
                        addEdge(point_j, current_index);
                        printf("Now P* has %d edges\n", point_j->edge_count);
                    }
                    else{
                        printf("something is off\n");
                    }
                }
            }
            // Go to another point
            i++;
        }
    }
    free(shuffled_point_indexes);
    

}
