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
    printf("All graph-points initiallized\n");


    printf("Freeing base vectors..\n");
    for (int i = 0; i < base_num_points; i++) {
        free(base_vectors[i]);
    }
    free(base_vectors);
    printf("Graph created\n");

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
    // printf("Adding edge from %d to %d | Edge count before adding: %d\n", point->index, toIndex, point->edge_count);
    point->edges[point->edge_count] = toIndex;
    point->edge_count += 1;
    // printf("Edge count after adding: %d\n", point->edge_count);

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
    printf("\n\nRobust Prune: STARTING\n");
    printf(" -> V_size: %d\n", V_size);

    Point *p = &graph->points[p_index];


    // ============== V <- V U Nout(p) ============== //
    V_size += p->edge_count;

    V = realloc(V, (V_size + 1) * sizeof(int));
    if (!V) {
        printf("Memory allocation of V failed!\n");
        exit(1);
    }

    for (int i = 0; i < p->edge_count; i++) {
        V[i + V_size - p->edge_count] = p->edges[i];

        // Set the edge to point -1 ( NULL )
        p->edges[i] = -1;
    }

    p->edge_count = 0;
    int m =0;
    int min_distance = INF;
    int min_index = -1;

    while(V_size > 0) {
        
        // ============== p* <- arg min p'ε V d(p, p') ============== //
        // 
        min_distance = INF;
        min_index = -1;
        // printf("m: %d | mint_index: %d\n", m, min_index);
        Point *p_star = NULL;
        for (int i = 0; i < V_size; i++) {
            p_star = &graph->points[V[i]];
            float distance = squared_euclidean_distance(graph->points[p_index].coordinates, p_star->coordinates, graph->num_dimensions);
            if (distance < min_distance && V[i] != p_index) {
                min_distance = distance;
                min_index = i;
            }
        }

        if (min_index == -1) {
            printf("Error: min_index not found\n");
            exit(1);
        }

        printf("m: %d min_index: V[%d] = %d\n", m, min_index, V[min_index]);
        // ============== Nout(p) <- Nout(p) U {p*} ============== //
        if(p->edge_count < R){
            addEdge(p, V[min_index]);
        }
        // ============== if |Nout(p) = R| then break ============== //
        if (p->edge_count == R) {
            break;
        }

        // ============== for p' ε V do (if a d(p*, p') < d(p, p') then remove p' from V) ============== //
        printf("\n");for (int i = 0; i < V_size; i++) {
            // printf("==Checking V[%d] = %d \n", i, V[i]);
            if ((a * squared_euclidean_distance(graph->points[V[i]].coordinates, p_star->coordinates, graph->num_dimensions)) <= min_distance) {
                // printf("===Removing V[%d] = %d\n", i, V[i]);
                V[i] = V[V_size - 1];
                V_size--;
                i--;
            }
        }

        V = realloc(V, (V_size + 1) * sizeof(int));
        if (!V) {
            printf("Memory allocation of V failed!\n");
            exit(1);
        }
        m++;
    }
    printf(" -> edges after pruning point %d: %d\n", p_index, p->edge_count);
    for(int i = 0; i < p->edge_count; i++) {
        printf("%d ", p->edges[i]);
    }
    printf("\n");

    printf("Robust Prune: FINISHED\n");
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
    printf("= Sorting array of size %d...\n", array_size);
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

    printf("= Done sorting array\n");
}


// Helper function to check if Lamda \ V is not empty
int exists_in_difference(int *Lamda, int Lamda_size, bool *Visited) {
    printf("Looking a point in Lamda that isnt visited in a list of size %d\n", Lamda_size);
    for (int i = 0; i < Lamda_size; i++) {
        int index = Lamda[i];
        if(!Visited[index]){
            printf("Found one!\n====Lamda[i] = %d=====\n", Lamda[i]);
            return Visited[index];
        }
    }
    printf("Didnt find any\n");
    return 0;
}

// Function to add an element to a dynamically resized array
void add_to_dynamic_array(int **array, int *size, int element) {
    
    if (array == NULL || size == NULL) {
        fprintf(stderr, "Invalid input: null pointer\n");
        return;
    }
    
    int *temp = realloc(*array, (*size + 1) * sizeof(int));
    if (temp == NULL) {
        perror("Failed to allocate memory");
        return;
    }

    *array = temp;
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

int *get_the_difference(int *Lamda, int Lamda_size, int *V, int V_size, int *Lamda_minus_V_size) {

    *Lamda_minus_V_size = 0;
    int *Lamda_minus_V = (int *)malloc(sizeof(int));
    for(int i = 0; i < Lamda_size; i++) {
        if(!is_in_array(V, V_size, Lamda[i])) {
            add_to_dynamic_array(&Lamda_minus_V, Lamda_minus_V_size, Lamda[i]);
        }
    }
    return Lamda_minus_V;
}

// Greedy search function
void greedy_search(Graph *graph, float *Xq, int start_index, int **V, int *V_size, int **Lamda, int *Lamda_size, int L, int k) {

    // Allocate initial space for V and Lamda
    *V = (int *)malloc(sizeof(int) * graph->num_points);
    *V_size = 0;
    *Lamda = (int *)malloc(sizeof(int) * graph->num_points);
    *Lamda_size = 0;

    // Initialize Lamda with the start point
    (*Lamda)[(*Lamda_size)++] = start_index;

    // Create an int array that contains points that are in Lamda but not in V
    int *Lamda_minus_V;
    int Lamda_minus_V_size = 0;
    Lamda_minus_V = get_the_difference(*Lamda, *Lamda_size, *V, *V_size, &Lamda_minus_V_size);

    Point *p_star = NULL;

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
        p_star = &graph->points[closest_index];

        // Add p* to V
        add_to_dynamic_array(V, V_size, p_star->index);

        // Add the neighbors of p* to Lamda
        for (int i = 0; i < p_star->edge_count; i++) {
            int neighbor_index = p_star->edges[i];
            add_to_dynamic_array(Lamda, Lamda_size, neighbor_index);
        }

        // Check if |Lamda| > L
        if (*Lamda_size > L) {
            // Update Lamda to contain the L closest points to Xq
            sort_array(graph, *Lamda, *Lamda_size, Xq);
            *Lamda = (int *)realloc(*Lamda, L * sizeof(int));
            *Lamda_size = L;
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
    printf("No duplicates found.\n");
}

void vamana_indexing(Graph *graph, int k, int L, float a, int R, FILE *outputfd) {

    printf("Starting Vamana Indexing\n");

    printf("k = %d\nL = %d\na = %f\nR = %d\n", k, L, a, R);
    
    int max = graph->num_points;                                                                    // max - Number of points in graph
    int num_sample_points = graph->num_points / 10;                                                 // num_sample_points - Number of sample points     


    //  =============== MEDOID CALCULATION ================ //
    int *sample_point_indexes = sample_points(max, num_sample_points);                              // sample_point_indexes - Indexes of sampled points
    int medoid_index = calculate_medoid(graph, sample_point_indexes, num_sample_points);            // medoid_index - Index of medoid
    printf("Medoid index: %d\n", medoid_index);




    // =============== FIRST INITIALLIZATIONS ================ //
    bool *shuffled_point_indexes = (bool*)calloc(max, sizeof(bool));                                // shuffled_point_indexes - Indexes of points in a random order
    
    int s_index = -1;         // Temporary index for random selection
    int *V;                         // V - Visited List
    int V_size = 0;                 // V_size - Size of V
    int *lamda;                         // L~ - k nearest neighbors list
    int lamda_size = 0;             // lamda_size - Size of L~

    int i = 0;

    // =============== VAMANA INDEXING FOR S(i) ================ // 
    while (i < max) {
        printf("Points remaining %d\n", max - i);
        // Randomly select an index from the graph
        s_index = rand() % max;
        // Ensure no duplicate indices are selected
        if (shuffled_point_indexes[s_index] == false) {

            // Mark the current index as visited
            shuffled_point_indexes[s_index] = true;
            printf("Indexing point %d\n", s_index);


            // =============== GREEDY SEARCH ================ //
            printf("Calling greedy search with s: %d, Xq: %d, V_size: %d, lamda_size: %d, L: %d, k: %d\n", medoid_index, graph->points[s_index].index, V_size, lamda_size, L, 1);
            greedy_search(graph, graph->points[s_index].coordinates, medoid_index, &V, &V_size, &lamda, &lamda_size, L, 1);

            printf("V size: %d\n", V_size);
            for(int h = 0; h < V_size; h++) {
                printf("V[%d] = %d\n", h, V[h]);
            }
            check_for_duplicates(V, V_size);

            printf("Closest point to Xq %d : %d \n", s_index, lamda[0]);
            

            exit(3);

            // =============== ROBUST PRUNE ================ //
            printf("Calling robust prune with s_index: %d and V_size: %d\n", s_index, V_size);
            robustPrune(graph, s_index, V, V_size, a, R);

            // =============== FOR J IN Nout(S(i)) ================ //
            Point *point_j = NULL;
    
            for (int j = 0; j < graph->points[s_index].edge_count; j++) {
                printf("Vamana: Checking Nout(%d) -> point %d\n", graph->points[s_index].index, graph->points[s_index].edges[j]);
                
                if(graph->points[s_index].edges[j] == -1) {
                    printf("Invalid edge index\n");
                    exit(1);   
                }
                point_j = &graph->points[graph->points[s_index].edges[j]];
                // ============= if |Nout(j) U {s(i)} | > R then ================ //
                if(point_j->edge_count + 1 > R) {

                    // =============== ROBUST PRUNE ( j, Nout(j) U {s(i)}, a, R)================ //

                    // Create the Nout(j) U {s(i)} list
                    int *temp_v = (int*)malloc((point_j->edge_count + 1) * sizeof(int));
                    int temp_v_size = point_j->edge_count + 1;

                    for(int m = 0; m < point_j->edge_count; m++) {
                        temp_v[m] = point_j->edges[m];
                    }
                    temp_v[point_j->edge_count] = s_index;

                    // Robust Prune
                    robustPrune(graph, point_j->index, temp_v, temp_v_size, a, R);
                    free(temp_v);
                }
                else{
                    // =============== UPDATE Nout(j) <- Nout(j) U {s(i)} ================ //
                    point_j->edge_count++;
                    point_j->edges[point_j->edge_count - 1] = s_index;
                    if(point_j->edge_count >=R) {
                        printf("Point %d has %d/%dedges\n", point_j->index, point_j->edge_count, R);
                    }
                }

            }


            i++;
            free(V);
            free(lamda);
        }
    }
    free(shuffled_point_indexes);
    

}
