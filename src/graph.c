#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "../headers/graph.h"


/**
 * Creates a random graph with the given number of vectors and dimensions
 * and random edges within the given maximum number of edges per point.
 * 
 * @param base_vectors The 2D array of vectors to use as points in the graph
 * @param base_num_dimensions The number of dimensions of each vector
 * @param max_edges The maximum number of edges each point can have
 * @param base_num_vectors The number of vectors in the graph
 * @return A pointer to the created graph
 */
Graph* create_random_graph(float **base_vectors, int base_num_dimensions, int max_edges, int base_num_vectors) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    
    if (!graph) {
        printf("Memory allocation failed!\n");
        exit(1);
    }

    graph->points = (Point*)malloc(base_num_vectors * sizeof(Point));

    if (!graph->points) {
        printf("Memory allocation failed!\n");
        exit(1);
    }

    // Create the graph
    printf("Creating random graph with %d vectors and %d dimensions\n", base_num_vectors, base_num_dimensions);

    // Initialize points with vector values
    for (int i = 0; i < base_num_vectors; i++) {
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
    add_random_edges(graph, max_edges, base_num_vectors);
    printf("Freeing base vectors..\n");
    for (int i = 0; i < base_num_vectors; i++) {
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
 * @param base_num_vectors The number of vectors in the graph
 */
void add_random_edges(Graph* graph, int max_edges, int base_num_vectors) {
    printf("Adding random edges..\n");
    for (int i = 0; i < base_num_vectors; i++) {
        for (int j = 0; j < max_edges; j++) {
            int random_index;
            do {
                // Generate a random index of a point to connect to
                random_index = rand() % base_num_vectors;
            } while (random_index == i); // Avoid self-loop
            // printf("Connecting point %d to point %d\n", i, random_index);

            // Add the edge if there is space
            if ((graph->points[i].edge_count + 1) < max_edges) {
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
 * @param base_num_vectors The number of vectors in the graph
 * @param base_num_coordinates The number of coordinates per vector
 * @param outputfd The file descriptor to write to
 */
void fprint_graph_coordinates(Graph* graph, int base_num_vectors, int base_num_coordinates, FILE *outputfd) {
    printf("Printing graph coordinates..\n");
    
    // Iterate over each point in the graph
    for (int i = 0; i < base_num_vectors; i++) {
        Point p = graph->points[i];
        
        // Print the coordinates of the current point
        fprintf(outputfd, "Point %d: ( ", p.index);
        for(int j = 0; j < base_num_coordinates; j++) {
            if (j == base_num_coordinates - 1) {
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
 * @param base_num_vectors The number of vectors in the graph
 * @param outputfd The file descriptor to write to
 */
void fprint_graph(Graph* graph, int base_num_vectors, FILE *outputfd) {
    printf("Printing graph..\n");
    for (int i = 0; i < base_num_vectors; i++) {
        Point p = graph->points[i];
        fprintf(outputfd, " Point %d edges: ", p.index);
        
        for (int j = 0; j < p.edge_count; j++) {
            fprintf(outputfd, "%d ", p.edges[j].to);
        }
        fprintf(outputfd, "\n");
    }
}

// Function to perform RobustPrune on a point in the graph
// void robust_prune(Graph* graph, int point_index, int *candidate_set, float a, int R, int base_num_coordinates)
void robust_prune(Graph* graph, int point_index, int *V, int *V_size, float a, int R, int base_num_dimensions) {
    
    printf("Performing RobustPrune on point %d\n", point_index);
    printf("Visited count: %d\n", *V_size);
    for (int i = 0; i < *V_size; i++) {
        printf("%d ", V[i]);
    }
    printf("\n");

    Point* p = &graph->points[point_index];
    
    // Clear the edges of the point (to be re-added during pruning)
    p->edge_count = 0;

    while (*(V_size) > 0 && p->edge_count < R) {
        // Find the closest point p* in the candidate set
        int best_index = -1;
        float min_distance = INF;

        // Find the closest point
        for (int i = 0; i < *(V_size); i++) {
            int candidate_index = V[i];
            Point* candidate_point = &graph->points[candidate_index];
            double distance = squared_euclidean_distance(p->coordinates, candidate_point->coordinates, base_num_dimensions);

            if (distance < min_distance) {
                min_distance = distance;
                best_index = i;
            }
        }

        if (best_index == -1) {
            break;
        }

        // Add p* to the neighbors of p
        int p_star_index = V[best_index];
        p->edges[p->edge_count].to = p_star_index;
        p->edge_count++;

        // Remove p* from the candidate set
        for (int i = best_index; i < *(V_size) - 1; i++) {
            V[i] = V[i + 1];
        }
        *(V_size)--;

        // Prune other points in the candidate set
        for (int i = 0; i < *(V_size); ) {
            int candidate_index = V[i];
            Point* candidate_point = &graph->points[candidate_index];
            double distance_p_star = squared_euclidean_distance(graph->points[p_star_index].coordinates, candidate_point->coordinates, base_num_dimensions);
            double distance_p = squared_euclidean_distance(p->coordinates, candidate_point->coordinates, base_num_dimensions);

            // IF:
            // After removing a candidate that doesn’t satisfy
            // the condition, you do not increment the loop counter (i stays the same).
            // This is because the next candidate has been shifted into
            // the current position in the array,
            // and you need to check that candidate next.
            //
            // ELSE: 
            // You increment the counter (i++) as usual,
            // moving on to check the next candidate in the set,
            // which is located to the same index
            if (a * distance_p_star <= distance_p) {
                // Remove the point if it doesn't satisfy the condition
                for (int j = i; j < *(V_size) - 1; j++) {
                    V[j] = V[j + 1];
                }
                *(V_size)--;
            } else {
                i++;
            }
        }
    }
}


/**
 * Inserts a new node into the sorted list l, keeping L closest nodes
 * 
 * @param l The sorted list of node indices
 * @param distances The corresponding distances of the nodes in l
 * @param new_node The index of the new node to insert
 * @param new_distance The distance of the new node
 * @param L The number of nodes to keep in the list
 */
void insert_l_closest(int *l, float *distances, int new_node, float new_distance, int L) {
    int i = L - 1;
    while (i >= 0 && distances[i] > new_distance) {
        if (i < L - 1) {
            distances[i + 1] = distances[i];
            l[i + 1] = l[i];
        }
        i--;
    }
    if (i < L - 1) {
        distances[i + 1] = new_distance;
        l[i + 1] = new_node;
    }
}

/**
 * After GreedySearch, find L closest nodes from visited V to Xq.
 *
 * @param V The set of visited nodes
 * @param V_size The size of the set
 * @param l The sorted list of node indices
 * @param L The number of nodes to keep in the list
 * @param Xq The query vector
 * @param graph The graph
 * @param dimensions The number of dimensions
 */
void find_L_closest(int *V, int V_size, int *l, int L, float *Xq, Graph *graph, int dimensions) {
    // Initialize distances with a large value
    float *distances = (float *)malloc(L * sizeof(float));
    for (int i = 0; i < L; i++) distances[i] = FLT_MAX;

    // Iterate over the visited nodes and find L closest nodes
    for (int i = 0; i < V_size; i++) {
        int node = V[i];
        float dist = (float)squared_euclidean_distance(graph->points[node].coordinates, Xq, dimensions);
        insert_l_closest(l, distances, node, dist, L);
    }

    free(distances);
}

// Greedy search algorithm with logging for debugging
/**
 * @brief Greedy search algorithm with logging for debugging
 * 
 * @param graph The graph to search
 * @param dimensions The number of dimensions in the graph
 * @param V The set of visited nodes
 * @param V_size The size of the set
 * @param l The sorted list of node indices
 * @param L The number of nodes to keep in the list
 * @param Xq The query vector
 * @param k The number of neighbors to consider
 * @param Xs The starting node
 */
void GreedySearch(Graph *graph, const int dimensions, int **V, int *V_size, int *l, const int L, float *Xq, const int k, const int Xs) {
    
    if (graph == NULL || V == NULL || *V == NULL || V_size == NULL || l == NULL || Xq == NULL) {
        printf("Error: One or more input parameters are NULL\n");
        return;
    }

    printf("Starting GreedySearch...\n");

    (*V)[*V_size] = Xs;  // Start search from Xs
    (*V_size)++;

    // Reallocate memory for V
    *V = realloc(*V, ((*V_size) + 1) * sizeof(int));
    if (*V == NULL) {
        printf("Error: Failed to allocate memory for V\n");
        exit(1);
    }

    int current_node = Xs;
    
    // Perform greedy search
    while (1) {
        float min_distance = FLT_MAX;
        int closest_node = -1;

        // Look at neighbors of the current node
        Point *p = &graph->points[current_node];

        for (int i = 0; i < p->edge_count; i++) {
            int neighbor = p->edges[i].to;
            if (!V_contains(*V, *V_size, neighbor)) {
                float dist = squared_euclidean_distance(graph->points[neighbor].coordinates, Xq, dimensions);
                if (dist < min_distance) {
                    min_distance = dist;
                    closest_node = neighbor;
                }
            }
        }
        printf("Closest neighbor of node %d: %d\n", current_node, closest_node);

        // If no more neighbors to explore, stop search
        if (closest_node == -1) {
            printf("No more neighbors to explore. Stopping search.\n");
            break;
        }

        // Add the closest neighbor to visited nodes V
        (*V)[*V_size] = closest_node;
        (*V_size)++;
        *V = realloc(*V, ((*V_size) + 1) * sizeof(int));
        if (*V == NULL) {
            printf("Error: Failed to allocate memory for V\n");
            exit(1);
        }
        current_node = closest_node;  // Move to the next node
    }

    // After search, find the L closest nodes from the visited ones
    find_L_closest(*V, *V_size, l, L, Xq, graph, dimensions);
    printf("GreedySearch completed successfully\n");
}

/**
 * Checks if a node is in the visited array V
 * 
 * @param V The visited array
 * @param V_size The size of the array
 * @param node The node to search for
 * @return 1 if the node is in the array, 0 otherwise
 */
int V_contains(int *V, int V_size, int node) {
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
