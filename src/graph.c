#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "../headers/graph.h"
#include "../headers/knnGraph.h"

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

// Function to randomly connect points with edges
void add_random_edges(Graph* graph, int max_edges, int base_num_vectors) {
    printf("Adding random edges..\n");
    for (int i = 0; i < base_num_vectors; i++) {
        int temp_edge_count = rand() % (max_edges + 1) + 1; // Random number of edges (0 to max_edges)
        // printf("Point %d will have %d edges\n", i, temp_edge_count);
        
        for (int j = 0; j < temp_edge_count; j++) {
            int random_index;
            do {
                random_index = rand() % base_num_vectors; // Random point to connect to
            } while (random_index == i); // Avoid self-loop
            // printf("Connecting point %d to point %d\n", i, random_index);

            // Add the edge if there is space
            if (graph->points[i].edge_count < temp_edge_count) {
                graph->points[i].edges[graph->points[i].edge_count].to = random_index;
                graph->points[i].edge_count++;
            }
            // printf("Point %d now has %d edges\n", i, graph->points[i].edge_count);
        }
        // printf("Point %d has edges\n", i);
    }
    printf("Edges added\n");
}

// Function to print the graph for debugging
void fprint_graph_coordinates(Graph* graph, int base_num_vectors, int base_num_coordinates, FILE *outputfd) {
    printf("Printing graph..\n");
    for (int i = 0; i < base_num_vectors; i++) {
        Point p = graph->points[i];
        fprintf(outputfd, "Point %d: ( ", p.index);
        for(int j = 0; j < base_num_coordinates; j++) {
            if (j == base_num_coordinates - 1) {
                fprintf(outputfd, "%f )\n", p.coordinates[j]);
                break;
            }
            fprintf(outputfd, "%.2f, ", p.coordinates[j]);
        }
        fprintf(outputfd, " Point %d edges: ", p.index);
        
        for (int j = 0; j < p.edge_count; j++) {
            fprintf(outputfd, "%d ", p.edges[j].to);
        }
        fprintf(outputfd, "\n");
    }
    printf("Printing done\n");
}

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
    printf("Printing done\n");
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


// Insert into sorted list l, keeping L closest nodes
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

// After GreedySearch, find L closest nodes from visited V to Xq
void find_L_closest(int *V, int V_size, int *l, int L, float *Xq, Graph *graph, int dimensions) {
    printf("Finding %d closest nodes out of %d..\n", L, V_size);
    float *distances = (float *)malloc(L * sizeof(float));
    for (int i = 0; i < L; i++) distances[i] = FLT_MAX; // Initialize distances with a large value

    for (int i = 0; i < V_size; i++) {
        int node = V[i];
        float dist = (float)squared_euclidean_distance(graph->points[node].coordinates, Xq, dimensions);
        insert_l_closest(l, distances, node, dist, L);
    }

    free(distances);
}

// Greedy search algorithm with logging for debugging
void GreedySearch(Graph *graph, const int dimensions, int **V, int *V_size, int *l, const int L, float *Xq, const int k, const int Xs) {
    
    if (graph == NULL || V == NULL || *V == NULL || V_size == NULL || l == NULL || Xq == NULL) {
        printf("Error: One or more input parameters are NULL\n");
        return;
    }

    printf("Starting GreedySearch...\n");

    (*V)[*V_size] = Xs;  // Start search from Xs
    (*V_size)++;
    printf("Added starting node Xs to visited nodes V\n");

    // Reallocate memory for V
    *V = realloc(*V, ((*V_size) + 1) * sizeof(int));
    if (*V == NULL) {
        printf("Error: Failed to allocate memory for V\n");
        exit(1);
    }

    int current_node = Xs;
    
    // Perform greedy search
    while (1) {
        printf("Current node: %d\n", current_node);
        float min_distance = FLT_MAX;
        int closest_node = -1;

        // Look at neighbors of the current node
        Point *p = &graph->points[current_node];

        // printf("Node %d has %d edges\n", current_node, p->edge_count);
        // for (int i = 0; i < p->edge_count; i++) {
        //     printf(" %d\n", p->edges[i].to);
        // }
        // printf("\n");
        // printf(" V_size: %d, V: ", *V_size);
        // for (int j = 0; j < *V_size; j++) {
        //     printf(" %d,", (*V)[j]);
        // }
        // printf("\n");
        for (int i = 0; i < p->edge_count; i++) {
            
            int neighbor = p->edges[i].to;
            printf("Checking neighbor %d\n", neighbor); 
            printf(" V_size: %d, V: ", *V_size);
            // for (int j = 0; j < *V_size; j++) {
            //     printf(" %d,", (*V)[j]);
            // }
            printf("\n");
            if (!V_contains(*V, *V_size, neighbor)) {
                printf("Calculating distance to neighbor %d\n", neighbor);
                float dist = squared_euclidean_distance(graph->points[neighbor].coordinates, Xq, dimensions);
                if (dist < min_distance) {
                    min_distance = dist;
                    closest_node = neighbor;
                }
            }
            printf("Checked neighbor %d\n", neighbor);
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
        printf("Added closest neighbor %d to visited nodes V\n", closest_node);
        current_node = closest_node;  // Move to the next node
    }

    // After search, find the L closest nodes from the visited ones
    find_L_closest(*V, *V_size, l, L, Xq, graph, dimensions);
    printf("GreedySearch completed successfully\n");
}

// Check if node is in visited array
int V_contains(int *V, int V_size, int node) {
    printf("Checking if %d is in V..\n", node);
    for (int i = 0; i < V_size; i++) {
        if (V[i] == node) {
            printf("Found %d in V\n", node);
            return 1;
        }
    }
    printf("%d not found in V\n", node);
    return 0;
}
