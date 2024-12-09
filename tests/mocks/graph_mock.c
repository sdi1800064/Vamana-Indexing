#include <stdlib.h>
#include "structs.h"

Graph* create_graph() {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->num_points = 8;
    graph->points = (Point *)malloc(sizeof(Point) * 8);
    int R = 5;

    for (int i = 0; i < 8; i++) {
        graph->points[i].index = i;
        graph->points[i].category = i;
        graph->points[i].coordinates = (float *)malloc(sizeof(float) * 100);
        for (int j = 0; j < 100; j++) {
            graph->points[i].coordinates[j] = j * i * 0.1;
        }
        graph->points[i].edge_count = R;
    }

    for(int i = 0; i < 8; i++) {
        graph->points[i].edges = (int *)malloc(sizeof(int) * R);
    }

    graph->points[0].edges[0] = 1;
    graph->points[0].edges[1] = 2;
    graph->points[0].edges[2] = 3;
    graph->points[0].edges[3] = 4;
    graph->points[0].edges[4] = 5;

    graph->points[1].edges[0] = 0;
    graph->points[1].edges[1] = 2;
    graph->points[1].edges[2] = 3;
    graph->points[1].edges[3] = 4;
    graph->points[1].edges[4] = 5;

    graph->points[2].edges[0] = 7;
    graph->points[2].edges[1] = 6;
    graph->points[2].edges[2] = 1;
    graph->points[2].edges[3] = 3;
    graph->points[2].edges[4] = 4;

    graph->points[3].edges[0] = 6;
    graph->points[3].edges[1] = 4;
    graph->points[3].edges[2] = 7;
    graph->points[3].edges[3] = 5;
    graph->points[3].edges[4] = 2;

    graph->points[4].edges[0] = 6;
    graph->points[4].edges[1] = 5;
    graph->points[4].edges[2] = 7;
    graph->points[4].edges[3] = 1;
    graph->points[4].edges[4] = 2;    

    graph->points[5].edges[0] = 6;
    graph->points[5].edges[1] = 4;
    graph->points[5].edges[2] = 7;
    graph->points[5].edges[3] = 3;
    graph->points[5].edges[4] = 2;

    graph->points[6].edges[0] = 7;
    graph->points[6].edges[1] = 5;
    graph->points[6].edges[2] = 4;
    graph->points[6].edges[3] = 3;
    graph->points[6].edges[4] = 2;

    graph->points[7].edges[0] = 0;
    graph->points[7].edges[1] = 1;
    graph->points[7].edges[2] = 2;
    graph->points[7].edges[3] = 3;
    graph->points[7].edges[4] = 4;

    return graph;
}