#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main() {
    char new_graph_file_name[100];
    int R = 20;
    int L = 200;
    float a = 1.2;
    int numOfPoints = 10001;
    snprintf(new_graph_file_name, sizeof(new_graph_file_name), "%s_R%d_L%d_a%.2f_#%d.bin", 
            "graphs/stitchedGraph", R, L, a, numOfPoints);
    FILE* output;
    output = fopen(new_graph_file_name, "w");
    if (output == NULL) {
        perror("Error opening output file");
        exit(EXIT_FAILURE);
    }
    printf("Found it!\n");
    fclose(output);

    return 0;
}