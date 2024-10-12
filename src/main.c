#include <stdio.h>
#include <stdlib.h>
#include "../headers/fvecs.h"

#define FILENAME "testSets/siftsmall/siftsmall_base.fvecs"

int main() {
    float** vectors;
    int num_vectors;
    int num_dimensions;

    // Read the vectors from the file
    read_fvecs(FILENAME, &vectors, &num_vectors, &num_dimensions);

    // Example: Print the first vector's components
    printf("First 10 dimentions:\n");
    for (int i = 0; i < num_vectors; i++) {
        printf("%d : ", i+1);
        for(int j=0; j < num_dimensions; j++){
                printf("%.1f, ", vectors[i][j]);
        }    
         printf("\n");   
    }
    printf("\n");

    // Free the allocated memory
    for (int i = 0; i < num_vectors; i++) {
        free(vectors[i]);
    }
    free(vectors);

    return 0;
}
