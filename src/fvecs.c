#include <stdio.h>
#include <stdlib.h>
#include "../headers/fvecs.h"

// Function to read an fvecs or an ivecs file
void read_fvecs(const char* filename, float*** vectors, int* num_vectors, int* dimension) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        printf("Error opening file.\n");
        exit(EXIT_FAILURE);
    }

    // Get the size of the file
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);

    // Read the first 4 bytes to get the dimensionality (d)
    fread(dimension, sizeof(int), 1, file);

    // Calculate the number of vectors
    *num_vectors = file_size / ((*dimension + 1) * sizeof(float));

    // Allocate memory for the vectors
    *vectors = (float**)malloc((*num_vectors) * sizeof(float*));
    for (int i = 0; i < *num_vectors; i++) {
        (*vectors)[i] = (float*)malloc((*dimension) * sizeof(float));
    }

    // Set file pointer back to the beginning
    fseek(file, 0, SEEK_SET);

    // Read all vectors from the file
    for (int i = 0; i < *num_vectors; i++) {
        int dim;
        fread(&dim, sizeof(int), 1, file);  // Read the dimension (should always match *dimension)

        if (dim != *dimension) {
            printf("Error: Dimensionality mismatch.\n");
            exit(EXIT_FAILURE);
        }

        // Read the vector components
        fread((*vectors)[i], sizeof(float), *dimension, file);
    }

    fclose(file);
}

void fprintFloatVectors(float** vectors, int num_vectors, int dimension, FILE *outputfd) {
    for (int i = 0; i < num_vectors; i++) {
        fprintf(outputfd, "Point %d: ", i);
        for (int j = 0; j < dimension; j++) {
            fprintf(outputfd, "%.1f ", vectors[i][j]);
        }
        fprintf(outputfd, "\n");
    }
}