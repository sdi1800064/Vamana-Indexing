// fvecs.h

#ifndef FVEC_H
#define FVEC_H

// Function prototype for reading an fvecs file
void read_fvecs(const char* filename, float*** vectors, int* num_vectors, int* dimension);
void fprintFloatVectors(float** vectors, int num_vectors, int dimension, FILE *outputfd);

#endif // FVEC_H
