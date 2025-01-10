#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main() {
    struct timespec start, end;
    long seconds, nanoseconds;
    double elapsed;

    // Get the start time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Code to measure
    for (volatile int i = 0; i < 100000000; i++); // Example workload

    // Get the end time
    clock_gettime(CLOCK_MONOTONIC, &end);

    // Calculate elapsed time
    seconds = end.tv_sec - start.tv_sec;
    nanoseconds = end.tv_nsec - start.tv_nsec;
    elapsed = seconds + nanoseconds*1e-9;

    printf("Elapsed time: %.9f seconds\n", elapsed);

    return 0;
