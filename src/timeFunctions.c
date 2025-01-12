#define _POSIX_C_SOURCE 199309L
#include <time.h>

double get_current_time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

double get_elapsed_time(double start_time) {
    return get_current_time() - start_time;
}

