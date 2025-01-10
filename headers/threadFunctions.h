#ifndef threadFunctions_h
#define threadFunctions_h

void* thread_function_vamana_indexing(void* args);
Graph* threadStitchedVamanaIndexing(DatasetInfo* dataset, int L_small, float a, int R_small, int numOfThreads);

#endif
