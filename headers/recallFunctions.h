#ifndef RECALLFUNCTIONS_H
#define RECALLFUNCTIONS_H

void calculateRecallStitched(DatasetInfo* dataset, QueryInfo* querySet, int** groundTruthSet, Graph* graphs, int stitchedGraphs_count, int L, int k, int numOfThreads);

#endif