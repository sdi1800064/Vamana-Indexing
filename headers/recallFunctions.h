#ifndef RECALLFUNCTIONS_H
#define RECALLFUNCTIONS_H

void calculateRecallStitched(DatasetInfo* dataSet, QueryInfo* querySet, int** groundTruthSet, Graph* stitchedGraphs, int stitchedGraphs_count, int L, int k, int numOfThreads, FILE* resultFile);

#endif