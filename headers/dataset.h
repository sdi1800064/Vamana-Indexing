// fvecs.h

#ifndef DATASET_H
#define DATASET_H

#include <stdio.h>
#include "graph.h"
#include "structs.h"

void add_or_increment_filter(filterInfo *filters, int point_category, int point_index);

DatasetInfo* read_dataset(const char *filename);
void free_dataset(DatasetInfo *dataset);
void print_dataset(DatasetInfo *dataset);
void cprint_dataset(DatasetInfo *dataset, int target_category);

QueryInfo* read_query_dataset(const char *filename);
void free_query_dataset(QueryInfo *dataset);
void print_query_dataset(QueryInfo *dataset);
void cprint_query_dataset(QueryInfo *dataset, int category);



#endif // FVEC_H
