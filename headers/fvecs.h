// fvecs.h

#ifndef FVEC_H
#define FVEC_H

#include <stdio.h>
#include "graph.h"

void add_or_increment_filter(filterInfo *filters, int category);

DatasetInfo* read_dataset(const char *filename, uint32_t *total_vectors, filterInfo *filters);
void free_dataset(DatasetInfo *dataset);
void print_dataset(DatasetInfo *dataset);
void cprint_dataset(DatasetInfo *dataset, int target_category);

QueryInfo* read_query_dataset(const char *filename, uint32_t *total_vectors);
void free_query_dataset(QueryInfo *dataset);
void print_query_dataset(QueryInfo *dataset);
void cprint_query_dataset(QueryInfo *dataset, int category);



#endif // FVEC_H
