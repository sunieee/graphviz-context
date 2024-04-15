// graphviz_functions.h
#ifndef GRAPHVIZ_FUNCTIONS_H
#define GRAPHVIZ_FUNCTIONS_H

#include <stdbool.h>
#include "types.h"  // 或其他包含 graph_t 和 node_t 的定义的文件

bool is_left_subgraph(node_t *n);
bool is_right_subgraph(node_t *n);
bool is_context(node_t *n);

#endif // GRAPHVIZ_FUNCTIONS_H
