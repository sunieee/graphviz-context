# Graphviz-context

An enhanced Graphviz DOT layout supporting context proxy nodes to present influence from outer part of layered graph to focused part of layered graph.  

[TOC]

## Originally build

Graphviz - Graph Drawing Programs from AT&T Research and Lucent Bell Labs

See doc/build.html for prerequisites and detailed build notes.



## Guide

输入DOT文件仅适用于布局包含context节点和边的层次图，需要遵守特定规则：

- cross_type 指定特定特定交叉类型
- 主体部分指定3个子图分别是 left, focus, right
- left和right子图包含所有层次的节点，并从前到后连接，设置边的权重为10000
- focus边的权重设置为 `focus_edge_weight`
- 添加虚拟边，用于连接零散的连通块



### Node Type

| Type ID | Type       | description                              | Cost/panelty | Name                                         |
| ------- | ---------- | ---------------------------------------- | ------------ | -------------------------------------------- |
| 0       | $$v_f$$    | Node of focus subgraph                   | $$v_f$$      | <id>                                         |
| 1       | $$v_{ff}$$ | Virtual node between $$v_f$$ and $$v_f$$ |              | <rank>\_<id>\_<id>                           |
| 2       | $$v_c$$    | Node of context                          | $$v_c$$      | l<year> r<year>                              |
| 3       | $$v_{cf}$$ | Virtual node between $$v_f$$ and $$v_c$$ |              | <rank>\_l<year>\_<id>, <rank>\_<id>\_r<year> |

### Edge Type

| Type ID | Type         | description                                        | Cost/panelty                  | Name                              |
| ------- | ------------ | -------------------------------------------------- | ----------------------------- | --------------------------------- |
| 0       | $$e_{ff}$$   | Edge between $$v_f$$and $$v_f$$(direct connection) | focus_edge_weight = 10        | <head node name>-<tail node name> |
| 1       | $$e_{fff}$$  | Edge between $$v_f$$and $$v_{ff}$$                 |                               |                                   |
| 2       | $$e_{ffff}$$ | Edge between $$v_{ff}$$and $$v_{ff}$$              |                               |                                   |
| 3       | $$e_{cf}$$   | Edge between $$v_c$$and $$v_f$$(direct connection) | weight (size of context edge) |                                   |
| 4       | $$e_{cff}$$  | Edge between $$v_f$$and $$v_{cf}$$                 |                               |                                   |
| 5       | $$e_{ccf}$$  | Edge between $$v_c$$and $$v_{cf}$$                 |                               |                                   |
| 6       | $$e_{cfcf}$$ | Edge between $$v_{cf}$$and $$v_{cf}$$              |                               |                                   |
| 7       | $$e_{cc}$$   | Edge between $$v_c$$and $$v_{cc}$$                 | 10000 (inf)                   |                                   |
| 8       | Undefined    |                                                    |                               |                                   |

### Crossing Type

| Type ID | Type           | description       | Cost/panelty            | Name                        |
| ------- | -------------- | ----------------- | ----------------------- | --------------------------- |
| 0       |                | All crossing      | A = I + E + H           | (<edge1 name>,<edge2 name>) |
| 1       |                | Weighted crossing | weight(e1) * weight(e2) |                             |
| 2       | $$(e_f, e_f)$$ | Internal crossing |                         |                             |
| 3       | $$(e_c, e_c)$$ | External crossing |                         |                             |
| 4       | $$(e_c, e_f)$$ | Hybrid crossing   |                         |                             |



## Change log

### flexible optimization(4.26)

These modifications aim to provide a more robust and flexible optimization process by enabling different types of crossing optimizations based on edge weights and other criteria.

**New Features**

- **Edge Weight Assignment**: Added functionality to assign weights to edges, enhancing the existing layout optimization.
- **Crossing Types Optimization**: Introduced four different types of crossings to optimize:
  - All Cross
  - Weighted Cross
  - Internal Cross
  - External Cross

**Modifications**

- Edge Penalty Calculation:
  - Replaced direct access to edge penalties with a new function `get_edge_penalty`, which internally calls `get_edge_penalty_with_type` depending on the crossing type set on the graph. This abstraction allows handling different types of edge weights based on the crossing optimization requirement.
- Cross Calculation Functions:
  - Updated `in_cross` and `out_cross` functions to utilize `get_edge_penalty` instead of directly accessing edge penalties.
  - Modified `local_cross` to accept a `type` parameter to specify the type of crossing calculation.
  - Introduced a new version of the `rcross` function that takes a `type` argument, allowing different handling based on the crossing type.

**Detailed Modifications**

- Type Determination Functionality:
  - Extended `get_type` function to include more types based on different criteria, simplifying the handling of node and edge types across different crossing optimizations.
- Documentation and Comments: Improved inline comments for better understanding of changes and functionalities added, particularly explaining the purpose of each major modification or addition.
- Error Handling: Added error logging in the `get_crossing_type` function to handle cases where the crossing type is not set or is invalid.

- Syntax and Semantic Corrections:
  - Corrected various semantic issues to ensure that the new functionalities integrate smoothly with the existing codebase without causing regressions in the optimization process.

### context constrain(4.16)

Support special context nodes in the DOT graph layout, ensuring that nodes prefixed with 'l' and 'r' are positioned on the extreme left and right of the layout, respectively. These changes enhance the DOT layout engine's capability to handle layout constraints dictated by node naming conventions.

**New Features**

- **Context Node Integration**: Implemented context nodes to enforce specific positioning constraints in the DOT layout. Nodes starting with 'l' are forced to the leftmost position, while nodes beginning with 'r' are forced to the rightmost position in the graph layout.
- **Constraint Handling in Layout**: Enhanced layout processing to respect the new context node constraints during various layout steps such as rank building and transposition.

**Modifications**

- Utility Functions:
  - Added new utility functions `is_left_subgraph`, `is_right_subgraph`, and `is_context` to check if a node qualifies as a context node based on its naming convention.
  - Modified `get_name` and `get_type` functions to handle context nodes and virtual nodes with custom logic for name derivation and type classification.
- Graph Processing Adjustments:
  - Adjusted the `examine_order` function to ensure context nodes are correctly positioned at the start or end of their ranks, depending on their type (left or right context).
  - Updated the `transpose_step` and `mincross_step` functions to skip context nodes when adjusting node positions, preserving their designated placements.

**Detailed Modifications**

- Error Handling: Implemented checks and assertions specifically for virtual and context nodes to ensure stability and correctness of the graph layout operations.

- Node Position Adjustments:
  - Integrated checks in the `build_ranks` and `reorder` functions to accommodate the special positioning requirements of context nodes, ensuring they are placed appropriately according to their designated constraints.
- Graph Layout Integrity:
  - Ensured that the introduction of context nodes does not disrupt the overall layout integrity by maintaining existing graph layout algorithms' functionality with additional checks and balances for context nodes.



###### 
