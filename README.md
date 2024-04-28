# Graphviz-context

An enhanced Graphviz DOT layout supporting context proxy nodes to present influence from outer part of layered graph to focused part of layered graph.  

[TOC]

## Originally build

Graphviz - Graph Drawing Programs from AT&T Research and Lucent Bell Labs

See doc/build.html for prerequisites and detailed build notes.



## Guide

Input DOT files are only applicable to layouts that contain context nodes and edges, and specific rules need to be followed.

- Graph parameters:
  - cross_type: Specify the type of crossover as the optimization target, see the Crossing Type table for details
  - Concentrate = true: Use edge bundling algorithm
  - concentrate_type: When `concentrate = true `, fine bundling control, see Change Log > Enhanced Edge Bundling
- The main part specifies three subgraphs: left, focus, and right.
  - The left and right subgraphs contain nodes of all levels and are connected from front to back, with the weight of the edge set to 10000.
  - The weight of the focus edge is set to `focus_edge_weight`
- The layer part is the level where each node is located, with the same rank as the year context node
- The context defines the connection between the conext node and the focus node
- Add virtual edges to connect scattered connected blocks. If they are not connected, an error will be reported.



### Node Type

| Type ID | Type       | description                              | Cost/panelty | Name                  |
| ------- | ---------- | ---------------------------------------- | ------------ | --------------------- |
| 0       | $$v_f$$    | Node of focus subgraph                   | $$v_f$$      | <id>                  |
| 1       | $$v_{ff}$$ | Virtual node between $$v_f$$ and $$v_f$$ |              | <rank>\_<id>\_<id>    |
| 2       | $$v_c$$    | Node of context                          | $$v_c$$      | l<year> r<year>       |
| 3       | $$v_{cf}$$ | Virtual node from $$v_c$$ to  $$v_f$$    |              | <rank>\_l<year>\_<id> |
| 4       | $$v_{fc}$$ | Virtual node from $$v_f$$ to  $$v_c$$    |              | <rank>\_<id>\_r<year> |

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

### Enhanced Edge Bundling(4.28)

The latest update to our graph visualization library introduces refined edge bundling capabilities with a new attribute `concentrate_type`. This attribute expands the traditional edge bundling method, allowing for more nuanced control over how edges are concentrated based on their contextual relationships. The change log details modifications in the `mergevirtual` and `dot_concentrate` functions, enabling specific bundling behaviors:

- **Default bundling (0)**: Engages the standard edge concentration without additional context.
- **Concentrate context edges (1)**: Focuses on bundling edges identified as context-specific.
- **Concentrate focus edges (2)**: Bundles edges identified as focus-specific.
- **Separate treatment of context and focus edges (3)**: Applies concentration to both context and focus edges independently.
- **Context edge concentration starting from focus nodes (4)**: Concentrates context edges initiating from focus nodes.
- **Exclude focus node initiation for context edges (5)**: Concentrates context edges without starting from focus nodes.

These enhancements are aimed at providing users with greater flexibility in visualizing complex graphs, particularly useful in applications requiring detailed analysis of network interactions and relationships.

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
- Error Handling: Added error logging in the `get_attr(g, "crossing_type")` function to handle cases where the crossing type is not set or is invalid.

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


