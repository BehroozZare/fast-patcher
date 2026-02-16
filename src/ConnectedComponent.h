#pragma once
#include <vector>

/**
 * @brief Connected-component detection for CSR graphs.
 *
 * Provides a BFS-based algorithm to label every vertex of an undirected
 * graph (stored in Compressed Sparse Row format) with its connected-component
 * identifier.
 *
 * @par Graph convention
 * | Symbol | Meaning                              |
 * |--------|--------------------------------------|
 * | Gp     | Row pointers (size G_N + 1)          |
 * | Gi     | Column indices (size Gp[G_N])        |
 * | G_N    | Number of vertices                   |
 */
namespace ConnectedComponent {

/**
 * @brief Compute connected components of a CSR graph via BFS.
 *
 * Walks every unvisited vertex, performs a BFS from it, and labels all
 * reachable vertices with the same connected-component id.
 *
 * @param Gp   Row pointers of the CSR adjacency (size G_N + 1).
 * @param Gi   Column indices of the CSR adjacency (size Gp[G_N]).
 * @param G_N  Number of vertices in the graph.
 * @param[out] vertex_to_cc  Resized to G_N; vertex_to_cc[v] is the
 *             zero-based connected-component id assigned to vertex v.
 * @return The total number of connected components found.
 */
int compute(
    const int* Gp, const int* Gi, int G_N,
    std::vector<int>& vertex_to_cc);

} // namespace ConnectedComponent
