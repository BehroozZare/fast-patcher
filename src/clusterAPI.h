#pragma once
#include <vector>
#include <tuple>
#include "Lloyd.h"

/**
 * @brief Public API for the fast-patcher graph clustering library.
 */
namespace Patcher {

/**
 * @brief Partition an undirected CSR graph into clusters of bounded size.
 *
 * This is the main entry point of the library.  The algorithm pipeline is:
 *   1. Detect connected components so that clustering never spans
 *      disconnected subgraphs.
 *   2. Select initial seeds (random or Morton-code based).
 *   3. Iteratively refine clusters with Lloyd-style BFS assignment
 *      followed by centre updates.
 *   4. Split any cluster that exceeds @p patch_size and repeat
 *      until all clusters are within the size bound.
 *
 * @param n                Number of vertices in the graph.
 * @param Gp               Row pointers of the CSR adjacency (size n + 1).
 * @param Gi               Column indices of the CSR adjacency (size Gp[n]).
 * @param patch_size        Maximum allowed number of vertices per cluster.
 * @param opt              Pointer to algorithm options, or @c nullptr to
 *                         use default settings.
 * @param[out] vertex_to_cluster  Resized to @p n on return;
 *             vertex_to_cluster[v] is the zero-based cluster id for vertex v.
 * @param vertex_positions Optional pointer to 3D vertex positions
 *                         (x, y, z) used for Morton-code seed selection.
 *                         Normalization and quantization are handled
 *                         internally.  If @c nullptr and
 *                         LloydOptions::MORTON_CODE is requested, the
 *                         algorithm falls back to random seeding.
 *
 * @see LloydOptions, Lloyd::multi_source_bfs, Lloyd::update_centers
 */
void create_clusters(
    int n,
    const int* Gp,
    const int* Gi,
    int patch_size,
    const LloydOptions* opt,
    std::vector<int>& vertex_to_cluster,
    const std::vector<std::tuple<double, double, double>>* vertex_positions = nullptr);

} // namespace Patcher
