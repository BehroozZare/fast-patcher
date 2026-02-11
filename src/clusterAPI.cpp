#include "clusterAPI.h"
#include "ConnectedComponent.h"
#include "BFS.h"
#include <cassert>

namespace Patcher {

void create_clusters(
    int n,
    const int* Gp,
    const int* Gi,
    int patch_size,
    const LloydOptions* opt,
    std::vector<int>& vertex_to_cluster)
{
    // ── Step 0: Load inputs & defaults ────────────────────────────────
    LloydOptions defaults;
    if (opt == nullptr)
        opt = &defaults;

    // ── Step 1: Connected components ──────────────────────────────────
    // Label every vertex with its CC id so we never cluster across
    // disconnected subgraphs.
    std::vector<int> vertex_to_cc;
    int num_cc = ConnectedComponent::compute(Gp, Gi, n, vertex_to_cc);

    //Count the number of vertices in each CC
    std::vector<int> num_vertices_per_cc(num_cc, 0);
    for (int v = 0; v < n; v++)
        num_vertices_per_cc[vertex_to_cc[v]]++;

    // Prepare the output
    vertex_to_cluster.assign(n, -1);

    // ── Step 2: Initialize seeds ──────────────────────────────────
    std::vector<int> seeds;
    Lloyd::initialize_seeds(vertex_to_cc, num_vertices_per_cc, patch_size, opt->random_seed, seeds);

    // ── Step 3: Lloyd-style iterative refinement ──────────────────
    std::vector<int> node_to_cluster, dist;
    int total_num_clusters = seeds.size();

    //Init node_to_cluster with the seeds
    node_to_cluster.assign(n, -1);
    for(int i = 0; i < seeds.size(); i++) {
        node_to_cluster[seeds[i]] = i;
    }

    // Outer loop: repeat rounds of Lloyd iterations, then check
    // whether any cluster is oversized.
    for(int iter = 0; iter < opt->lloyd_max_iterations; iter+=opt->lloyd_iters_to_add_seed) {
        std::vector<int> cluster_sizes(total_num_clusters, 0);
        // Run a fixed number of Lloyd iterations.
        for (int seed_iter = 0; seed_iter < opt->lloyd_iters_to_add_seed; seed_iter++) {
            // Step 4.1 + 4.2: Multi-source BFS from seeds, assign vertices.
            node_to_cluster.assign(n, -1);
            dist.assign(n, -1);
            for (int i = 0; i < static_cast<int>(seeds.size()); i++) {
                node_to_cluster[seeds[i]] = i;
                dist[seeds[i]] = 0;
            }
            BFS::multi_source_bfs(Gp, Gi, n, node_to_cluster, dist);

            //Compute the maximum size of the clusters
            int max_cluster_size = 0;
            cluster_sizes.clear();
            cluster_sizes.resize(total_num_clusters, 0);
            for(int i = 0; i < n; i++) {
                cluster_sizes[node_to_cluster[i]]++;
            }
            max_cluster_size = *std::max_element(cluster_sizes.begin(), cluster_sizes.end());
            if(max_cluster_size <= patch_size) {
                break;
            }
            // Step 4.3: Update cluster centers.
            int num_clusters = static_cast<int>(seeds.size());
            Lloyd::update_centers(
                Gp, Gi, n, node_to_cluster,
                num_clusters, opt->random_seed, seeds);
        }

        // Step 4.4: Check for oversized clusters and split them.
        int added = Lloyd::split_oversized(
            Gp, Gi, n, node_to_cluster,
            cluster_sizes,
            total_num_clusters, patch_size, opt->random_seed, seeds);
        if(added == 0) {
            break;
        }
        total_num_clusters += added;
    }

    // ── Step 5: Final labeling ────────────────────────────────────────
    // Copy the cluster assignments to the output vector.
    vertex_to_cluster = node_to_cluster;
}

} // namespace Patcher
