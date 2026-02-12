#include "clusterAPI.h"
#include "ConnectedComponent.h"
#include "Lloyd.h"
#include <cassert>
#include <algorithm>
#include <spdlog/spdlog.h>
#include <chrono>

namespace Patcher {

void create_clusters(
    int n,
    const int* Gp,
    const int* Gi,
    int patch_size,
    const LloydOptions* opt,
    std::vector<int>& vertex_to_cluster)
{
    // Prepare the output
    vertex_to_cluster.assign(n, -1);
    // ── Step 0: Load inputs & defaults ────────────────────────────────
    LloydOptions defaults;
    if (opt == nullptr)
        opt = &defaults;

    // ── Step 1: Connected components ──────────────────────────────────
    // Label every vertex with its CC id so we never cluster across
    // disconnected subgraphs.
    spdlog::info("Step 1: Connected components");
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<int> vertex_to_cc;
    int num_cc = ConnectedComponent::compute(Gp, Gi, n, vertex_to_cc);
    auto end = std::chrono::high_resolution_clock::now();
    double step_1_runtime = std::chrono::duration<double, std::milli>(end - start).count();
    spdlog::info("Connected components runtime: {:.2f} ms", step_1_runtime);


    // ── Step 2: Initialize seeds ──────────────────────────────────
    spdlog::info("Step 2: Initialize seeds");
    start = std::chrono::high_resolution_clock::now();
    //Count the number of vertices in each CC
    std::vector<int> num_vertices_per_cc(num_cc, 0);
    for (int v = 0; v < n; v++)
        num_vertices_per_cc[vertex_to_cc[v]]++;


    std::vector<int> seeds;
    Lloyd::initialize_seeds(vertex_to_cc, num_vertices_per_cc, patch_size, opt->random_seed, seeds);
    end = std::chrono::high_resolution_clock::now();
    double step_2_runtime = std::chrono::duration<double, std::milli>(end - start).count();
    spdlog::info("Initialize seeds runtime: {:.2f} ms", step_2_runtime);

    // ── Step 3: Lloyd-style iterative refinement ──────────────────
    spdlog::info("Step 3: Lloyd-style iterative refinement");
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> node_to_cluster, dist;

    //Init node_to_cluster with the seeds
    node_to_cluster.assign(n, -1);
    for(int i = 0; i < seeds.size(); i++) {
        node_to_cluster[seeds[i]] = i;
    }

    // Outer loop: repeat rounds of Lloyd iterations, then check
    // whether any cluster is oversized.
    for(int iter = 0; iter < opt->lloyd_max_iterations; iter+=opt->lloyd_iters_to_add_seed) {
        std::vector<int> cluster_sizes(seeds.size(), 0);
        // Run a fixed number of Lloyd iterations.
        for (int seed_iter = 0; seed_iter < opt->lloyd_iters_to_add_seed; seed_iter++) {
            // Step 4.1 + 4.2: Multi-source BFS from seeds, assign vertices.
            node_to_cluster.assign(n, -1);
            dist.assign(n, -1);
            for (int i = 0; i < static_cast<int>(seeds.size()); i++) {
                node_to_cluster[seeds[i]] = i;
                dist[seeds[i]] = 0;
            }
            Lloyd::multi_source_bfs(Gp, Gi, n, seeds, node_to_cluster, dist);

            //Compute the maximum size of the clusters
            cluster_sizes.clear();
            cluster_sizes.resize(seeds.size(), 0);
            for(int i = 0; i < n; i++) {
                cluster_sizes[node_to_cluster[i]]++;
            }
            int max_cluster_size = *std::max_element(cluster_sizes.begin(), cluster_sizes.end());
            int min_cluster_size = *std::min_element(cluster_sizes.begin(), cluster_sizes.end());
            if(max_cluster_size <= patch_size) {
                break;
            }
// #ifndef NDEBUG
            spdlog::info("Lloyd iteration {}, max cluster {}, min cluster {}", iter, max_cluster_size, min_cluster_size);
// #endif
            // Step 4.3: Update cluster centers.
            int num_clusters = static_cast<int>(seeds.size());
            Lloyd::update_centers(
                Gp, Gi, n, node_to_cluster,
                num_clusters, opt->random_seed, seeds);
        }

        // Step 4.4: Merge and Split clusters
        Lloyd::merge_and_split_clusters(
            Gp, Gi, n, node_to_cluster,
            cluster_sizes, seeds.size(), patch_size, opt->random_seed, seeds);
    }
    end = std::chrono::high_resolution_clock::now();
    double step_3_runtime = std::chrono::duration<double, std::milli>(end - start).count();
    spdlog::info("Lloyd-style iterative refinement runtime: {:.2f} ms", step_3_runtime);
    // ── Step 5: Final labeling ────────────────────────────────────────
    // Copy the cluster assignments to the output vector.
    vertex_to_cluster = node_to_cluster;
}

} // namespace Patcher
