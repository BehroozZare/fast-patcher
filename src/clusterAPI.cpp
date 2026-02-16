#include "clusterAPI.h"
#include "ConnectedComponent.h"
#include "Lloyd.h"
#include "bit_array.h"
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
    std::vector<int>& vertex_to_cluster,
    const std::vector<std::tuple<double, double, double>>* vertex_positions)
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
    if (vertex_positions != nullptr &&
        opt->seed_selection_method == LloydOptions::SeedSelectionMethod::MORTON_CODE) {
        spdlog::info("Seed selection: Morton code");
        Lloyd::initialize_seeds_morton_code(
            vertex_to_cc, num_vertices_per_cc, *vertex_positions,
            patch_size, *opt, seeds);
    } else {
        if (opt->seed_selection_method == LloydOptions::SeedSelectionMethod::MORTON_CODE
            && vertex_positions == nullptr) {
            spdlog::warn("Morton code requested but no vertex positions provided; falling back to random");
        }
        spdlog::info("Seed selection: Random");
        Lloyd::initialize_seeds_random(
            vertex_to_cc, num_vertices_per_cc, patch_size, *opt, seeds);
    }
    end = std::chrono::high_resolution_clock::now();
    double step_2_runtime = std::chrono::duration<double, std::milli>(end - start).count();
    spdlog::info("Initialize seeds runtime: {:.2f} ms", step_2_runtime);

    // ── Step 3: Lloyd-style iterative refinement ──────────────────
    spdlog::info("Step 3: Lloyd-style iterative refinement");
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> node_to_cluster, dist;
    dist.assign(n, -1);

    //Init node_to_cluster with the seeds
    node_to_cluster.assign(n, -1);
    for(int i = 0; i < seeds.size(); i++) {
        node_to_cluster[seeds[i]] = i;
    }

    // Active mask: 1 = vertex participates in BFS, 0 = frozen.
    // Initially all vertices are active (no clustering yet).
    Lloyd::BitArray active_mask(n);
    active_mask.set_all();

    // Outer loop: repeat rounds of Lloyd iterations, then check
    // whether any cluster is oversized.
    for(int iter = 0; iter < opt->lloyd_max_iterations; iter+=opt->lloyd_iters_to_add_seed) {
        std::vector<int> cluster_sizes(seeds.size(), 0);
        int max_cluster_size = 0;
        int min_cluster_size = 0;
        // Run a fixed number of Lloyd iterations.
        for (int seed_iter = 0; seed_iter < opt->lloyd_iters_to_add_seed; seed_iter++) {
            auto iter_start = std::chrono::high_resolution_clock::now();

            // Step 4.1 + 4.2: Multi-source BFS from seeds, assign vertices.
            // Only reset active vertices; frozen ones keep their assignments.
            for (int v = 0; v < n; v++) {
                if (active_mask.get(v)) {
                    node_to_cluster[v] = -1;
                    dist[v] = -1;
                }
            }
            for (int i = 0; i < static_cast<int>(seeds.size()); i++) {
                if (active_mask.get(seeds[i])) {
                    node_to_cluster[seeds[i]] = i;
                    dist[seeds[i]] = 0;
                }
            }
            Lloyd::multi_source_bfs(Gp, Gi, n, seeds, node_to_cluster, dist, &active_mask);

            auto bfs_end = std::chrono::high_resolution_clock::now();
            double bfs_ms = std::chrono::duration<double, std::milli>(bfs_end - iter_start).count();

            //Compute the maximum size of the clusters
            cluster_sizes.clear();
            cluster_sizes.resize(seeds.size(), 0);
            for(int i = 0; i < n; i++) {
                cluster_sizes[node_to_cluster[i]]++;
            }
            max_cluster_size = *std::max_element(cluster_sizes.begin(), cluster_sizes.end());
            min_cluster_size = *std::min_element(cluster_sizes.begin(), cluster_sizes.end());
            if(max_cluster_size <= patch_size) {
                spdlog::info("Iteration [{},{}]: converged (max cluster {} <= patch_size {}), BFS {:.2f} ms",
                             iter, seed_iter, max_cluster_size, patch_size, bfs_ms);
                // break;
            }

            // Step 4.3: Update cluster centers.
            int num_clusters = static_cast<int>(seeds.size());
            Lloyd::update_centers(
                Gp, Gi, n, node_to_cluster,
                num_clusters, *opt, seeds, &active_mask);

            auto iter_end = std::chrono::high_resolution_clock::now();
            double center_ms = std::chrono::duration<double, std::milli>(iter_end - bfs_end).count();
            double iter_ms = std::chrono::duration<double, std::milli>(iter_end - iter_start).count();

            spdlog::info("Iteration [{},{}]: #seeds={}, max_cluster={}, min_cluster={}, BFS {:.2f} ms, center_update {:.2f} ms, total {:.2f} ms",
                         iter, seed_iter, seeds.size(), max_cluster_size, min_cluster_size, bfs_ms, center_ms, iter_ms);
        }
        if(max_cluster_size <= patch_size) {
            break;
        }

        // Step 4.4: Merge and Split clusters
        auto split_start = std::chrono::high_resolution_clock::now();
        Lloyd::merge_and_split_clusters(
            Gp, Gi, n, node_to_cluster,
            cluster_sizes, seeds.size(), patch_size, *opt, seeds);
        auto split_end = std::chrono::high_resolution_clock::now();
        double split_ms = std::chrono::duration<double, std::milli>(split_end - split_start).count();
        spdlog::info("Outer iteration {}: merge_and_split {:.2f} ms, #seeds after split={}",
                     iter, split_ms, seeds.size());

        // Rebuild the active mask: only vertices in oversized clusters
        // are active for the next round of Lloyd iterations.
        active_mask.clear_all();
        for (int v = 0; v < n; v++) {
            if (cluster_sizes[node_to_cluster[v]] >= patch_size)
                active_mask.set(v);
        }
    }
    end = std::chrono::high_resolution_clock::now();
    double step_3_runtime = std::chrono::duration<double, std::milli>(end - start).count();
    spdlog::info("Lloyd-style iterative refinement runtime: {:.2f} ms", step_3_runtime);
    // ── Step 5: Final labeling ────────────────────────────────────────
    // Copy the cluster assignments to the output vector.
    vertex_to_cluster = node_to_cluster;
}

} // namespace Patcher
