//
// Morton-code seed visualization on a procedural uniform grid.
// Creates an NxN grid, applies Morton-code seeding, runs a single
// multi-source BFS, and visualizes the initial Voronoi partition
// together with the seed positions and the Morton Z-order curve.
//
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>
#include <cmath>

#include "ConnectedComponent.h"
#include "Lloyd.h"
#include "polyscope_visualizer.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

// ── CLI ─────────────────────────────────────────────────────────────────
struct CLIArgs
{
    int grid_size   = 20;
    int num_patches = 16;
    int bits        = 10;

    CLIArgs(int argc, char* argv[])
    {
        CLI::App app{"Morton-Code Seed Visualization on Uniform Grid"};
        app.add_option("-n,--grid_size", grid_size,
                       "grid dimension N (creates NxN grid, default: 20)");
        app.add_option("-p,--num_patches", num_patches,
                       "number of clusters / patches (default: 16)");
        app.add_option("-b,--bits", bits,
                       "Morton code bits per dimension (default: 10)");

        try {
            app.parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            exit(app.exit(e));
        }
    }
};

// ── Build a 4-connected CSR adjacency for an NxN grid ──────────────────
// Vertex (r, c) has index r * N + c.
// Neighbors: up, down, left, right (when inside bounds).
static void build_grid_csr(int N,
                           std::vector<int>& Gp,
                           std::vector<int>& Gi)
{
    int num_vertices = N * N;
    Gp.resize(num_vertices + 1, 0);

    // First pass: count neighbors per vertex
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            int v = r * N + c;
            int deg = 0;
            if (r > 0)     deg++;  // up
            if (r < N - 1) deg++;  // down
            if (c > 0)     deg++;  // left
            if (c < N - 1) deg++;  // right
            Gp[v + 1] = deg;
        }
    }
    // Prefix sum
    for (int i = 1; i <= num_vertices; i++)
        Gp[i] += Gp[i - 1];

    Gi.resize(Gp[num_vertices]);

    // Second pass: fill column indices (sorted ascending per row)
    std::vector<int> offset(num_vertices, 0);
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            int v = r * N + c;
            // Neighbors added in ascending index order
            if (r > 0)     Gi[Gp[v] + offset[v]++] = (r - 1) * N + c;
            if (c > 0)     Gi[Gp[v] + offset[v]++] = r * N + (c - 1);
            if (c < N - 1) Gi[Gp[v] + offset[v]++] = r * N + (c + 1);
            if (r < N - 1) Gi[Gp[v] + offset[v]++] = (r + 1) * N + c;
        }
    }
}

// ── main ────────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    CLIArgs args(argc, argv);

    const int N = args.grid_size;
    const int num_vertices = N * N;
    const int patch_size = std::max(1, num_vertices / args.num_patches);

    spdlog::info("Grid: {}x{} = {} vertices", N, N, num_vertices);
    spdlog::info("Requested patches: {}, target patch size: {}", args.num_patches, patch_size);

    // ── Generate vertex positions ────────────────────────────────────
    Eigen::MatrixXd V(num_vertices, 3);
    std::vector<std::tuple<double,double,double>> vertex_positions(num_vertices);
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            int v = r * N + c;
            double x = static_cast<double>(c);
            double y = static_cast<double>(r);
            double z = 0.0;
            V(v, 0) = x;
            V(v, 1) = y;
            V(v, 2) = z;
            vertex_positions[v] = {x, y, z};
        }
    }

    // ── Triangulate the grid (2 tris per quad cell) ──────────────────
    int num_cells = (N - 1) * (N - 1);
    Eigen::MatrixXi F(2 * num_cells, 3);
    int fi = 0;
    for (int r = 0; r < N - 1; r++) {
        for (int c = 0; c < N - 1; c++) {
            int v00 = r * N + c;
            int v10 = (r + 1) * N + c;
            int v01 = r * N + (c + 1);
            int v11 = (r + 1) * N + (c + 1);
            // lower-left triangle
            F(fi, 0) = v00; F(fi, 1) = v10; F(fi, 2) = v01;
            fi++;
            // upper-right triangle
            F(fi, 0) = v10; F(fi, 1) = v11; F(fi, 2) = v01;
            fi++;
        }
    }

    // ── Build CSR adjacency ──────────────────────────────────────────
    std::vector<int> Gp, Gi;
    build_grid_csr(N, Gp, Gi);
    spdlog::info("Graph: {} vertices, {} directed edges", num_vertices, (int)Gi.size());

    // ── Connected components (trivially 1 for connected grid) ────────
    std::vector<int> vertex_to_cc;
    int num_cc = ConnectedComponent::compute(Gp.data(), Gi.data(), num_vertices, vertex_to_cc);
    spdlog::info("Connected components: {}", num_cc);

    std::vector<int> num_vertices_per_cc(num_cc, 0);
    for (int v = 0; v < num_vertices; v++)
        num_vertices_per_cc[vertex_to_cc[v]]++;

    // ── Morton-code seeding ──────────────────────────────────────────
    LloydOptions opts;
    opts.seed_selection_method = LloydOptions::SeedSelectionMethod::MORTON_CODE;
    opts.morton_code_bits_per_dimension = args.bits;

    std::vector<int> morton_seeds;
    std::vector<int> sorted_indices;
    auto t0 = std::chrono::high_resolution_clock::now();
    Lloyd::initialize_seeds_morton_code(
        vertex_to_cc, num_vertices_per_cc, vertex_positions,
        patch_size, opts, morton_seeds, sorted_indices);
    auto t1 = std::chrono::high_resolution_clock::now();
    double morton_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    spdlog::info("Morton seeds: {} seeds in {:.2f} ms", morton_seeds.size(), morton_ms);

    // ── Multi-source BFS from Morton seeds ───────────────────────────
    std::vector<int> node_to_cluster(num_vertices, -1);
    std::vector<int> dist(num_vertices, -1);
    for (int i = 0; i < static_cast<int>(morton_seeds.size()); i++) {
        node_to_cluster[morton_seeds[i]] = i;
        dist[morton_seeds[i]] = 0;
    }
    Lloyd::multi_source_bfs(Gp.data(), Gi.data(), num_vertices,
                            morton_seeds, node_to_cluster, dist);

    int num_clusters = static_cast<int>(morton_seeds.size());
    spdlog::info("Initial BFS partition: {} clusters", num_clusters);

    // ── Cluster size stats ───────────────────────────────────────────
    std::vector<int> cluster_sizes(num_clusters, 0);
    for (int label : node_to_cluster)
        if (label >= 0) cluster_sizes[label]++;
    int min_sz = *std::min_element(cluster_sizes.begin(), cluster_sizes.end());
    int max_sz = *std::max_element(cluster_sizes.begin(), cluster_sizes.end());
    double avg_sz = static_cast<double>(num_vertices) / num_clusters;
    spdlog::info("  Min: {}, Max: {}, Avg: {:.1f}", min_sz, max_sz, avg_sz);

    // ── Build Morton Z-order curve (edges in sorted order) ───────────
    // sorted_indices maps local CC index -> position in Morton order.
    // For a single CC, cc_vertices[0] holds the global vertex ids.
    // We reconstruct the global Morton-ordered vertex list.
    std::vector<int> cc_vertices(num_vertices);
    std::iota(cc_vertices.begin(), cc_vertices.end(), 0);

    Eigen::MatrixXd curve_nodes(num_vertices, 3);
    Eigen::MatrixXi curve_edges(num_vertices - 1, 2);
    for (int i = 0; i < num_vertices; i++) {
        int global_v = cc_vertices[sorted_indices[i]];
        curve_nodes.row(i) = V.row(global_v);
    }
    for (int i = 0; i < num_vertices - 1; i++) {
        curve_edges(i, 0) = i;
        curve_edges(i, 1) = i + 1;
    }

    // ── Polyscope visualization ──────────────────────────────────────
    polyscope::options::programName = "Morton Seeds on Grid";
    polyscope::init();

    // Colored surface mesh showing initial partition
    Eigen::MatrixXd cluster_colors = PatchVisualizer::generate_cluster_colors(num_clusters);
    Eigen::MatrixXd vertex_colors(num_vertices, 3);
    for (int i = 0; i < num_vertices; i++)
        vertex_colors.row(i) = cluster_colors.row(node_to_cluster[i]);

    auto* sm = polyscope::registerSurfaceMesh("grid", V, F);
    sm->addVertexColorQuantity("initial partition", vertex_colors)->setEnabled(true);
    sm->setEdgeWidth(0.8);

    // Seed positions as red point cloud
    Eigen::MatrixXd seed_pos(num_clusters, 3);
    for (int i = 0; i < num_clusters; i++)
        seed_pos.row(i) = V.row(morton_seeds[i]);

    auto* pc = polyscope::registerPointCloud("morton seeds", seed_pos);
    pc->setPointRadius(0.008);
    Eigen::MatrixXd seed_colors_mat(num_clusters, 3);
    for (int i = 0; i < num_clusters; i++)
        seed_colors_mat.row(i) << 1.0, 0.0, 0.0;
    pc->addColorQuantity("seed color", seed_colors_mat)->setEnabled(true);

    // Morton Z-order curve
    auto* cn = polyscope::registerCurveNetwork("morton curve", curve_nodes, curve_edges);
    cn->setRadius(0.001);
    cn->setColor({0.2, 0.2, 0.2});
    cn->setEnabled(false);  // off by default, toggle in UI

    // ── ImGui stats panel ────────────────────────────────────────────
    static std::vector<float> s_bin_counts;
    static int s_num_bins, s_bin_width, s_patch_size;
    static int s_min_size, s_max_size, s_num_clusters, s_grid_size;
    static float s_avg_size, s_max_bin;

    s_num_clusters = num_clusters;
    s_min_size = min_sz;
    s_max_size = max_sz;
    s_avg_size = static_cast<float>(avg_sz);
    s_patch_size = patch_size;
    s_grid_size = N;
    s_num_bins = 20;
    s_bin_width = std::max(1, s_patch_size / s_num_bins);

    s_bin_counts.assign(s_num_bins, 0.0f);
    for (int sz : cluster_sizes) {
        int bin = sz / s_bin_width;
        if (bin >= s_num_bins) bin = s_num_bins - 1;
        s_bin_counts[bin] += 1.0f;
    }
    s_max_bin = *std::max_element(s_bin_counts.begin(), s_bin_counts.end());

    polyscope::state::userCallback = []() {
        ImGui::Begin("Morton Seeds on Grid -- Stats");
        ImGui::Text("Grid: %dx%d = %d vertices",
                     s_grid_size, s_grid_size, s_grid_size * s_grid_size);
        ImGui::Text("Seeds (clusters): %d", s_num_clusters);
        ImGui::Text("Min: %d  |  Max: %d  |  Avg: %.1f  |  Target: %d",
                     s_min_size, s_max_size, s_avg_size, s_patch_size);
        ImGui::Separator();
        ImGui::PlotHistogram("##sizes", s_bin_counts.data(), s_num_bins,
                             0, nullptr, 0.0f, s_max_bin * 1.1f,
                             ImVec2(400, 200));
        ImGui::Text("Bin width: %d  |  Range: [0, %d]  |  Bins: %d",
                     s_bin_width, s_patch_size, s_num_bins);
        ImGui::End();
    };

    polyscope::show();
    return 0;
}
