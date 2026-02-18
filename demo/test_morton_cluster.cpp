//
// Morton-code seed visualization demo.
// Loads a mesh, computes initial seeds via Morton code (and optionally
// random), runs a single BFS from those seeds, and visualizes the
// initial Voronoi partition together with the seed positions.
//
#include <igl/read_triangle_mesh.h>
#include <igl/readMESH.h>
#include <igl/cotmatrix.h>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <chrono>
#include <unsupported/Eigen/SparseExtra>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <numeric>
#include <tuple>

#include "remove_diagonal.h"
#include "ConnectedComponent.h"
#include "Lloyd.h"
#include "polyscope_visualizer.h"
#include "polyscope/point_cloud.h"

struct CLIArgs
{
    std::string input_mesh;
    int num_patches = 512;

    CLIArgs(int argc, char* argv[])
    {
        CLI::App app{"Morton-Code Seed Visualization"};
        app.add_option("-i,--input", input_mesh,
                       "input mesh file (.mesh for tet, .obj/.off/.ply/.stl for triangle)")
            ->required();
        app.add_option("-n,--num_patches", num_patches, "number of patches (default: 512)");

        try {
            app.parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            exit(app.exit(e));
        }
    }
};

std::string get_extension(const std::string& path) {
    std::filesystem::path p(path);
    std::string ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

int main(int argc, char* argv[])
{
    CLIArgs args(argc, argv);

    std::string ext = get_extension(args.input_mesh);
    spdlog::info("Loading mesh from: {}", args.input_mesh);
    spdlog::info("Number of patches: {}", args.num_patches);

    // ── Load mesh ────────────────────────────────────────────────────
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    bool is_tet_mesh = (ext == ".mesh");

    if (is_tet_mesh) {
        Eigen::MatrixXi OF;
        if (!igl::readMESH(args.input_mesh, V, F, OF)) {
            spdlog::error("Failed to read tet mesh: {}", args.input_mesh);
            return 1;
        }
        spdlog::info("Loaded tet mesh: {} vertices, {} tetrahedra", V.rows(), F.rows());
    } else if (ext == ".obj" || ext == ".off" || ext == ".ply" || ext == ".stl") {
        if (!igl::read_triangle_mesh(args.input_mesh, V, F)) {
            spdlog::error("Failed to read triangle mesh: {}", args.input_mesh);
            return 1;
        }
        spdlog::info("Loaded triangle mesh: {} vertices, {} triangles", V.rows(), F.rows());
    } else {
        spdlog::error("Unsupported format: {}", ext);
        return 1;
    }

    int num_vertices = static_cast<int>(V.rows());
    int patch_size = num_vertices / args.num_patches;

    // ── Build CSR adjacency from cotangent Laplacian ─────────────────
    spdlog::info("Computing cotangent matrix...");
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    std::vector<int> Gp, Gi;
    Patcher::remove_diagonal(
        L.rows(), L.outerIndexPtr(), L.innerIndexPtr(), Gp, Gi);
    spdlog::info("Graph: {} vertices, {} edges", num_vertices, Gi.size());

    // ── Connected components ─────────────────────────────────────────
    std::vector<int> vertex_to_cc;
    int num_cc = ConnectedComponent::compute(Gp.data(), Gi.data(), num_vertices, vertex_to_cc);
    spdlog::info("Connected components: {}", num_cc);

    std::vector<int> num_vertices_per_cc(num_cc, 0);
    for (int v = 0; v < num_vertices; v++)
        num_vertices_per_cc[vertex_to_cc[v]]++;

    // ── Vertex positions for Morton code ─────────────────────────────
    std::vector<std::tuple<double,double,double>> vertex_positions(num_vertices);
    for (int i = 0; i < num_vertices; i++)
        vertex_positions[i] = {V(i,0), V(i,1), V(i,2)};

    // ── Initialize seeds via Morton code ─────────────────────────────
    LloydOptions opts;
    opts.seed_selection_method = LloydOptions::SeedSelectionMethod::MORTON_CODE;

    std::vector<int> morton_seeds;
    std::vector<int> sorted_indices;
    auto t0 = std::chrono::high_resolution_clock::now();
    Lloyd::initialize_seeds_morton_code(
        vertex_to_cc, num_vertices_per_cc, vertex_positions,
        patch_size, opts, morton_seeds, sorted_indices);
    auto t1 = std::chrono::high_resolution_clock::now();
    double morton_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    spdlog::info("Morton seeds: {} seeds in {:.2f} ms", morton_seeds.size(), morton_ms);

    // ── Single BFS from Morton seeds to get initial partition ────────
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
        cluster_sizes[label]++;
    int min_sz = *std::min_element(cluster_sizes.begin(), cluster_sizes.end());
    int max_sz = *std::max_element(cluster_sizes.begin(), cluster_sizes.end());
    double avg_sz = static_cast<double>(num_vertices) / num_clusters;
    spdlog::info("  Min: {}, Max: {}, Avg: {:.1f}", min_sz, max_sz, avg_sz);

    // ── Polyscope visualization ──────────────────────────────────────
    polyscope::options::programName = "Morton Seed Viewer";
    polyscope::init();

    // Colored mesh showing the initial Voronoi partition
    Eigen::MatrixXd cluster_colors = PatchVisualizer::generate_cluster_colors(num_clusters);
    Eigen::MatrixXd vertex_colors(num_vertices, 3);
    for (int i = 0; i < num_vertices; i++)
        vertex_colors.row(i) = cluster_colors.row(node_to_cluster[i]);

    if (is_tet_mesh) {
        auto* vm = polyscope::registerTetMesh("mesh", V, F);
        vm->addVertexColorQuantity("initial partition", vertex_colors)->setEnabled(true);
        vm->setEdgeWidth(0.8);
    } else {
        auto* sm = polyscope::registerSurfaceMesh("mesh", V, F);
        sm->addVertexColorQuantity("initial partition", vertex_colors)->setEnabled(true);
        sm->setEdgeWidth(0.8);
    }

    // Seed positions as a point cloud
    Eigen::MatrixXd seed_pos(num_clusters, 3);
    for (int i = 0; i < num_clusters; i++) {
        int v = morton_seeds[i];
        seed_pos.row(i) = V.row(v);
    }
    auto* pc = polyscope::registerPointCloud("morton seeds", seed_pos);
    pc->setPointRadius(0.005);
    Eigen::MatrixXd seed_colors(num_clusters, 3);
    for (int i = 0; i < num_clusters; i++)
        seed_colors.row(i) << 1.0, 0.0, 0.0;
    pc->addColorQuantity("seed color", seed_colors)->setEnabled(true);

    // ImGui panel with cluster-size histogram
    static std::vector<float> s_bin_counts;
    static int s_num_bins, s_bin_width, s_patch_size;
    static int s_min_size, s_max_size, s_num_clusters;
    static float s_avg_size, s_max_bin;

    s_num_clusters = num_clusters;
    s_min_size = min_sz;
    s_max_size = max_sz;
    s_avg_size = static_cast<float>(avg_sz);
    s_patch_size = patch_size;
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
        ImGui::Begin("Morton Seed -- Initial Partition Stats");
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
