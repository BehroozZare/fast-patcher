#pragma once

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <array>

#include <Eigen/Core>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "imgui.h"

namespace PatchVisualizer {

// Convert HSL to RGB (h in [0,360), s,l in [0,1]) → RGB in [0,1]
inline std::array<double, 3> hsl_to_rgb(double h, double s, double l) {
    double c = (1.0 - std::abs(2.0 * l - 1.0)) * s;
    double hp = h / 60.0;
    double x = c * (1.0 - std::abs(std::fmod(hp, 2.0) - 1.0));
    double r1 = 0, g1 = 0, b1 = 0;

    if      (hp < 1) { r1 = c; g1 = x; }
    else if (hp < 2) { r1 = x; g1 = c; }
    else if (hp < 3) { g1 = c; b1 = x; }
    else if (hp < 4) { g1 = x; b1 = c; }
    else if (hp < 5) { r1 = x; b1 = c; }
    else             { r1 = c; b1 = x; }

    double m = l - c / 2.0;
    return {r1 + m, g1 + m, b1 + m};
}

// Generate `n` visually distinct colors using the golden-angle hue spacing.
// The golden angle (~137.5 deg) maximises hue separation between consecutive
// indices, which is important because BFS tends to assign nearby cluster IDs
// to spatially adjacent patches.  Cycling through several saturation/lightness
// bands adds extra contrast.
inline Eigen::MatrixXd generate_cluster_colors(int n, int seed = 42) {
    Eigen::MatrixXd colors(n, 3);
    if (n == 0) return colors;
    if (n == 1) {
        colors.row(0) << 0.267, 0.467, 0.667; // pleasant blue
        return colors;
    }

    constexpr double golden_angle = 137.50776405; // degrees
    constexpr int    num_bands    = 4;
    constexpr double sat_bands[num_bands] = {0.72, 0.58, 0.85, 0.65};
    constexpr double lit_bands[num_bands] = {0.55, 0.42, 0.63, 0.48};

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> hue_jitter(-5.0, 5.0);

    for (int i = 0; i < n; ++i) {
        double hue = std::fmod(i * golden_angle + hue_jitter(rng) + 360.0, 360.0);
        int    band = i % num_bands;
        auto   rgb  = hsl_to_rgb(hue, sat_bands[band], lit_bands[band]);
        colors.row(i) << rgb[0], rgb[1], rgb[2];
    }
    return colors;
}

// Visualize mesh patches using Polyscope, with a patch size histogram.
//
//   V                - Nx3 vertex positions
//   F                - element connectivity (Fx3 triangles or Fx4 tetrahedra)
//   vertex_to_cluster - per-vertex cluster label (size N)
//   is_tet_mesh      - true for volumetric (tet) meshes
//   patch_size       - target patch size (used for histogram binning)
//
inline void visualize_patches(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const std::vector<int>& vertex_to_cluster,
    bool is_tet_mesh,
    int patch_size)
{
    // --- Polyscope init ---
    polyscope::options::programName = "Fast Patcher – Patch Viewer";
    polyscope::init();

    int num_clusters = *std::max_element(vertex_to_cluster.begin(),
                                         vertex_to_cluster.end()) + 1;

    // Build per-vertex color matrix from cluster labels
    Eigen::MatrixXd cluster_colors = generate_cluster_colors(num_clusters);
    int num_verts = static_cast<int>(V.rows());
    Eigen::MatrixXd vertex_colors(num_verts, 3);
    for (int i = 0; i < num_verts; ++i) {
        vertex_colors.row(i) = cluster_colors.row(vertex_to_cluster[i]);
    }

    if (is_tet_mesh) {
        auto* vm = polyscope::registerTetMesh("mesh", V, F);
        vm->addVertexColorQuantity("patch", vertex_colors)->setEnabled(true);
    } else {
        auto* sm = polyscope::registerSurfaceMesh("mesh", V, F);
        sm->addVertexColorQuantity("patch", vertex_colors)->setEnabled(true);
    }

    // --- Compute cluster sizes and bin into histogram ---
    static std::vector<float> s_bin_counts;
    static int s_num_bins;
    static float s_max_bin;
    static int s_min_size, s_max_size, s_num_clusters, s_bin_width, s_patch_size;
    static float s_avg_size;

    std::vector<int> cluster_sizes(num_clusters, 0);
    for (int label : vertex_to_cluster)
        cluster_sizes[label]++;

    s_num_clusters = num_clusters;
    s_min_size = *std::min_element(cluster_sizes.begin(), cluster_sizes.end());
    s_max_size = *std::max_element(cluster_sizes.begin(), cluster_sizes.end());
    s_avg_size = static_cast<float>(std::accumulate(cluster_sizes.begin(), cluster_sizes.end(), 0))
                 / num_clusters;
    s_patch_size = patch_size;
    s_num_bins = 20;
    s_bin_width = std::max(1, patch_size / s_num_bins);

    s_bin_counts.assign(s_num_bins, 0.0f);
    for (int sz : cluster_sizes) {
        int bin = sz / s_bin_width;
        if (bin >= s_num_bins) bin = s_num_bins - 1;
        s_bin_counts[bin] += 1.0f;
    }
    s_max_bin = *std::max_element(s_bin_counts.begin(), s_bin_counts.end());

    polyscope::state::userCallback = []() {
        ImGui::Begin("Patch Size Distribution");
        ImGui::Text("Clusters: %d  |  Min: %d  |  Max: %d  |  Avg: %.1f",
                     s_num_clusters, s_min_size, s_max_size, s_avg_size);
        ImGui::Separator();
        ImGui::PlotHistogram("##sizes", s_bin_counts.data(), s_num_bins,
                             0, nullptr, 0.0f, s_max_bin * 1.1f,
                             ImVec2(400, 200));
        ImGui::Text("Bin width: %d  |  Range: [0, %d]  |  Bins: %d",
                     s_bin_width, s_patch_size, s_num_bins);
        ImGui::End();
    };

    polyscope::show();
}

} // namespace PatchVisualizer
