#include <vector>
#include <queue>
#include <cmath>
#include <cstddef>

#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TObject.h>
#include <TString.h>
#include <TBox.h>
#include <TMarker.h>
#include <TF1.h>
#include <TLatex.h>

#include "H2GCROC_Common.hxx"


struct ClusterResult {
    std::vector<int> xs;
    std::vector<int> ys;

    double sumE = 0.0;          // sum of KEPT cells
    double sumE_full = 0.0;     // sum of ALL cells in connected component (optional)

    double x_cm_cell = 0.0;
    double y_cm_cell = 0.0;
    double x_cm_cm = 0.0;
    double y_cm_cm = 0.0;
};

inline bool in_range(int x, int y, int NX, int NY) {
    return (x >= 0 && x < NX && y >= 0 && y < NY);
}


struct CellHit {
    int x;
    int y;
    double e;
};


std::vector<ClusterResult> cluster_uniform_grid_topcells(
    const std::vector<int>& channel_x_positions,
    const std::vector<int>& channel_y_positions,
    const std::vector<double>& channel_value_positions,
    int NX = 16,
    int NY = 12,
    double pitch_cm = 1.0,
    double seed_thr = 80.0,          // start cluster only if E >= seed_thr
    double neighbor_thr = 30.0,      // grow cluster if neighbor E >= neighbor_thr
    bool use_8_neighbors = true,
    double min_cluster_sum = 0.0,    // filter by KEPT sumE
    int top_cells = 30,              // keep top K cells by energy per cluster
    int max_clusters = -1            // keep top N clusters by sumE; -1 = no limit
) {
    // 1) Build dense energy grid
    std::vector<double> E(NX * NY, 0.0);
    const size_t N = channel_value_positions.size();
    if (channel_x_positions.size() != N || channel_y_positions.size() != N) return {};

    for (size_t i = 0; i < N; ++i) {
        int x = channel_x_positions[i];
        int y = channel_y_positions[i];
        if (!in_range(x, y, NX, NY)) continue;
        E[y * NX + x] += channel_value_positions[i]; // accumulate duplicates
    }

    // 2) Visited map (each cell belongs to at most one cluster)
    std::vector<unsigned char> vis(NX * NY, 0);

    // Neighbor offsets
    static const int dx4[4] = {+1, -1,  0,  0};
    static const int dy4[4] = { 0,  0, +1, -1};
    static const int dx8[8] = {+1, -1,  0,  0, +1, +1, -1, -1};
    static const int dy8[8] = { 0,  0, +1, -1, +1, -1, +1, -1};

    const int* dx = use_8_neighbors ? dx8 : dx4;
    const int* dy = use_8_neighbors ? dy8 : dy4;
    const int nnb = use_8_neighbors ? 8 : 4;

    std::vector<ClusterResult> clusters;
    clusters.reserve(32);

    // 3) Find connected components (dual threshold)
    for (int y0 = 0; y0 < NY; ++y0) {
        for (int x0 = 0; x0 < NX; ++x0) {
            const int idx0 = y0 * NX + x0;
            if (vis[idx0]) continue;
            if (E[idx0] < seed_thr) continue; // seed condition

            std::queue<std::pair<int,int>> q;
            q.push({x0, y0});
            vis[idx0] = 1;

            std::vector<CellHit> hits;
            hits.reserve(64);

            double sum_full = 0.0;

            while (!q.empty()) {
                auto [x, y] = q.front();
                q.pop();

                const int idx = y * NX + x;
                const double w = E[idx];

                hits.push_back({x, y, w});
                sum_full += w;

                for (int k = 0; k < nnb; ++k) {
                    int xn = x + dx[k];
                    int yn = y + dy[k];
                    if (!in_range(xn, yn, NX, NY)) continue;
                    const int idn = yn * NX + xn;
                    if (vis[idn]) continue;
                    if (E[idn] < neighbor_thr) continue; // grow condition

                    vis[idn] = 1;
                    q.push({xn, yn});
                }
            }

            if (hits.empty()) continue;

            // 4) Keep top K cells by energy within this connected component
            int K = top_cells;
            if (K <= 0) K = (int)hits.size(); // <=0 means keep all
            if ((int)hits.size() > K) {
                // nth_element: O(n) average, very fast
                std::nth_element(
                    hits.begin(), hits.begin() + K, hits.end(),
                    [](const CellHit& a, const CellHit& b){ return a.e > b.e; }
                );
                hits.resize(K);
                // optional: sort kept cells (nice for debug/plots)
                std::sort(
                    hits.begin(), hits.end(),
                    [](const CellHit& a, const CellHit& b){ return a.e > b.e; }
                );
            } else {
                std::sort(
                    hits.begin(), hits.end(),
                    [](const CellHit& a, const CellHit& b){ return a.e > b.e; }
                );
            }

            // 5) Build ClusterResult from KEPT cells and recompute COM
            ClusterResult c;
            c.sumE_full = sum_full;

            double sw = 0.0;
            double xw = 0.0, yw = 0.0;

            c.xs.reserve(hits.size());
            c.ys.reserve(hits.size());

            for (const auto& h : hits) {
                c.xs.push_back(h.x);
                c.ys.push_back(h.y);
                c.sumE += h.e;

                sw += h.e;
                xw += h.e * double(h.x);
                yw += h.e * double(h.y);
            }

            // filter by KEPT sumE (more meaningful when using top-cells)
            if (c.sumE < min_cluster_sum) continue;

            if (sw > 0.0) {
                c.x_cm_cell = xw / sw;
                c.y_cm_cell = yw / sw;
            } else {
                c.x_cm_cell = x0;
                c.y_cm_cell = y0;
            }

            c.x_cm_cm = (c.x_cm_cell + 0.5) * pitch_cm;
            c.y_cm_cm = (c.y_cm_cell + 0.5) * pitch_cm;

            clusters.push_back(std::move(c));
        }
    }

    // 6) Limit total number of clusters (keep the most energetic clusters)
    if (!clusters.empty()) {
        std::sort(
            clusters.begin(), clusters.end(),
            [](const ClusterResult& a, const ClusterResult& b){
                return a.sumE > b.sumE; // compare by KEPT sumE
            }
        );
        if (max_clusters > 0 && (int)clusters.size() > max_clusters) {
            clusters.resize(max_clusters);
        }
    }

    return clusters;
}