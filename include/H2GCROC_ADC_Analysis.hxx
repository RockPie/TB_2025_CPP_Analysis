#include <vector>
#include <cmath>
#include <limits>
#include <utility>  // std::swap
#include <functional>
#include <array>
#include <cassert>
#include <string>

#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TObject.h>
#include <TString.h>
#include "H2GCROC_Common.hxx"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

// the average of the smallest 2 samples among the first 3 samples
inline void pedestal_subtraction_2minoffirst3(std::vector<int>& data, double& pedestal) {
    if (data.size() < 3) {
        pedestal = 0.0;
        return;
    }
    std::vector<int> first3(data.begin(), data.begin() + 3);
    std::sort(first3.begin(), first3.end());
    pedestal = (double(first3[0]) + double(first3[1])) / 2.0;
    for (auto& val : data) {
        val -= pedestal;
        if (val < 0) {
            val = 0;
        }
    }
}

inline double pedestal_median_of_first3(std::vector<int>& data) {
    if (data.size() < 3) {
        return 0.0;
    }
    std::vector<int> first3(data.begin(), data.begin() + 3);
    std::sort(first3.begin(), first3.end());
    return double(first3[1]);
}

inline double pedestal_average_of_first3(std::vector<int>& data) {
    if (data.size() < 3) {
        return 0.0;
    }
    return (double(data[0]) + double(data[1]) + double(data[2])) / 3.0;
}

// sigma clipping method to calculate pedestal sigma
// k: clipping threshold in unit of sigma
// n_iter_max: maximum number of iterations
inline double pedestal_sigma_clipping(const std::vector<int>& data, double k, int n_iter_max) {
    using std::size_t;
    if (data.empty() || k <= 0.0 || n_iter_max <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Working buffers (reuse capacity to avoid reallocations)
    std::vector<double> cur;  cur.reserve(data.size());
    std::vector<double> tmp;  tmp.reserve(data.size());
    for (int v : data) cur.push_back(static_cast<double>(v));

    auto mean_and_sigma = [](const std::vector<double>& x) -> std::pair<double,double> {
        // Welford's algorithm
        double mean = 0.0;
        double M2   = 0.0;
        size_t n = 0;
        for (double v : x) {
            ++n;
            double delta = v - mean;
            mean += delta / static_cast<double>(n);
            double delta2 = v - mean;
            M2 += delta * delta2;
        }
        if (n < 2) return { (n ? mean : std::numeric_limits<double>::quiet_NaN()), 0.0 };
        double var = M2 / static_cast<double>(n); // population variance for pedestal noise
        return { mean, std::sqrt(var) };
    };

    for (int it = 0; it < n_iter_max; ++it) {
        auto [mu, sigma] = mean_and_sigma(cur);
        if (!std::isfinite(mu)) return mu;          // handle NaN (n==0)
        if (sigma == 0.0)        return mu;          // all equal ⇒ converged

        const double lo = mu - k * sigma;
        const double hi = mu + k * sigma;

        tmp.clear();
        tmp.reserve(cur.size()); // ensure capacity, no reallocation
        for (double v : cur) {
            if (v >= lo && v <= hi) tmp.push_back(v);
        }

        if (tmp.size() == cur.size() || tmp.empty()) {
            // No change or all rejected ⇒ stop; use current mean
            return mu;
        }

        cur.swap(tmp); // keep filtered set and reuse buffers
    }

    // Recompute final mean after last iteration
    return mean_and_sigma(cur).first;
}

inline void draw_on_pad(TPad* pad, TObject* obj, bool minimalist_axis, bool th2_logz)
{
    if (!pad || !obj) return;
    pad->cd();
    pad->SetMargin(0, 0, 0, 0);
    pad->SetBorderMode(0);
    pad->SetFrameBorderMode(0);
    pad->SetFillStyle(0);

    if (auto* h2 = dynamic_cast<TH2*>(obj)) {
        h2->SetStats(0);
        if (minimalist_axis) {
            h2->GetXaxis()->SetLabelSize(0);
            h2->GetYaxis()->SetLabelSize(0);
            h2->GetXaxis()->SetTitleSize(0);
            h2->GetYaxis()->SetTitleSize(0);
            h2->GetXaxis()->SetTickLength(0);
            h2->GetYaxis()->SetTickLength(0);
        }
        pad->SetLogz(th2_logz ? 1 : 0);
        h2->Draw("COLZ");
        return;
    }
    if (auto* h1 = dynamic_cast<TH1*>(obj)) {
        h1->SetStats(0);
        if (minimalist_axis) {
            h1->GetXaxis()->SetLabelSize(0);
            h1->GetYaxis()->SetLabelSize(0);
            h1->GetXaxis()->SetTitleSize(0);
            h1->GetYaxis()->SetTitleSize(0);
            h1->GetXaxis()->SetTickLength(0);
            h1->GetYaxis()->SetTickLength(0);
        }
        pad->SetLogz(0);
        h1->Draw("HIST");
        return;
    }
    if (auto* gr = dynamic_cast<TGraphErrors*>(obj)) {
        if (minimalist_axis) {
            gr->GetXaxis()->SetLabelSize(0);
            gr->GetYaxis()->SetLabelSize(0);
            gr->GetXaxis()->SetTitleSize(0);
            gr->GetYaxis()->SetTitleSize(0);
            gr->GetXaxis()->SetTickLength(0);
            gr->GetYaxis()->SetTickLength(0);
        }
        pad->SetLogz(0);
        gr->Draw("APL");
        return;
    }
    pad->SetLogz(0);
    obj->Draw();
}

// ============ 无缝网格 ============
inline void build_gapless_grid(TCanvas& canvas, int NX, int NY,
                               std::vector<std::vector<TPad*>>& pads)
{
    assert(NX > 0 && NY > 0);
    pads.assign(NY, std::vector<TPad*>(NX, nullptr));

    // 统一外观（注意需要 #include <TStyle.h>）
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);

    const double dx  = 1.0 / NX;
    const double dy  = 1.0 / NY;
    const double eps = 1e-6;

    for (int col = 0; col < NX; ++col) {
        for (int row = 0; row < NY; ++row) {
            const double x1 = col * dx;
            const double x2 = (col + 1) * dx;
            const double y1 = 1.0 - (row + 1) * dy;
            const double y2 = 1.0 - row * dy;

            const double xl = (col == 0)    ? 0.0 : x1 - eps;
            const double xr = (col == NX-1) ? 1.0 : x2 + eps;
            const double yb = (row == NY-1) ? 0.0 : y1 - eps;
            const double yt = (row == 0)    ? 1.0 : y2 + eps;

            TString nm = Form("pad_r%d_c%d", row, col);
            auto* p = new TPad(nm, nm, xl, yb, xr, yt);
            p->SetMargin(0, 0, 0, 0);
            p->SetBorderMode(0);
            p->SetFrameBorderMode(0);
            p->SetFillStyle(0);
            p->Draw();
            pads[row][col] = p;
        }
    }
}

// ============ 拓扑封装（固定顺序） ============
// 用 chan2pad 把 “(vldb, channel) 的线性索引” 映射到 “pad 的线性索引”。
// 线性索引定义：
//   chan_linear = vldb * channels_per_vldb + channel
//   pad_linear  = row * NX + col      （row-major）
//
// 这样 items[i] 就永远画到 pad_linear = chan2pad[i] 的 pad 上。
struct MosaicTopology {
    int NX{16};
    int NY{12};
    int vldb_number{0};
    int channels_per_vldb{0};
    bool reverse_row{true};     // 是否最后做行反转（与你当前图片一致）
    bool minimalist_axis{true};
    bool th2_logz{true};

    // 固定顺序 LUT：长度应为 vldb_number * channels_per_vldb。
    // 若某通道不想画，置为 -1。
    std::vector<int> chan2pad;  // pad 的线性索引（0..NX*NY-1）或 -1

    inline int pad_linear_of(int col, int row) const {
        if (reverse_row) row = (NY - 1) - row;
        return row * NX + col;
    }

    // 安全性检查
    bool valid() const {
        const int need = vldb_number * channels_per_vldb;
        if (NX <= 0 || NY <= 0 || vldb_number <= 0 || channels_per_vldb <= 0) return false;
        if ((int)chan2pad.size() != need) return false;
        return true;
    }
};

// ============ 绘制：严格按拓扑顺序 ============
inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TObject*>& items,
                              const MosaicTopology& topo)
{
    if (!topo.valid()) return;

    std::vector<std::vector<TPad*>> pads;
    build_gapless_grid(canvas, topo.NX, topo.NY, pads);

    const int pad_count = topo.NX * topo.NY;

    // 预构建 pad* 的线性数组，O(1) 访问
    std::vector<TPad*> pad_linear(pad_count, nullptr);
    for (int r = 0; r < topo.NY; ++r) {
        for (int c = 0; c < topo.NX; ++c) {
            int pr = r;
            if (topo.reverse_row) pr = (topo.NY - 1) - r; // 因为 pad_linear 用的就是最终的行序
            pad_linear[pr * topo.NX + c] = pads[r][c];
        }
    }

    const int n_items_expected = topo.vldb_number * topo.channels_per_vldb;
    const int n_items = std::min<int>( (int)items.size(), n_items_expected );

    for (int i = 0; i < n_items; ++i) {
        TObject* obj = items[i];
        if (!obj) continue;

        int pad_lnr = topo.chan2pad[i];
        if (pad_lnr < 0 || pad_lnr >= pad_count) continue;

        TPad* pad = pad_linear[pad_lnr];
        if (!pad) continue;

        draw_on_pad(pad, obj, topo.minimalist_axis, topo.th2_logz);
    }

    canvas.Update();
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TH2D*>& h2,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(h2.size());
    for (auto* p : h2) objs.push_back(static_cast<TObject*>(p));
    draw_mosaic_fixed(canvas, objs, topo);
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TH1D*>& h1,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(h1.size());
    for (auto* p : h1) objs.push_back(static_cast<TObject*>(p));
    // 对 TH1*，如果你不想 LogZ，在构建 topo 时把 th2_logz=false 即可（只影响 TH2）
    draw_mosaic_fixed(canvas, objs, topo);
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TGraphErrors*>& gr,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(gr.size());
    for (auto* p : gr) objs.push_back(static_cast<TObject*>(p));
    draw_mosaic_fixed(canvas, objs, topo);
}

std::vector<int> build_chan2pad_LUT(
    int vldb_number,
    int channels_per_vldb,
    int NX, int NY,
    int board_cols, int board_rows,
    const json& sipm_board,          // 通道到(row,col)
    const json& board_loc,           // 板到(row,col)
    const json& board_rotation,      // 板是否旋转
    const json& board_flip           // 板是否翻转（可不使用）
) {
    std::vector<int> lut(vldb_number * channels_per_vldb, -1);

    for (int vldb_id = 0; vldb_id < vldb_number; ++vldb_id) {
        for (int ch = 0; ch < channels_per_vldb; ++ch) {

            // ---------------- 通道拓扑逻辑 ----------------
            int channel_index_in_h2g       = ch % 76;
            int channel_index_in_h2g_half  = channel_index_in_h2g % 38;
            int asic_id                    = ch / 76;
            int half_id                    = channel_index_in_h2g / 38;

            if (channel_index_in_h2g_half == 19) continue;
            if (channel_index_in_h2g_half == 0)  continue;

            int channel_index_no_CM_Calib = channel_index_in_h2g_half;
            if (channel_index_in_h2g_half > 19)      channel_index_no_CM_Calib -= 2;
            else if (channel_index_in_h2g_half > 0)  channel_index_no_CM_Calib -= 1;

            std::string chn_key = std::to_string(channel_index_no_CM_Calib);
            int col = -1, row = -1;
            if (sipm_board.contains(chn_key)) {
                col = sipm_board.at(chn_key).at("col").get<int>();
                row = sipm_board.at(chn_key).at("row").get<int>();
            } else {
                continue;
            }

            std::string board_key = std::to_string(half_id + asic_id * 2 + vldb_id * 4);
            int board_col = -1, board_row = -1;
            int board_rotated = 0;
            int board_flipped = 0;

            if (board_loc.contains(board_key)) {
                board_col = board_loc.at(board_key).at("col").get<int>();
                board_row = board_loc.at(board_key).at("row").get<int>();
            }
            if (board_rotation.contains(board_key)) {
                board_rotated = board_rotation.at(board_key).get<int>();
            }
            if (board_flip.contains(board_key)) {
                board_flipped = board_flip.at(board_key).get<int>();
            }

            int uni_col = -1, uni_row = -1;
            if (board_rotated == 0) {
                uni_col = board_col * board_cols + col;
                uni_row = board_row * board_rows + row;
            } else {
                // 旋转180°
                uni_col = board_col * board_cols + (board_cols - 1 - col);
                uni_row = board_row * board_rows + (board_rows - 1 - row);
            }

            if (uni_col < 0 || uni_row < 0) continue;
            if (uni_col >= NX || uni_row >= NY) continue;

            int pad_linear = uni_row * NX + uni_col;
            int chan_linear = vldb_id * channels_per_vldb + ch;

            lut[chan_linear] = pad_linear;
        }
    }
    return lut;
}