#include <bits/stdc++.h>
#include <Eigen/Dense>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

inline int from_total_channel_to_valid_channel(int _total_channel){
    // remove channels 76-151
    if (_total_channel >= 152){
        return _total_channel - 76;
    }
    if (_total_channel >= 76){
        return -1;
    } else {
        return _total_channel;
    }
}

inline int from_valid_channel_to_total_channel(int _valid_channel){
    // add back channels 76-151
    if (_valid_channel >= 0 && _valid_channel < 76){
        return _valid_channel;
    } else if (_valid_channel >= 76){
        return _valid_channel + 76;
    } else {
        return -1;
    }
}

struct AlignConfig {
    int    NCH            = 196;     // 通道数
    double DT_MIN         = -50.0;   // Δt 直方图范围（ns）
    double DT_MAX         = +50.0;
    int    DT_NBINS       = 800;     // 直方图 bin
    int    MIN_ENTRIES    = 50;      // 每个通道对参与解的最小样本数
    double Tbc_ns         = 0.0;     // 束交周期（设 0 关闭展开；例如 25ns）
    int    REF_CH         = 0;       // 参考通道（每个连通分量内部会自动选一个）
    double HUBER_C_SCALE  = 1.345;   // IRLS Huber 缩放
    bool   RUN_IRLS       = true;    // 是否做一轮 IRLS
    TFile* out_root = nullptr; // 输出 ROOT 文件指针（nullptr 则新建）
    // std::string out_root  = "toa_align.root";
    // std::string out_txt   = "toa_offsets.txt";
};

inline int pair_index(int i, int j, int N) {
    // i<j
    return i*(2*N - i - 1)/2 + (j - i - 1);
}

inline double unwrap_bc(double dt, double Tbc) {
    if (Tbc <= 0) return dt;
    dt = std::fmod(dt + 0.5*Tbc, Tbc);
    if (dt < 0) dt += Tbc;
    return dt - 0.5*Tbc;
}

struct PairStat {
    int i, j;
    TH1D* h = nullptr;
    double mu = 0;
    double sigma = 0;
    double w = 0;
};

class ToAAligner {
public:
    AlignConfig cfg;
    std::vector<PairStat> pairs; // size = N*(N-1)/2
    int N;
    // 诊断：校正前后 Δt 均值热图
    TH2D* hMeanBefore = nullptr;
    TH2D* hMeanAfter  = nullptr;

    ToAAligner(const AlignConfig& c): cfg(c), N(c.NCH) {
        makePairs();
        hMeanBefore = new TH2D("hMeanBefore","Mean dt before;ch i;ch j",
                               N, -0.5, N-0.5, N, -0.5, N-0.5);
        hMeanAfter  = new TH2D("hMeanAfter","Mean dt after;ch i;ch j",
                               N, -0.5, N-0.5, N, -0.5, N-0.5);
        hMeanBefore->SetDirectory(nullptr);
        hMeanAfter->SetDirectory(nullptr);
    }

    // 事件输入：传入本 event 的 (channel, toa_ns) 命中
    void addEvent(const std::vector<std::pair<int,double>>& hits) {
        if (hits.size() < 2) return;
        // 两两配对
        for (size_t a=0; a<hits.size(); ++a) {
            int i = hits[a].first; if (i<0 || i>=N) continue;
            double ti = hits[a].second;
            for (size_t b=a+1; b<hits.size(); ++b) {
                int j = hits[b].first; if (j<0 || j>=N) continue;
                double tj = hits[b].second;
                int ii = std::min(i,j), jj = std::max(i,j);
                double dt = (i<j) ? (ti - tj) : (tj - ti);
                dt = unwrap_bc(dt, cfg.Tbc_ns);
                pairs[pair_index(ii,jj,N)].h->Fill(dt);
            }
        }
    }

    // 统计量：为每个 pair 估计 mu/sigma/weight
    // void finalizePairStats() {
    //     for (auto& p : pairs) {
    //         const int n = (int)p.h->GetEntries();
    //         if (n < cfg.MIN_ENTRIES) { p.w = 0; continue; }
    //         // 分位数（稳健）
    //         double qs[3] = {0.16, 0.50, 0.84}; double qv[3];
    //         p.h->GetQuantiles(3, qv, qs);
    //         p.mu    = qv[1];
    //         p.sigma = 0.5*(qv[2]-qv[0]);
    //         if (p.sigma <= 1e-12) {
    //             // 退化保护
    //             double rms = p.h->GetRMS();
    //             p.sigma = (rms>1e-12) ? rms : 1.0;
    //         }
    //         p.w = n / (p.sigma*p.sigma);
    //     }
    // }
    void finalize_pairstats(std::vector<PairStat>& pairs, double min_entries=50) {
        for (auto& p : pairs) {
            const int n = (int)p.h->GetEntries();
            if (n < min_entries) { p.w = 0; continue; }

            // 1) 先找众数（最高 bin）
            int    imax = p.h->GetMaximumBin();
            double mode = p.h->GetBinCenter(imax);

            // 2) 以众数为中心取一个自适应窗口（先用直方图 RMS 估个尺度）
            double rms = p.h->GetRMS();
            if (rms < 1e-6) rms = (p.h->GetXaxis()->GetXmax() - p.h->GetXaxis()->GetXmin()) / 100.0;
            double win = std::clamp(2.0*rms, 1.0, 10.0); // 1~10 ns 之间
            double lo = mode - win, hi = mode + win;

            // 3) 统计窗口内的样本，做中位/分位估计（用 GetQuantiles 的话先临时复制子直方图，
            //    这里用更快的近似：在窗口内做加权均值/方差，再回退到分位）
            //    简化做法：再次用分位，但只对 [lo,hi] 范围。
            std::unique_ptr<TH1D> hwin((TH1D*)p.h->Clone("tmp"));
            hwin->SetDirectory(nullptr);
            for (int b=1; b<=hwin->GetNbinsX(); ++b) {
                double x = hwin->GetBinCenter(b);
                if (x < lo || x > hi) hwin->SetBinContent(b, 0.0);
            }

            int nw = (int)hwin->GetEntries();
            if (nw < min_entries/2) { // 窗太窄，退回全量分位
                double qs[3] = {0.16, 0.50, 0.84}, qv[3];
                p.h->GetQuantiles(3, qv, qs);
                p.mu    = qv[1];
                p.sigma = 0.5*(qv[2]-qv[0]);
            } else {
                double qs[3] = {0.16, 0.50, 0.84}, qv[3];
                hwin->GetQuantiles(3, qv, qs);
                p.mu    = qv[1];
                p.sigma = 0.5*(qv[2]-qv[0]);
            }

            if (p.sigma <= 1e-6) {
                double r = p.h->GetRMS();
                p.sigma = (r>1e-6) ? r : 1.0;
            }
            p.w = std::max(1.0, (double)n) / (p.sigma*p.sigma);
        }
    }

    // 根据有权边构建邻接表，找连通分量
    std::vector<std::vector<int>> connectedComponents() const {
        std::vector<std::vector<int>> adj(N);
        for (auto& p : pairs) {
            if (p.w > 0) {
                adj[p.i].push_back(p.j);
                adj[p.j].push_back(p.i);
            }
        }
        std::vector<int> vis(N,0);
        std::vector<std::vector<int>> comps;
        for (int s=0; s<N; ++s) if (!vis[s]) {
            std::vector<int> comp;
            std::vector<int> stk = {s}; vis[s]=1;
            while (!stk.empty()) {
                int u = stk.back(); stk.pop_back(); comp.push_back(u);
                for (int v : adj[u]) if (!vis[v]) { vis[v]=1; stk.push_back(v); }
            }
            comps.push_back(std::move(comp));
        }
        return comps;
    }

    // 在指定分量里解 offsets（ref 为分量内参考通道索引—not channel id）
    Eigen::VectorXd solveComponent(const std::vector<int>& comp, int ref_idx_in_comp = 0) {
        const int M = (int)comp.size();
        // map: ch -> local index
        std::unordered_map<int,int> loc;
        loc.reserve(M);
        for (int k=0;k<M;++k) loc[comp[k]] = k;

        Eigen::MatrixXd L = Eigen::MatrixXd::Zero(M,M);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(M);

        for (auto& p : pairs) {
            if (p.w <= 0) continue;
            auto it_i = loc.find(p.i);
            auto it_j = loc.find(p.j);
            if (it_i==loc.end() || it_j==loc.end()) continue;
            int i = it_i->second, j = it_j->second;
            double w = p.w, mu = p.mu;
            L(i,i) += w; L(j,j) += w;
            L(i,j) -= w; L(j,i) -= w;
            b(i)   += w*mu; b(j) -= w*mu;
        }

        // 删掉参考行列（o_ref = 0）
        const int K = M-1;
        Eigen::MatrixXd Lr(K,K);
        Eigen::VectorXd br(K);
        for (int r=0, rr=0; r<M; ++r) if (r!=ref_idx_in_comp) {
            for (int c=0, cc=0; c<M; ++c) if (c!=ref_idx_in_comp) {
                Lr(rr,cc) = L(r,c); ++cc;
            }
            br(rr) = b(r); ++rr;
        }
        Eigen::VectorXd orr = Lr.ldlt().solve(br);

        Eigen::VectorXd o = Eigen::VectorXd::Zero(M);
        for (int r=0, rr=0; r<M; ++r) if (r!=ref_idx_in_comp) o(r) = orr(rr++);
        return o;
    }

    // 计算残差并进行一轮 Huber IRLS 降权
    void irlsHuberUpdate(const std::vector<double>& o_global) {
        for (auto& p : pairs) {
            if (p.w <= 0 || p.sigma <= 0) continue;
            double r = (o_global[p.i] - o_global[p.j]) - p.mu;
            double z = r / (cfg.HUBER_C_SCALE * p.sigma);
            double psi_over_z = 1.0;
            if (std::abs(z) > 1.0) psi_over_z = 1.0 / std::abs(z);
            double neww = p.w * psi_over_z;
            p.w = std::max(neww, 1e-12);
        }
    }

    // 运行完整求解：finalize -> 分量解 -> 可选 IRLS -> 再解
    std::vector<double> solveAll(bool with_irls=true) {
        finalize_pairstats(pairs, cfg.MIN_ENTRIES);

        // 先粗略填充“校正前”均值热图（μ）
        for (auto& p : pairs) if (p.w>0) {
            hMeanBefore->SetBinContent(p.i+1, p.j+1, p.mu);
            hMeanBefore->SetBinContent(p.j+1, p.i+1, -p.mu);
        }

        std::vector<double> offsets(N, 0.0);
        auto comps = connectedComponents();

        // 第一次解：各分量独立求
        for (auto& comp : comps) {
            if (comp.size()==1) { offsets[comp[0]] = 0.0; continue; }
            // 参考取分量内 index 0
            Eigen::VectorXd o = solveComponent(comp, 0);
            for (int k=0;k<(int)comp.size();++k) offsets[comp[k]] = o(k);
        }

        if (with_irls && cfg.RUN_IRLS) {
            irlsHuberUpdate(offsets);
            // 第二次解（再跑一次）
            for (auto& comp : comps) {
                if (comp.size()==1) { offsets[comp[0]] = 0.0; continue; }
                Eigen::VectorXd o = solveComponent(comp, 0);
                for (int k=0;k<(int)comp.size();++k) offsets[comp[k]] = o(k);
            }
        }

        // “校正后”均值热图（根据 offsets 计算残差的均值：近似 -μ + (o_i - o_j)）
        for (auto& p : pairs) if (p.w>0) {
            double rmean = (offsets[p.i] - offsets[p.j]) - p.mu;
            hMeanAfter->SetBinContent(p.i+1, p.j+1, rmean);
            hMeanAfter->SetBinContent(p.j+1, p.i+1, -rmean);
        }

        return offsets;
    }

    // 保存结果
    void saveResults(const std::vector<double>& offsets) {
        // {
        //     std::ofstream ofs(cfg.out_txt);
        //     ofs << "# channel  offset_ns\n";
        //     for (int i=0;i<N;++i) ofs << i << " " << std::setprecision(9) << offsets[i] << "\n";
        // }
        // ROOT
        if (cfg.out_root == nullptr) {
            cfg.out_root = TFile::Open("toa_align.root", "RECREATE");
        }
            // TFile* fout = TFile::Open(cfg.out_root.c_str(), "RECREATE");
        // save th1ds in a directory
        cfg.out_root->mkdir("pair_histograms");
        cfg.out_root->cd("pair_histograms");
        for (auto& p : pairs) if (p.h) p.h->Write();
        cfg.out_root->cd();
        hMeanBefore->Write();
        hMeanAfter->Write();

        gStyle->SetOptStat(0);
        TCanvas c1("c_before","before",1000,800);
        hMeanBefore->Draw("COLZ");
        c1.Write("c_before");
        TCanvas c2("c_after","after",1000,800);
        hMeanAfter->Draw("COLZ");
        c2.Write("c_after");

        cfg.out_root->mkdir("ToA_Offsets");
        cfg.out_root->cd("ToA_Offsets");
        for (size_t chn = 0; chn < offsets.size(); chn++) {
            TParameter<double> *p_offset = new TParameter<double>(("toa_offset_ch" + std::to_string(chn)).c_str(), offsets[chn]);
            p_offset->Write();
        }
        cfg.out_root->cd();

        // TDirectory* toa_offset_dir = output_root->mkdir("ToA_Offsets");
        // toa_offset_dir->cd();
        // spdlog::info("Calculated ToA offsets for each channel:");
        // for (size_t chn = 0; chn < toa_offsets.size(); chn++) {
        //     spdlog::info("Channel {}: ToA Offset = {:.3f} ns", chn, toa_offsets[chn]);
        //     TParameter<double> *p_offset = new TParameter<double>(("toa_offset_ch" + std::to_string(chn)).c_str(), toa_offsets[chn]);
        //     p_offset->Write();
        // }
    }

private:
    void makePairs() {
        pairs.reserve(N*(N-1)/2);
        for (int i=0;i<N;i++) for (int j=i+1;j<N;j++) {
            auto name = Form("h_dt_%03d_%03d", i, j);
            auto title= Form("dt = t[%d]-t[%d]", i, j);
            TH1D* h = new TH1D(name, title, cfg.DT_NBINS, cfg.DT_MIN, cfg.DT_MAX);
            h->SetDirectory(nullptr);
            pairs.push_back(PairStat{i,j,h});
        }
    }
};
