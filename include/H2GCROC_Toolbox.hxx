#ifndef H2GCROC_TOOLBOX_HXX
#define H2GCROC_TOOLBOX_HXX

#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include <Eigen/Sparse>

inline double crystallball_left(double *x, double *par) {
    // par[0] = A
    // par[1] = a
    // par[2] = n
    // par[3] = mean
    // par[4] = sigma
    const double N      = par[0];
    const double alpha  = par[1];
    const double n      = par[2];
    const double mean   = par[3];
    const double sigma  = par[4];

    double z = (x[0] - mean) / sigma;
    if (z > -alpha) {
        return N * exp(-0.5 * z * z);
    } else {
        double A = pow(n / fabs(alpha), n) * exp(-0.5 * alpha * alpha);
        double B = n / fabs(alpha) - fabs(alpha);
        return N * A * pow(B - z, -n);
    }
}

struct Est { double avg, err_stat, err_sys; };

static Est combine_one(const std::vector<double>& x, const std::vector<double>& dx) {
    const size_t N = x.size();
    if (N == 0) return {0, 0, 0};
    std::vector<double> w(N);
    for (size_t i=0;i<N;++i) {
        double d = dx[i];
        if (!(d>0)) d = 1e30;
        w[i] = 1.0/(d*d);
    }
    const double Sw  = std::accumulate(w.begin(), w.end(), 0.0);
    const double Sw2 = std::accumulate(w.begin(), w.end(), 0.0,
                         [&](double s,double wi){ return s + wi*wi; });

    // 加权均值（fixed-effect）
    double xbar = 0.0;
    if (Sw > 0) {
        for (size_t i=0;i<N;++i) xbar += w[i]*x[i];
        xbar /= Sw;
    } else {
        xbar = std::accumulate(x.begin(), x.end(), 0.0)/double(N);
    }

    double Q = 0.0;
    for (size_t i=0;i<N;++i) Q += w[i]*(x[i]-xbar)*(x[i]-xbar);
    const double dof = (N>1)? (N-1) : 0.0;

    double tau2 = 0.0;
    double denom = Sw - (Sw2 / (Sw>0 ? Sw : 1.0));
    if (denom > 0.0) tau2 = std::max(0.0, (Q - dof) / denom);

    double err_stat = (Sw>0)? std::sqrt(1.0/Sw) : 0.0;
    double err_sys  = std::sqrt(tau2);

    return {xbar, err_stat, err_sys};
}

void mean_sigma_list_calculator(
    const std::vector<double>& means,
    const std::vector<double>& mean_errs,
    const std::vector<double>& sigmas,
    const std::vector<double>& sigma_errs,
    double& mean_avg,
    double& mean_err_sys,
    double& mean_err_stat,
    double& sigma_avg,
    double& sigma_err_sys,
    double& sigma_err_stat
) {
    Est m = combine_one(means,  mean_errs);
    Est s = combine_one(sigmas, sigma_errs);
    mean_avg      = m.avg;
    mean_err_stat = m.err_stat;
    mean_err_sys  = m.err_sys;
    sigma_avg      = s.avg;
    sigma_err_stat = s.err_stat;
    sigma_err_sys  = s.err_sys;
}

struct FastSmooth {
  double tv;          // total variation
  double tv_norm;     // tv / sum
  double chi2_adj;    // sum (dy^2 / (y_i+y_{i+1}+eps))
  double chi2_adj_ndf;// chi2_adj / (N-1)
  int    n;           // bins used
  double sum;         // sum of counts
};

inline FastSmooth FastSmoothMetrics(const TH1D* h) {
  const int n = h->GetNbinsX();
  FastSmooth m{}; if (n < 2) return m;

  const double* arr = h->GetArray();
  double sum = 0.0, tv = 0.0, chi2 = 0.0;
  const double eps = 1e-12;

  double prev = arr[1];
  sum += prev;

  for (int i = 2; i <= n; ++i) {
    double yi = arr[i];
    sum += yi;
    double dy = yi - prev;
    tv += std::abs(dy);
    chi2 += (dy*dy) / (yi + prev + eps);
    prev = yi;
  }

  m.n = n;
  m.sum = sum;
  m.tv = tv;
  m.tv_norm = (sum > 0.0) ? tv / sum : 0.0;
  m.chi2_adj = chi2;
  m.chi2_adj_ndf = chi2 / (n - 1);
  return m;
}

inline void TH2D_Slice_X_Gaussian_Fit(TH2D& h2d,
                                      TGraphErrors& gr,
                                      std::vector<TH1D*>* out_slices,   // may be nullptr
                                      std::vector<TF1*>*  out_fits,     // may be nullptr
                                      int xbin_step,
                                      int xbin_min = 0,
                                      int xbin_max = 0,
                                      int min_entries = 100)
{
    if (xbin_step <= 0) xbin_step = 1;

    const int nbinsx = h2d.GetNbinsX();
    if (xbin_min <= 0) xbin_min = 1;
    if (xbin_max <= 0 || xbin_max > nbinsx) xbin_max = nbinsx;

    // Avoid auto-registering temporary histos in gDirectory
    const bool old_adddir = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);

    for (int ix = xbin_min; ix <= xbin_max; ix += xbin_step) {
        // If you want to aggregate multiple X bins per step, use:
        // const int jx = std::min(ix + xbin_step - 1, xbin_max);
        // auto* raw = h2d.ProjectionY(Form("%s_py_%d_%d", h2d.GetName(), ix, jx), ix, jx);
        auto* raw = h2d.ProjectionY(Form("%s_py_ix%d", h2d.GetName(), ix), ix, ix);

        // Manage lifetime: start wrapped, detach from any directory
        std::unique_ptr<TH1D> hslice(raw);
        hslice->SetDirectory(nullptr);

        if (!hslice || hslice->GetEntries() < min_entries) continue;

        // Robust seeds via quantiles
        double q[3] = {0.16, 0.50, 0.84};
        double v[3] = {0., 0., 0.};
        hslice->GetQuantiles(3, v, q);
        double seed_mean  = v[1];
        double seed_sigma = 0.5 * (v[2] - v[0]);
        if (!(seed_sigma > 0)) { seed_mean = hslice->GetMean(); seed_sigma = hslice->GetRMS(); }
        if (!(seed_sigma > 0)) continue;

        // Fit window and function
        const double ylo = seed_mean - 2.5 * seed_sigma;
        const double yhi = seed_mean + 2.5 * seed_sigma;
        TF1 fgaus("fgaus", "gaus", ylo, yhi);
        fgaus.SetParameters(hslice->GetMaximum(), seed_mean, seed_sigma);

        const int fitStatus = hslice->Fit(&fgaus, "QNR"); // quiet, no draw, respect range

        double y    = seed_mean;
        double yerr = seed_sigma / std::sqrt(std::max(1.0, hslice->GetEntries()*1.0));
        if (fitStatus == 0) {
            y    = fgaus.GetParameter(1);
            yerr = fgaus.GetParError(1);
        }

        // X value and error
        const double x    = h2d.GetXaxis()->GetBinCenter(ix);
        const double xerr = 0.5 * h2d.GetXaxis()->GetBinWidth(ix);
        const int    ip   = gr.GetN();
        gr.SetPoint(ip, x, y);
        gr.SetPointError(ip, xerr, yerr);

        // Hand off ownership if requested
        if (out_slices) {
            // Give caller a heap-owned copy with a unique name
            auto* kept = static_cast<TH1D*>(hslice->Clone(Form("%s_kept", hslice->GetName())));
            kept->SetDirectory(nullptr);
            out_slices->push_back(kept);
        }

        if (out_fits) {
            // Clone the fitted function with a unique name so it outlives the stack
            auto* fcopy = static_cast<TF1*>(fgaus.Clone(Form("fgaus_%s_ix%d", h2d.GetName(), ix)));
            out_fits->push_back(fcopy);
        }

        // If not exporting, unique_ptr cleans the slice here.
    }

    TH1::AddDirectory(old_adddir);
}

struct EdgeAccum { double sum_wd = 0.0; double sum_w = 0.0; int samples = 0; };

inline uint64_t undirected_key(int a, int b){
    int p = (a < b) ? a : b;
    int q = (a < b) ? b : a;
    return ( (uint64_t)(uint32_t)p << 32 ) | (uint64_t)(uint32_t)q;
}

// Solve per-pixel offsets using pairwise diffs, but read ToA/ADC from per-channel matrices via pid->channel mapping.
inline void diff_toa_offset_calculator_mapped(
    const std::vector<std::vector<double>>& event_chn_toa_matrix, // [events][channels]
    const std::vector<std::vector<double>>& event_chn_adc_matrix, // [events][channels]
    std::vector<double>& pixel_offset_list,                       // [pixels]; will resize to num_pixels
    const std::vector<std::array<int,4>>& neigh,                  // [pixels][4], -1 if none
    const std::vector<int>& deg,                                  // [pixels]
    const std::vector<int>& pid2ch,                               // [pixels] -> channel col (or -1)
    double min_edge_weight = 0.0,                                 // drop edges with total weight < this
    int    min_edge_samples = 1,                                  // drop edges seen in < this many events
    bool   verbose = true
){
    const int num_events = (int)event_chn_toa_matrix.size();
    if (num_events == 0) return;

    const int num_pixels = (int)neigh.size();
    if (num_pixels == 0) return;

    if ((int)deg.size() != num_pixels || (int)pid2ch.size() != num_pixels)
        throw std::runtime_error("mapped solver: size mismatch (deg/neigh/pid2ch).");

    const int ncols0 = (int)event_chn_toa_matrix[0].size();
    for (int e=0; e<num_events; ++e){
        if ((int)event_chn_toa_matrix[e].size() != ncols0 ||
            (int)event_chn_adc_matrix[e].size() != ncols0)
            throw std::runtime_error("mapped solver: event matrices have inconsistent column counts.");
    }

    if (pixel_offset_list.size() != (size_t)num_pixels)
        pixel_offset_list.assign(num_pixels, 0.0);

    std::unordered_map<uint64_t, EdgeAccum> accum;
    accum.reserve((size_t)num_pixels*2);

    long long total_pairs = 0, used_pairs = 0;

    // Accumulate differences over events using pid->channel to read matrix columns.
    for (int e=0; e<num_events; ++e){
        const auto& toa = event_chn_toa_matrix[e];
        const auto& adc = event_chn_adc_matrix[e];

        for (int p=0; p<num_pixels; ++p){
            int ch_p = pid2ch[p];
            if (ch_p < 0 || ch_p >= ncols0) continue;
            double tp = toa[ch_p];
            if (!std::isfinite(tp)) continue;

            double Ap = adc[ch_p];
            for (int k=0; k<deg[p]; ++k){
                int q = neigh[p][k];
                if (q < 0 || q >= num_pixels) continue;
                if (p >= q) continue; // undirected once

                ++total_pairs;

                int ch_q = pid2ch[q];
                if (ch_q < 0 || ch_q >= ncols0) continue;
                double tq = toa[ch_q];
                if (!std::isfinite(tq)) continue;

                // difference (wrap here if needed)
                double d = tp - tq;

                // weight ~ sqrt(min(Ap,Aq)) (fallback 1)
                double Aq = adc[ch_q];
                double w = 1.0;
                if (std::isfinite(Ap) && std::isfinite(Aq)) {
                    double mA = std::max(0.0, std::min(Ap, Aq));
                    if (mA > 0.0) w = std::sqrt(mA);
                }

                auto& ac = accum[undirected_key(p,q)];
                ac.sum_wd += w * d;
                ac.sum_w  += w;
                ac.samples += 1;
                ++used_pairs;
            }
        }
    }

    // Build adjacency from kept edges
    std::vector<std::vector<int>> adj(num_pixels);
    int kept_edges = 0;
    for (const auto& kv : accum){
        const auto& ac = kv.second;
        if (ac.sum_w < min_edge_weight) continue;
        if (ac.samples < min_edge_samples) continue;
        int p = (int)(kv.first >> 32);
        int q = (int)(kv.first & 0xffffffffu);
        adj[p].push_back(q);
        adj[q].push_back(p);
        ++kept_edges;
    }

    if (verbose){
        std::fprintf(stderr, "[diff_toa_offset_mapped] events=%d, pixels=%d, cols=%d, total_pairs=%lld, used_pairs=%lld, kept_edges=%d\n",
                     num_events, num_pixels, ncols0, total_pairs, used_pairs, kept_edges);
    }

    if (kept_edges == 0){
        std::fill(pixel_offset_list.begin(), pixel_offset_list.end(), 0.0);
        return;
    }

    using SpMat   = Eigen::SparseMatrix<double>;
    using Triplet = Eigen::Triplet<double>;

    // Connected components solve
    std::vector<int> comp(num_pixels, -1);
    int cid = 0;

    auto solve_component = [&](const std::vector<int>& nodes){
        // local indexing
        const int Mfull = (int)nodes.size();
        std::vector<int> lidx(num_pixels, -1);
        for (int i=0;i<Mfull;++i) lidx[nodes[i]] = i;

        std::vector<Triplet> trips;
        trips.reserve(Mfull*5);
        Eigen::VectorXd b_full = Eigen::VectorXd::Zero(Mfull);

        for (const auto& kv : accum){
            const auto& ac = kv.second;
            if (ac.sum_w < min_edge_weight || ac.samples < min_edge_samples) continue;
            int gp = (int)(kv.first >> 32);
            int gq = (int)(kv.first & 0xffffffffu);
            int p = lidx[gp], q = lidx[gq];
            if (p < 0 || q < 0) continue;
            double W = ac.sum_w;
            double dbar = ac.sum_wd / std::max(1e-30, ac.sum_w);

            trips.emplace_back(p,p,  W);
            trips.emplace_back(q,q,  W);
            trips.emplace_back(p,q, -W);
            trips.emplace_back(q,p, -W);

            b_full[p] += W * dbar;
            b_full[q] -= W * dbar;
        }

        if (Mfull == 1){ pixel_offset_list[nodes[0]] = 0.0; return; }

        // choose a reference with degree>0 if possible
        int ref_local = 0;
        for (int i=0;i<Mfull;++i){ if (!adj[nodes[i]].empty()){ ref_local = i; break; } }

        SpMat A_full(Mfull, Mfull);
        A_full.setFromTriplets(trips.begin(), trips.end());

        // build reduced system
        const int M = Mfull - 1;
        SpMat A_red(M, M);
        std::vector<Triplet> trips_red; trips_red.reserve(trips.size());
        std::vector<int> mapIdx(Mfull, -1);
        for (int i=0,r=0;i<Mfull;++i) if (i != ref_local) mapIdx[i] = r++;

        for (int k=0;k<A_full.outerSize();++k){
            for (SpMat::InnerIterator it(A_full,k); it; ++it){
                int i = it.row(), j = it.col();
                if (i == ref_local || j == ref_local) continue;
                trips_red.emplace_back(mapIdx[i], mapIdx[j], it.value());
            }
        }
        A_red.setFromTriplets(trips_red.begin(), trips_red.end());

        Eigen::VectorXd b_red(M);
        for (int i=0,r=0;i<Mfull;++i) if (i != ref_local) b_red[r++] = b_full[i];

        // Solve: LDLT; fallback CG; fallback CG with tiny Tikhonov
        Eigen::VectorXd x_red;
        bool ok = false;
        {
            Eigen::SimplicialLDLT<SpMat> ldlt;
            ldlt.compute(A_red);
            if (ldlt.info() == Eigen::Success) {
                x_red = ldlt.solve(b_red);
                if (ldlt.info() == Eigen::Success) ok = true;
            }
        }
        if (!ok){
            Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> cg;
            cg.setMaxIterations(std::max(1000, 10*M));
            cg.setTolerance(1e-10);
            cg.compute(A_red);
            x_red = cg.solve(b_red);
            if (cg.info() == Eigen::Success) ok = true;
        }
        if (!ok){
            const double lambda = 1e-6;
            SpMat A_reg = A_red;
            for (int i=0;i<M;++i) A_reg.coeffRef(i,i) += lambda;
            Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> cg;
            cg.setMaxIterations(std::max(1500, 10*M));
            cg.setTolerance(1e-10);
            cg.compute(A_reg);
            x_red = cg.solve(b_red);
            if (cg.info() != Eigen::Success)
                throw std::runtime_error("mapped solver: regularized CG failed.");
        }

        for (int i=0,r=0;i<Mfull;++i){
            int g = nodes[i];
            if (i == ref_local) pixel_offset_list[g] = 0.0;
            else                pixel_offset_list[g] = x_red[r++];
        }
    };

    for (int i=0;i<num_pixels;++i){
        if (comp[i] != -1) continue;
        if (adj[i].empty()) { comp[i] = -2; continue; }
        std::vector<int> nodes; nodes.reserve(64);
        std::queue<int> q; q.push(i); comp[i]=cid; nodes.push_back(i);
        while(!q.empty()){
            int u=q.front(); q.pop();
            for (int v: adj[u]) if (comp[v]==-1){ comp[v]=cid; q.push(v); nodes.push_back(v); }
        }
        solve_component(nodes);
        ++cid;
    }

    // median-center over used pixels
    std::vector<char> used(num_pixels, 0);
    for (const auto& kv : accum){
        const auto& ac = kv.second;
        if (ac.sum_w < min_edge_weight || ac.samples < min_edge_samples) continue;
        int p = (int)(kv.first >> 32);
        int q = (int)(kv.first & 0xffffffffu);
        used[p]=1; used[q]=1;
    }
    std::vector<double> vals; vals.reserve(num_pixels);
    for (int i=0;i<num_pixels;++i) if (used[i]) vals.push_back(pixel_offset_list[i]);
    if (!vals.empty()){
        std::nth_element(vals.begin(), vals.begin()+vals.size()/2, vals.end());
        double med = vals[vals.size()/2];
        for (int i=0;i<num_pixels;++i) if (used[i]) pixel_offset_list[i] -= med;
    }

    if (verbose){
        int used_nodes = 0, iso_nodes = 0;
        for (int i=0;i<num_pixels;++i){
            if (!adj[i].empty()) ++used_nodes;
            else ++iso_nodes;
        }
        std::fprintf(stderr,"[diff_toa_offset_mapped] components=%d used_pixels=%d iso_pixels=%d\n",
                     cid, used_nodes, iso_nodes);
    }
}

inline void diff_toa_offset_calculator_2(
    std::vector<std::vector<double>>& event_chn_toa_matrix,
    std::vector<std::vector<double>>& event_chn_adc_matrix,
    std::vector<double>& toa_offset_list,
    std::vector<std::array<int,4>>& neigh,
    std::vector<int>& deg
) {
    (void)neigh; // not used in this method
    (void)deg;   // not used in this method

    const int num_events = static_cast<int>(event_chn_toa_matrix.size());
    if (num_events == 0) return;

    // Basic shape checks
    const int num_channels = static_cast<int>(event_chn_toa_matrix[0].size());
    if (num_channels == 0) return;
    if (static_cast<int>(event_chn_adc_matrix.size()) != num_events)
        throw std::runtime_error("diff_toa_offset_calculator: adc matrix row count mismatch");
    for (int e = 0; e < num_events; ++e) {
        if ((int)event_chn_toa_matrix[e].size() != num_channels ||
            (int)event_chn_adc_matrix[e].size() != num_channels) {
            throw std::runtime_error("diff_toa_offset_calculator: matrices must have same number of columns per event");
        }
    }

    if (toa_offset_list.empty()) toa_offset_list.assign(num_channels, 0.0);
    else if ((int)toa_offset_list.size() != num_channels) toa_offset_list.resize(num_channels, 0.0);

    // Tunables for the reference-time estimator
    const int   K_earliest = 4;        // number of earliest hits to form t_ref (3–5 is typical)
    const double alpha     = 0.5;      // weight exponent: w ~ ADC^alpha
    const double eps_adc   = 1e-9;     // avoid zero/negative ADC in pow
    const bool  use_trim   = true;     // optional: trim outliers when computing residuals
    const double trim_ns   = 12.5;     // discard residuals with |r| > trim_ns (set <=0 to disable)

    // Per-channel residual pools
    std::vector<std::vector<double>> residuals(num_channels);
    for (auto& v : residuals) v.reserve(std::min(128, std::max(8, num_events / 8)));

    // Temporary container used per event
    struct Hit { double t; double w; int ch; };
    std::vector<Hit> hits; hits.reserve(256);

    for (int e = 0; e < num_events; ++e) {
        const auto& toa_row = event_chn_toa_matrix[e];
        const auto& adc_row = event_chn_adc_matrix[e];

        // 1) Gather valid hits (finite ToA)
        hits.clear();
        hits.shrink_to_fit(); // keep capacity but allow growth
        for (int ch = 0; ch < num_channels; ++ch) {
            const double t  = toa_row[ch];
            if (!std::isfinite(t)) continue;
            const double A  = adc_row[ch];
            const double Aw = std::isfinite(A) ? std::max(A, 0.0) : 0.0;
            const double w  = std::pow(std::max(Aw, eps_adc), alpha);
            hits.push_back({t, w, ch});
        }
        if (hits.size() < 2) continue; // need at least 2 hits to form a sensible t_ref

        // 2) Select the earliest K hits (partial selection)
        const int K = std::min(K_earliest, static_cast<int>(hits.size()));
        std::nth_element(hits.begin(), hits.begin() + (K - 1), hits.end(),
                         [](const Hit& a, const Hit& b){ return a.t < b.t; });

        // Sort those K by time for weighted-median computation
        std::sort(hits.begin(), hits.begin() + K, [](const Hit& a, const Hit& b){ return a.t < b.t; });

        // Weighted median over the K earliest
        double wsum = 0.0; for (int i = 0; i < K; ++i) wsum += hits[i].w;
        double half = 0.5 * (wsum > 0.0 ? wsum : static_cast<double>(K));
        double acc  = 0.0;
        double t_ref = hits[0].t; // fallback
        if (wsum > 0.0) {
            for (int i = 0; i < K; ++i) {
                acc += hits[i].w;
                if (acc >= half) { t_ref = hits[i].t; break; }
            }
        } else {
            // All weights zero → plain median of the K earliest
            t_ref = (K % 2) ? hits[K/2].t : 0.5 * (hits[K/2 - 1].t + hits[K/2].t);
        }

        // 3) Accumulate residuals for all valid channels in this event
        for (int ch = 0; ch < num_channels; ++ch) {
            const double t = toa_row[ch];
            if (!std::isfinite(t)) continue;
            double r = t - t_ref;

            if (use_trim && trim_ns > 0.0) {
                if (!(std::abs(r) <= trim_ns)) continue;
            }

            residuals[ch].push_back(r);
        }
    } // end loop over events

    // Per-channel median of residuals
    auto median_of = [](std::vector<double>& v) -> double {
        if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
        const size_t n = v.size();
        const size_t mid = n / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        double med = v[mid];
        if ((n % 2) == 0) {
            // need the lower mid
            auto it = std::max_element(v.begin(), v.begin() + mid);
            const double lower = *it;
            med = 0.5 * (lower + v[mid]);
        }
        return med;
    };

    std::vector<double> used_offsets; used_offsets.reserve(num_channels);
    for (int ch = 0; ch < num_channels; ++ch) {
        double med = median_of(residuals[ch]);
        toa_offset_list[ch] = med; // may be NaN if insufficient stats
        if (std::isfinite(med)) used_offsets.push_back(med);
    }

    // Global centering: subtract overall median across channels with valid offsets
    if (!used_offsets.empty()) {
        const size_t mid = used_offsets.size() / 2;
        std::nth_element(used_offsets.begin(), used_offsets.begin() + mid, used_offsets.end());
        double global_med = used_offsets[mid];
        if ((used_offsets.size() % 2) == 0) {
            auto it = std::max_element(used_offsets.begin(), used_offsets.begin() + mid);
            global_med = 0.5 * (*it + used_offsets[mid]);
        }
        for (double& x : toa_offset_list) {
            if (std::isfinite(x)) x -= global_med;
            else x = 0.0; // optional: replace NaN by 0 so downstream math is safe
        }
    } else {
        // No valid offsets at all → zeros
        std::fill(toa_offset_list.begin(), toa_offset_list.end(), 0.0);
    }
}

inline double decode_toa_value_ns_with_dnl(UInt_t val2, std::vector<double>& coarse_dnl_lookup, std::vector<double>& fine_dnl_lookup) {
    constexpr double scale0 = 0.025;
    constexpr double scale1 = 0.2;
    constexpr double scale2 = 6.25;

    unsigned int part0 = val2 & 0x07;
    unsigned int part1 = (val2 >> 3) & 0x1F;
    unsigned int part2 = (val2 >> 8) & 0x03;

    return scale0 * fine_dnl_lookup[part0] * 8.0 +
           scale1 * coarse_dnl_lookup[part1] * 32.0 +
           scale2 * static_cast<double>(part2);
}

#endif