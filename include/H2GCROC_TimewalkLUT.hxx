// RobustTimewalkLUT.hxx
#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <utility>

#include "TH2D.h"
#include "TH1D.h"
#include "TGraph.h"

struct TimewalkLUT {
  // Bin centers on ADC axis and the corresponding robust central tendency of ToA residual
  std::vector<double> adc_centers;
  std::vector<double> dt_median;     // correction value (to subtract)
  std::vector<double> dt_sigma;      // ~1.4826*MAD as robustness estimate

  // Evaluate by piecewise-linear interpolation; clamp at edges.
  double eval(double adc) const {
    if (adc_centers.empty()) return 0.0;
    if (adc <= adc_centers.front()) return dt_median.front();
    if (adc >= adc_centers.back())  return dt_median.back();
    auto it = std::upper_bound(adc_centers.begin(), adc_centers.end(), adc);
    size_t i = std::max<size_t>(1, std::distance(adc_centers.begin(), it)) - 1;
    double x0 = adc_centers[i], x1 = adc_centers[i+1];
    double y0 = dt_median[i],   y1 = dt_median[i+1];
    double t = (adc - x0) / (x1 - x0 + 1e-12);
    return y0 + t * (y1 - y0);
  }

  // Export as TGraph for quick plotting
  TGraph* AsGraph() const {
    auto *g = new TGraph((int)adc_centers.size());
    for (int i = 0; i < (int)adc_centers.size(); ++i)
      g->SetPoint(i, adc_centers[i], dt_median[i]);
    return g;
  }
};

// ----- helpers -----

// Robust median from a vector (no copies if you pass a local buffer)
template <class T>
static inline double median_inplace(std::vector<T>& v) {
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  size_t n = v.size();
  size_t k = n / 2;
  std::nth_element(v.begin(), v.begin()+k, v.end());
  double m = v[k];
  if ((n & 1) == 0) { // even
    auto max_it = std::max_element(v.begin(), v.begin()+k);
    m = 0.5*(m + *max_it);
  }
  return m;
}

// Return ~1.4826*MAD (robust sigma) around given median
template <class T>
static inline double mad_sigma(const std::vector<T>& v, double med) {
  if (v.empty()) return 0.0;
  std::vector<double> dev; dev.reserve(v.size());
  for (auto &x : v) dev.push_back(std::abs((double)x - med));
  // median of deviations
  size_t n = dev.size(), k = n/2;
  std::nth_element(dev.begin(), dev.begin()+k, dev.end());
  double mdev = dev[k];
  if ((n & 1) == 0) {
    auto max_it = std::max_element(dev.begin(), dev.begin()+k);
    mdev = 0.5*(mdev + *max_it);
  }
  return 1.4826 * mdev;
}

// Pool Adjacent Violators (PAV) for isotonic regression (monotone decreasing).
// Minimizes sum w_i*(y_i - f_i)^2 subject to f_1 >= f_2 >= ... >= f_n
static inline void isotonic_decreasing(std::vector<double>& y, const std::vector<double>& w) {
  const int n = (int)y.size();
  if (n==0) return;
  // We implement decreasing by flipping sign: apply increasing on (-y)
  std::vector<double> yy(n), ww = w;
  for (int i=0;i<n;++i) yy[i] = -y[i];
  // PAV for increasing
  std::vector<double> v = yy, lambda = ww;
  int m = 0; // number of blocks-1
  std::vector<int> block_start(n), block_end(n);
  for (int i=0;i<n;++i) {
    block_start[m] = i; block_end[m] = i;
    while (m>0 && v[m-1] > v[m]) {
      // merge blocks m-1 and m with weights
      double num = lambda[m-1]*v[m-1] + lambda[m]*v[m];
      double den = lambda[m-1] + lambda[m];
      v[m-1] = num/den;
      lambda[m-1] = den;
      block_end[m-1] = block_end[m];
      --m;
    }
    ++m;
  }
  // Expand blocks back to full solution
  std::vector<double> sol(n);
  for (int b=0;b<m;++b)
    for (int i=block_start[b]; i<=block_end[b]; ++i)
      sol[i] = v[b];
  // Flip sign back to decreasing
  for (int i=0;i<n;++i) y[i] = -sol[i];
}

static inline void ReanchorLUTAtADC(TimewalkLUT& lut, double anchor_adc) {
  if (lut.adc_centers.empty()) return;
  const double offset = lut.eval(anchor_adc);   // f(A*)
  for (auto& y : lut.dt_median) y -= offset;    // y_i ← y_i - f(A*)
}

static void WeightedMovingMedian(std::vector<double>& y,
                                 const std::vector<double>& sigma,
                                 int win = 2) {
  const size_t n = y.size();
  std::vector<double> out(n);
  auto weight = [&](size_t i){ return 1.0 / (sigma[i] + 1e-9); };
  for (size_t i=0;i<n;++i){
    size_t L = (i>=(size_t)win? i-win : 0);
    size_t R = std::min(i+win, n-1);
    std::vector<std::pair<double,double>> buf; buf.reserve(R-L+1);
    for(size_t j=L;j<=R;++j) buf.emplace_back(y[j], weight(j));
    std::sort(buf.begin(), buf.end(),
              [](auto&a, auto&b){ return a.first < b.first; });
    double tot=0, acc=0; for(auto& p:buf) tot+=p.second;
    double med = buf.empty()? y[i] : buf.back().first;
    for(auto& p:buf){ acc+=p.second; if(acc>=0.5*tot){ med=p.first; break; } }
    out[i]=med;
  }
  y.swap(out);
}

static void TVDenoise1D(std::vector<double>& y,
                        const std::vector<double>& sigma,
                        double lambda = 0.002,   // 单位同 y（ns）；可调：0.001–0.01
                        int iters = 10) {
  const size_t n = y.size();
  if (n<3) return;
  std::vector<double> f = y, fnew(n);
  std::vector<double> w(n);
  for (size_t i=0;i<n;++i) w[i]=1.0/(sigma[i]+1e-9); // 稳健权

  for (int it=0; it<iters; ++it){
    fnew[0]=f[0]; fnew[n-1]=f[n-1];
    for (size_t i=1;i<n-1;++i){
      // 二次项拟合 + TV 正则近似：三点平滑，抑制折线
      double data = w[i]*y[i];
      double reg  = lambda*( (f[i-1]-2*f[i]+f[i+1]) );
      double denom= w[i] + 2*lambda + 1e-12;
      fnew[i] = (data + lambda*(f[i-1]+f[i+1]) - reg) / denom;
    }
    f.swap(fnew);
  }
  y.swap(f);
}

// 带容差的单调递减 PAV：仅当 v[m-1] > v[m] + eps 合并
static inline void isotonic_decreasing_tol(std::vector<double>& y,
                                           const std::vector<double>& w,
                                           double eps_const_ns = 0.0,
                                           const std::vector<double>* sigma = nullptr,
                                           double k_sigma = 0.0) {
  const int n = (int)y.size();
  if (n==0) return;

  // 转成递增问题：-y
  std::vector<double> v(n), lambda = w;
  for (int i=0;i<n;++i) v[i] = -y[i];

  auto local_eps = [&](int a, int b)->double {
    double e = eps_const_ns;
    if (sigma && k_sigma>0.0) {
      double sa = std::max(1e-9, (*sigma)[a]);
      double sb = std::max(1e-9, (*sigma)[b]);
      // 依据不确定度放宽容差（可改成别的形式）
      e = std::max(e, k_sigma * std::sqrt(sa*sa + sb*sb));
    }
    return e;
  };

  int m = 0;
  std::vector<int> bs(n), be(n);
  for (int i=0;i<n;++i) {
    bs[m] = be[m] = i;
    // 注意这里的“容差”体现在合并判据上：
    while (m>0 && v[m-1] > v[m] + local_eps(be[m-1], bs[m])) {
      double num = lambda[m-1]*v[m-1] + lambda[m]*v[m];
      double den = lambda[m-1] + lambda[m];
      v[m-1] = num/den;
      lambda[m-1] = den;
      be[m-1] = be[m];
      --m;
    }
    ++m;
  }
  // 展开块
  std::vector<double> sol(n);
  for (int b=0;b<m;++b)
    for (int i=bs[b]; i<=be[b]; ++i)
      sol[i] = v[b];

  // 还原为递减
  for (int i=0;i<n;++i) y[i] = -sol[i];
}

// ----- main builder -----
//
// Build a robust time-walk LUT from a TH2D (x=ADC, y=ToA residual).
// - For每个ADC列，做自适应并箱（左右扩展到 +/- expand bins），直到样本数≥Nmin，计算加权中位数与MAD。
// - 可选沿ADC轴做单调递减的等序回归（isotonic）。
// - 返回 piecewise-linear 可插值的 LUT。
//
inline TimewalkLUT BuildFitFreeTimewalkLUT(
    const TH2D* h2,
    int    Nmin              = 300,   // 每个ADC列期望的最少样本
    int    max_expand_bins   = 6,     // 最多左右扩展的列数
    bool   enforce_monotonic = true   // 沿ADC单调（通常 Δt 随ADC增大而减小）
) {
  if (!h2) throw std::runtime_error("BuildFitFreeTimewalkLUT: null TH2D");

  const int nx = h2->GetNbinsX();
  const int ny = h2->GetNbinsY();

  TimewalkLUT lut;
  lut.adc_centers.reserve(nx);
  lut.dt_median.reserve(nx);
  lut.dt_sigma.reserve(nx);

  // 预取轴信息
  std::vector<double> xcent(nx), ycent(ny);
  for (int ix=1; ix<=nx; ++ix) xcent[ix-1] = h2->GetXaxis()->GetBinCenter(ix);
  for (int iy=1; iy<=ny; ++iy) ycent[iy-1] = h2->GetYaxis()->GetBinCenter(iy);

  // 为了鲁棒中位数：把某个ADC邻域内的 (y 值) 展开成样本向量，用计数作重复或权重。
  // 这里采用“权重展开”为重复次数可能非常大；因此我们改用“加权中位数”实现：
  auto weighted_median = [&](int ix_center, int expand)->std::pair<double,double> {
    // 收集该 ADC 并箱窗口内所有 (ybin, weight)
    std::vector<double> ys; ys.reserve(ny);
    std::vector<double> ws; ws.reserve(ny);
    int ix_lo = std::max(1, ix_center - expand);
    int ix_hi = std::min(nx, ix_center + expand);
    double total = 0.0;

    for (int iy=1; iy<=ny; ++iy) {
      double wsum = 0.0;
      for (int ix=ix_lo; ix<=ix_hi; ++ix)
        wsum += h2->GetBinContent(ix, iy);
      if (wsum > 0) {
        ys.push_back(ycent[iy-1]);
        ws.push_back(wsum);
        total += wsum;
      }
    }
    if (total < 1.0) return std::make_pair(std::numeric_limits<double>::quiet_NaN(), 0.0);

    // 计算加权中位数（按 y 升序的累计权重过半位置）
    std::vector<int> ord(ys.size());
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return ys[a] < ys[b]; });

    double cum = 0.0;
    double med = ys[ord.back()];
    for (size_t k=0; k<ord.size(); ++k) {
      cum += ws[ord[k]];
      if (cum >= 0.5*total) { med = ys[ord[k]]; break; }
    }

    // 估算加权 MAD（简化：近似为未加权 MAD，用 ws 作为重复的近似）
    // 为效率，这里以“邻域 68% 分位宽度/1.349”近似 sigma，也很稳健：
    // 取 16% 与 84% 权重分位作为 IQR68。
    auto weighted_quantile = [&](double q)->double{
      double target = q*total, acc=0.0;
      for (size_t k=0; k<ord.size(); ++k) {
        acc += ws[ord[k]];
        if (acc >= target) return ys[ord[k]];
      }
      return ys[ord.back()];
    };
    double y16 = weighted_quantile(0.16);
    double y84 = weighted_quantile(0.84);
    double sigma = (y84 - y16) / 2.0; // ~ 1 sigma for near-symmetric

    return std::make_pair(med, sigma);
  };

  // 主循环：每个 ADC 列做自适应并箱，直到样本数达到 Nmin（或达到最大扩展）
  for (int ix=1; ix<=nx; ++ix) {
    int expand = 0;
    double sumN = 0.0;
    auto column_counts = [&](int e)->double{
      int lo = std::max(1, ix - e), hi = std::min(nx, ix + e);
      double s = 0.0;
      for (int j=lo; j<=hi; ++j) s += h2->Integral(j,j, 1, ny);
      return s;
    };
    sumN = column_counts(0);
    while (sumN < Nmin && expand < max_expand_bins) {
      ++expand;
      sumN = column_counts(expand);
    }

    auto [med, sigma] = weighted_median(ix, expand);
    if (std::isnan(med)) continue; // 没有统计量，跳过

    lut.adc_centers.push_back(xcent[ix-1]);
    lut.dt_median.push_back(med);
    lut.dt_sigma.push_back(sigma);
  }

  // // 可选：沿 ADC 轴做单调性约束（通常 Δt 随 ADC 增大而下降；如相反请先取反）
  // if (enforce_monotonic && lut.adc_centers.size() >= 3) {
  //   // 若你期望 "随 ADC 增大，校正 Δt 减小"，则这是“单调递减”约束
  //   std::vector<double> w(lut.dt_sigma.size(), 1.0);
  //   // 用 1/sigma^2 当权更合理，但为稳健防止小sigma过大权重，这里用 1/(sigma+eps)
  //   for (size_t i=0;i<w.size();++i) w[i] = 1.0 / (lut.dt_sigma[i] + 1e-6);
  //   isotonic_decreasing(lut.dt_median, w);
  // }
  // 可选：沿 ADC 轴做“带容差”的单调性约束（Δt 随 ADC 增大而下降）
// 仅当相邻点的上扬幅度显著超过容差时才合并，避免被压成一条直线
  if (enforce_monotonic && lut.adc_centers.size() >= 3) {
    std::vector<double> w(lut.dt_sigma.size());
    for (size_t i = 0; i < w.size(); ++i)
      w[i] = 1.0 / (lut.dt_sigma[i] + 1e-6);   // 稳健权重

    // 可调参数：
    const double eps_const_ns = 0.02;   // 常数容差（单位与你的Δt一致，示例：25ps）
    const double k_sigma      = 1.0;   // 按不确定度自适应放宽的系数（0~1 之间比较常用）

    // 说明：两种容差取较大者生效：max(eps_const_ns, k_sigma*sqrt(sigma_i^2+sigma_{i+1}^2))
    isotonic_decreasing_tol(lut.dt_median, w, eps_const_ns, &lut.dt_sigma, k_sigma);
  }

  WeightedMovingMedian(lut.dt_median, lut.dt_sigma, 2);
  TVDenoise1D(lut.dt_median, lut.dt_sigma, 0.003, 10);
  ReanchorLUTAtADC(lut, /*anchor_adc=*/900.0);

  return lut;
}