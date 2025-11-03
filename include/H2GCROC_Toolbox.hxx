#ifndef H2GCROC_TOOLBOX_HXX
#define H2GCROC_TOOLBOX_HXX

#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"

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

#endif