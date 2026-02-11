#ifndef COMMON_HPP
#define COMMON_HPP

#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include <string>
#include <regex>
#include <filesystem>
#include <sstream>
#include <cstdlib>
#include <system_error>
#include <numeric>
#include <span>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>
#include <cstdint>

#include "TCanvas.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TH2.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TVirtualFitter.h"
#include "TParameter.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadExecutor.hxx>

#include <nlohmann/json.hpp>

#include "cxxopts.hpp"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#ifndef FPGA_CHANNEL_NUMBER
#define FPGA_CHANNEL_NUMBER 152
#endif

#ifndef FPGA_CHANNEL_NUMBER_VALID
#define FPGA_CHANNEL_NUMBER_VALID 144
#endif

using json = nlohmann::json;

namespace tb::logging {

class Stream {
public:
    explicit Stream(spdlog::level::level_enum level) : level_(level) {}
    ~Stream() {
        spdlog::log(level_, buffer_.str());
    }

    template <typename T>
    Stream& operator<<(const T& value) {
        buffer_ << value;
        return *this;
    }

    Stream& operator<<(std::ostream& (*manip)(std::ostream&)) {
        buffer_ << manip;
        return *this;
    }

private:
    spdlog::level::level_enum level_;
    std::ostringstream buffer_;
};

} // namespace tb::logging

#define TB_LOG_LEVEL_INFO spdlog::level::info
#define TB_LOG_LEVEL_WARNING spdlog::level::warn
#define TB_LOG_LEVEL_ERROR spdlog::level::err
#define TB_LOG_LEVEL_DEBUG spdlog::level::debug

#define CANVAS_TITLE "FoCal-H Prototype 3"
#define TESTBEAM_TITLE "SPS H2 2025 October Testbeam"

#ifndef LOG
#define LOG(level) ::tb::logging::Stream(TB_LOG_LEVEL_##level)
#endif

inline void configure_logger(bool verbose) {
    static std::shared_ptr<spdlog::logger> logger = []() {
        auto sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        sink->set_pattern("%H:%M:%S [%^%l%$] %v");
        auto instance = std::make_shared<spdlog::logger>("tb_logger", sink);
        spdlog::register_logger(instance);
        spdlog::set_default_logger(instance);
        return instance;
    }();

    logger->set_level(verbose ? spdlog::level::debug : spdlog::level::info);
    spdlog::set_level(logger->level());
    spdlog::flush_on(spdlog::level::info);
}

struct ScriptOptions {
    std::string input_file;
    std::string output_file;
    std::string output_folder;
    int n_events;
    bool verbose;
    bool focal;
    bool timewalk;
    std::string script_name;
    std::string script_version;
    std::string pedestal_file;
    std::string csv_file;
    std::string timewalk_file;
    std::string fitting_file;
};

namespace tb::detail {

inline std::string script_basename(const char* arg0) {
    std::string name = arg0 ? std::string(arg0) : std::string("program");
    auto pos = name.find_last_of("/\\");
    if (pos != std::string::npos) {
        name = name.substr(pos + 1);
    }
    auto dot = name.find_last_of('.');
    if (dot != std::string::npos) {
        name = name.substr(0, dot);
    }
    return name;
}

inline void ensure_output_folder(ScriptOptions& opts, const std::string& override_folder) {
    if (!override_folder.empty()) {
        opts.output_folder = override_folder;
    } else if (opts.output_folder.empty()) {
        std::filesystem::path output_path(opts.output_file);
        if (output_path.has_parent_path()) {
            opts.output_folder = output_path.parent_path().string();
        }
        if (opts.output_folder.empty()) {
            opts.output_folder = "./dump/" + opts.script_name;
        }
    }

    if (!opts.output_folder.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(opts.output_folder, ec);
        if (ec) {
            spdlog::error("Failed to create output folder {}: {}", opts.output_folder, ec.message());
            std::exit(1);
        }
    }
}

inline void validate_existing_file(const std::string& path, const std::string& label) {
    if (path.empty()) {
        spdlog::error("{} must be provided.", label);
        std::exit(1);
    }
    if (::access(path.c_str(), F_OK) != 0) {
        spdlog::error("{} {} does not exist.", label, path);
        std::exit(1);
    }
}

inline void warn_existing_output(const std::string& path) {
    if (!path.empty() && ::access(path.c_str(), F_OK) == 0) {
        spdlog::warn("Output file {} already exists and will be overwritten.", path);
    }
}

} // namespace tb::detail

inline ScriptOptions parse_arguments_single_root(int argc, char **argv, const std::string& version = "0.1") {
    ScriptOptions opts{};
    opts.script_version = version;
    opts.script_name = tb::detail::script_basename(argc > 0 ? argv[0] : "program");
    opts.n_events = -1;
    opts.verbose = false;
    opts.focal = false;
    opts.timewalk = false;

    configure_logger(false);

    cxxopts::Options options(opts.script_name, "Process a single ROOT input file");
    options.add_options()
        ("i,input", "Input ROOT file", cxxopts::value<std::string>())
        ("o,output", "Output ROOT file", cxxopts::value<std::string>())
        ("n,events", "Number of events to process", cxxopts::value<int>()->default_value("-1"))
        ("f,focal", "Enable FoCal mode", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("t,timewalk", "Enable timewalk correction", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("p,pedestal", "Pedestal ROOT file", cxxopts::value<std::string>()->default_value(""))
        ("output-folder", "Override output folder", cxxopts::value<std::string>()->default_value(""))
        ("verbose", "Verbose logging", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("help", "Print help");

    cxxopts::ParseResult result;
    try {
        result = options.parse(argc, argv);
    } catch (const cxxopts::exceptions::exception& err) {
        spdlog::error("{}", err.what());
        std::cout << options.help() << std::endl;
        std::exit(1);
    }

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        std::exit(0);
    }

    if (!result.count("input") || !result.count("output")) {
        spdlog::error("Both --input and --output must be specified.");
        std::cout << options.help() << std::endl;
        std::exit(1);
    }

    opts.input_file  = result["input"].as<std::string>();
    opts.output_file = result["output"].as<std::string>();
    opts.pedestal_file = result["pedestal"].as<std::string>();
    opts.n_events    = result["events"].as<int>();
    opts.focal       = result["focal"].as<bool>();
    opts.timewalk    = result["timewalk"].as<bool>();
    opts.verbose     = result["verbose"].as<bool>();

    tb::detail::validate_existing_file(opts.input_file, "Input file");
    if (!opts.pedestal_file.empty()) {
        tb::detail::validate_existing_file(opts.pedestal_file, "Pedestal file");
    }

    tb::detail::warn_existing_output(opts.output_file);
    tb::detail::ensure_output_folder(opts, result["output-folder"].as<std::string>());

    configure_logger(opts.verbose);
    spdlog::info("{} v{}", opts.script_name, opts.script_version);
    spdlog::info("Input file: {}", opts.input_file);
    spdlog::info("Output file: {}", opts.output_file);

    return opts;
}

inline ScriptOptions parse_arguments_single_json(int argc, char **argv, const std::string& version = "0.1") {
    ScriptOptions opts{};
    opts.script_version = version;
    opts.script_name = tb::detail::script_basename(argc > 0 ? argv[0] : "program");
    opts.n_events = -1;
    opts.verbose = false;
    opts.focal = false;
    opts.timewalk = false;

    configure_logger(false);

    cxxopts::Options options(opts.script_name, "Process a JSON configuration alongside ROOT files");
    options.add_options()
        ("i,input", "Input JSON configuration file", cxxopts::value<std::string>())
        ("o,output", "Output ROOT file", cxxopts::value<std::string>())
        ("p,pedestal", "Pedestal ROOT file", cxxopts::value<std::string>())
        ("c,csv", "CSV summary output", cxxopts::value<std::string>()->default_value(""))
        ("timewalk-file", "Timewalk calibration file", cxxopts::value<std::string>()->default_value(""))
        ("fitting-file", "Fitting configuration file", cxxopts::value<std::string>()->default_value(""))
        ("n,events", "Number of events to process", cxxopts::value<int>()->default_value("-1"))
        ("f,focal", "Enable FoCal mode", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("t,timewalk", "Enable timewalk correction", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("output-folder", "Override output folder", cxxopts::value<std::string>()->default_value(""))
        ("verbose", "Verbose logging", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("help", "Print help");

    cxxopts::ParseResult result;
    try {
        result = options.parse(argc, argv);
    } catch (const cxxopts::exceptions::exception& err) {
        spdlog::error("{}", err.what());
        std::cout << options.help() << std::endl;
        std::exit(1);
    }

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        std::exit(0);
    }

    if (!result.count("input") || !result.count("output") || !result.count("pedestal")) {
        spdlog::error("--input, --output, and --pedestal must be specified.");
        std::cout << options.help() << std::endl;
        std::exit(1);
    }

    opts.input_file     = result["input"].as<std::string>();
    opts.output_file    = result["output"].as<std::string>();
    opts.pedestal_file  = result["pedestal"].as<std::string>();
    opts.csv_file       = result["csv"].as<std::string>();
    opts.timewalk_file  = result["timewalk-file"].as<std::string>();
    opts.fitting_file   = result["fitting-file"].as<std::string>();
    opts.n_events       = result["events"].as<int>();
    opts.focal          = result["focal"].as<bool>();
    opts.timewalk       = result["timewalk"].as<bool>();
    opts.verbose        = result["verbose"].as<bool>();

    tb::detail::validate_existing_file(opts.input_file, "JSON configuration");
    tb::detail::validate_existing_file(opts.pedestal_file, "Pedestal file");
    if (!opts.timewalk_file.empty()) {
        tb::detail::validate_existing_file(opts.timewalk_file, "Timewalk file");
    }
    if (!opts.fitting_file.empty()) {
        tb::detail::validate_existing_file(opts.fitting_file, "Fitting file");
    }

    tb::detail::warn_existing_output(opts.output_file);
    tb::detail::ensure_output_folder(opts, result["output-folder"].as<std::string>());

    configure_logger(opts.verbose);
    spdlog::info("{} v{}", opts.script_name, opts.script_version);
    spdlog::info("Configuration file: {}", opts.input_file);
    spdlog::info("Output file: {}", opts.output_file);

    return opts;
}

inline ScriptOptions parse_arguments_single_root_single_csv(int argc, char **argv, const std::string& version = "0.1") {
    ScriptOptions opts{};
    opts.script_version = version;
    opts.script_name = tb::detail::script_basename(argc > 0 ? argv[0] : "program");
    opts.n_events = -1;
    opts.verbose = false;
    opts.focal = false;
    opts.timewalk = false;

    configure_logger(false);

    cxxopts::Options options(opts.script_name, "Process a ROOT file and emit CSV summaries");
    options.add_options()
        ("i,input", "Input ROOT file", cxxopts::value<std::string>())
        ("o,output", "Output ROOT file", cxxopts::value<std::string>()->default_value(""))
        ("c,csv", "CSV output file", cxxopts::value<std::string>())
        ("n,events", "Number of events", cxxopts::value<int>()->default_value("-1"))
        ("output-folder", "Override output folder", cxxopts::value<std::string>()->default_value(""))
        ("verbose", "Verbose logging", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("help", "Print help");

    cxxopts::ParseResult result;
    try {
        result = options.parse(argc, argv);
    } catch (const cxxopts::exceptions::exception& err) {
        spdlog::error("{}", err.what());
        std::cout << options.help() << std::endl;
        std::exit(1);
    }

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        std::exit(0);
    }

    if (!result.count("input") || !result.count("csv")) {
        spdlog::error("--input and --csv must be specified.");
        std::cout << options.help() << std::endl;
        std::exit(1);
    }

    opts.input_file  = result["input"].as<std::string>();
    opts.csv_file    = result["csv"].as<std::string>();
    opts.output_file = result["output"].as<std::string>();
    opts.n_events    = result["events"].as<int>();
    opts.verbose     = result["verbose"].as<bool>();

    tb::detail::validate_existing_file(opts.input_file, "Input file");
    tb::detail::warn_existing_output(opts.output_file);
    tb::detail::ensure_output_folder(opts, result["output-folder"].as<std::string>());

    configure_logger(opts.verbose);
    spdlog::info("{} v{}", opts.script_name, opts.script_version);
    spdlog::info("Input file: {}", opts.input_file);
    if (!opts.output_file.empty()) {
        spdlog::info("Output file: {}", opts.output_file);
    }
    spdlog::info("CSV file: {}", opts.csv_file);

    return opts;
}

int get_valid_fpga_channel(int fpga_channel);
int get_total_fpga_channel(int fpga_channel);
inline int get_unified_fpga_channel(int fpga_id, int fpga_channel){
    return fpga_id * FPGA_CHANNEL_NUMBER + fpga_channel;
}
inline int get_unified_valid_fpga_channel(int fpga_id, int fpga_channel){
    return fpga_id * FPGA_CHANNEL_NUMBER_VALID + fpga_channel;
}

inline UInt_t decode_tot_value(UInt_t val1) {
    UInt_t mask = -(val1 >= 512);  // mask is 0xFFFFFFFF if val1 >= 512, else 0
    return (val1 & ~mask) | ((val1 - 512) * 8 & mask);
}

inline double decode_toa_value_ns(UInt_t val2) {
    constexpr double scale0 = 0.025;
    constexpr double scale1 = 0.2;
    constexpr double scale2 = 6.25;

    UInt_t part0 = val2 & 0x07;
    UInt_t part1 = (val2 >> 3) & 0x1F;
    UInt_t part2 = (val2 >> 8) & 0x1F;

    return part0 * scale0 + part1 * scale1 + part2 * scale2;
}

inline void toa_manual_correction(int vldb_id, int asic_id, int half_id, UInt_t& toa_max, double& toa_value_ns) {
    // Apply manual correction to ToA value based on known issues
    // This is a placeholder for actual correction logic
    if (vldb_id == 1 && asic_id == 1) {
        if (toa_max < 260) {
            toa_value_ns += 25.0;  // correct for the known offset issue
        }
    }
    if (vldb_id == 1 && asic_id == 0) {
        if (toa_max > 768) {
            toa_value_ns -= 25.0;  // correct for the known offset issue
        }
    }
    if (vldb_id == 0 && half_id == 0) {
        if (toa_max > 544) {
            toa_value_ns -= 25.0;  // correct for the known offset issue
        }
    }
    if (vldb_id == 0 && half_id == 1) {
        if (toa_max > 780) {
            toa_value_ns -= 25.0;  // correct for the known offset issue
        }
    }
}

inline double find_closest_toa(const json& timewalk_json, double adc_max) {
    // Fast checks and one-time extraction
    auto it_adc = timewalk_json.find("adc_max_minus_pedestal");
    auto it_toa = timewalk_json.find("toa_ns");
    if (it_adc == timewalk_json.end() || it_toa == timewalk_json.end()
        || !it_adc->is_array() || !it_toa->is_array()) {
        return 0.0;
    }

    // One conversion (avoid converting per element)
    // NOTE: This copies once. If you call this many times, pull these out once upstream and
    // pass the vectors by reference instead.
    std::vector<double> adc = it_adc->get<std::vector<double>>();
    std::vector<double> toa = it_toa->get<std::vector<double>>();

    const size_t n = adc.size();
    if (n == 0 || toa.size() != n) return 0.0;

    // If ADC is sorted (common in lookup tables), do O(log N) search
    if (std::is_sorted(adc.begin(), adc.end())) {
        auto it = std::lower_bound(adc.begin(), adc.end(), adc_max);
        size_t idx;
        if (it == adc.begin()) {
            idx = 0;
        } else if (it == adc.end()) {
            idx = n - 1;
        } else {
            // Pick the closer of neighbors
            const double hi = *it;
            const double lo = *(it - 1);
            idx = (adc_max - lo <= hi - adc_max) ? size_t((it - 1) - adc.begin())
                                                 : size_t(it - adc.begin());
        }
        return toa[idx];
    }

    // Otherwise, tight linear scan (no extra JSON lookups)
    size_t best = 0;
    double best_diff = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < n; ++i) {
        double d = std::abs(adc[i] - adc_max);
        if (d < best_diff) { best_diff = d; best = i; }
    }
    return toa[best];
}

class GlobalChannelPainter {
public:
    GlobalChannelPainter(const std::string& mapping_file);
    GlobalChannelPainter(const std::string& mapping_file, const std::string& channel_mapping_file);
    ~GlobalChannelPainter();

    TCanvas* get_canvas() { return painter_canvas; }
    void draw_global_channel_hists2D(std::vector <TH2D*> hists, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title);
    void draw_global_channel_hists1D(std::vector <TH1D*> hists, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title);
    void draw_global_channel_hists1D_group(std::vector <std::vector <TH1D*>> hists_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, std::vector <EColor> colors, std::vector <std::string> legend_labels);
    void draw_global_channel_hists1D_run_group(std::vector <std::vector <TH1D*>> hists_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, std::vector <EColor> colors, std::vector <std::string> legend_labels);
    void draw_global_channel_canvas(std::vector <TCanvas*> canvas_list,  std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title);
    void draw_canvas_components(TCanvas* source_canvas, TPad* target_pad);

    TCanvas* draw_module_channel_canvas(std::vector <TCanvas*> canvas_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, int module_index);

private:
    void clear_canvas();

private:
    TCanvas *painter_canvas;
    json mapping_json;

    json focal_mapping_json;
    json focal_channel_mapping_json;

    bool is_EEEMCal_mapping;
    bool is_FoCal_mapping;

    std::vector <int> module_fpga_list;
    std::vector <int> module_asic_list;
    std::vector <int> module_connector_list;
    std::vector <std::vector <int>> connector_list_list;

    std::vector <std::vector <int>> focal_module_board_list;
    std::vector <std::vector <int>> focal_module_channel_list;
    std::unordered_map <int, int> focal_channel_map;
    std::unordered_map <int, int> focal_fpga_map;

    std::vector <TCanvas*> sub_canvas_list;
};


inline int64_t posmod(int64_t a, int64_t m) {
    int64_t r = a % m;
    return (r < 0) ? (r + m) : r; // [0, m)
}

inline bool is_0_or_half_mod(int64_t value, int64_t period, int64_t tol) {
    // distance to 0 (i.e., k*period)
    int64_t r = posmod(value, period);
    int64_t dist0 = std::min(r, period - r);

    // distance to half-period (i.e., k*period + period/2)
    int64_t half = period / 2; // assume even period (40)
    int64_t d = std::llabs(r - half);
    int64_t distHalf = std::min(d, period - d);

    return (dist0 <= tol) || (distHalf <= tol);
}

// Accept normal, half-period, and ±1-bit slip (×2 or ÷2) cases
inline bool is_aligned_with_bit_slip(int64_t diff, int64_t period, int64_t tol) {
    if (is_0_or_half_mod(diff, period, tol)) return true;

    // Consider a left-shift (current diff ≈ 2 * correct) → undo with /2
    int64_t halfDiff = diff / 2; // truncates toward 0; tol should cover 1-tick rounding
    if (is_0_or_half_mod(halfDiff, period, tol)) return true;

    // Consider a right-shift (current diff ≈ correct / 2) → undo with *2 (safely)
    if (std::llabs(diff) <= (std::numeric_limits<int64_t>::max() / 2)) {
        int64_t doubleDiff = diff * 2;
        if (is_0_or_half_mod(doubleDiff, period, tol)) return true;
    }

    return false;
}

static std::string format_decimal(double value, int precision = 2)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

double quantile_from_tf1(TF1* f, double q, double xmin, double xmax)
{
    const double I_tot = f->Integral(xmin, xmax);
    double lo = xmin, hi = xmax;
    for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (lo + hi);
        double I   = f->Integral(xmin, mid);
        if (I / I_tot < q)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

struct MeanRMS90 {
    // central values
    double mean90;      
    double rms90;       
    double sigma_gauss; // = rms90 / k90
    double res_sigma_over_mean90; // (sigma_gauss / mean90)

    // Gaussian-approx stat errors
    double err_mean90;  
    double err_rms90;   
    double err_sigma_gauss;
    double err_res_sigma_over_mean90;

    // window info
    int    binLo, binHi;
    int    Neff;        // effective events inside the 90% window
};

MeanRMS90 computeMeanRMS90(const TH1* h, 
                           double frac = 0.90, 
                           bool useBinWidth = false,
                           long long nEventsOverride = -1)
{
    const int nb = h->GetNbinsX();
    MeanRMS90 out{};
    if (nb <= 0) { out.mean90 = out.rms90 = out.sigma_gauss = out.res_sigma_over_mean90 = NAN; return out; }

    // 预取边界、中心、权重（用于均值/RMS 计算）
    std::vector<double> edge(nb+1), x(nb+1), w(nb+1);
    for (int i=1; i<=nb; ++i) {
        edge[i-1] = h->GetXaxis()->GetBinLowEdge(i);
        x[i]      = h->GetXaxis()->GetBinCenter(i);
        double bw = h->GetXaxis()->GetBinWidth(i);
        double cnt= h->GetBinContent(i);
        w[i]      = useBinWidth ? (cnt * bw) : cnt; // 仅用于均值/RMS 计算
    }
    edge[nb] = h->GetXaxis()->GetBinUpEdge(nb);

    // 总计权重（若是计数直方图，就是总事件数）
    double totalW = 0.0;
    for (int i=1; i<=nb; ++i) totalW += w[i];
    if (totalW <= 0) { out.mean90 = out.rms90 = out.sigma_gauss = out.res_sigma_over_mean90 = NAN; return out; }

    const double targetW = frac * totalW;

    // 前缀和
    std::vector<double> csum(nb+1, 0.0);
    for (int i=1; i<=nb; ++i) csum[i] = csum[i-1] + w[i];

    // 双指针找“包含 frac 的最窄区间”
    int bestL = 1, bestR = nb;
    double bestWidth = edge[nb] - edge[0];

    int R = 1;
    for (int L = 1; L <= nb; ++L) {
        if (R < L) R = L;
        while (R <= nb && (csum[R] - csum[L-1]) < targetW) ++R;
        if (R > nb) break;
        double width = edge[R] - edge[L-1];
        if (width < bestWidth) { bestWidth = width; bestL = L; bestR = R; }
    }

    // 计算 Mean90 / RMS90
    double W = 0.0, Sx = 0.0, Sx2 = 0.0;
    for (int i=bestL; i<=bestR; ++i) {
        W   += w[i];
        Sx  += w[i] * x[i];
        Sx2 += w[i] * x[i] * x[i];
    }
    double mean90 = (W > 0) ? (Sx / W) : NAN;
    double var90  = (W > 0) ? (Sx2 / W - mean90*mean90) : NAN;
    double rms90  = (var90 > 0) ? std::sqrt(var90) : 0.0;

    // 高斯等效
    constexpr double k90 = 0.7893132990990863; // precise coefficient for central 90%
    double sigma_gauss = rms90 / k90;
    double res = (!std::isnan(mean90) && mean90 != 0.0) ? (sigma_gauss / mean90) : NAN;

    // -------- 高斯近似统计误差 ----------
    // 需要 N_eff：90% 区间内的“事件数”。对于计数直方图，W≈N_eff。
    long long Neff = -1;
    if (!useBinWidth && nEventsOverride < 0) {
        // 计数直方图：直接用该窗口内的计数总和
        double sumCounts = 0.0;
        for (int i=bestL; i<=bestR; ++i) sumCounts += h->GetBinContent(i);
        Neff = static_cast<long long>(std::llround(sumCounts));
    } else if (nEventsOverride >= 0) {
        // 提供了总事件数：用 0.9*N 作为 Neff 的近似
        Neff = static_cast<long long>(std::llround(frac * nEventsOverride));
    } else {
        // 无法确定 Neff：给出 NAN 提醒
        Neff = -1;
    }

    double err_rms90 = NAN, err_mean90 = NAN, err_sigma = NAN, err_res = NAN;
    if (Neff > 1 && std::isfinite(rms90)) {
        // ΔRMS90/RMS90 ≈ 1/sqrt(2*(Neff-1))
        double rel = 1.0 / std::sqrt(2.0 * (double)(Neff - 1));
        err_rms90  = rms90 * rel;

        // ΔMean90 ≈ RMS90 / sqrt(Neff)
        err_mean90 = rms90 / std::sqrt((double)Neff);

        // Δσ_Gauss = ΔRMS90 / k90
        err_sigma  = err_rms90 / k90;

        // Δ(R=σ/Mean) 误差传播（忽略两者相关性的近似）
        if (std::isfinite(res) && mean90 != 0.0 && sigma_gauss != 0.0) {
            double rel_sigma = err_sigma / sigma_gauss;
            double rel_mean  = err_mean90 / mean90;
            err_res = res * std::sqrt(rel_sigma*rel_sigma + rel_mean*rel_mean);
        }
    }

    out.mean90  = mean90;
    out.rms90   = rms90;
    out.sigma_gauss = sigma_gauss;
    out.res_sigma_over_mean90 = res;

    out.err_mean90 = err_mean90;
    out.err_rms90  = err_rms90;
    out.err_sigma_gauss = err_sigma;
    out.err_res_sigma_over_mean90 = err_res;

    out.binLo = bestL;
    out.binHi = bestR;
    out.Neff  = (int)Neff;
    return out;
}

#endif // COMMON_HPP