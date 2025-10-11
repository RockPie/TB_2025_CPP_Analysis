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
#define TESTBEAM_TITLE "CERN SPS H2 2025 October Testbeam"

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

#endif // COMMON_HPP