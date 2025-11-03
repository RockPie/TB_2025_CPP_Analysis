#include <cmath>
#include <fstream>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"

double resolution_func(double *x, double *par) {
    // par[0] = p0 (constant term)
    // par[1] = p1 (stochastic term)
    // par[2] = p2 (noise term)
    double E = x[0];
    if (E == 0) return 0;
    double res = std::sqrt( std::pow(par[0], 2) + std::pow(par[1]/std::sqrt(E), 2) + std::pow(par[2]/E, 2) );
    return res;
}

int main(int argc, char **argv) {
    gROOT->SetBatch(kTRUE);
    std::string script_input_file, script_output_file;
    int script_n_events = -1;
    bool script_verbose = false;
    std::string script_name = __FILE__;
    const std::string script_version = "0.1";
    std::string script_output_folder;

    script_name = script_name.substr(script_name.find_last_of("/\\") + 1).substr(0, script_name.find_last_of("."));

    configure_logger(false);
    
    cxxopts::Options options(script_name, "Generate heatmaps from machine gun data");

    options.add_options()
        ("f,file", "Input .json config file", cxxopts::value<std::string>())
        ("o,output", "Output .root file", cxxopts::value<std::string>())
        ("e,events", "Number of events to process", cxxopts::value<int>()->default_value("-1"))
        ("v,verbose", "Verbose mode", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("h,help", "Print help");

    cxxopts::ParseResult parsed;
    try {
        parsed = options.parse(argc, argv);
    } catch (const cxxopts::exceptions::exception& err) {
        spdlog::error("{}", err.what());
        std::cout << options.help() << std::endl;
        return 1;
    }

    if (parsed.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (!parsed.count("file") || !parsed.count("output")) {
        spdlog::error("Both --file and --output must be specified");
        std::cout << options.help() << std::endl;
        return 1;
    }

    // config/106_scan_0.json
    script_input_file  = parsed["file"].as<std::string>();
    json config_json;
    std::ifstream config_ifs(script_input_file);
    if (!config_ifs.is_open()) {
        spdlog::error("Failed to open config json file {}", script_input_file);
        return 1;
    }
    config_ifs >> config_json;
    config_ifs.close();
    script_output_file = parsed["output"].as<std::string>();
    script_n_events    = parsed["events"].as<int>();
    script_verbose     = parsed["verbose"].as<bool>();

    configure_logger(script_verbose);

    // load configuration from json
    std::vector<int> run_numbers;
    // it's either run_energies or run_labels giving the information
    // look for run_energies first
    // if not found, look for run_labels
    std::vector<int> run_energies;
    std::vector<std::string> run_labels;
    std::string beam_type;
    std::string config_description;
    bool use_peak_integral = false;
    bool enable_gaussian_fit = false;
    int target_event_number = -1;
    bool found_run_energies = false;
    try {
        run_numbers = config_json.at("run_numbers").get<std::vector<int>>();
        if (config_json.contains("run_energies")) {
            run_energies = config_json.at("run_energies").get<std::vector<int>>();
            found_run_energies = true;
        } else if (config_json.contains("run_labels")) {
            run_labels = config_json.at("run_labels").get<std::vector<std::string>>();
        } else {
            spdlog::error("Config json file {} must contain either run_energies or run_labels!", script_input_file);
            return 1;
        }
        if (config_json.contains("enable_gaussian_fit")) {
            enable_gaussian_fit = config_json.at("enable_gaussian_fit").get<bool>();
        }
        if (config_json.contains("use_peak_integral")) {
            use_peak_integral = config_json.at("use_peak_integral").get<bool>();
        }
        beam_type = config_json.at("beam_type").get<std::string>();
        config_description = config_json.at("description").get<std::string>();
        if (config_json.contains("n_events")) {
            target_event_number = config_json.at("n_events").get<int>();
        }
    } catch (const json::exception& e) {
        spdlog::error("Failed to parse config json file {}: {}", script_input_file, e.what());
        return 1;
    }

    std::vector<std::string> run_files;

    const std::string run_file_folder = "dump/304_RawADC/beamtests/";
    // try to find the input root file from run number
    if (run_numbers.empty()) {
        spdlog::error("No run numbers specified in config json file {}", script_input_file);
        return 1;
    }
    // run files always 4 digits
    for (const auto& run_number : run_numbers) {
        std::string run_file = run_file_folder + "Run" + fmt::format("{:04d}", run_number) + ".root";
        if (access(run_file.c_str(), F_OK) == -1) {
            spdlog::error("Run file {} does not exist!", run_file);
            return 1;
        } else {
            spdlog::info("Found run file: {}", run_file);
        }
        run_files.push_back(run_file);
    }

    // sort the run files by run energy descending, if run energies are provided
    if (found_run_energies){
        std::vector<size_t> indices(run_files.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&run_energies](size_t i1, size_t i2) {
            return run_energies[i1] > run_energies[i2];
        });
        std::vector<std::string> sorted_run_files;
        std::vector<int> sorted_run_energies;
        std::vector<int> sorted_run_numbers;
        for (const auto& index : indices) {
            sorted_run_files.push_back(run_files[index]);
            sorted_run_energies.push_back(run_energies[index]);
            sorted_run_numbers.push_back(run_numbers[index]);
        }
        run_files = sorted_run_files;
        run_energies = sorted_run_energies;
        run_numbers = sorted_run_numbers;
    }
    

    // print configuration
    spdlog::info("Script version: {}", script_version);
    spdlog::info("Run numbers: {}", fmt::join(run_numbers, ", "));
    if (use_peak_integral)
        spdlog::info("Using peak integral for ADC sum");
    spdlog::info("Beam type: {}", beam_type);
    if (found_run_energies)
        spdlog::info("Run energies: {}", fmt::join(run_energies, ", "));
    else
        spdlog::info("Run labels: {}", fmt::join(run_labels, ", "));

    if (enable_gaussian_fit)
        spdlog::info("Gaussian fit enabled");
    else
        spdlog::info("Gaussian fit disabled");

    spdlog::info("Input config file: {}", script_input_file);

    if (script_output_file.substr(script_output_file.find_last_of(".") + 1) != "root") {
        spdlog::error("Output file {} should end with .root!", script_output_file);
        return 1;
    }

    script_output_folder = script_output_file.substr(0, script_output_file.find_last_of("/\\"));
    if (script_output_folder.empty()) {
        script_output_folder = "./dump/" + script_name;
    }
    if (access(script_output_folder.c_str(), F_OK) == -1) {
        spdlog::info("Creating output folder {}", script_output_folder);
        if (mkdir(script_output_folder.c_str(), 0777) == -1) {
            spdlog::error("Failed to create output folder {}", script_output_folder);
            return 1;
        }
    }
    if (access(script_output_file.c_str(), F_OK) != -1) {
        spdlog::warn("Output file {} already exists!", script_output_file);
    }

    spdlog::info("Script name: {}", script_name);
    spdlog::info("Input config file: {}", script_input_file);
    spdlog::info("Output file: {} in {}", script_output_file, script_output_folder);
    spdlog::info("Number of events: {}", script_n_events);

    // std::string mapping_json_file = "config/mapping_Oct2025.json";
    // std::ifstream mapping_json_ifs(mapping_json_file);
    // if (!mapping_json_ifs.is_open()) {
    //     spdlog::error("Failed to open mapping json file {}", mapping_json_file);
    //     return 1;
    // }
    // json mapping_json;
    // mapping_json_ifs >> mapping_json;
    // const auto& sipm_board = mapping_json.at("SiPM_Board");
    // const auto& board_loc = mapping_json.at("Board_Loc");
    // const auto& board_rotation = mapping_json.at("Board_Rotation");
    // const auto& board_flip = mapping_json.at("Board_Flip");

    // * --- Read input file ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    const int machine_gun_samples = 16;
    const int vldb_number = 2;

    // first go through all run files to get total number of events
    std::vector<int> run_file_event_counts;
    run_file_event_counts.reserve(run_files.size());
    for (const auto& run_file : run_files) {
        TFile *input_root = new TFile(run_file.c_str(), "READ");
        if (input_root->IsZombie()) {
            spdlog::error("Failed to open input file {}", run_file);
            return 1;
        }
        input_root->Close();
    }
    // TParameter<double> param_gaus_fit_mean("gaus_fit_mean", mean_avg); 
    // TParameter<double> param_gaus_fit_sigma("gaus_fit_sigma", sigma_avg);
    // TParameter<double> param_gaus_fit_resolution("gaus_fit_resolution", resolution);
    // TParameter<double> param_gaus_fit_mean_err_stat("gaus_fit_mean_err_stat", mean_err_stat);
    // TParameter<double> param_gaus_fit_mean_err_sys("gaus_fit_mean_err_sys", mean_err_sys);
    // TParameter<double> param_gaus_fit_sigma_err_stat("gaus_fit_sigma_err_stat", sigma_err);
    // TParameter<double> param_gaus_fit_sigma_err_sys("gaus_fit_sigma_err_sys", sigma_err);
    // TParameter<double> param_gaus_fit_resolution_err("gaus_fit_resolution_err", resolution_err);

    // TParameter<double> param_cb_fit_mean("cb_fit_mean", cb_mean_avg);
    // TParameter<double> param_cb_fit_sigma("cb_fit_sigma", cb_sigma_avg);
    // TParameter<double> param_cb_fit_effective_sigma("cb_fit_effective_sigma", cb_effective_sigma_avg);
    // TParameter<double> param_cb_fit_resolution("cb_fit_resolution", cb_resolution);
    // TParameter<double> param_cb_fit_mean_err_stat("cb_fit_mean_err_stat", cb_mean_err_stat);
    // TParameter<double> param_cb_fit_mean_err_sys("cb_fit_mean_err_sys", cb_mean_err_sys);
    // TParameter<double> param_cb_fit_sigma_err_stat("cb_fit_sigma_err_stat", cb_sigma_err);
    // TParameter<double> param_cb_fit_sigma_err_sys("cb_fit_sigma_err_sys", cb_sigma_err);
    // TParameter<double> param_cb_fit_effective_sigma_err_stat("cb_fit_effective_sigma_err_stat", cb_effective_sigma_err_stat);
    // TParameter<double> param_cb_fit_effective_sigma_err_sys("cb_fit_effective_sigma_err_sys", cb_effective_sigma_err_sys);
    // TParameter<double> param_cb_fit_resolution_err("cb_fit_resolution_err", cb_resolution_err);

    std::vector<double> list_param_gaus_fit_mean;
    std::vector<double> list_param_gaus_fit_sigma;
    std::vector<double> list_param_gaus_fit_resolution;
    std::vector<double> list_param_gaus_fit_mean_err_stat;
    std::vector<double> list_param_gaus_fit_mean_err_sys;
    std::vector<double> list_param_gaus_fit_sigma_err_stat;
    std::vector<double> list_param_gaus_fit_sigma_err_sys;
    std::vector<double> list_param_gaus_fit_resolution_err;

    std::vector<double> list_param_cb_fit_mean;;
    std::vector<double> list_param_cb_fit_sigma;
    std::vector<double> list_param_cb_fit_effective_sigma;
    std::vector<double> list_param_cb_fit_resolution;
    std::vector<double> list_param_cb_fit_effective_resolution;
    std::vector<double> list_param_cb_fit_mean_err_stat;
    std::vector<double> list_param_cb_fit_mean_err_sys;
    std::vector<double> list_param_cb_fit_sigma_err_stat;
    std::vector<double> list_param_cb_fit_sigma_err_sys;
    std::vector<double> list_param_cb_fit_effective_sigma_err_stat;
    std::vector<double> list_param_cb_fit_effective_sigma_err_sys;
    std::vector<double> list_param_cb_fit_resolution_err;
    std::vector<double> list_param_cb_fit_effective_resolution_err;

    auto getParam = [](TDirectory* dir, const char* name, bool required = true) -> double {
        TParameter<double>* p = nullptr;
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
        p = dir->Get<TParameter<double>>(name);
    #else
        dir->GetObject(name, p);
    #endif
        if (!p) {
            if (required) {
                spdlog::error("Missing TParameter<double> key '{}' in directory '{}'",
                            name, dir->GetPath());
                throw std::runtime_error(std::string("Missing key: ") + name);
            } else {
                spdlog::warn("Missing optional key '{}' in directory '{}'; will use NaN.",
                            name, dir->GetPath());
                return std::numeric_limits<double>::quiet_NaN();
            }
        }
        return p->GetVal();
    };

    for (int _file_index = 0; _file_index < static_cast<int>(run_files.size()); ++_file_index) {
        const auto& run_file = run_files[_file_index];
        spdlog::info("Processing run file: {}", run_file);

        TFile* input_root = TFile::Open(run_file.c_str(), "READ");
        if (!input_root || input_root->IsZombie()) {
            spdlog::error("Failed to open input file {}", run_file);
            if (input_root) { input_root->Close(); delete input_root; }
            continue; // 或者 return 1; 视你整体流程而定
        }

        // 如果参数写在子目录里，比如 "Results"；否则就用根目录
        TDirectory* dir = input_root; 
        // TDirectory* dir = input_root->GetDirectory("Results"); // 若你写在子目录，请改用这行
        if (!dir) {
            spdlog::error("Directory with parameters is missing in file {}", run_file);
            input_root->Close(); delete input_root;
            continue;
        }

        try {
            double v_gaus_fit_mean             = getParam(dir, "gaus_fit_mean",             /*required=*/true);
            double v_gaus_fit_sigma            = getParam(dir, "gaus_fit_sigma",            /*required=*/true);
            double v_gaus_fit_resolution       = getParam(dir, "gaus_fit_resolution",       /*required=*/false);
            double v_gaus_fit_mean_err_stat    = getParam(dir, "gaus_fit_mean_err_stat",    /*required=*/false);
            double v_gaus_fit_mean_err_sys     = getParam(dir, "gaus_fit_mean_err_sys",     /*required=*/false);
            double v_gaus_fit_sigma_err_stat   = getParam(dir, "gaus_fit_sigma_err_stat",   /*required=*/false);
            double v_gaus_fit_sigma_err_sys    = getParam(dir, "gaus_fit_sigma_err_sys",    /*required=*/false);
            double v_gaus_fit_resolution_err   = getParam(dir, "gaus_fit_resolution_err",   /*required=*/false);

            list_param_gaus_fit_mean.push_back(v_gaus_fit_mean);
            list_param_gaus_fit_sigma.push_back(v_gaus_fit_sigma);
            list_param_gaus_fit_resolution.push_back(v_gaus_fit_resolution);
            list_param_gaus_fit_mean_err_stat.push_back(v_gaus_fit_mean_err_stat);
            list_param_gaus_fit_mean_err_sys.push_back(v_gaus_fit_mean_err_sys);
            list_param_gaus_fit_sigma_err_stat.push_back(v_gaus_fit_sigma_err_stat);
            list_param_gaus_fit_sigma_err_sys.push_back(v_gaus_fit_sigma_err_sys);
            list_param_gaus_fit_resolution_err.push_back(v_gaus_fit_resolution_err);

            double v_cb_fit_mean                   = getParam(dir, "cb_fit_mean",                   /*required=*/true);
            double v_cb_fit_sigma                  = getParam(dir, "cb_fit_sigma",                  /*required=*/true);
            double v_cb_fit_effective_sigma        = getParam(dir, "cb_fit_effective_sigma",        /*required=*/false);
            double v_cb_fit_resolution             = getParam(dir, "cb_fit_resolution",             /*required=*/false);
            double v_cb_fit_effective_resolution    = getParam(dir, "cb_fit_effective_resolution",    /*required=*/false);
            double v_cb_fit_mean_err_stat          = getParam(dir, "cb_fit_mean_err_stat",          /*required=*/false);
            double v_cb_fit_mean_err_sys           = getParam(dir, "cb_fit_mean_err_sys",           /*required=*/false);
            double v_cb_fit_sigma_err_stat         = getParam(dir, "cb_fit_sigma_err_stat",         /*required=*/false);
            double v_cb_fit_sigma_err_sys          = getParam(dir, "cb_fit_sigma_err_sys",          /*required=*/false);
            double v_cb_fit_effective_sigma_err_stat = getParam(dir, "cb_fit_effective_sigma_err_stat", /*required=*/false);
            double v_cb_fit_effective_sigma_err_sys  = getParam(dir, "cb_fit_effective_sigma_err_sys",  /*required=*/false);
            double v_cb_fit_resolution_err         = getParam(dir, "cb_fit_resolution_err",         /*required=*/false);
            double v_cb_fit_effective_resolution_err = getParam(dir, "cb_fit_effective_resolution_err", /*required=*/false);

            list_param_cb_fit_mean.push_back(v_cb_fit_mean);
            list_param_cb_fit_sigma.push_back(v_cb_fit_sigma);
            list_param_cb_fit_effective_sigma.push_back(v_cb_fit_effective_sigma);
            list_param_cb_fit_resolution.push_back(v_cb_fit_resolution);
            list_param_cb_fit_effective_resolution.push_back(v_cb_fit_effective_resolution);
            list_param_cb_fit_mean_err_stat.push_back(v_cb_fit_mean_err_stat);
            list_param_cb_fit_mean_err_sys.push_back(v_cb_fit_mean_err_sys);
            list_param_cb_fit_sigma_err_stat.push_back(v_cb_fit_sigma_err_stat);
            list_param_cb_fit_sigma_err_sys.push_back(v_cb_fit_sigma_err_sys);
            list_param_cb_fit_effective_sigma_err_stat.push_back(v_cb_fit_effective_sigma_err_stat);
            list_param_cb_fit_effective_sigma_err_sys.push_back(v_cb_fit_effective_sigma_err_sys);
            list_param_cb_fit_resolution_err.push_back(v_cb_fit_resolution_err);
            list_param_cb_fit_effective_resolution_err.push_back(v_cb_fit_effective_resolution_err);

        } catch (const std::exception& e) {
            spdlog::error("While reading parameters from {}: {}", run_file, e.what());
            input_root->Close(); delete input_root;
            continue; // 或者 return 1;
        }

        input_root->Close();
        delete input_root;
    }

    // make a color list
    std::vector<int> color_list = {kRed, kBlue, kGreen+2, kMagenta+2, kCyan+2, kOrange+7, kViolet+7, kSpring+7, kTeal+7, kAzure+7, kPink+7};

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    
    auto output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }
    output_root->cd();
    const std::string out_pdf = script_output_file.substr(0, script_output_file.find_last_of(".")) + ".pdf";

    double energy_axis_min = 40.0;
    double energy_axis_max = 400.0;

    double axis_title_size = 0.04;
    double axis_label_size = 0.035;
    double legend_text_size = 0.025;

    // ! === draw beam energy vs gaussian fit mean ====================
    TCanvas* canvas_linearity = new TCanvas("canvas_linearity", "Gaussian Fit Mean vs Beam Energy", 800, 600);
    // leave some margin for legend
    canvas_linearity->SetLeftMargin(0.15);
    canvas_linearity->SetBottomMargin(0.10);
    TGraphErrors* graph_linearity_gaus = new TGraphErrors();
    TGraphErrors* graph_linearity_cb = new TGraphErrors();
    for (size_t i = 0; i < list_param_gaus_fit_mean.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y_gaus = list_param_gaus_fit_mean[i];
        double y_cb = list_param_cb_fit_mean[i];
        double yerr_gaus = list_param_gaus_fit_mean_err_stat[i];
        double yerr_cb = list_param_cb_fit_mean_err_stat[i];
        graph_linearity_gaus->SetPoint(i, x, y_gaus);
        graph_linearity_gaus->SetPointError(i, 0.03 * x, yerr_gaus);
        graph_linearity_cb->SetPoint(i, x, y_cb);
        graph_linearity_cb->SetPointError(i, 0.03 * x, yerr_cb);
    }
    graph_linearity_gaus->SetMarkerStyle(20);
    graph_linearity_gaus->SetMarkerColor(kCyan+2);
    graph_linearity_gaus->SetLineColor(kCyan+2);
    graph_linearity_gaus->SetTitle(";Beam Energy [GeV];Gaussian Fit Mean [ADC counts]");
    graph_linearity_gaus->GetXaxis()->SetTitleSize(axis_title_size);
    graph_linearity_gaus->GetYaxis()->SetTitleSize(axis_title_size);
    graph_linearity_gaus->GetXaxis()->SetLabelSize(axis_label_size);
    graph_linearity_gaus->GetYaxis()->SetLabelSize(axis_label_size);
    graph_linearity_gaus->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);

    // linear fit
    TF1* fit_linearity_gaus = new TF1("fit_linearity_gaus", "pol1", energy_axis_min, energy_axis_max);
    fit_linearity_gaus->SetLineColor(kCyan+2);
    // set initial parameters from the average of the data points
    double init_a = (list_param_gaus_fit_mean.back() - list_param_gaus_fit_mean.front()) /
                    ( (found_run_energies ? static_cast<double>(run_energies.back()) : static_cast<double>(list_param_gaus_fit_mean.size()-1)) -
                      (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0) );
    double init_b = list_param_gaus_fit_mean.front() - init_a * (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0);
    fit_linearity_gaus->SetParameters(init_b, init_a);
    graph_linearity_gaus->Fit(fit_linearity_gaus, "R");
    double gaus_fit_a = fit_linearity_gaus->GetParameter(1);;
    double gaus_fit_b = fit_linearity_gaus->GetParameter(0);
    graph_linearity_gaus->Draw("AP");
    fit_linearity_gaus->Draw("SAME");

    graph_linearity_cb->SetMarkerStyle(21);
    graph_linearity_cb->SetMarkerColor(kMagenta+2);
    graph_linearity_cb->SetLineColor(kMagenta+2);

    // linear fit
    TF1* fit_linearity_cb = new TF1("fit_linearity_cb", "pol1", energy_axis_min, energy_axis_max);
    fit_linearity_cb->SetLineColor(kMagenta+2);
    // set initial parameters from the average of the data points
    init_a = (list_param_cb_fit_mean.back() - list_param_cb_fit_mean.front()) /
                    ( (found_run_energies ? static_cast<double>(run_energies.back()) : static_cast<double>(list_param_cb_fit_mean.size()-1)) -
                      (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0) );
    init_b = list_param_cb_fit_mean.front() - init_a * (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0);
    fit_linearity_cb->SetParameters(init_b, init_a);
    graph_linearity_cb->Fit(fit_linearity_cb, "R");
    double cb_fit_a = fit_linearity_cb->GetParameter(1);;
    double cb_fit_b = fit_linearity_cb->GetParameter(0);
    graph_linearity_cb->Draw("P SAME");
    fit_linearity_cb->Draw("SAME");
    // bottom right legend
    TLegend* legend_linearity = new TLegend(0.55, 0.15, 0.90, 0.45);
    legend_linearity->SetBorderSize(0);
    legend_linearity->SetFillStyle(0);
    legend_linearity->SetTextSize(legend_text_size);
    legend_linearity->AddEntry(graph_linearity_gaus, "Gaussian Fit Mean", "lp");
    // add fit function info to legend
    legend_linearity->AddEntry(fit_linearity_gaus, fmt::format("Fit: ADC = {:.2f}E + {:.2f}", gaus_fit_a, gaus_fit_b).c_str(), "l");
    legend_linearity->AddEntry(graph_linearity_cb, "Crystal Ball Fit Mean", "lp");
    legend_linearity->AddEntry(fit_linearity_cb, fmt::format("Fit: ADC = {:.2f}E + {:.2f}", cb_fit_a, cb_fit_b).c_str(), "l");
    legend_linearity->Draw();


    const double latex_x_start = 0.18;
    const double latex_y_start = 0.85;
    const double latex_y_step  = 0.05;
    TLatex latex_linearity;
    latex_linearity.SetNDC();
    latex_linearity.SetTextSize(0.05);
    latex_linearity.SetTextFont(62);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    latex_linearity.SetTextSize(0.035);
    latex_linearity.SetTextFont(42);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step, annotation_testbeam_title.c_str());
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, config_description.c_str());
    // write date
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    char date_str[100];
    if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y", std::localtime(&now_c))) {
        latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    }

    canvas_linearity->Print((out_pdf + "(").c_str());
    canvas_linearity->Write();
    canvas_linearity->Close();

    // ! === calculate and draw non-linearity ====================
    TCanvas* canvas_nonlinearity = new TCanvas("canvas_nonlinearity", "Non-Linearity vs Beam Energy", 800, 300);
    canvas_nonlinearity->SetLeftMargin(0.15);
    canvas_nonlinearity->SetBottomMargin(0.23);
    TGraphErrors* graph_nonlinearity_gaus = new TGraphErrors();
    TGraphErrors* graph_nonlinearity_cb = new TGraphErrors();
    for (size_t i = 0; i < list_param_gaus_fit_mean.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y_gaus = list_param_gaus_fit_mean[i];
        double y_cb = list_param_cb_fit_mean[i];
        double yerr_gaus = list_param_gaus_fit_mean_err_stat[i];
        double yerr_cb = list_param_cb_fit_mean_err_stat[i];
        double y_gaus_fit = gaus_fit_a * x + gaus_fit_b;
        double y_cb_fit = cb_fit_a * x + cb_fit_b;
        double nonlinearity_gaus = (y_gaus - y_gaus_fit) / y_gaus_fit * 100.0; // in percentage
        double nonlinearity_cb = (y_cb - y_cb_fit) / y_cb_fit * 100.0; // in percentage
        double nonlinearity_gaus_err = (yerr_gaus / y_gaus_fit) * 100.0;
        double nonlinearity_cb_err = (yerr_cb / y_cb_fit) * 100.0;
        graph_nonlinearity_gaus->SetPoint(i, x, nonlinearity_gaus);
        graph_nonlinearity_gaus->SetPointError(i, 0.03 * x, nonlinearity_gaus_err);
        graph_nonlinearity_cb->SetPoint(i, x, nonlinearity_cb);
        graph_nonlinearity_cb->SetPointError(i, 0.03 * x, nonlinearity_cb_err);
    }
    graph_nonlinearity_gaus->SetMarkerStyle(20);
    graph_nonlinearity_gaus->SetMarkerColor(kCyan+2);
    graph_nonlinearity_gaus->SetLineColor(kCyan+2);
    graph_nonlinearity_gaus->SetTitle(";Beam Energy [GeV];Non-Linearity [%]");
    graph_nonlinearity_gaus->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);
    graph_nonlinearity_gaus->GetYaxis()->SetRangeUser(-8.0, 8.0);
    graph_nonlinearity_gaus->GetXaxis()->SetTitleSize(axis_title_size*2.0);
    graph_nonlinearity_gaus->GetYaxis()->SetTitleSize(axis_title_size*2.0);
    graph_nonlinearity_gaus->GetXaxis()->SetLabelSize(axis_label_size*2.0);
    graph_nonlinearity_gaus->GetYaxis()->SetLabelSize(axis_label_size*2.0);
    graph_nonlinearity_gaus->Draw("AP");

    graph_nonlinearity_cb->SetMarkerStyle(21);
    graph_nonlinearity_cb->SetMarkerColor(kMagenta+2);
    graph_nonlinearity_cb->SetLineColor(kMagenta+2);
    graph_nonlinearity_cb->Draw("P SAME");

    // draw horizontal line at y=0
    TLine* line_zero = new TLine(energy_axis_min, 0.0, energy_axis_max, 0.0);
    line_zero->SetLineStyle(2);
    line_zero->SetLineColor(kBlack);
    line_zero->Draw("SAME");

    // legend
    // TLegend* legend_nonlinearity = new TLegend(0.55, 0.65, 0.90, 0.89);
    // legend_nonlinearity->SetBorderSize(0);
    // legend_nonlinearity->SetFillStyle(0);
    // legend_nonlinearity->SetTextSize(legend_text_size*2.0);
    // legend_nonlinearity->AddEntry(graph_nonlinearity_gaus, "Gaussian Fit Mean", "lp");
    // legend_nonlinearity->AddEntry(graph_nonlinearity_cb, "Crystal Ball Fit Mean", "lp");
    // legend_nonlinearity->Draw();

    // latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    // latex_linearity.SetTextSize(0.035);
    // latex_linearity.SetTextFont(42);
    // latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step, annotation_testbeam_title.c_str());
    // latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, config_description.c_str());
    // // write date
    // if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y ", std::localtime(&now_c))) {
    //     latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    // }   

    canvas_nonlinearity->Print(out_pdf.c_str());
    canvas_nonlinearity->Write();
    canvas_nonlinearity->Close();

    // ! === Resolution vs Beam Energy ====================
    TCanvas* canvas_resolution = new TCanvas("canvas_resolution", "Energy Resolution vs Beam Energy", 800, 600);
    canvas_resolution->SetLeftMargin(0.15);
    canvas_resolution->SetBottomMargin(0.10);
    TGraphErrors* graph_resolution_gaus = new TGraphErrors();
    TGraphErrors* graph_resolution_cb = new TGraphErrors();
    TGraphErrors* graph_resolution_cb_eff = new TGraphErrors();
    for (size_t i = 0; i < list_param_gaus_fit_resolution.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y_gaus = list_param_gaus_fit_resolution[i];
        double y_cb = list_param_cb_fit_resolution[i];
        double y_cb_eff = list_param_cb_fit_effective_resolution[i];
        double yerr_gaus = list_param_gaus_fit_resolution_err[i];
        double yerr_cb = list_param_cb_fit_resolution_err[i];
        double yerr_cb_eff = list_param_cb_fit_effective_resolution_err[i];
        graph_resolution_gaus->SetPoint(i, x, y_gaus);
        graph_resolution_gaus->SetPointError(i, 0.03 * x, yerr_gaus);
        graph_resolution_cb->SetPoint(i, x, y_cb);
        graph_resolution_cb->SetPointError(i, 0.03 * x, yerr_cb);
        graph_resolution_cb_eff->SetPoint(i, x, y_cb_eff);
        graph_resolution_cb_eff->SetPointError(i, 0.03 * x, yerr_cb_eff);
    }
    graph_resolution_gaus->SetMarkerStyle(20);
    graph_resolution_gaus->SetMarkerColor(kCyan+2);
    graph_resolution_gaus->SetLineColor(kCyan+2);
    graph_resolution_gaus->SetTitle(";Beam Energy [GeV];Energy Resolution [%]");
    graph_resolution_gaus->GetXaxis()->SetTitleSize(axis_title_size);
    graph_resolution_gaus->GetYaxis()->SetTitleSize(axis_title_size);
    graph_resolution_gaus->GetXaxis()->SetLabelSize(axis_label_size);
    graph_resolution_gaus->GetYaxis()->SetLabelSize(axis_label_size);
    graph_resolution_gaus->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);
    graph_resolution_gaus->GetYaxis()->SetRangeUser(6.0, 40.0);

    // fit with function: sqrt( (a/sqrt(E))^2 + b^2 + (c/E)^2 )
    TF1* fit_resolution_gaus = new TF1("fit_resolution_gaus", resolution_func, energy_axis_min, energy_axis_max, 3);
    // guess c parameter to be the smallest resolution point
    double min_resolution = *std::min_element(list_param_gaus_fit_resolution.begin(), list_param_gaus_fit_resolution.end());
    fit_resolution_gaus->SetParameters(min_resolution, 100.0, 0.0);
    fit_resolution_gaus->SetLineColor(kCyan+2);
    graph_resolution_gaus->Fit(fit_resolution_gaus, "RQ");
    double gaus_res_fit_a = fabs(fit_resolution_gaus->GetParameter(1)); // stochastic term (already in %)
    double gaus_res_fit_b = fabs(fit_resolution_gaus->GetParameter(2)); // noise term (already in %)
    double gaus_res_fit_c = fabs(fit_resolution_gaus->GetParameter(0)); // constant term (already in %)
    graph_resolution_gaus->Draw("AP");
    fit_resolution_gaus->Draw("SAME");

    graph_resolution_cb->SetMarkerStyle(21);
    graph_resolution_cb->SetMarkerColor(kMagenta+2);
    graph_resolution_cb->SetLineColor(kMagenta+2);
    TF1* fit_resolution_cb = new TF1("fit_resolution_cb", resolution_func, energy_axis_min, energy_axis_max, 3);
    min_resolution = *std::min_element(list_param_cb_fit_resolution.begin(), list_param_cb_fit_resolution.end());
    fit_resolution_cb->SetParameters(min_resolution, 100.0, 0.0);
    fit_resolution_cb->SetLineColor(kMagenta+2);
    graph_resolution_cb->Fit(fit_resolution_cb, "RQ");
    double cb_res_fit_a = fabs(fit_resolution_cb->GetParameter(1)); // stochastic term (already in %)
    double cb_res_fit_b = fabs(fit_resolution_cb->GetParameter(2)); // noise term (already in %)
    double cb_res_fit_c = fabs(fit_resolution_cb->GetParameter(0)); // constant term (already in %)
    graph_resolution_cb->Draw("P SAME");
    fit_resolution_cb->Draw("SAME");

    graph_resolution_cb_eff->SetMarkerStyle(22);
    graph_resolution_cb_eff->SetMarkerColor(kGreen+2);
    graph_resolution_cb_eff->SetLineColor(kGreen+2);
    graph_resolution_cb_eff->Draw("P SAME");

    // legend
    TLegend* legend_resolution = new TLegend(0.55, 0.55, 0.90, 0.89);
    legend_resolution->SetBorderSize(0);
    legend_resolution->SetFillStyle(0);
    legend_resolution->SetTextSize(legend_text_size);
    legend_resolution->AddEntry(graph_resolution_gaus, "Gaussian Fit", "lp");
    legend_resolution->AddEntry(
        fit_resolution_gaus,
        fmt::format("#sigma/E = {:.2f}/#sqrt{{E}} #oplus {:.2f}/E #oplus {:.2f}",
                    gaus_res_fit_a, gaus_res_fit_b, gaus_res_fit_c).c_str(), "l");
    legend_resolution->AddEntry(graph_resolution_cb, "Crystal Ball Core", "lp");
    legend_resolution->AddEntry(fit_resolution_cb,
        fmt::format("#sigma/E = {:.2f}/#sqrt{{E}} #oplus {:.2f}/E #oplus {:.2f}",
                    cb_res_fit_a, cb_res_fit_b, cb_res_fit_c).c_str(), "l");
    legend_resolution->AddEntry(graph_resolution_cb_eff, "68% Effective", "lp");
    legend_resolution->Draw();

    latex_linearity.SetTextSize(0.05);
    latex_linearity.SetTextFont(62);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    latex_linearity.SetTextSize(0.035);
    latex_linearity.SetTextFont(42);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step, annotation_testbeam_title.c_str());
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, config_description.c_str());
    // write date
    if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y ", std::localtime(&now_c))) {
        latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    }

    canvas_resolution->Print(out_pdf.c_str());
    canvas_resolution->Write();
    canvas_resolution->Close();



    // write dummy canvas
    TCanvas* canvas_dummy = new TCanvas("canvas_dummy", "Dummy Canvas", 800, 600);
    canvas_dummy->Print((out_pdf + ")").c_str());

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}