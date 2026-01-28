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

    const std::string run_file_folder = "dump/316_ADC_ToT/beamtests/";
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

    // std::vector<double> list_param_cb_fit_mean;
    // std::vector<double> list_param_cb_fit_sigma;
    // std::vector<double> list_param_cb_fit_effective_sigma;
    // std::vector<double> list_param_cb_fit_resolution;
    // std::vector<double> list_param_cb_fit_effective_resolution;
    // std::vector<double> list_param_cb_fit_mean_err_stat;
    // std::vector<double> list_param_cb_fit_mean_err_sys;
    // std::vector<double> list_param_cb_fit_sigma_err_stat;
    // std::vector<double> list_param_cb_fit_sigma_err_sys;
    // std::vector<double> list_param_cb_fit_effective_sigma_err_stat;
    // std::vector<double> list_param_cb_fit_effective_sigma_err_sys;
    // std::vector<double> list_param_cb_fit_resolution_err;
    // std::vector<double> list_param_cb_fit_effective_resolution_err;

    std::vector<double> list_param_tot_fit_mean;
    std::vector<double> list_param_tot_fit_sigma;
    std::vector<double> list_param_tot_fit_effective_sigma;
    std::vector<double> list_param_tot_fit_resolution;
    std::vector<double> list_param_tot_fit_effective_resolution;
    std::vector<double> list_param_tot_fit_mean_err_stat;
    std::vector<double> list_param_tot_fit_mean_err_sys;
    std::vector<double> list_param_tot_fit_sigma_err_stat;
    std::vector<double> list_param_tot_fit_sigma_err_sys;
    std::vector<double> list_param_tot_fit_effective_sigma_err_stat;
    std::vector<double> list_param_tot_fit_effective_sigma_err_sys;
    std::vector<double> list_param_tot_fit_resolution_err;
    std::vector<double> list_param_tot_fit_effective_resolution_err;

    std::vector<TH1D*> adc_tot_histograms;

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
            continue;
        }

        // TH1D* h1_adc_tot_combined = new TH1D("h1_adc_tot_combined", "ADC in ToT Channels Distribution;ADC Sum in ToT Channels;Counts", n_bins, x_hist_min_adc_in_tot, x_hist_max_adc_in_tot);
        auto h1_adc_tot_combined = dynamic_cast<TH1D*>(input_root->Get("h1_adc_tot_combined"));
        if (!h1_adc_tot_combined) {
            spdlog::error("Failed to get histogram 'h1_adc_tot_combined' from file {}", run_file);
            input_root->Close();
            delete input_root;
            continue;
        }
        h1_adc_tot_combined->SetDirectory(nullptr); // Detach from file
        spdlog::info("Successfully read histogram 'h1_adc_tot_combined' from file {}", run_file);
        adc_tot_histograms.push_back( h1_adc_tot_combined );

        TDirectory* dir = input_root;
        if (!dir) {
            spdlog::error("Directory with parameters is missing in file {}", run_file);
            input_root->Close(); delete input_root;
            continue;
        }

        try {
            double v_tot_cb_fit_mean                   = getParam(dir, "tot_cb_fit_mean",                   /*required=*/true);
            double v_tot_cb_fit_sigma                  = getParam(dir, "tot_cb_fit_sigma",                  /*required=*/true);
            double v_tot_cb_fit_effective_sigma        = getParam(dir, "tot_cb_fit_effective_sigma",        /*required=*/false);
            double v_tot_cb_fit_resolution             = getParam(dir, "tot_cb_fit_resolution",             /*required=*/false);
            double v_tot_cb_fit_effective_resolution    = getParam(dir, "tot_cb_fit_effective_resolution",    /*required=*/false);
            double v_tot_cb_fit_mean_err_stat          = getParam(dir, "tot_cb_fit_mean_err_stat",          /*required=*/false);
            double v_tot_cb_fit_mean_err_sys           = getParam(dir, "tot_cb_fit_mean_err_sys",           /*required=*/false);
            double v_tot_cb_fit_sigma_err_stat         = getParam(dir, "tot_cb_fit_sigma_err_stat",         /*required=*/false);
            double v_tot_cb_fit_sigma_err_sys          = getParam(dir, "tot_cb_fit_sigma_err_sys",          /*required=*/false);
            double v_tot_cb_fit_effective_sigma_err_stat = getParam(dir, "tot_cb_fit_effective_sigma_err_stat", /*required=*/false);
            double v_tot_cb_fit_effective_sigma_err_sys  = getParam(dir, "tot_cb_fit_effective_sigma_err_sys",  /*required=*/false);
            double v_tot_cb_fit_resolution_err         = getParam(dir, "tot_cb_fit_resolution_err",         /*required=*/false);
            double v_tot_cb_fit_effective_resolution_err = getParam(dir, "tot_cb_fit_effective_resolution_err", /*required=*/false);

            list_param_tot_fit_mean.push_back(v_tot_cb_fit_mean);
            list_param_tot_fit_sigma.push_back(v_tot_cb_fit_sigma);
            list_param_tot_fit_effective_sigma.push_back(v_tot_cb_fit_effective_sigma);
            list_param_tot_fit_resolution.push_back(v_tot_cb_fit_resolution);
            list_param_tot_fit_effective_resolution.push_back(v_tot_cb_fit_effective_resolution);
            list_param_tot_fit_mean_err_stat.push_back(v_tot_cb_fit_mean_err_stat);
            list_param_tot_fit_mean_err_sys.push_back(v_tot_cb_fit_mean_err_sys);
            list_param_tot_fit_sigma_err_stat.push_back(v_tot_cb_fit_sigma_err_stat);
            list_param_tot_fit_sigma_err_sys.push_back(v_tot_cb_fit_sigma_err_sys);
            list_param_tot_fit_effective_sigma_err_stat.push_back(v_tot_cb_fit_effective_sigma_err_stat);
            list_param_tot_fit_effective_sigma_err_sys.push_back(v_tot_cb_fit_effective_sigma_err_sys);
            list_param_tot_fit_resolution_err.push_back(v_tot_cb_fit_resolution_err);
            list_param_tot_fit_effective_resolution_err.push_back(v_tot_cb_fit_effective_resolution_err);

        } catch (const std::exception& e) {
            spdlog::error("While reading parameters from {}: {}", run_file, e.what());
            input_root->Close(); delete input_root;
            continue;
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
    TGraphErrors* graph_linearity_cb = new TGraphErrors();
    for (size_t i = 0; i < list_param_tot_fit_mean.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y_cb = list_param_tot_fit_mean[i];
        double yerr_cb = list_param_tot_fit_mean_err_stat[i];
        graph_linearity_cb->SetPoint(i, x, y_cb);
        graph_linearity_cb->SetPointError(i, 0.03 * x, yerr_cb);
    }
    graph_linearity_cb->SetMarkerStyle(21);
    graph_linearity_cb->SetMarkerColor(kMagenta+2);
    graph_linearity_cb->SetLineColor(kMagenta+2);
    graph_linearity_cb->SetTitle(";Beam Energy [GeV];Fit Mean [ADC counts]");
    graph_linearity_cb->GetXaxis()->SetTitleSize(axis_title_size);
    graph_linearity_cb->GetYaxis()->SetTitleSize(axis_title_size);
    graph_linearity_cb->GetXaxis()->SetLabelSize(axis_label_size);
    graph_linearity_cb->GetYaxis()->SetLabelSize(axis_label_size);
    graph_linearity_cb->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);

    // linear fit
    TF1* fit_linearity_cb = new TF1("fit_linearity_cb", "pol1", energy_axis_min, energy_axis_max);
    fit_linearity_cb->SetLineColor(kMagenta+2);
    double init_a = (list_param_tot_fit_mean.back() - list_param_tot_fit_mean.front()) /
                    ( (found_run_energies ? static_cast<double>(run_energies.back()) : static_cast<double>(list_param_tot_fit_mean.size()-1)) -
                      (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0) );
    double init_b = list_param_tot_fit_mean.front() - init_a * (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0);
    fit_linearity_cb->SetParameters(init_b, init_a);
    graph_linearity_cb->Fit(fit_linearity_cb, "R");
    double cb_fit_a = fit_linearity_cb->GetParameter(1);;
    double cb_fit_b = fit_linearity_cb->GetParameter(0);
    graph_linearity_cb->Draw("AP SAME");
    fit_linearity_cb->Draw("SAME");
    TLegend* legend_linearity = new TLegend(0.55, 0.15, 0.90, 0.35);
    legend_linearity->SetBorderSize(0);
    legend_linearity->SetFillStyle(0);
    legend_linearity->SetTextSize(legend_text_size);
    legend_linearity->AddEntry(graph_linearity_cb, "Crystal Ball Fit Mean", "lp");
    // add fit function info to legend
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
    std::string out_png = out_pdf.substr(0, out_pdf.find_last_of(".")) + "_pic0.pdf";
    canvas_linearity->SaveAs(out_png.c_str());
    canvas_linearity->Write();
    canvas_linearity->Close();


    // ! === draw ADC sum distributions ====================
    TCanvas* canvas_adc_sum_distribution = new TCanvas("canvas_adc_sum_distribution", "ADC Sum Distributions", 800, 600);
    double adc_sum_90pct_max = 0.0;
    for (const auto& hist_adc_sum : adc_tot_histograms) {
        double current_90pct_max = 0;
        // go through the histogram to find the 90% max bin
        int n_bins = hist_adc_sum->GetNbinsX();
        double total_counts = hist_adc_sum->Integral();
        double cumulative_counts = 0.0;
        for (int bin = 1; bin <= n_bins; ++bin) {
            cumulative_counts += hist_adc_sum->GetBinContent(bin);
            if (cumulative_counts >= 0.9 * total_counts) {
                current_90pct_max = hist_adc_sum->GetBinCenter(bin);
                break;
            }
        }
        if (current_90pct_max > adc_sum_90pct_max) {
            adc_sum_90pct_max = current_90pct_max;
        }
    }

    double adc_tot_sum_y_max = 0.0;
    for (auto& hist_adc_sum : adc_tot_histograms) {
        hist_adc_sum->GetXaxis()->SetRangeUser(0, adc_sum_90pct_max * 1.3);
        // normalize the histogram
        double integral = hist_adc_sum->Integral();
        if (integral > 0) {
            hist_adc_sum->Scale(1.0 / integral);
            double current_y_max = hist_adc_sum->GetMaximum();
            if (current_y_max > adc_tot_sum_y_max) {
                adc_tot_sum_y_max = current_y_max;
            }
        }
    }

    for (size_t i = 0; i < adc_tot_histograms.size(); ++i) {
        TH1D* hist_adc_sum = adc_tot_histograms[i];
        // set y axis range
        hist_adc_sum->GetYaxis()->SetRangeUser(0, adc_tot_sum_y_max * 1.3);
        hist_adc_sum->SetLineColor( color_list[i % color_list.size()] );
        hist_adc_sum->SetLineWidth(2);
        // print some info
        spdlog::info("Drawing ADC sum distribution for run {}", i);
        spdlog::info("  Mean: {:.2f}, RMS: {:.2f}", hist_adc_sum->GetMean(), hist_adc_sum->GetRMS());
        spdlog::info("  Entries: {}", hist_adc_sum->GetEntries());
        if (i == 0) {
            format_1d_hist_canvas(canvas_adc_sum_distribution, hist_adc_sum,
                color_list[i % color_list.size()], annotation_canvas_title, annotation_testbeam_title, "ADC+ToT Sum Distributions");
        } else {
            hist_adc_sum->Draw("HIST SAME");
        }
    }

    // top right legend
    TLegend* legend_adc_sum = new TLegend(0.75, 0.55, 0.95, 0.89);
    legend_adc_sum->SetBorderSize(0);
    legend_adc_sum->SetFillStyle(0);
    legend_adc_sum->SetTextSize(legend_text_size);
    // for (size_t i = 0; i < hist_input_adc_sums.size(); ++i) {
    //     const auto& hist_adc_sum = hist_input_adc_sums[i];
    //     std::string label;
    //     if (found_run_energies) {
    //         label = fmt::format("{} GeV", run_energies[i]);
    //     } else {
    //         label = run_labels[i];
    //     }
    //     legend_adc_sum->AddEntry(hist_adc_sum, label.c_str(), "l");
    // }
    // legend_adc_sum->Draw("SAME");   
    for (size_t i = 0; i < adc_tot_histograms.size(); ++i) {
        const auto& hist_adc_sum = adc_tot_histograms[i];
        std::string label;
        if (found_run_energies) {
            label = fmt::format("{} GeV", run_energies[i]);
        } else {
            label = run_labels[i];
        }
        legend_adc_sum->AddEntry(hist_adc_sum, label.c_str(), "l");
    }
    legend_adc_sum->Draw("SAME");


    canvas_adc_sum_distribution->Modified();
    canvas_adc_sum_distribution->Update();
    canvas_adc_sum_distribution->Write();
    canvas_adc_sum_distribution->Print(out_pdf.c_str());
    // save as a png file
    out_png = out_pdf.substr(0, out_pdf.find_last_of(".")) + "_pic1.pdf";
    canvas_adc_sum_distribution->SaveAs(out_png.c_str());
    canvas_adc_sum_distribution->Close();

    // ! === calculate and draw non-linearity ====================
    TCanvas* canvas_nonlinearity = new TCanvas("canvas_nonlinearity", "Non-Linearity vs Beam Energy", 800, 300);
    canvas_nonlinearity->SetLeftMargin(0.15);
    canvas_nonlinearity->SetBottomMargin(0.23);
    // TGraphErrors* graph_nonlinearity_gaus = new TGraphErrors();
    TGraphErrors* graph_nonlinearity_cb = new TGraphErrors();
    // for (size_t i = 0; i < list_param_gaus_fit_mean.size(); ++i) {
    //     double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
    //     double y_cb_fit = cb_fit_a * x + cb_fit_b;
    //     double nonlinearity_gaus = (y_gaus - y_gaus_fit) / y_gaus_fit * 100.0; // in percentage
    //     double nonlinearity_cb = (y_cb - y_cb_fit) / y_cb_fit * 100.0; // in percentage
    //     double nonlinearity_gaus_err = (yerr_gaus / y_gaus_fit) * 100.0;
    //     double nonlinearity_cb_err = (yerr_cb / y_cb_fit) * 100.0;
    //     graph_nonlinearity_cb->SetPoint(i, x, nonlinearity_cb);
    //     graph_nonlinearity_cb->SetPointError(i, 0.03 * x, nonlinearity_cb_err);
    // }

    // graph_nonlinearity_cb->SetMarkerStyle(21);
    // graph_nonlinearity_cb->SetMarkerColor(kMagenta+2);
    // graph_nonlinearity_cb->SetLineColor(kMagenta+2);
    // graph_nonlinearity_cb->Draw("P SAME");
    for (size_t i = 0; i < list_param_tot_fit_mean.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y_cb_fit = cb_fit_a * x + cb_fit_b;
        double nonlinearity_cb = (list_param_tot_fit_mean[i] - y_cb_fit) / y_cb_fit * 100.0; // in percentage
        double nonlinearity_cb_err = (list_param_tot_fit_mean_err_stat[i] / y_cb_fit) * 100.0;
        graph_nonlinearity_cb->SetPoint(i, x, nonlinearity_cb);
        graph_nonlinearity_cb->SetPointError(i, 0.03 * x, nonlinearity_cb_err);
    }

    graph_nonlinearity_cb->SetMarkerStyle(21);
    graph_nonlinearity_cb->SetMarkerColor(kMagenta+2);
    graph_nonlinearity_cb->SetLineColor(kMagenta+2);
    graph_nonlinearity_cb->SetTitle(";Beam Energy [GeV];Non-Linearity [%]");
    double title_size_ratio = 600.0 / 300.0; // canvas height ratio
    graph_nonlinearity_cb->GetXaxis()->SetTitleSize(axis_title_size * title_size_ratio);
    graph_nonlinearity_cb->GetYaxis()->SetTitleSize(axis_title_size * title_size_ratio);
    graph_nonlinearity_cb->GetXaxis()->SetLabelSize(axis_label_size * title_size_ratio);
    graph_nonlinearity_cb->GetYaxis()->SetLabelSize(axis_label_size * title_size_ratio);
    graph_nonlinearity_cb->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);
    graph_nonlinearity_cb->GetYaxis()->SetRangeUser(-8.0, 8.0);
    graph_nonlinearity_cb->Draw("AP");

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
    out_png = out_pdf.substr(0, out_pdf.find_last_of(".")) + "_pic2.pdf";
    canvas_nonlinearity->SaveAs(out_png.c_str());
    canvas_nonlinearity->Write();
    canvas_nonlinearity->Close();

    // ! === Resolution vs Beam Energy ====================
    TCanvas* canvas_resolution = new TCanvas("canvas_resolution", "Energy Resolution vs Beam Energy", 800, 600);
    canvas_resolution->SetLeftMargin(0.15);
    canvas_resolution->SetBottomMargin(0.10);
    // TGraphErrors* graph_resolution_gaus = new TGraphErrors();
    TGraphErrors* graph_resolution_cb = new TGraphErrors();
    TGraphErrors* graph_resolution_cb_eff = new TGraphErrors();
    // for (size_t i = 0; i < list_param_gaus_fit_resolution.size(); ++i) {
    //     double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
    //     double y_gaus = list_param_gaus_fit_resolution[i];
    //     double y_cb = list_param_cb_fit_resolution[i];
    //     double y_cb_eff = list_param_cb_fit_effective_resolution[i];
    //     double yerr_gaus = list_param_gaus_fit_resolution_err[i];
    //     double yerr_cb = list_param_cb_fit_resolution_err[i];
    //     double yerr_cb_eff = list_param_cb_fit_effective_resolution_err[i];
    //     graph_resolution_gaus->SetPoint(i, x, y_gaus);
    //     graph_resolution_gaus->SetPointError(i, 0.03 * x, yerr_gaus);
    //     graph_resolution_cb->SetPoint(i, x, y_cb);
    //     graph_resolution_cb->SetPointError(i, 0.03 * x, yerr_cb);
    //     graph_resolution_cb_eff->SetPoint(i, x, y_cb_eff);
    //     graph_resolution_cb_eff->SetPointError(i, 0.03 * x, yerr_cb_eff);
    // }
    // graph_resolution_gaus->SetMarkerStyle(20);
    // graph_resolution_gaus->SetMarkerColor(kCyan+2);
    // graph_resolution_gaus->SetLineColor(kCyan+2);
    // graph_resolution_gaus->SetTitle(";Beam Energy [GeV];Energy Resolution [%]");
    // graph_resolution_gaus->GetXaxis()->SetTitleSize(axis_title_size);
    // graph_resolution_gaus->GetYaxis()->SetTitleSize(axis_title_size);
    // graph_resolution_gaus->GetXaxis()->SetLabelSize(axis_label_size);
    // graph_resolution_gaus->GetYaxis()->SetLabelSize(axis_label_size);
    // graph_resolution_gaus->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);
    // graph_resolution_gaus->GetYaxis()->SetRangeUser(6.0, 40.0);

    for (size_t i = 0; i < list_param_tot_fit_resolution.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y_cb = list_param_tot_fit_resolution[i];
        double y_cb_eff = list_param_tot_fit_effective_resolution[i];
        double yerr_cb = list_param_tot_fit_resolution_err[i];
        double yerr_cb_eff = list_param_tot_fit_effective_resolution_err[i];
        graph_resolution_cb->SetPoint(i, x, y_cb);
        graph_resolution_cb->SetPointError(i, 0.03 * x, yerr_cb);
        graph_resolution_cb_eff->SetPoint(i, x, y_cb_eff);
        graph_resolution_cb_eff->SetPointError(i, 0.03 * x, yerr_cb_eff);
    }
    graph_resolution_cb->SetMarkerStyle(21);
    graph_resolution_cb->SetMarkerColor(kMagenta+2);
    graph_resolution_cb->SetLineColor(kMagenta+2);
    graph_resolution_cb->SetTitle(";Beam Energy [GeV];Energy Resolution [%]");
    graph_resolution_cb->GetXaxis()->SetTitleSize(axis_title_size);
    graph_resolution_cb->GetYaxis()->SetTitleSize(axis_title_size);
    graph_resolution_cb->GetXaxis()->SetLabelSize(axis_label_size); 
    graph_resolution_cb->GetYaxis()->SetLabelSize(axis_label_size);
    graph_resolution_cb->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);
    graph_resolution_cb->GetYaxis()->SetRangeUser(6.0, 40.0);

    // fit with function: sqrt( (a/sqrt(E))^2 + b^2 + (c/E)^2 )
    // TF1* fit_resolution_gaus = new TF1("fit_resolution_gaus", resolution_func, energy_axis_min, energy_axis_max, 3);
    // // guess c parameter to be the smallest resolution point
    // double min_resolution = *std::min_element(list_param_gaus_fit_resolution.begin(), list_param_gaus_fit_resolution.end());
    // fit_resolution_gaus->SetParameters(min_resolution, 100.0, 0.0);
    // fit_resolution_gaus->SetLineColor(kCyan+2);
    // graph_resolution_gaus->Fit(fit_resolution_gaus, "RQ");
    // double gaus_res_fit_a = fabs(fit_resolution_gaus->GetParameter(1)); // stochastic term (already in %)
    // double gaus_res_fit_b = fabs(fit_resolution_gaus->GetParameter(2)); // noise term (already in %)
    // double gaus_res_fit_c = fabs(fit_resolution_gaus->GetParameter(0)); // constant term (already in %)
    // graph_resolution_gaus->Draw("AP");
    // fit_resolution_gaus->Draw("SAME");

    // graph_resolution_cb->SetMarkerStyle(21);
    // graph_resolution_cb->SetMarkerColor(kMagenta+2);
    // graph_resolution_cb->SetLineColor(kMagenta+2);
    // TF1* fit_resolution_cb = new TF1("fit_resolution_cb", resolution_func, energy_axis_min, energy_axis_max, 3);
    // min_resolution = *std::min_element(list_param_cb_fit_resolution.begin(), list_param_cb_fit_resolution.end());
    // fit_resolution_cb->SetParameters(min_resolution, 100.0, 0.0);
    // fit_resolution_cb->SetLineColor(kMagenta+2);
    // graph_resolution_cb->Fit(fit_resolution_cb, "RQ");
    // double cb_res_fit_a = fabs(fit_resolution_cb->GetParameter(1)); // stochastic term (already in %)
    // double cb_res_fit_b = fabs(fit_resolution_cb->GetParameter(2)); // noise term (already in %)
    // double cb_res_fit_c = fabs(fit_resolution_cb->GetParameter(0)); // constant term (already in %)
    // graph_resolution_cb->Draw("P SAME");
    // fit_resolution_cb->Draw("SAME");

    // graph_resolution_cb_eff->SetMarkerStyle(22);
    // graph_resolution_cb_eff->SetMarkerColor(kGreen+2);
    // graph_resolution_cb_eff->SetLineColor(kGreen+2);
    // graph_resolution_cb_eff->Draw("P SAME");

    TF1* fit_resolution_cb = new TF1("fit_resolution_cb", resolution_func, energy_axis_min, energy_axis_max, 3);
    double min_resolution = *std::min_element(list_param_tot_fit_resolution.begin(), list_param_tot_fit_resolution.end());
    fit_resolution_cb->SetParameters(min_resolution, 100.0, 0.0);
    fit_resolution_cb->SetLineColor(kMagenta+2);
    graph_resolution_cb->Fit(fit_resolution_cb, "RQ");
    double cb_res_fit_a = fabs(fit_resolution_cb->GetParameter(1)); // stochastic term (already in %)
    double cb_res_fit_b = fabs(fit_resolution_cb->GetParameter(2)); // noise term (already in %)
    double cb_res_fit_c = fabs(fit_resolution_cb->GetParameter(0)); // constant term (already in %)
    graph_resolution_cb->Draw("AP");
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
    out_png = out_pdf.substr(0, out_pdf.find_last_of(".")) + "_pic3.pdf";
    canvas_resolution->SaveAs(out_png.c_str());
    canvas_resolution->Write();
    canvas_resolution->Close();

    // draw resolution compare with previous testbeam results if available
    TCanvas* canvas_resolution_to_previous = new TCanvas("canvas_resolution_to_previous", "Energy Resolution Comparison", 800, 600);
    canvas_resolution_to_previous->SetLeftMargin(0.15);
    canvas_resolution_to_previous->SetBottomMargin(0.10);
    TGraphErrors* graph_resolution_previous = new TGraphErrors();
    TGraphErrors* graph_resolution_MC = new TGraphErrors();
    std::vector<double> mc_energies = {60, 80, 100, 150, 200, 250, 300, 350};
    std::vector<double> mc_resolutions = {8.75, 7.1, 6.5, 5.9, 5.2, 4.3, 4.9, 5.0}; // in %
    std::vector<double> mc_resolution_errors = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<double> previous_energies = {60, 80, 100, 150, 200, 250, 300};
    std::vector<double> previous_resolutions = {0.202903,0.193606,0.180118,0.157828,0.143362,0.135984,0.128683,0.124429};
    std::vector<double> previous_resolution_errors = {0.0223685347754385,
                0.0276764080400618,
                0.0256984028686609,
                0.0235984254347615,
                0.016989855708628,
                0.0151618517338747,
                0.0109004483394033,
                0.0102402795860269};
    for (size_t i = 0; i < mc_energies.size(); ++i) {
        graph_resolution_MC->SetPoint(i, mc_energies[i], mc_resolutions[i]);
        graph_resolution_MC->SetPointError(i, mc_energies[i]*0.03, mc_resolution_errors[i]);
    }
    for (size_t i = 0; i < previous_energies.size(); ++i) {
        graph_resolution_previous->SetPoint(i, previous_energies[i], previous_resolutions[i]*100.0); // convert to %
        graph_resolution_previous->SetPointError(i, previous_energies[i]*0.03, previous_resolution_errors[i]*100.0);
    }
    graph_resolution_previous->SetMarkerStyle(25);
    graph_resolution_previous->SetMarkerColor(kRed);
    graph_resolution_previous->SetLineColor(kRed);
    graph_resolution_previous->SetTitle(";Beam Energy [GeV];Energy Resolution [%]");
    graph_resolution_previous->GetXaxis()->SetTitleSize(axis_title_size);
    graph_resolution_previous->GetYaxis()->SetTitleSize(axis_title_size);
    graph_resolution_previous->GetXaxis()->SetLabelSize(axis_label_size);
    graph_resolution_previous->GetYaxis()->SetLabelSize(axis_label_size);
    graph_resolution_previous->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);
    graph_resolution_previous->GetYaxis()->SetRangeUser(0.0, 40.0);
    graph_resolution_previous->Draw("AP");
    // fit with function: sqrt( (a/sqrt(E))^2 + b^2 + (c/E)^2 )
    TF1* fit_resolution_previous = new TF1("fit_resolution_previous", resolution_func, energy_axis_min, energy_axis_max, 3);
    double min_resolution_previous = *std::min_element(previous_resolutions.begin(), previous_resolutions.end());
    fit_resolution_previous->SetParameters(min_resolution_previous, 100.0, 0.0);
    fit_resolution_previous->SetLineColor(kRed);
    graph_resolution_previous->Fit(fit_resolution_previous, "RQ");
    double previous_res_fit_a = fabs(fit_resolution_previous->GetParameter(1)); // stochastic term (already in %)
    double previous_res_fit_b = fabs(fit_resolution_previous->GetParameter(2)); // noise term (already in %)
    double previous_res_fit_c = fabs(fit_resolution_previous->GetParameter(0)); // constant term (already in %)
    fit_resolution_previous->Draw("SAME");

    // then draw only the crystal ball core from this analysis
    graph_resolution_cb->Draw("P SAME");

    graph_resolution_MC->SetMarkerStyle(24);
    graph_resolution_MC->SetMarkerColor(kBlue);
    graph_resolution_MC->SetLineColor(kBlue);
    graph_resolution_MC->Draw("P SAME");
    // fit with function: sqrt( (a/sqrt(E))^2 + b^2 + (c/E)^2 )
    TF1* fit_resolution_MC = new TF1("fit_resolution_MC", resolution_func, energy_axis_min, energy_axis_max, 3);
    double min_resolution_MC = *std::min_element(mc_resolutions.begin(), mc_resolutions.end());
    fit_resolution_MC->SetParameters(min_resolution_MC, 100.0, 0.0);
    fit_resolution_MC->SetLineColor(kBlue);
    graph_resolution_MC->Fit(fit_resolution_MC, "RQ");
    double mc_res_fit_a = fabs(fit_resolution_MC->GetParameter(1)); // stochastic term (already in %)
    double mc_res_fit_b = fabs(fit_resolution_MC->GetParameter(2)); // noise term (already in %)
    double mc_res_fit_c = fabs(fit_resolution_MC->GetParameter(0)); // constant term (already in %)
    fit_resolution_MC->Draw("SAME");

    // legend
    TLegend* legend_resolution_compare = new TLegend(0.55, 0.55, 0.90, 0.89);
    legend_resolution_compare->SetBorderSize(0);
    legend_resolution_compare->SetFillStyle(0);
    legend_resolution_compare->SetTextSize(legend_text_size);
    legend_resolution_compare->AddEntry(graph_resolution_previous, "FoCal-H Prototype 2", "lp");
    legend_resolution_compare->AddEntry(graph_resolution_previous, fmt::format("#sigma/E = {:.2f}/#sqrt{{E}} #oplus {:.2f}/E #oplus {:.2f}", previous_res_fit_a, previous_res_fit_b, previous_res_fit_c).c_str(), "l");
    legend_resolution_compare->AddEntry(graph_resolution_cb, "FoCal-H Prototype 3", "lp");
    legend_resolution_compare->AddEntry(graph_resolution_cb, fmt::format("#sigma/E = {:.2f}/#sqrt{{E}} #oplus {:.2f}/E #oplus {:.2f}", cb_res_fit_a, cb_res_fit_b, cb_res_fit_c).c_str(), "l");
    legend_resolution_compare->AddEntry(graph_resolution_MC, "MC for Prototype 3", "lp");
    legend_resolution_compare->AddEntry(graph_resolution_MC, fmt::format("#sigma/E = {:.2f}/#sqrt{{E}} #oplus {:.2f}/E #oplus {:.2f}", mc_res_fit_a, mc_res_fit_b, mc_res_fit_c).c_str(), "l");
    legend_resolution_compare->Draw();

    latex_linearity.SetTextSize(0.05);
    latex_linearity.SetTextFont(62);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    latex_linearity.SetTextSize(0.035);
    latex_linearity.SetTextFont(42);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step, annotation_testbeam_title.c_str());
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, "Comparison with Previous Results");
    // write date
    if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y ", std::localtime(&now_c))) {
        latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    }
    canvas_resolution_to_previous->Print(out_pdf.c_str());
    out_png = out_pdf.substr(0, out_pdf.find_last_of(".")) + "_pic4.pdf";
    canvas_resolution_to_previous->SaveAs(out_png.c_str());
    canvas_resolution_to_previous->Write();
    canvas_resolution_to_previous->Close();



    // write dummy canvas
    TCanvas* canvas_dummy = new TCanvas("canvas_dummy", "Dummy Canvas", 800, 600);
    canvas_dummy->Print((out_pdf + ")").c_str());

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}