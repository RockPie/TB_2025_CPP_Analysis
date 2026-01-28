#include <cmath>
#include <fstream>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"

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

    const std::string run_file_folder = "dump/311_RawToT/beamtests/";
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

    std::vector<TH1D*> hist_input_tot_sums;

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

    std::vector<double> tot_adc_slope_factors;
    // from 0.5 to 3.0, step 0.25
    // for (double slope = 0.5; slope <= 8.0; slope += 0.05) {
    for (double slope = 1.8; slope <= 2.4; slope += 0.1) {
        tot_adc_slope_factors.push_back(slope);
    }
    std::vector<double> tot_adc_offset_factors = {-200.0, -100.0, 0.0, 100.0, 200.0};

    int number_of_combine_factors = tot_adc_slope_factors.size() * tot_adc_offset_factors.size();
    std::vector<std::vector<std::vector<double>>> combined_value_matrix(run_files.size(), std::vector<std::vector<double>>(number_of_combine_factors));
    std::vector<std::vector<double>> raw_value_matrix(run_files.size());

    std::vector<std::vector<double>> example_value_matrix(run_files.size());

    double example_tot_adc_slope = 4.0;
    // double example_tot_adc_offset = 900.0 - example_tot_adc_slope * 512.0;
    double example_tot_adc_offset = -1000.0;
    spdlog::info("Example TOT ADC Slope: {}, Offset: {}", example_tot_adc_slope, example_tot_adc_offset);

    for (int _file_index = 0; _file_index < static_cast<int>(run_files.size()); ++_file_index) {
        const auto& run_file = run_files[_file_index];
        spdlog::info("Processing run file: {}", run_file);

        TFile* input_root = TFile::Open(run_file.c_str(), "READ");
        if (!input_root || input_root->IsZombie()) {
            spdlog::error("Failed to open input file {}", run_file);
            if (input_root) { input_root->Close(); delete input_root; }
            continue;
        }

        // TCanvas tot_sum_distribution_canvas("tot_sum_distribution_canvas", "ToT Sum Distribution Canvas", 800, 600);
        // find the TCanvas
        auto tot_sum_distribution_canvas_ptr = dynamic_cast<TCanvas*>(input_root->Get("tot_sum_distribution_canvas"));
        if (tot_sum_distribution_canvas_ptr) {
            TH1D* h1_tot_sum_distribution = dynamic_cast<TH1D*>(tot_sum_distribution_canvas_ptr->GetPrimitive("h1_tot_sum_distribution"));
            if (h1_tot_sum_distribution) {
                h1_tot_sum_distribution->SetDirectory(nullptr); // Detach from file
                hist_input_tot_sums.push_back( h1_tot_sum_distribution );
                spdlog::info("Successfully read histogram 'h1_tot_sum_distribution' from canvas in file {}", run_file);
            } else {
                spdlog::error("Failed to get histogram 'h1_tot_sum_distribution' from canvas in file {}", run_file);
                continue;
            }
        }

        TTree* input_tree = input_root->Get<TTree>("tot_tree");
        if (!input_tree) {
            spdlog::error("Input tree 'tot_Tree' not found in file {}", run_file);
            input_root->Close(); delete input_root;
            continue; // 或者 return 1;
        }
        auto tree_entry_count = input_tree->GetEntries();
        spdlog::info("Input tree 'tot_Tree' has {} entries", tree_entry_count);
        auto& run_combined_values = combined_value_matrix[_file_index];
        auto& run_raw_values = raw_value_matrix[_file_index];
        auto& run_example_values = example_value_matrix[_file_index];
        for (int slope_idx = 0; slope_idx < static_cast<int>(tot_adc_slope_factors.size()); ++slope_idx) {
            double slope = tot_adc_slope_factors[slope_idx];
            for (int offset_idx = 0; offset_idx < static_cast<int>(tot_adc_offset_factors.size()); ++offset_idx) {
                double offset = tot_adc_offset_factors[offset_idx];
                auto & combined_values = run_combined_values[slope_idx * tot_adc_offset_factors.size() + offset_idx];
                combined_values.reserve(tree_entry_count);
            }
        }
        auto *out_branch_adc_sum = new double;
        auto *out_branch_tot_sum = new double;
        auto *out_branch_adc_in_tot_channels = new double;
        input_tree->SetBranchAddress("adc_sum", out_branch_adc_sum);
        input_tree->SetBranchAddress("tot_sum", out_branch_tot_sum);
        input_tree->SetBranchAddress("adc_in_tot_channels", out_branch_adc_in_tot_channels);
        for (Long64_t entry = 0; entry < tree_entry_count; ++entry) {
            input_tree->GetEntry(entry);
            double adc_sum = *out_branch_adc_sum;
            double tot_sum = *out_branch_tot_sum;
            double adc_in_tot_channels = *out_branch_adc_in_tot_channels;
            for (int slope_idx = 0; slope_idx < static_cast<int>(tot_adc_slope_factors.size()); ++slope_idx) {
                double slope = tot_adc_slope_factors[slope_idx];
                for (int offset_idx = 0; offset_idx < static_cast<int>(tot_adc_offset_factors.size()); ++offset_idx) {
                    // double offset = tot_adc_offset_factors[offset_idx];
                    // TODO: remove this temporary fix
                    // double offset = 900.0 - 512.0 * slope;
                    double offset = 0;
                    double combined_value = adc_sum - adc_in_tot_channels + (tot_sum * slope + offset);
                    auto & combined_values = run_combined_values[slope_idx * tot_adc_offset_factors.size() + offset_idx];
                    combined_values.push_back(combined_value);
                    if (slope_idx == 0 && offset_idx == 0) {
                        run_raw_values.push_back(adc_sum);
                        run_example_values.push_back(
                            adc_sum - adc_in_tot_channels + (tot_sum * example_tot_adc_slope + example_tot_adc_offset)
                        );
                    }
                }
            }
        }

        TDirectory* dir = input_root; 
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

    // ! === draw raw ToT sum comparison ====================
    TCanvas* canvas_tot_sum_comparison = new TCanvas("canvas_tot_sum_comparison", "Raw ToT Sum Comparison", 800, 600);
    double tot_sum_90pct_max = 0.0;
    for (const auto& hist : hist_input_tot_sums) {
        double current_90pct_max = 0;
        int n_bins = hist->GetNbinsX();
        double total_entries = hist->GetEntries();
        double cumulative_entries = 0.0;
        for (int bin = 1; bin <= n_bins; ++bin) {
            cumulative_entries += hist->GetBinContent(bin);
            if (cumulative_entries >= 0.9 * total_entries) {
                current_90pct_max = hist->GetBinCenter(bin);
                break;
            }
        }
        if (current_90pct_max > tot_sum_90pct_max) {
            tot_sum_90pct_max = current_90pct_max;
        }
    }
    // reset the binning
    double tot_sum_y_max = 0.0;
    for (auto & hist : hist_input_tot_sums) {
        hist->GetXaxis()->SetRangeUser(0.0, tot_sum_90pct_max * 1.3);
        // normalize
        hist->Scale(1.0 / hist->Integral());
        if (hist->GetMaximum() > tot_sum_y_max) {
            tot_sum_y_max = hist->GetMaximum();
        }
    }
    // set y axis max
    tot_sum_y_max *= 1.2;
    for (auto & hist : hist_input_tot_sums) {
        hist->GetYaxis()->SetRangeUser(0.0, tot_sum_y_max);
    }
    // draw all histograms
    for (size_t i = 0; i < hist_input_tot_sums.size(); ++i) {
        auto & hist = hist_input_tot_sums[i];
        hist->SetLineColor(color_list[i % color_list.size()]);
        hist->GetYaxis()->SetTitle("Normalized Entries");
        if (i == 0) {
            format_1d_hist_canvas(canvas_tot_sum_comparison, hist,
                                color_list[0],annotation_canvas_title,
                                annotation_testbeam_title, "ToT Sum Distributions");
        } else {
            hist->Draw("HIST SAME");
        }
    }
    // legend at top right
    TLegend* legend_tot_sum = new TLegend(0.75, 0.55, 0.90, 0.89);
    legend_tot_sum->SetBorderSize(0);
    legend_tot_sum->SetFillStyle(0);
    legend_tot_sum->SetTextSize(legend_text_size);
    for (size_t i = 0; i < hist_input_tot_sums.size(); ++i) {
        std::string label;
        if (found_run_energies) {
            label = fmt::format("{} GeV", run_energies[i]);
        } else if (i < run_labels.size()) {
            label = run_labels[i];
        } else {
            label = fmt::format("Run {}", run_numbers[i]);
        }
        legend_tot_sum->AddEntry(hist_input_tot_sums[i], label.c_str(), "l");
    }
    legend_tot_sum->Draw("SAME");

    canvas_tot_sum_comparison->Modified();
    canvas_tot_sum_comparison->Update();
    canvas_tot_sum_comparison->Print(out_pdf.c_str());
    canvas_tot_sum_comparison->Write();
    canvas_tot_sum_comparison->Close();

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

    // ! === for each run in the matrix, draw in one canvas ===
    std::vector<double> fit_range_sigmas = {2.0, 2.25, 2.5, 2.75, 3.0}; // in sigma
    std::vector<double> fit_range_offsets = {-0.2, 0.0, 0.2, 0.4}; // in sigma

    std::vector<double> raw_cb_mean_avgs;
    std::vector<double> raw_cb_mean_errs;
    std::vector<double> raw_cb_resolution_avgs;
    std::vector<double> raw_cb_resolution_errs;
    std::vector<double> example_cb_mean_avgs;
    std::vector<double> example_cb_mean_errs;
    std::vector<double> example_cb_resolution_avgs;
    std::vector<double> example_cb_resolution_errs;
    std::vector<std::vector<double>> combined_cb_mean_avgs(tot_adc_slope_factors.size() * tot_adc_offset_factors.size());
    std::vector<std::vector<double>> combined_cb_mean_errs(tot_adc_slope_factors.size() * tot_adc_offset_factors.size());

    double max_fit_mean = 0.0;
    for (size_t run_idx = 0; run_idx < run_files.size(); ++run_idx) {
        const auto& run_file = run_files[run_idx];
        auto& run_raw_values = raw_value_matrix[run_idx];
        auto& run_example_values = example_value_matrix[run_idx];
        auto& run_combined_values = combined_value_matrix[run_idx];
        TCanvas* canvas_combined = new TCanvas(
            fmt::format("canvas_combined_run{}", run_numbers[run_idx]).c_str(),
            fmt::format("Combined ADC Values - Run {}", run_numbers[run_idx]).c_str(),
            1200, 800);

        std::vector<double> sorted_raw_values = run_raw_values;
        std::sort(sorted_raw_values.begin(), sorted_raw_values.end());
        size_t index_90_raw = static_cast<size_t>(0.9 * sorted_raw_values.size());
        if (index_90_raw >= sorted_raw_values.size()) {
            index_90_raw = sorted_raw_values.size() - 1;
        }
        double max_90_percent_x = sorted_raw_values[index_90_raw];
        for (int slope_idx = 0; slope_idx < static_cast<int>(tot_adc_slope_factors.size()); ++slope_idx) {
            double slope = tot_adc_slope_factors[slope_idx];
            for (int offset_idx = 0; offset_idx < static_cast<int>(tot_adc_offset_factors.size()); ++offset_idx) {
                double offset = tot_adc_offset_factors[offset_idx];
                auto & combined_values = run_combined_values[slope_idx * tot_adc_offset_factors.size() + offset_idx];

                // find 90th percentile
                std::vector<double> sorted_values = combined_values;
                std::sort(sorted_values.begin(), sorted_values.end());
                size_t index_90 = static_cast<size_t>(0.9 * sorted_values.size());
                if (index_90 >= sorted_values.size()) {
                    index_90 = sorted_values.size() - 1;
                }
                double value_90 = sorted_values[index_90];
                if (value_90 > max_90_percent_x) {
                    max_90_percent_x = value_90;
                }
            }
        }
        TLegend* legend_combined = new TLegend(0.50, 0.25, 0.89, 0.89);
        legend_combined->SetBorderSize(0);
        legend_combined->SetFillStyle(0);

        // draw raw ADC sum histogram first
        TH1D* hist_raw = new TH1D(
            fmt::format("hist_raw_run{}", run_numbers[run_idx]).c_str(),
            fmt::format("Raw ADC Sum - Run {}", run_numbers[run_idx]).c_str(),
            256, 0.0, max_90_percent_x * 1.6);
        // no stat box or title
        hist_raw->SetStats(kFALSE);
        hist_raw->SetTitle(""); 
        for (const auto& val : run_raw_values) {
            hist_raw->Fill(val);
        }
        hist_raw->GetXaxis()->SetTitle("ADC Sum");
        hist_raw->GetYaxis()->SetTitle("Counts");
        hist_raw->SetLineColor(kBlack);
        hist_raw->Draw();
        legend_combined->AddEntry(hist_raw, "Raw ADC Sum", "l");

        double raw_cb_mean_avg=0.0;
        double raw_cb_mean_err_sys=0.0;
        double raw_cb_mean_err_stat=0.0;
        double raw_cb_sigma_avg=0.0;
        double raw_cb_sigma_err_sys=0.0;
        double raw_cb_sigma_err_stat=0.0;
        double raw_cb_effective_sigma_avg=0.0;
        double raw_cb_effective_sigma_err_sys=0.0;
        double raw_cb_effective_sigma_err_stat=0.0;
        double raw_cb_resolution_avg=0.0;
        double raw_cb_resolution_err=0.0;

        crystalball_fit_th1d(*canvas_combined, *hist_raw, fit_range_sigmas, fit_range_offsets, kBlack, raw_cb_mean_avg, raw_cb_mean_err_sys, raw_cb_mean_err_stat, raw_cb_sigma_avg, raw_cb_sigma_err_sys, raw_cb_sigma_err_stat, raw_cb_resolution_avg, raw_cb_resolution_err);

        raw_cb_mean_avgs.push_back(raw_cb_mean_avg);
        raw_cb_mean_errs.push_back(std::sqrt(raw_cb_mean_err_sys * raw_cb_mean_err_sys + raw_cb_mean_err_stat * raw_cb_mean_err_stat));
        raw_cb_resolution_avgs.push_back(raw_cb_resolution_avg);
        raw_cb_resolution_errs.push_back(raw_cb_resolution_err);
        if (raw_cb_mean_avg > max_fit_mean) {
            max_fit_mean = raw_cb_mean_avg;
        }

        // draw example ADC sum histogram next
        TH1D* hist_example = new TH1D(
            fmt::format("hist_example_run{}", run_numbers[run_idx]).c_str(),
            fmt::format("Example ADC Values - Run {}", run_numbers[run_idx]).c_str(),
            256, 0.0, max_90_percent_x * 1.6);
        hist_example->SetStats(kFALSE);
        hist_example->SetTitle("");
        for (const auto& val : run_example_values) {
            hist_example->Fill(val);
        }
        hist_example->SetLineColor(kGray+2);
        hist_example->Draw("SAME");

        double example_cb_mean_avg=0.0;
        double example_cb_mean_err_sys=0.0;
        double example_cb_mean_err_stat=0.0;
        double example_cb_sigma_avg=0.0;
        double example_cb_sigma_err_sys=0.0;
        double example_cb_sigma_err_stat=0.0;
        double example_cb_effective_sigma_avg=0.0;
        double example_cb_effective_sigma_err_sys=0.0;
        double example_cb_effective_sigma_err_stat=0.0;
        double example_cb_resolution_avg=0.0;
        double example_cb_resolution_err=0.0;

        crystalball_fit_th1d(*canvas_combined, *hist_example, fit_range_sigmas, fit_range_offsets, kGray+2, example_cb_mean_avg, example_cb_mean_err_sys, example_cb_mean_err_stat, example_cb_sigma_avg, example_cb_sigma_err_sys, example_cb_sigma_err_stat, example_cb_resolution_avg, example_cb_resolution_err);

        example_cb_mean_avgs.push_back(example_cb_mean_avg);
        example_cb_mean_errs.push_back(std::sqrt(example_cb_mean_err_sys * example_cb_mean_err_sys + example_cb_mean_err_stat * example_cb_mean_err_stat));
        example_cb_resolution_avgs.push_back(example_cb_resolution_avg);
        example_cb_resolution_errs.push_back(example_cb_resolution_err);
        if (example_cb_mean_avg > max_fit_mean) {
            max_fit_mean = example_cb_mean_avg;
        }

        for (int slope_idx = 0; slope_idx < static_cast<int>(tot_adc_slope_factors.size()); ++slope_idx) {
            double slope = tot_adc_slope_factors[slope_idx];
            for (int offset_idx = 0; offset_idx < static_cast<int>(tot_adc_offset_factors.size()); ++offset_idx) {
                double offset = tot_adc_offset_factors[offset_idx];
                auto & combined_values = run_combined_values[slope_idx * tot_adc_offset_factors.size() + offset_idx];
                auto hist_color = color_list[(slope_idx * tot_adc_offset_factors.size() + offset_idx) % color_list.size()];

                TH1D* hist_combined = new TH1D(
                    fmt::format("hist_combined_run{}_slope{}_offset{}", run_numbers[run_idx], slope, offset).c_str(),
                    fmt::format("Combined ADC Values - Run {} (slope={}, offset={})", run_numbers[run_idx], slope, offset).c_str(),
                    256, 0.0, max_90_percent_x * 1.6);
                hist_combined->SetStats(kFALSE);
                hist_combined->SetTitle("");
                for (const auto& val : combined_values) {
                    hist_combined->Fill(val);
                }
                hist_combined->SetLineColor(hist_color);
                hist_combined->Draw("SAME");

                legend_combined->AddEntry(
                    hist_combined,
                    fmt::format("Combined (slope={:.1f}, offset={:.1f})", slope, offset).c_str(),
                    "l");

                double cb_mean_avg=0.0;
                double cb_mean_err_sys=0.0;
                double cb_mean_err_stat=0.0;
                double cb_sigma_avg=0.0;
                double cb_sigma_err_sys=0.0;
                double cb_sigma_err_stat=0.0;
                double cb_effective_sigma_avg=0.0;
                double cb_effective_sigma_err_sys=0.0;
                double cb_effective_sigma_err_stat=0.0;
                double cb_resolution_avg=0.0;
                double cb_resolution_err=0.0;

                crystalball_fit_th1d(*canvas_combined, *hist_combined, fit_range_sigmas, fit_range_offsets, hist_color, cb_mean_avg, cb_mean_err_sys, cb_mean_err_stat, cb_sigma_avg, cb_sigma_err_sys, cb_sigma_err_stat, cb_resolution_avg, cb_resolution_err);

                combined_cb_mean_avgs[slope_idx * tot_adc_offset_factors.size() + offset_idx].push_back(cb_mean_avg);
                combined_cb_mean_errs[slope_idx * tot_adc_offset_factors.size() + offset_idx].push_back(std::sqrt(cb_mean_err_sys * cb_mean_err_sys + cb_mean_err_stat * cb_mean_err_stat));

                if (cb_mean_avg > max_fit_mean) {
                    max_fit_mean = cb_mean_avg;
                }
            }
        }
        legend_combined->Draw();
        canvas_combined->Print(out_pdf.c_str());
        canvas_combined->Write();
        canvas_combined->Close();
    }

    // ! === draw canvas with raw and combined linearity ====================
    std::vector<std::vector<double>> nonlinearity_combined_cb_mean_avgs(tot_adc_offset_factors.size(), std::vector<double>(tot_adc_slope_factors.size(), 0.0));
    double nonlinearity_raw_cb_mean_max = 0.0;
    double max_nonlinearity = 0.0;

    std::vector<std::vector<double>> fit_offset_cmbined_cb_mean_avgs(tot_adc_offset_factors.size());
    std::vector<std::vector<double>> fit_slope_cmbined_cb_mean_avgs(tot_adc_slope_factors.size());
    std::vector<double> fit_offset_raw_cb_mean_avgs;
    std::vector<double> fit_slope_raw_cb_mean_avgs;
    double max_fit_offset = 0.0;
    double max_fit_slope = 0.0;

    TCanvas* canvas_final_linearity = new TCanvas("canvas_final_linearity", "Final Linearity Comparison", 800, 600);
    canvas_final_linearity->SetLeftMargin(0.15);
    canvas_final_linearity->SetBottomMargin(0.10);
    TLegend* legend_final_linearity = new TLegend(0.52, 0.45, 0.89, 0.89);
    legend_final_linearity->SetBorderSize(0);
    legend_final_linearity->SetFillStyle(0);
    TGraphErrors* graph_final_linearity_raw = new TGraphErrors();
    for (size_t i = 0; i < raw_cb_mean_avgs.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y = raw_cb_mean_avgs[i];
        double yerr = raw_cb_mean_errs[i];
        graph_final_linearity_raw->SetPoint(i, x, y);
        graph_final_linearity_raw->SetPointError(i, 0.03 * x, yerr);
    }
    graph_final_linearity_raw->SetMarkerStyle(20);
    graph_final_linearity_raw->SetMarkerColor(kBlack);
    graph_final_linearity_raw->SetLineColor(kBlack);
    graph_final_linearity_raw->SetTitle(";Beam Energy [GeV];Crystal Ball Fit Mean [ADC counts]");
    graph_final_linearity_raw->GetXaxis()->SetTitleSize(axis_title_size);
    graph_final_linearity_raw->GetYaxis()->SetTitleSize(axis_title_size);
    graph_final_linearity_raw->GetXaxis()->SetLabelSize(axis_label_size);
    graph_final_linearity_raw->GetYaxis()->SetLabelSize(axis_label_size);
    graph_final_linearity_raw->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);
    graph_final_linearity_raw->GetYaxis()->SetRangeUser(0.0, max_fit_mean * 1.2);
    graph_final_linearity_raw->Draw("AP");

    // do linear fit for raw
    TF1* fit_final_linearity_raw = new TF1("fit_final_linearity_raw", "pol1", energy_axis_min, energy_axis_max);
    fit_final_linearity_raw->SetLineColor(kBlack);
    // set initial parameters from the average of the data points
    init_a = (raw_cb_mean_avgs.back() - raw_cb_mean_avgs.front()) /
                    ( (found_run_energies ? static_cast<double>(run_energies.back()) : static_cast<double>(raw_cb_mean_avgs.size()-1)) -
                      (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0) );
    init_b = raw_cb_mean_avgs.front() - init_a * (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0);
    fit_final_linearity_raw->SetParameters(init_b, init_a);
    graph_final_linearity_raw->Fit(fit_final_linearity_raw, "R");
    graph_final_linearity_raw->Draw("P SAME");

    double nonlinearity_raw_max = 0.0;
    for (size_t i = 0; i < raw_cb_mean_avgs.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y = raw_cb_mean_avgs[i];
        double y_fit = fit_final_linearity_raw->Eval(x);
        double nonlinearity = fabs((y - y_fit) / y_fit) * 100.0; // in percentage
        if (nonlinearity > nonlinearity_raw_max) {
            nonlinearity_raw_max = nonlinearity;
        }
    }
    nonlinearity_raw_cb_mean_max = nonlinearity_raw_max;
    if (nonlinearity_raw_cb_mean_max > max_nonlinearity) {
        max_nonlinearity = nonlinearity_raw_cb_mean_max;
    }

    double fit_offset_raw = fit_final_linearity_raw->GetParameter(0);
    double fit_slope_raw = fit_final_linearity_raw->GetParameter(1);
    fit_offset_raw_cb_mean_avgs.push_back(fit_offset_raw);
    fit_slope_raw_cb_mean_avgs.push_back(fit_slope_raw);
    if (fabs(fit_offset_raw) > max_fit_offset) {
        max_fit_offset = fabs(fit_offset_raw);
    }
    if (fabs(fit_slope_raw) > max_fit_slope) {
        max_fit_slope = fabs(fit_slope_raw);
    }

    legend_final_linearity->AddEntry(graph_final_linearity_raw, Form("Raw ADC Sum Fit: ADC = %.2fE + %.2f (%.2f%%)", fit_final_linearity_raw->GetParameter(1), fit_final_linearity_raw->GetParameter(0), nonlinearity_raw_max), "l");

    for (size_t combo_idx = 0; combo_idx < combined_cb_mean_avgs.size(); ++combo_idx) {
        TGraphErrors* graph_final_linearity_combined = new TGraphErrors();
        for (size_t i = 0; i < combined_cb_mean_avgs[combo_idx].size(); ++i) {
            double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
            double y = combined_cb_mean_avgs[combo_idx][i];
            double yerr = combined_cb_mean_errs[combo_idx][i];
            graph_final_linearity_combined->SetPoint(i, x, y);
            graph_final_linearity_combined->SetPointError(i, 0.03 * x, yerr);
        }
        int color = color_list[combo_idx % color_list.size()];
        graph_final_linearity_combined->SetMarkerStyle(21 + combo_idx % 4);
        graph_final_linearity_combined->SetMarkerColor(color);
        graph_final_linearity_combined->SetLineColor(color);
        graph_final_linearity_combined->Draw("P SAME");

        // do linear fit for combined
        TF1* fit_final_linearity_combined = new TF1(
            fmt::format("fit_final_linearity_combined_{}", combo_idx).c_str(),
            "pol1", energy_axis_min, energy_axis_max);
        fit_final_linearity_combined->SetLineColor(color);
        // set initial parameters from the average of the data points
        init_a = (combined_cb_mean_avgs[combo_idx].back() - combined_cb_mean_avgs[combo_idx].front()) /
                        ( (found_run_energies ? static_cast<double>(run_energies.back()) : static_cast<double>(combined_cb_mean_avgs[combo_idx].size()-1)) -
                          (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0) );
        init_b = combined_cb_mean_avgs[combo_idx].front() - init_a * (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0);
        fit_final_linearity_combined->SetParameters(init_b, init_a);
        graph_final_linearity_combined->Fit(fit_final_linearity_combined, "R");
        graph_final_linearity_combined->Draw("P SAME");

        double nonlinearity_combined_max = 0.0;
        for (size_t i = 0; i < combined_cb_mean_avgs[combo_idx].size(); ++i) {
            double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
            double y = combined_cb_mean_avgs[combo_idx][i];
            double y_fit = fit_final_linearity_combined->Eval(x);
            double nonlinearity = fabs((y - y_fit) / y_fit) * 100.0; // in percentage
            if (nonlinearity > nonlinearity_combined_max) {
                nonlinearity_combined_max = nonlinearity;
            }
        }

        int offset_idx = combo_idx % tot_adc_offset_factors.size();
        int slope_idx = combo_idx / tot_adc_offset_factors.size();
        nonlinearity_combined_cb_mean_avgs[offset_idx][slope_idx] = nonlinearity_combined_max;
        if (nonlinearity_combined_max > max_nonlinearity) {
            max_nonlinearity = nonlinearity_combined_max;
        }

        double fit_offset_combined = fit_final_linearity_combined->GetParameter(0);
        double fit_slope_combined = fit_final_linearity_combined->GetParameter(1);
        fit_offset_cmbined_cb_mean_avgs[offset_idx].push_back(fit_offset_combined);
        fit_slope_cmbined_cb_mean_avgs[slope_idx].push_back(fit_slope_combined);
        if (fabs(fit_offset_combined) > max_fit_offset) {
            max_fit_offset = fabs(fit_offset_combined);
        }
        if (fabs(fit_slope_combined) > max_fit_slope) {
            max_fit_slope = fabs(fit_slope_combined);
        }

        legend_final_linearity->AddEntry(
            graph_final_linearity_combined,
            fmt::format("Combined Fit (slope={:.1f}, offset={:.1f}): ADC = {:.2f}E + {:.2f} (max nonlinearity = {:.2f}%)",
                        tot_adc_slope_factors[combo_idx / tot_adc_offset_factors.size()],
                        tot_adc_offset_factors[combo_idx % tot_adc_offset_factors.size()],
                        fit_final_linearity_combined->GetParameter(1),
                        fit_final_linearity_combined->GetParameter(0),
                        nonlinearity_combined_max).c_str(),
            "l");
    }

    latex_linearity.SetTextSize(0.05);
    latex_linearity.SetTextFont(62);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    latex_linearity.SetTextSize(0.035);
    latex_linearity.SetTextFont(42);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step , annotation_testbeam_title.c_str());
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, config_description.c_str());
    // write date
    if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y ", std::localtime(&now_c))) {
        latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    }

    legend_final_linearity->Draw();
    canvas_final_linearity->Print(out_pdf.c_str());
    canvas_final_linearity->Write();
    canvas_final_linearity->Close();


    // ! === draw nonlinearity graph, each offset is one line ===
    max_nonlinearity = 10.0;
    TCanvas* canvas_final_nonlinearity = new TCanvas("canvas_final_nonlinearity", "Final Non-Linearity Comparison", 800, 600);
    canvas_final_nonlinearity->SetLeftMargin(0.15);
    canvas_final_nonlinearity->SetBottomMargin(0.10);
    TLegend* legend_final_nonlinearity = new TLegend(0.55, 0.65, 0.89, 0.89);
    legend_final_nonlinearity->SetBorderSize(0);
    legend_final_nonlinearity->SetFillStyle(0);

    for (size_t offset_idx = 0; offset_idx < tot_adc_offset_factors.size(); ++offset_idx) {
        TGraph* graph_final_nonlinearity = new TGraph();
        for (size_t slope_idx = 0; slope_idx < tot_adc_slope_factors.size(); ++slope_idx) {
            double x = tot_adc_slope_factors[slope_idx];
            double y = nonlinearity_combined_cb_mean_avgs[offset_idx][slope_idx];
            graph_final_nonlinearity->SetPoint(slope_idx, x, y);
        }
        int color = color_list[offset_idx % color_list.size()];
        graph_final_nonlinearity->SetMarkerStyle(21 + offset_idx % 4);
        graph_final_nonlinearity->SetMarkerColor(color);
        graph_final_nonlinearity->SetLineColor(color);
        graph_final_nonlinearity->SetTitle(";TOT ADC Slope Factor;Max Non-Linearity [%]");
        graph_final_nonlinearity->GetXaxis()->SetTitleSize(axis_title_size);
        graph_final_nonlinearity->GetYaxis()->SetTitleSize(axis_title_size);
        graph_final_nonlinearity->GetXaxis()->SetLabelSize(axis_label_size);
        graph_final_nonlinearity->GetYaxis()->SetLabelSize(axis_label_size);
        graph_final_nonlinearity->GetXaxis()->SetLimits(
            tot_adc_slope_factors.front() - 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()),
            tot_adc_slope_factors.back() + 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()));
        graph_final_nonlinearity->GetYaxis()->SetRangeUser(0.0, max_nonlinearity * 1.3);
        if (offset_idx == 0) {
            graph_final_nonlinearity->Draw("AP");
        } else {
            graph_final_nonlinearity->Draw("P SAME");
        }
        legend_final_nonlinearity->AddEntry(
            graph_final_nonlinearity,
            fmt::format("Offset={:.1f}", tot_adc_offset_factors[offset_idx]).c_str(),
            "lp");
    }

    // draw raw non-linearity as horizontal line
    TLine* line_raw_nonlinearity = new TLine(
        tot_adc_slope_factors.front() - 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()),
        nonlinearity_raw_cb_mean_max,
        tot_adc_slope_factors.back() + 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()),
        nonlinearity_raw_cb_mean_max);
    line_raw_nonlinearity->SetLineStyle(2);
    line_raw_nonlinearity->SetLineColor(kBlack);
    line_raw_nonlinearity->Draw("SAME");
    legend_final_nonlinearity->AddEntry(
        line_raw_nonlinearity,
        fmt::format("Raw ADC Sum Non-Linearity: {:.2f}%", nonlinearity_raw_cb_mean_max).c_str(),
        "l");

    latex_linearity.SetTextSize(0.05);
    latex_linearity.SetTextFont(62);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    latex_linearity.SetTextSize(0.035);
    latex_linearity.SetTextFont(42);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step , annotation_testbeam_title.c_str());
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, config_description.c_str());
    // write date
    if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y ", std::localtime(&now_c))) {
        latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    }
    
    legend_final_nonlinearity->Draw();
    canvas_final_nonlinearity->Print(out_pdf.c_str());
    canvas_final_nonlinearity->Write();
    canvas_final_nonlinearity->Close();

    // draw the offset in the same style
    TCanvas* canvas_final_fit_offset = new TCanvas("canvas_final_fit_offset", "Final Fit Offset Comparison", 800, 600);
    canvas_final_fit_offset->SetLeftMargin(0.15);
    canvas_final_fit_offset->SetBottomMargin(0.10);
    TLegend* legend_final_fit_offset = new TLegend(0.55, 0.65, 0.89, 0.89);
    legend_final_fit_offset->SetBorderSize(0);
    legend_final_fit_offset->SetFillStyle(0);

    for (size_t offset_idx = 0; offset_idx < tot_adc_offset_factors.size(); ++offset_idx) {
        TGraph* graph_final_fit_offset = new TGraph();
        for (size_t slope_idx = 0; slope_idx < tot_adc_slope_factors.size(); ++slope_idx) {
            double x = tot_adc_slope_factors[slope_idx];
            double y = fit_offset_cmbined_cb_mean_avgs[offset_idx][slope_idx];
            graph_final_fit_offset->SetPoint(slope_idx, x, y);
        }
        int color = color_list[offset_idx % color_list.size()];
        graph_final_fit_offset->SetMarkerStyle(21 + offset_idx % 4);
        graph_final_fit_offset->SetMarkerColor(color);
        graph_final_fit_offset->SetLineColor(color);
        
        if (offset_idx == 0) {
            graph_final_fit_offset->SetTitle(";TOT ADC Slope Factor;Fit Offset [ADC counts]");
            graph_final_fit_offset->Draw("AP");
            
            // Set axis properties after first draw
            graph_final_fit_offset->GetXaxis()->SetTitleSize(axis_title_size);
            graph_final_fit_offset->GetYaxis()->SetTitleSize(axis_title_size);
            graph_final_fit_offset->GetXaxis()->SetLabelSize(axis_label_size);
            graph_final_fit_offset->GetYaxis()->SetLabelSize(axis_label_size);
            graph_final_fit_offset->GetXaxis()->SetLimits(
                tot_adc_slope_factors.front() - 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()),
                tot_adc_slope_factors.back() + 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()));
            graph_final_fit_offset->GetYaxis()->SetRangeUser(-max_fit_offset * 1.3, max_fit_offset * 1.3);
        } else {
            graph_final_fit_offset->Draw("P SAME");
        }
        legend_final_fit_offset->AddEntry(
            graph_final_fit_offset,
            fmt::format("Offset={:.1f}", tot_adc_offset_factors[offset_idx]).c_str(),
            "lp");
    }

    // draw raw fit offset as horizontal line
    TLine* line_raw_fit_offset = new TLine(
        tot_adc_slope_factors.front() - 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()),
        fit_offset_raw_cb_mean_avgs[0],
        tot_adc_slope_factors.back() + 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()),
        fit_offset_raw_cb_mean_avgs[0]);
    line_raw_fit_offset->SetLineStyle(2);
    line_raw_fit_offset->SetLineColor(kBlack);
    line_raw_fit_offset->Draw("SAME");
    latex_linearity.SetTextSize(0.05);
    latex_linearity.SetTextFont(62);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    latex_linearity.SetTextSize(0.035);
    latex_linearity.SetTextFont(42);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step , annotation_testbeam_title.c_str());
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, config_description.c_str());
    // write date
    if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y ", std::localtime(&now_c))) {
        latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    }
    legend_final_fit_offset->Draw();
    canvas_final_fit_offset->Print(out_pdf.c_str());
    canvas_final_fit_offset->Write();
    canvas_final_fit_offset->Close();

    TCanvas* canvas_final_fit_slope = new TCanvas("canvas_final_fit_slope", "Final Fit Slope Comparison", 800, 600);
    canvas_final_fit_slope->SetLeftMargin(0.15);
    canvas_final_fit_slope->SetBottomMargin(0.10);
    TLegend* legend_final_fit_slope = new TLegend(0.55, 0.65, 0.89, 0.89);
    legend_final_fit_slope->SetBorderSize(0);
    legend_final_fit_slope->SetFillStyle(0);

    for (size_t offset_idx = 0; offset_idx < tot_adc_offset_factors.size(); ++offset_idx) {
        TGraph* graph_final_fit_slope = new TGraph();
        for (size_t slope_idx = 0; slope_idx < tot_adc_slope_factors.size(); ++slope_idx) {
            double x = tot_adc_slope_factors[slope_idx];
            double y = fit_slope_cmbined_cb_mean_avgs[slope_idx][offset_idx];
            graph_final_fit_slope->SetPoint(slope_idx, x, y);
        }
        int color = color_list[offset_idx % color_list.size()];
        graph_final_fit_slope->SetMarkerStyle(21 + offset_idx % 4);
        graph_final_fit_slope->SetMarkerColor(color);
        graph_final_fit_slope->SetLineColor(color);
        
        if (offset_idx == 0) {
            graph_final_fit_slope->SetTitle(";TOT ADC Slope Factor;Fit Slope [ADC counts/GeV]");
            graph_final_fit_slope->Draw("AP");
            
            // Set axis properties after first draw
            graph_final_fit_slope->GetXaxis()->SetTitleSize(axis_title_size);
            graph_final_fit_slope->GetYaxis()->SetTitleSize(axis_title_size);
            graph_final_fit_slope->GetXaxis()->SetLabelSize(axis_label_size);
            graph_final_fit_slope->GetYaxis()->SetLabelSize(axis_label_size);
            graph_final_fit_slope->GetXaxis()->SetLimits(
                tot_adc_slope_factors.front() - 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()),
                tot_adc_slope_factors.back() + 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()));
            graph_final_fit_slope->GetYaxis()->SetRangeUser(0, max_fit_slope * 1.5);
        } else {
            graph_final_fit_slope->Draw("P SAME");
        }
        legend_final_fit_slope->AddEntry(
            graph_final_fit_slope,
            fmt::format("Offset={:.1f}", tot_adc_offset_factors[offset_idx]).c_str(),
            "lp");
    }
    // draw raw fit slope as horizontal line
    TLine* line_raw_fit_slope = new TLine(
        tot_adc_slope_factors.front() - 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()),
        fit_slope_raw_cb_mean_avgs[0],
        tot_adc_slope_factors.back() + 0.1 * (tot_adc_slope_factors.back() - tot_adc_slope_factors.front()),
        fit_slope_raw_cb_mean_avgs[0]);
    line_raw_fit_slope->SetLineStyle(2);
    line_raw_fit_slope->SetLineColor(kBlack);
    line_raw_fit_slope->Draw("SAME");
    latex_linearity.SetTextSize(0.05);
    latex_linearity.SetTextFont(62);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    latex_linearity.SetTextSize(0.035);
    latex_linearity.SetTextFont(42);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step , annotation_testbeam_title.c_str());
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, config_description.c_str());
    // write date
    if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y ", std::localtime(&now_c))) {
        latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    }
    legend_final_fit_slope->Draw();
    canvas_final_fit_slope->Print(out_pdf.c_str());
    canvas_final_fit_slope->Write();
    canvas_final_fit_slope->Close();

    // draw the raw linearity and example linearity in the same canvas for comparison
    TCanvas* canvas_linearity_comparison = new TCanvas("canvas_linearity_comparison", "Linearity Comparison", 800, 600);
    canvas_linearity_comparison->SetLeftMargin(0.15);
    canvas_linearity_comparison->SetBottomMargin(0.10);
    TLegend* legend_linearity_comparison = new TLegend(0.52, 0.45, 0.89, 0.89);
    legend_linearity_comparison->SetBorderSize(0);
    legend_linearity_comparison->SetFillStyle(0);
    graph_final_linearity_raw->Draw("AP");
    legend_linearity_comparison->AddEntry(graph_final_linearity_raw, Form("Raw: ADC = %.2fE + %.2f (%.2f%%)", fit_final_linearity_raw->GetParameter(1), fit_final_linearity_raw->GetParameter(0), nonlinearity_raw_cb_mean_max), "l");
    TGraphErrors* graph_final_linearity_example = new TGraphErrors();
    for (size_t i = 0; i < example_cb_mean_avgs.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y = example_cb_mean_avgs[i];
        double yerr = example_cb_mean_errs[i];
        graph_final_linearity_example->SetPoint(i, x, y);
        graph_final_linearity_example->SetPointError(i, 0.03 * x, yerr);
    }
    graph_final_linearity_example->SetMarkerStyle(21);
    graph_final_linearity_example->SetMarkerColor(kCyan+2);
    graph_final_linearity_example->SetLineColor(kCyan+2);
    graph_final_linearity_example->Draw("P SAME");

    // linear fit for example
    TF1* fit_final_linearity_example = new TF1("fit_final_linearity_example", "pol1", energy_axis_min, energy_axis_max);
    fit_final_linearity_example->SetLineColor(kCyan+2);
    // set initial parameters from the average of the data points
    init_a = (example_cb_mean_avgs.back() - example_cb_mean_avgs.front()) /
                    ( (found_run_energies ? static_cast<double>(run_energies.back()) : static_cast<double>(example_cb_mean_avgs.size()-1)) -
                      (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0) );
    init_b = example_cb_mean_avgs.front() - init_a * (found_run_energies ? static_cast<double>(run_energies.front()) : 0.0);
    fit_final_linearity_example->SetParameters(init_b, init_a);
    graph_final_linearity_example->Fit(fit_final_linearity_example, "R");
    graph_final_linearity_example->Draw("P SAME");
    legend_linearity_comparison->AddEntry(graph_final_linearity_example, Form("Example: ADC = %.2fE + %.2f (%.2f%%)", fit_final_linearity_example->GetParameter(1), fit_final_linearity_example->GetParameter(0), [] (const TGraphErrors* graph, const TF1* fit) {
        double max_nonlinearity = 0.0;
        for (int i = 0; i < graph->GetN(); ++i) {
            double x, y;
            graph->GetPoint(i, x, y);
            double y_fit = fit->Eval(x);
            double nonlinearity = fabs((y - y_fit) / y_fit) * 100.0;
            if (nonlinearity > max_nonlinearity) {
                max_nonlinearity = nonlinearity;
            }
        }
        return max_nonlinearity;
    }(graph_final_linearity_example, fit_final_linearity_example)), "l");
    latex_linearity.SetTextSize(0.05);
    latex_linearity.SetTextFont(62);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    latex_linearity.SetTextSize(0.035);
    latex_linearity.SetTextFont(42);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step , annotation_testbeam_title.c_str());
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, config_description.c_str());
    // write date
    if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y ", std::localtime(&now_c))) {
        latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    }
    legend_linearity_comparison->Draw();
    canvas_linearity_comparison->Print(out_pdf.c_str());
    canvas_linearity_comparison->Write();
    canvas_linearity_comparison->Close();

    // draw the resolution graph of raw and example
    TCanvas* canvas_resolution_comparison = new TCanvas("canvas_resolution_comparison", "Resolution Comparison", 800, 600);
    canvas_resolution_comparison->SetLeftMargin(0.15);
    canvas_resolution_comparison->SetBottomMargin(0.10);
    TLegend* legend_resolution_comparison = new TLegend(0.52, 0.45, 0.89, 0.89);
    legend_resolution_comparison->SetBorderSize(0);
    legend_resolution_comparison->SetFillStyle(0);
    TGraphErrors* graph_resolution_comparison_raw = new TGraphErrors();
    for (size_t i = 0; i < raw_cb_resolution_avgs.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y = raw_cb_resolution_avgs[i];
        double yerr = raw_cb_resolution_errs[i];
        graph_resolution_comparison_raw->SetPoint(i, x, y);
        graph_resolution_comparison_raw->SetPointError(i, 0.03 * x, yerr);
    }
    graph_resolution_comparison_raw->SetMarkerStyle(20);
    graph_resolution_comparison_raw->SetMarkerColor(kBlack);
    graph_resolution_comparison_raw->SetLineColor(kBlack);
    graph_resolution_comparison_raw->SetTitle(";Beam Energy [GeV];Crystal Ball Fit Resolution [%]");
    graph_resolution_comparison_raw->GetXaxis()->SetTitleSize(axis_title_size);
    graph_resolution_comparison_raw->GetYaxis()->SetTitleSize(axis_title_size);
    graph_resolution_comparison_raw->GetXaxis()->SetLabelSize(axis_label_size);
    graph_resolution_comparison_raw->GetYaxis()->SetLabelSize(axis_label_size);
    graph_resolution_comparison_raw->GetXaxis()->SetLimits(energy_axis_min, energy_axis_max);
    graph_resolution_comparison_raw->GetYaxis()->SetRangeUser(6.0, 40.0);
    graph_resolution_comparison_raw->Draw("AP");

    legend_resolution_comparison->AddEntry(graph_resolution_comparison_raw, "Raw ADC Values", "l");
    TGraphErrors* graph_resolution_comparison_example = new TGraphErrors();
    for (size_t i = 0; i < example_cb_resolution_avgs.size(); ++i) {
        double x = found_run_energies ? static_cast<double>(run_energies[i]) : static_cast<double>(i);
        double y = example_cb_resolution_avgs[i];
        double yerr = example_cb_resolution_errs[i];
        graph_resolution_comparison_example->SetPoint(i, x, y);
        graph_resolution_comparison_example->SetPointError(i, 0.03 * x, yerr);
    }
    graph_resolution_comparison_example->SetMarkerStyle(21);
    graph_resolution_comparison_example->SetMarkerColor(kCyan+2);
    graph_resolution_comparison_example->SetLineColor(kCyan+2);
    graph_resolution_comparison_example->Draw("P SAME");

    legend_resolution_comparison->AddEntry(graph_resolution_comparison_example, "Example ADC Values", "l");
    latex_linearity.SetTextSize(0.05);
    latex_linearity.SetTextFont(62);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start, annotation_canvas_title.c_str());
    latex_linearity.SetTextSize(0.035);
    latex_linearity.SetTextFont(42);
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - latex_y_step , annotation_testbeam_title.c_str());
    latex_linearity.DrawLatex(latex_x_start, latex_y_start - 2 * latex_y_step, config_description.c_str());
    // write date
    if (std::strftime(date_str, sizeof(date_str), "%d-%m-%Y ", std::localtime(&now_c))) {
        latex_linearity.DrawLatex(latex_x_start, latex_y_start - 3 * latex_y_step, date_str);
    }
    legend_resolution_comparison->Draw();
    canvas_resolution_comparison->Print(out_pdf.c_str());
    canvas_resolution_comparison->Write();
    canvas_resolution_comparison->Close();


    // write dummy canvas
    TCanvas* canvas_dummy = new TCanvas("canvas_dummy", "Dummy Canvas", 800, 600);
    canvas_dummy->Print((out_pdf + ")").c_str());

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}