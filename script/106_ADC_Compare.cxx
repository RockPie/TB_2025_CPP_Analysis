#include <fstream>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"

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

    const std::string run_file_folder = "dump/102_EventMatch/beamtests/";
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
        TTree *input_tree = (TTree*) input_root->Get("data_tree");
        if (input_tree == nullptr) {
            spdlog::error("Failed to get data tree from input file {}", run_file);
            return 1;
        }
        run_file_event_counts.push_back(input_tree->GetEntries());
        input_root->Close();
    }
    // find the smallest number of events among all run files
    int min_run_file_event_count = *std::min_element(run_file_event_counts.begin(), run_file_event_counts.end());
    if (min_run_file_event_count < 1) {
        spdlog::error("No events to process");
        // create empty output file
        auto output_root = new TFile(script_output_file.c_str(), "RECREATE");
        if (output_root->IsZombie()) {
            spdlog::error("Failed to create output file {}", script_output_file);
            return 1;
        }
        output_root->Close();
        return 1;
    }
    spdlog::info("Minimum number of events among all run files: {}", min_run_file_event_count);
    if (target_event_number > 0 && target_event_number < min_run_file_event_count) {
        min_run_file_event_count = target_event_number;
        spdlog::info("Using target number of events from config: {}", min_run_file_event_count);
    }

    std::vector<std::vector<Double_t>*> files_adc_sums;
    for (int _file_index = 0; _file_index < run_files.size(); _file_index++) {
        auto* adc_sums = new std::vector<Double_t>;
        adc_sums->reserve(min_run_file_event_count);
        files_adc_sums.push_back(adc_sums);
    }

    std::vector<Double_t> adc_sums_sorting;
    adc_sums_sorting.reserve(min_run_file_event_count * run_files.size());

    // then second round to actually read the data
    for (int _file_index = 0; _file_index < run_files.size(); _file_index++) {
        const auto& run_file = run_files[_file_index];
        spdlog::info("Processing run file: {}", run_file);
        
        auto adc_sums = files_adc_sums[_file_index];

        TFile *input_root = new TFile(run_file.c_str(), "READ");
        if (input_root->IsZombie()) {
            spdlog::error("Failed to open input file {}", run_file);
            return 1;
        }
        TTree *input_tree = (TTree*) input_root->Get("data_tree");
        if (input_tree == nullptr) {
            spdlog::error("Failed to get data tree from input file {}", run_file);
            return 1;
        }

        std::vector<std::vector<ULong64_t>> timestamp_pools(vldb_number);
        std::vector<std::vector<ULong64_t*>> timestamp_pools_original(vldb_number);
        std::vector<std::vector<UInt_t*>> daqh_list_pools(vldb_number);
        std::vector<std::vector<Bool_t*>> tc_list_pools(vldb_number);
        std::vector<std::vector<Bool_t*>> tp_list_pools(vldb_number);
        std::vector<std::vector<UInt_t*>> val0_list_pools(vldb_number);
        std::vector<std::vector<UInt_t*>> val1_list_pools(vldb_number);
        std::vector<std::vector<UInt_t*>> val2_list_pools(vldb_number);
        std::vector<std::vector<UInt_t*>> crc32_list_pools(vldb_number);
        std::vector<std::vector<UInt_t*>> last_heartbeat_pools(vldb_number);

        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            auto *branch_timestamps = new ULong64_t[machine_gun_samples];  // 64 bits
            auto *branch_daqh_list  = new UInt_t[4 * machine_gun_samples];    // 32 bits
            auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
            auto *branch_crc32_list = new UInt_t[4 * machine_gun_samples];
            auto *branch_last_heartbeat = new UInt_t[machine_gun_samples];
            input_tree->SetBranchAddress(("timestamps_"    + std::to_string(vldb_id)).c_str(), branch_timestamps);
            input_tree->SetBranchAddress(("daqh_list_"     + std::to_string(vldb_id)).c_str(), branch_daqh_list);
            input_tree->SetBranchAddress(("tc_list_"       + std::to_string(vldb_id)).c_str(), branch_tc_list);
            input_tree->SetBranchAddress(("tp_list_"       + std::to_string(vldb_id)).c_str(), branch_tp_list);
            input_tree->SetBranchAddress(("val0_list_"     + std::to_string(vldb_id)).c_str(), branch_val0_list);
            input_tree->SetBranchAddress(("val1_list_"     + std::to_string(vldb_id)).c_str(), branch_val1_list);
            input_tree->SetBranchAddress(("val2_list_"     + std::to_string(vldb_id)).c_str(), branch_val2_list);
            input_tree->SetBranchAddress(("crc32_list_"    + std::to_string(vldb_id)).c_str(), branch_crc32_list);
            input_tree->SetBranchAddress(("last_heartbeat_"+ std::to_string(vldb_id)).c_str(), branch_last_heartbeat);
            timestamp_pools_original[vldb_id].push_back(branch_timestamps);
            daqh_list_pools[vldb_id].push_back(branch_daqh_list);
            tc_list_pools[vldb_id].push_back(branch_tc_list);
            tp_list_pools[vldb_id].push_back(branch_tp_list);
            val0_list_pools[vldb_id].push_back(branch_val0_list);
            val1_list_pools[vldb_id].push_back(branch_val1_list);
            val2_list_pools[vldb_id].push_back(branch_val2_list);
            crc32_list_pools[vldb_id].push_back(branch_crc32_list);
            last_heartbeat_pools[vldb_id].push_back(branch_last_heartbeat);
        }

        for (int _entry = 0; _entry < min_run_file_event_count; _entry++) {
            input_tree->GetEntry(_entry);
            // process data here
            Double_t adc_sum = 0.0;
            for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
                for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                    double adc_max = 0;
                    double adc_pedestal = 0;
                    std::vector<int> adc_samples;
                    adc_samples.reserve(machine_gun_samples);
                    for (int sample = 0; sample < machine_gun_samples; sample++) {
                        auto adc_value = val0_list_pools[vldb_id][0][FPGA_CHANNEL_NUMBER * sample + channel];
                        adc_samples.push_back(adc_value);
                        auto adc_value_double = double(adc_value);
                        if (adc_value_double > adc_max) {
                            adc_max = adc_value_double;
                        }
                    } // end of sample loop
                    // calculate pedestal as average of the smallest 2 samples among first 3 samples
                    pedestal_subtraction_2minoffirst3(adc_samples, adc_pedestal);
                    if (adc_max <= adc_pedestal) {
                        continue;
                    }
                    // spdlog::debug("Event {}: VLDB {} Channel {}: Pedestal = {}, Peak = {}", _entry, vldb_id, channel, adc_pedestal, adc_max);
                    adc_sum += Double_t(adc_max - adc_pedestal);
                } // end of channel loop
            } // end of vldb_id loop
            // spdlog::info("Event {}: ADC sum = {}", _entry, adc_sum);
            adc_sums->push_back(adc_sum);
            adc_sums_sorting.push_back(adc_sum);
        } // end of entry loop
        input_root->Close();
    }

    // make a color list
    std::vector<int> color_list = {kRed, kBlue, kGreen+2, kMagenta+2, kCyan+2, kOrange+7, kViolet+7, kSpring+7, kTeal+7, kAzure+7, kPink+7};

    // * calculate the 90% largest ADC sum
    std::sort(adc_sums_sorting.begin(), adc_sums_sorting.end());
    Double_t adc_sum_90pct = adc_sums_sorting[adc_sums_sorting.size() * 9 / 10];
    spdlog::info("90% largest ADC sum among all events: {}", adc_sum_90pct);
    const double hist_max_adc_sum = adc_sum_90pct * 1.3;

    // * create array of hist1d for each run file
    double hist_max_bin_value = 0;
    std::vector<TH1D*> run_adc_sum_hist1d_list;
    const int hist_bin_number = 256;
    for (int _file_index = 0; _file_index < run_files.size();_file_index++) {
        auto& adc_sums = files_adc_sums[_file_index];
        if (adc_sums->empty()) {
            spdlog::warn("No ADC sums calculated for run file {}, skipping histogram creation", run_files[_file_index]);
            continue;
        }
        TH1D* run_adc_sum_hist1d = new TH1D(
            fmt::format("run_adc_sum_{}", _file_index).c_str(),
            fmt::format("Run ADC Sum {};ADC Sum;Events", _file_index).c_str(),
            hist_bin_number,
            0, hist_max_adc_sum
        );
        run_adc_sum_hist1d->SetStats(0);
        run_adc_sum_hist1d->SetTitle("");
        for (const auto& adc_sum : *adc_sums) {
            run_adc_sum_hist1d->Fill(adc_sum);
        }
        run_adc_sum_hist1d_list.push_back(run_adc_sum_hist1d);
        if (run_adc_sum_hist1d->GetMaximum() > hist_max_bin_value) {
            hist_max_bin_value = run_adc_sum_hist1d->GetMaximum();
        }
    }

    std::vector<std::vector<TF1*>> run_adc_sum_hist1d_gauss_fit_list;
    std::vector<TF1*> run_adc_sum_hist1d_pre_fit_list;

    std::vector<double> gauss_fit_sigma_ranges_left  = {0.3, 0.4, 0.5};
    std::vector<double> gauss_fit_sigma_ranges_right = {2.0, 2.5, 3.0};
    if (enable_gaussian_fit) {
        for (int _hist_index = 0; _hist_index < run_adc_sum_hist1d_list.size(); _hist_index++) {
            auto* hist = run_adc_sum_hist1d_list[_hist_index];
            // pre-fit to find the peak position
            double hist_mean = hist->GetMean();
            double hist_rms = hist->GetRMS();
            double fit_min = hist_mean - 1.5 * hist_rms;
            if (fit_min < 0) fit_min = 0;
            double fit_max = hist_mean + 1.5 * hist_rms;
            TF1* pre_fit = new TF1(
                fmt::format("pre_fit_{}", _hist_index).c_str(),
                "gaus",
                fit_min, fit_max
            );
            pre_fit->SetParameters(hist->GetMaximum(), hist_mean, hist_rms);
            hist->Fit(pre_fit, "RQ0");
            run_adc_sum_hist1d_pre_fit_list.push_back(pre_fit);
            spdlog::info("Run {}: Pre-fit mean = {}, sigma = {}", _hist_index, pre_fit->GetParameter(1), pre_fit->GetParameter(2));
            // gaussian fit around the peak
            for (auto& _fit_left_sigma : gauss_fit_sigma_ranges_left) {
                for (auto& _fit_right_sigma : gauss_fit_sigma_ranges_right) {
                    double fit_min = pre_fit->GetParameter(1) - _fit_left_sigma * pre_fit->GetParameter(2);
                    if (fit_min < 0) fit_min = 0;
                    double fit_max = pre_fit->GetParameter(1) + _fit_right_sigma * pre_fit->GetParameter(2);
                    TF1* gauss_fit = new TF1(
                        fmt::format("gauss_fit_{}_left{}_right{}", _hist_index, _fit_left_sigma, _fit_right_sigma).c_str(),
                        "gaus",
                        fit_min, fit_max
                    );
                    gauss_fit->SetParameters(pre_fit->GetParameter(0), pre_fit->GetParameter(1), pre_fit->GetParameter(2));
                    hist->Fit(gauss_fit, "RQ0");
                    spdlog::info("Run {}: Gaussian fit (left {} sigma, right {} sigma): mean = {}, sigma = {}", 
                        _hist_index, _fit_left_sigma, _fit_right_sigma,
                        gauss_fit->GetParameter(1), gauss_fit->GetParameter(2)
                    );
                    if (run_adc_sum_hist1d_gauss_fit_list.size() <= _hist_index) {
                        run_adc_sum_hist1d_gauss_fit_list.emplace_back();
                    }
                    run_adc_sum_hist1d_gauss_fit_list[_hist_index].push_back(gauss_fit);
                } // end of fit right sigma loop
            } // end of fit left sigma loop
        } // end of run_adc_sum_hist1d_list loop    
    } // end of gaussian fit

    // * set y axis range
    for (auto* hist : run_adc_sum_hist1d_list) {
        hist->GetYaxis()->SetRangeUser(0, hist_max_bin_value * 1.3);
    }

    // * Draw them in the same canvas
    TCanvas* canvas_adc_sum = new TCanvas("canvas_adc_sum", "Canvas ADC Sum", 800, 600);
    TLegend* legend_adc_sum = new TLegend(0.7, 0.7, 0.89, 0.89);
    for (int _file_index = 0; _file_index < run_adc_sum_hist1d_list.size(); _file_index++) {
        auto* hist = run_adc_sum_hist1d_list[_file_index];
        hist->SetLineColor(color_list[_file_index % color_list.size()]);
        hist->SetLineWidth(2);
        std::string legend_entry;
        if (found_run_energies) {
            legend_entry = fmt::format("Run {} ({} GeV)", run_numbers[_file_index], run_energies[_file_index]);
        } else {
            legend_entry = fmt::format("Run {} ({})", run_numbers[_file_index], run_labels[_file_index]);
        }
        legend_adc_sum->AddEntry(hist, legend_entry.c_str(), "l");
        if (_file_index == 0) {
            hist->Draw("HIST");
        } else {
            hist->Draw("HIST SAME");
        }
        // draw all gaussian fits if enabled
        if (enable_gaussian_fit) {
            for (auto* gauss_fit : run_adc_sum_hist1d_gauss_fit_list[_file_index]) {
                gauss_fit->SetLineColorAlpha(color_list[_file_index % color_list.size()], 0.2);
                gauss_fit->Draw("SAME");
            }
        }
    }
    legend_adc_sum->SetBorderSize(0);
    legend_adc_sum->Draw();

    TCanvas* canvas_mean_adc_sum_vs_energy = new TCanvas("canvas_mean_adc_sum_vs_energy", "Canvas Mean ADC Sum vs Beam Energy", 800, 600);
    TLegend* legend_mean_adc_sum_vs_energy = new TLegend(0.7, 0.7, 0.89, 0.89);
    legend_mean_adc_sum_vs_energy->SetBorderSize(0);
    legend_mean_adc_sum_vs_energy->SetTextSize(0.03);

    // * If gaussian fit enabled, and beam energies provided, make a graph of mean ADC sum vs beam energy
    if (enable_gaussian_fit && found_run_energies) {
        TGraphErrors* graph_mean_adc_sum_vs_energy = new TGraphErrors();
        graph_mean_adc_sum_vs_energy->SetName("graph_mean_adc_sum_vs_energy");
        graph_mean_adc_sum_vs_energy->SetTitle("Mean ADC Sum vs Beam Energy;Beam Energy (GeV);Mean ADC Sum");
        for (int _file_index = 0; _file_index < run_adc_sum_hist1d_list.size(); _file_index++) {
            auto* hist = run_adc_sum_hist1d_list[_file_index];
            auto* pre_fit = run_adc_sum_hist1d_pre_fit_list[_file_index];
            double mean_adc_sum = pre_fit->GetParameter(1);
            double sigma_adc_sum = pre_fit->GetParameter(2);
            double beam_energy = run_energies[_file_index];
            graph_mean_adc_sum_vs_energy->SetPoint(_file_index, beam_energy, mean_adc_sum);
            graph_mean_adc_sum_vs_energy->SetPointError(_file_index, 0, sigma_adc_sum);
        }
        // sort the graph by x value (beam energy)
        graph_mean_adc_sum_vs_energy->Sort();
        // fit with a linear function
        TF1* linear_fit = new TF1("linear_fit", "[0] + [1]*x", 0, run_energies[0] * 1.2);
        linear_fit->SetParameters(0, hist_max_adc_sum / run_energies[0]);
        graph_mean_adc_sum_vs_energy->Fit(linear_fit, "RQ0");
        spdlog::info("Mean ADC Sum vs Beam Energy: Linear fit: offset = {}, slope = {}", linear_fit->GetParameter(0), linear_fit->GetParameter(1));
        
        graph_mean_adc_sum_vs_energy->SetMarkerStyle(20);
        graph_mean_adc_sum_vs_energy->SetMarkerSize(1);
        graph_mean_adc_sum_vs_energy->SetLineColor(kBlue);
        graph_mean_adc_sum_vs_energy->SetLineWidth(2);
        graph_mean_adc_sum_vs_energy->Draw("AP");
        linear_fit->SetLineColor(kRed);
        linear_fit->SetLineWidth(2);
        linear_fit->Draw("SAME");

        // print fit result on the legend
        legend_mean_adc_sum_vs_energy->AddEntry(linear_fit, fmt::format("Fit : offset = {:.1f}, slope = {:.1f}", linear_fit->GetParameter(0), linear_fit->GetParameter(1)).c_str(), "l");
        legend_mean_adc_sum_vs_energy->Draw();
    }

    // * Draw resolution vs beam energy if gaussian fit enabled and beam energies provided
    TCanvas* canvas_resolution_vs_energy = new TCanvas("canvas_resolution_vs_energy", "Canvas Resolution vs Beam Energy", 800, 600);
    TLegend *legend_resolution_vs_energy = new TLegend(0.7, 0.7, 0.89, 0.89);
    legend_resolution_vs_energy->SetBorderSize(0);
    legend_resolution_vs_energy->SetTextSize(0.03);
    if (enable_gaussian_fit && found_run_energies) {
        TGraphErrors* graph_resolution_vs_energy = new TGraphErrors();
        graph_resolution_vs_energy->SetName("graph_resolution_vs_energy");
        graph_resolution_vs_energy->SetTitle("Resolution vs Beam Energy;Beam Energy (GeV);Resolution");
        for (int _file_index = 0; _file_index < run_adc_sum_hist1d_list.size(); _file_index++) {
            auto* hist = run_adc_sum_hist1d_list[_file_index];
            auto* pre_fit = run_adc_sum_hist1d_pre_fit_list[_file_index];
            double mean_adc_sum = pre_fit->GetParameter(1);
            double sigma_adc_sum = pre_fit->GetParameter(2);
            double resolution = sigma_adc_sum / mean_adc_sum;
            double beam_energy = run_energies[_file_index];
            graph_resolution_vs_energy->SetPoint(_file_index, beam_energy, resolution);
            // assume 1% error on beam energy
            graph_resolution_vs_energy->SetPointError(_file_index, beam_energy * 0.01, resolution * 0.1);
        }
        // sort the graph by x value (beam energy)
        graph_resolution_vs_energy->Sort();
        // fit with a function sqrt(a^2/E + b^2/E^2 + c^2)
        TF1* resolution_fit = new TF1("resolution_fit", "sqrt([0]*[0]/x + [1]*[1]/(x*x) + [2]*[2])", 0, run_energies[0] * 1.2);
        resolution_fit->SetParameters(0.5, 0.03);
        graph_resolution_vs_energy->Fit(resolution_fit, "RQ0");
        spdlog::info("Resolution vs Beam Energy: Fit: a = {}, b = {}", resolution_fit->GetParameter(0), resolution_fit->GetParameter(1));   
        graph_resolution_vs_energy->SetMarkerStyle(21);
        graph_resolution_vs_energy->SetMarkerSize(1);
        graph_resolution_vs_energy->SetLineColor(kBlue);
        graph_resolution_vs_energy->SetLineWidth(2);
        graph_resolution_vs_energy->Draw("AP");
        resolution_fit->SetLineColor(kRed);
        resolution_fit->SetLineWidth(2);
        resolution_fit->Draw("SAME");

        // print fit result on the legend
        legend_resolution_vs_energy->AddEntry(resolution_fit, fmt::format("Fit: a = {:.3f}, b = {:.3f}, c = {:.3f}", resolution_fit->GetParameter(0), resolution_fit->GetParameter(1), resolution_fit->GetParameter(2)).c_str(), "l");
        legend_resolution_vs_energy->Draw();

        // fix y-axis range to 0.24 to 0.06
        graph_resolution_vs_energy->GetYaxis()->SetRangeUser(0.06, 0.24);
    }

    // * Save output file


    // TH2D* pedestal_distribution_th2d = new TH2D("pedestal_distribution", "Pedestal Distribution;Channel;ADC Value", FPGA_CHANNEL_NUMBER_VALID*vldb_number, 0, FPGA_CHANNEL_NUMBER_VALID*vldb_number, 512, 0, 512);
    // pedestal_distribution_th2d->SetStats(0);
    // pedestal_distribution_th2d->SetTitle("");

    // TH2D* peak_distribution_th2d = new TH2D("peak_distribution", "Peak Distribution;Channel;ADC Value", FPGA_CHANNEL_NUMBER_VALID*vldb_number, 0, FPGA_CHANNEL_NUMBER_VALID*vldb_number, 512, 0, 1024);
    // peak_distribution_th2d->SetStats(0);
    // peak_distribution_th2d->SetTitle("");

    // std::vector<double> event_adc_sum_list;
    // event_adc_sum_list.reserve(entry_max);
    // std::vector<double> event_adc_sum_list_copy;
    // event_adc_sum_list.reserve(entry_max);
    // // start event loop
    // for (int entry = 0; entry < entry_max; entry++) {
    //     input_tree->GetEntry(entry);
    //     double adc_sum_all = 0;
    //     for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
    //         // channel loop
    //         for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
    //             int channel_index_in_h2g = channel % 76;
    //             int channel_index_in_h2g_half = channel_index_in_h2g % 38;
    //             int asic_id = channel / 76;
    //             int half_id = channel_index_in_h2g / 38;

    //             if (vldb_id == 0 && asic_id == 1)
    //                 continue;
                
    //             if (channel_index_in_h2g_half == 19)
    //                 continue;
    //             if (channel_index_in_h2g_half == 0)
    //                 continue;

    //             int channel_index_no_CM_Calib = channel_index_in_h2g_half;
                
    //             if (channel_index_in_h2g_half > 19)
    //                 channel_index_no_CM_Calib -= 2;
    //             else if (channel_index_in_h2g_half > 0)
    //                 channel_index_no_CM_Calib -= 1;

    //             int channel_number_valid = asic_id * 72 + half_id * 36 + (channel_index_no_CM_Calib);

    //             // pedestal to be the average of the smallest 2 of the first 4 samples
    //             UInt_t adc_pedestal = 1024;
    //             std::vector<UInt_t> adc_samples;
    //             UInt_t adc_peak = 0;
    //             for (int sample = 0; sample < machine_gun_samples; sample++) {
    //                 auto &adc = val0_list_pools[vldb_id][0][channel + sample * FPGA_CHANNEL_NUMBER];
    //                 adc_samples.push_back(adc);
    //                 if (adc > adc_peak) {
    //                     adc_peak = adc;
    //                 }
    //             }
    //             std::sort(adc_samples.begin(), adc_samples.end());
    //             adc_pedestal = (adc_samples[0] + adc_samples[1]) / 2;

    //             auto adc_peak_subtracted = int(adc_peak) - int(adc_pedestal);
    //             if (adc_peak_subtracted < 0) {
    //                 adc_peak_subtracted = 0;
    //             }

    //             pedestal_distribution_th2d->Fill(vldb_id * FPGA_CHANNEL_NUMBER_VALID + channel_number_valid, adc_pedestal);
    //             peak_distribution_th2d->Fill(vldb_id * FPGA_CHANNEL_NUMBER_VALID + channel_number_valid, adc_peak_subtracted);

    //             adc_sum_all += adc_peak_subtracted;
    //         } // end of channel loop
    //     } // end of vldb loop
    //     event_adc_sum_list.push_back(adc_sum_all);
    //     event_adc_sum_list_copy.push_back(adc_sum_all);
    // } // end of event loop

    // // calculate the 90% max of all adc sums
    // // sort the adc sums
    // if (event_adc_sum_list.empty()) {
    //     spdlog::warn("No events processed, using default ADC sum range.");
    //     event_adc_sum_list.push_back(0.0);
    // }
    // std::sort(event_adc_sum_list.begin(), event_adc_sum_list.end());
    // // get the 90% max

    // double adc_sum_max = event_adc_sum_list[std::min(static_cast<size_t>(event_adc_sum_list.size() * 0.9), event_adc_sum_list.size() - 1)];
    // adc_sum_max = adc_sum_max * 1.3; // add 30% margin
    // // make 10 histograms each for 1/10 of the events
    // std::vector<TH1D*> run_adc_sum_partial_hist1d_list;

    // spdlog::info("ADC sum max: {}", adc_sum_max);

    // TH1D* run_adc_sum_all_hist1d = new TH1D("run_adc_sum_all", "Run ADC Sum All;ADC Sum;Counts", 256, 0, adc_sum_max);
    // run_adc_sum_all_hist1d->SetStats(0);
    // run_adc_sum_all_hist1d->SetTitle("");

    // double pre_fit_mean = run_adc_sum_all_hist1d->GetMean();
    // double pre_fit_sigma = run_adc_sum_all_hist1d->GetRMS();
    // double fit_min = std::max(0.0, pre_fit_mean - 2 * pre_fit_sigma);
    // double fit_max = std::min(adc_sum_max, pre_fit_mean + 2 * pre_fit_sigma);

    // TF1* gaus_fit_pre = new TF1("gaus_fit_pre", "gaus", fit_min, fit_max);
    // gaus_fit_pre->SetLineColorAlpha(kWhite, 0.0);
    // run_adc_sum_all_hist1d->Fit(gaus_fit_pre, "RQ");

    // double fit_mean = gaus_fit_pre->GetParameter(1);
    // double fit_sigma = gaus_fit_pre->GetParameter(2);

    // std::vector<double> fit_sigmas = {2.0, 2.5, 3.0};
    // std::vector<double> fit_offsets = {0.2, 0.4, 0.6};
    // std::vector<double> fit_results_means;
    // std::vector<double> fit_results_sigmas;
    // for (auto fit_sigma_multiplier : fit_sigmas) {
    //     for (auto fit_offset_multiplier : fit_offsets) {
    //         double fit_min_final = std::max(0.0, fit_mean - fit_sigma_multiplier * fit_sigma + fit_offset_multiplier * fit_sigma);
    //         double fit_max_final = std::min(adc_sum_max, fit_mean + fit_sigma_multiplier * fit_sigma + fit_offset_multiplier * fit_sigma);
    //         TF1* gaus_fit_final = new TF1(("gaus_fit_final_" + std::to_string(static_cast<int>(fit_sigma_multiplier * 10)) + "_" + std::to_string(static_cast<int>(fit_offset_multiplier * 10))).c_str(), "gaus", fit_min_final, fit_max_final);
    //         // set transparency
    //         gaus_fit_final->SetLineColorAlpha(kPink, 0.2);
    //         run_adc_sum_all_hist1d->Fit(gaus_fit_final, "RQ+");
    //         double final_fit_mean = gaus_fit_final->GetParameter(1);
    //         double final_fit_sigma = gaus_fit_final->GetParameter(2);
    //         fit_results_means.push_back(final_fit_mean);
    //         fit_results_sigmas.push_back(final_fit_sigma);
    //         spdlog::info("Final Fit with sigma range {} and offset {}: mean = {}, sigma = {}", fit_sigma_multiplier, fit_offset_multiplier, final_fit_mean, final_fit_sigma);
    //     }
    // }

    // // calculate the average of the fit results
    // double final_mean = 0.0;
    // double final_sigma = 0.0;
    // double final_mean_error = 0.0;
    // double final_sigma_error = 0.0;
    // if (!fit_results_means.empty()) {
    //     final_mean = std::accumulate(fit_results_means.begin(), fit_results_means.end(), 0.0) / fit_results_means.size();
    //     final_sigma = std::accumulate(fit_results_sigmas.begin(), fit_results_sigmas.end(), 0.0) / fit_results_sigmas.size();
    // } else {
    //     final_mean = fit_mean;
    //     final_sigma = fit_sigma;
    // }
    // if (fit_results_means.size() > 1) {
    //     double mean_variance = 0.0;
    //     double sigma_variance = 0.0;
    //     for (size_t i = 0; i < fit_results_means.size(); i++) {
    //         mean_variance += (fit_results_means[i] - final_mean) * (fit_results_means[i] - final_mean);
    //         sigma_variance += (fit_results_sigmas[i] - final_sigma) * (fit_results_sigmas[i] - final_sigma);
    //     }
    //     mean_variance /= (fit_results_means.size() - 1);
    //     sigma_variance /= (fit_results_sigmas.size() - 1);
    //     final_mean_error = std::sqrt(mean_variance / fit_results_means.size());
    //     final_sigma_error = std::sqrt(sigma_variance / fit_results_sigmas.size());
    // } else {
    //     final_mean_error = gaus_fit_pre->GetParError(1);
    //     final_sigma_error = gaus_fit_pre->GetParError(2);
    // }

    auto output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    output_root->cd();

    // save the adc sum canvas
    canvas_adc_sum->cd();
    TLatex latex_adc_sum;
    latex_adc_sum.SetTextColor(kGray+2);
    latex_adc_sum.SetNDC();
    latex_adc_sum.SetTextSize(0.05);
    latex_adc_sum.SetTextFont(62);
    latex_adc_sum.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());
    latex_adc_sum.SetTextSize(0.035);
    latex_adc_sum.SetTextFont(42);
    latex_adc_sum.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
    // write run number, date time
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream date_stream;
    date_stream << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    latex_adc_sum.DrawLatex(0.12, 0.76, config_description.c_str());
    // hardon or electron beam
    latex_adc_sum.DrawLatex(0.12, 0.72, ("Beam: " + beam_type).c_str());
    latex_adc_sum.DrawLatex(0.12, 0.68, date_stream.str().c_str());
    canvas_adc_sum->Write();
    canvas_adc_sum->Close();

    if (enable_gaussian_fit && found_run_energies) {
        canvas_mean_adc_sum_vs_energy->cd();
        TLatex mean_adc_sum_vs_energy_latex;
        mean_adc_sum_vs_energy_latex.SetTextColor(kGray+2);
        mean_adc_sum_vs_energy_latex.SetNDC();
        mean_adc_sum_vs_energy_latex.SetTextSize(0.05);
        mean_adc_sum_vs_energy_latex.SetTextFont(62);
        mean_adc_sum_vs_energy_latex.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());
        mean_adc_sum_vs_energy_latex.SetTextSize(0.035);
        mean_adc_sum_vs_energy_latex.SetTextFont(42);
        mean_adc_sum_vs_energy_latex.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
        mean_adc_sum_vs_energy_latex.DrawLatex(0.12, 0.76, ("Beam: " + beam_type).c_str());
        mean_adc_sum_vs_energy_latex.DrawLatex(0.12, 0.72, date_stream.str().c_str());
        canvas_mean_adc_sum_vs_energy->Write();
        canvas_mean_adc_sum_vs_energy->Close();
    }

    if (enable_gaussian_fit && found_run_energies) {
        canvas_resolution_vs_energy->cd();
        TLatex resolution_vs_energy_latex;
        resolution_vs_energy_latex.SetTextColor(kGray+2);
        resolution_vs_energy_latex.SetNDC();
        resolution_vs_energy_latex.SetTextSize(0.05);
        resolution_vs_energy_latex.SetTextFont(62);
        resolution_vs_energy_latex.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());
        resolution_vs_energy_latex.SetTextSize(0.035);
        resolution_vs_energy_latex.SetTextFont(42);
        resolution_vs_energy_latex.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
        resolution_vs_energy_latex.DrawLatex(0.12, 0.76, ("Beam: " + beam_type).c_str());
        resolution_vs_energy_latex.DrawLatex(0.12, 0.72, date_stream.str().c_str());
        canvas_resolution_vs_energy->Write();
        canvas_resolution_vs_energy->Close();
    }

    // TCanvas *pedestal_distribution_th2d_canvas = new TCanvas("pedestal_distribution_canvas", "pedestal_distribution_canvas", 1200, 600);
    // pedestal_distribution_th2d->Draw("COLZ");
    // pedestal_distribution_th2d_canvas->SetLogz();

    // TLatex latex_pedestal;
    // latex_pedestal.SetTextColor(kGray+2);
    // latex_pedestal.SetNDC();
    // latex_pedestal.SetTextSize(0.05);
    // latex_pedestal.SetTextFont(62);
    // latex_pedestal.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());
    // latex_pedestal.SetTextSize(0.035);
    // latex_pedestal.SetTextFont(42);
    // latex_pedestal.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
    // // write run number, date time
    // auto t = std::time(nullptr);
    // auto tm = *std::localtime(&t);
    // std::ostringstream date_stream;
    // date_stream << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    // latex_pedestal.DrawLatex(0.12, 0.76, (std::string("Run ") + script_input_run_number).c_str());
    // latex_pedestal.DrawLatex(0.12, 0.72, "Pedestal Distribution");
    // latex_pedestal.DrawLatex(0.12, 0.68, date_stream.str().c_str());    
    // pedestal_distribution_th2d_canvas->Write();
    // pedestal_distribution_th2d_canvas->Close();

    // TCanvas *peak_distribution_th2d_canvas = new TCanvas("peak_distribution_canvas", "peak_distribution_canvas", 1200, 600);
    // peak_distribution_th2d->Draw("COLZ");
    // peak_distribution_th2d_canvas->SetLogz();
    // TLatex latex_peak;
    // latex_peak.SetTextColor(kGray+2);
    // latex_peak.SetNDC();
    // latex_peak.SetTextSize(0.05);
    // latex_peak.SetTextFont(62);
    // latex_peak.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());
    // latex_peak.SetTextSize(0.035);
    // latex_peak.SetTextFont(42);
    // latex_peak.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
    // // write run number, date time
    // latex_peak.DrawLatex(0.12, 0.76, (std::string("Run ") + script_input_run_number).c_str());
    // latex_peak.DrawLatex(0.12, 0.72, "Peak Distribution");
    // latex_peak.DrawLatex(0.12, 0.68, date_stream.str().c_str());    
    // peak_distribution_th2d_canvas->Write();
    // peak_distribution_th2d_canvas->Close();

    // TCanvas *run_adc_sum_canvas = new TCanvas("run_adc_sum_all_canvas", "run_adc_sum_all_canvas", 800, 600);
    // // draw all the partial histograms in one canvas
    // run_adc_sum_all_hist1d->SetLineColor(kBlack);
    // run_adc_sum_all_hist1d->SetLineWidth(2);
    // run_adc_sum_all_hist1d->Draw("HIST");
    // TLegend *run_adc_sum_legend = new TLegend(0.80, 0.65, 0.89, 0.89);
    // run_adc_sum_legend->SetBorderSize(0);
    // run_adc_sum_legend->SetFillStyle(0);
    // run_adc_sum_legend->SetTextFont(42);
    // run_adc_sum_legend->SetTextSize(0.02);
    // run_adc_sum_legend->AddEntry(run_adc_sum_all_hist1d, "All Events", "l");

    // for (int i = 0; i < n_partial_hist; i++) {
    //     auto hist = run_adc_sum_partial_hist1d_list[i];
    //     hist->SetLineColor(kBlue + i * 2);
    //     hist->SetLineWidth(1);
    //     hist->Draw("HIST SAME");
    //     run_adc_sum_legend->AddEntry(hist, ("Partial " + std::to_string(i)).c_str(), "l");
    // }
    // run_adc_sum_legend->Draw();
    // TLatex latex_fit;
    // latex_fit.SetTextColor(kRed);
    // latex_fit.SetNDC();
    // latex_fit.SetTextSize(0.03);
    // latex_fit.SetTextFont(42);
    // latex_fit.DrawLatex(0.60, 0.60, ("Fit Mean: " + std::to_string(final_mean)).c_str());
    // latex_fit.DrawLatex(
    //     0.60, 0.55, 
    //     ("Fit Sigma: " + std::to_string(final_sigma) + 
    //     " (" + std::to_string(static_cast<int>(final_sigma / final_mean * 10000) / 100.0) + "%)").c_str()
    // );
    // // write run number
    // TLatex latex_run_adc_sum;
    // latex_run_adc_sum.SetTextColor(kGray+2);
    // latex_run_adc_sum.SetNDC();
    // latex_run_adc_sum.SetTextSize(0.05);
    // latex_run_adc_sum.SetTextFont(62);
    // latex_run_adc_sum.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());
    // latex_run_adc_sum.SetTextSize(0.035);
    // latex_run_adc_sum.SetTextFont(42);
    // latex_run_adc_sum.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
    // // write run number, date time
    // latex_run_adc_sum.DrawLatex(0.12, 0.76, (std::string("Run ") + script_input_run_number).c_str());
    // latex_run_adc_sum.DrawLatex(0.12, 0.72, "Run ADC Sum");
    // latex_run_adc_sum.DrawLatex(0.12, 0.68, date_stream.str().c_str());    
    // run_adc_sum_canvas->Write();
    // run_adc_sum_canvas->Close();

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}