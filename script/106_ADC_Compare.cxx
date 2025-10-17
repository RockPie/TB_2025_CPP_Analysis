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
    // if the target_event_number is -1, go through all events and normalize by number of events
    if (target_event_number == -1) {
        spdlog::info("No target number of events specified, using all events");
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

        auto target_event_count = min_run_file_event_count;
        if (script_n_events > 0 && script_n_events < min_run_file_event_count
            && script_n_events < run_file_event_counts[_file_index]) {
            target_event_count = script_n_events;
            spdlog::info("Using target number of events from command line: {}", target_event_count);
        }
        if (script_n_events > 0 && script_n_events >= min_run_file_event_count) {
            target_event_count = input_tree->GetEntries();
        }
        for (int _entry = 0; _entry < target_event_count; _entry++) {
            input_tree->GetEntry(_entry);
            // process data here
            Double_t adc_sum = 0.0;
            if (use_peak_integral) {
                for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
                    for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                        double adc_integral = 0;
                        double adc_pedestal = 0;
                        std::vector<int> adc_pedestal_samples;
                        for (int sample = 0; sample < machine_gun_samples; sample++) {
                            auto adc_value = val0_list_pools[vldb_id][0][FPGA_CHANNEL_NUMBER * sample + channel];
                            if (sample < 3) {
                                adc_pedestal_samples.push_back(adc_value);
                            }
                            adc_integral += double(adc_value);
                        } // end of sample loop
                        // calculate pedestal as average of the smallest 2 samples among first 3 samples
                        auto pedestal = pedestal_median_of_first3(adc_pedestal_samples);
                        adc_integral -= pedestal * machine_gun_samples;
                        if (adc_integral <= 0) {
                            continue;
                        }
                        // spdlog::debug("Event {}: VLDB {} Channel {}: Pedestal = {}, Integral = {}", _entry, vldb_id, channel, pedestal, adc_integral);
                        adc_sum += Double_t(adc_integral);
                    }
                }
            } else {
                for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
                    for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                        double adc_max = 0;
                        double adc_pedestal = 0;
                        std::vector<int> adc_pedestal_samples;
                        for (int sample = 0; sample < machine_gun_samples; sample++) {
                            auto adc_value = val0_list_pools[vldb_id][0][FPGA_CHANNEL_NUMBER * sample + channel];
                            if (sample < 3) {
                                adc_pedestal_samples.push_back(adc_value);
                            }
                            auto adc_value_double = double(adc_value);
                            if (adc_value_double > adc_max) {
                                adc_max = adc_value_double;
                            }
                        } // end of sample loop
                        // calculate pedestal as average of the smallest 2 samples among first 3 samples
                        auto pedestal = pedestal_median_of_first3(adc_pedestal_samples);
                        if (adc_max <= adc_pedestal) {
                            continue;
                        }
                        
                        // spdlog::debug("Event {}: VLDB {} Channel {}: Pedestal = {}, Peak = {}", _entry, vldb_id, channel, adc_pedestal, adc_max);
                        adc_sum += Double_t(adc_max) - Double_t(pedestal);
                    }
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
        
        // if script_n_events is -1, normalize by number of events
        if (target_event_number == -1) {
            run_adc_sum_hist1d->Scale(1.0 / adc_sums->size());
        }
        if (run_adc_sum_hist1d->GetMaximum() > hist_max_bin_value) {
            hist_max_bin_value = run_adc_sum_hist1d->GetMaximum();
        }
    }

    std::vector<std::vector<TF1*>> run_adc_sum_hist1d_gauss_fit_list;
    std::vector<TF1*> run_adc_sum_hist1d_pre_fit_list;

    std::vector<double> gauss_fit_sigma_ranges_left  = {0.15, 0.20, 0.25};
    std::vector<double> low_energy_gauss_fit_sigma_ranges_left  = {1.0, 1.5, 2.0};
    std::vector<double> gauss_fit_sigma_ranges_right = {2.0, 2.5, 3.0};
    if (enable_gaussian_fit) {
        for (int _hist_index = 0; _hist_index < run_adc_sum_hist1d_list.size(); _hist_index++) {
            auto* hist = run_adc_sum_hist1d_list[_hist_index];
            // pre-fit to find the peak position
            double hist_mean = hist->GetMean();
            double hist_rms  = hist->GetRMS();
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
            // choose different fit ranges for low energy runs
            auto _fit_left_sigma_ranges = gauss_fit_sigma_ranges_left;
            if (found_run_energies && run_energies[_hist_index] <= 100) {
                _fit_left_sigma_ranges = low_energy_gauss_fit_sigma_ranges_left;
            }
            // gaussian fit around the peak
            for (auto& _fit_left_sigma : _fit_left_sigma_ranges) {
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

    // * If gaussian fit enabled, and beam energies provided, make a graph of mean ADC sum vs beam energy
    if (enable_gaussian_fit && found_run_energies) {
        const double xmin = 40;
        const double xmax = 400;
        const Int_t kFont = 42;
        const float labelSize = 0.045;   // 刻度文字
        const float titleSize = 0.050;   // 轴标题
        const float yTitleOffsetTop = 1.15;
        const float yTitleOffsetBot = 1.10;
        const float xTitleOffsetBot = 1.10;
        const float bottomZoomRatio = 0.5;
        TGraphErrors* graph_mean_adc_sum_vs_energy = new TGraphErrors();
        graph_mean_adc_sum_vs_energy->SetName("graph_mean_adc_sum_vs_energy");
        graph_mean_adc_sum_vs_energy->SetTitle("");
        for (int _file_index = 0; _file_index < run_adc_sum_hist1d_list.size(); _file_index++) {
            auto* hist = run_adc_sum_hist1d_list[_file_index];
            // use all the fitting results to calculate the average mean ADC sum
            double mean_adc_sum = 0;
            double mean_adc_sum_err = 0;
            std::vector<double> mean_adc_sums;
            int n_fits = run_adc_sum_hist1d_gauss_fit_list[_file_index].size();
            for (auto* gauss_fit : run_adc_sum_hist1d_gauss_fit_list[_file_index]) {
                mean_adc_sum += gauss_fit->GetParameter(1);
                mean_adc_sums.push_back(gauss_fit->GetParameter(1));
            }
            if (n_fits > 0) {
                mean_adc_sum /= n_fits;
                // calculate standard deviation of the mean_adc_sums
                double sum_sq_diff = 0;
                for (const auto& mas : mean_adc_sums) {
                    sum_sq_diff += (mas - mean_adc_sum) * (mas - mean_adc_sum);
                }
                mean_adc_sum_err = std::sqrt(sum_sq_diff / n_fits);
            } else {
                mean_adc_sum = 0;
                mean_adc_sum_err = 0;
            }
            // assume 5% error on beam energy
            graph_mean_adc_sum_vs_energy->SetPoint(_file_index, run_energies[_file_index], mean_adc_sum);
            graph_mean_adc_sum_vs_energy->SetPointError(_file_index, run_energies[_file_index] * 0.05, mean_adc_sum_err);
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
        if (run_energies.size() >= 2) {
            graph_mean_adc_sum_vs_energy->GetXaxis()->SetRangeUser(40, 400);
        }
        
        linear_fit->SetLineColor(kRed);
        linear_fit->SetLineWidth(2);

        // create two plots in the same canvas
        // ratio 9:1
        canvas_mean_adc_sum_vs_energy->Divide(1, 2, 0, 0);
        auto* pad_mean = canvas_mean_adc_sum_vs_energy->cd(1);
        pad_mean->SetPad(0.0, 0.35, 1.0, 1.0);
        pad_mean->SetTopMargin(0.08);
        pad_mean->SetBottomMargin(0);
        pad_mean->SetLeftMargin(0.12);
        pad_mean->SetRightMargin(0.04);
        pad_mean->SetBorderMode(0);
        pad_mean->SetFrameBorderMode(0);
        pad_mean->SetTicks(1,1);

        auto* pad_dev = canvas_mean_adc_sum_vs_energy->cd(2);
        pad_dev->SetPad(0.0, 0.0, 1.0, 0.35);
        pad_dev->SetTopMargin(0);
        pad_dev->SetBottomMargin(0.35);
        pad_dev->SetLeftMargin(0.12);
        pad_dev->SetRightMargin(0.04);
        pad_dev->SetBorderMode(0);
        pad_dev->SetFrameBorderMode(0);
        pad_dev->SetTicks(1,1);

        pad_mean->cd();

        graph_mean_adc_sum_vs_energy->SetMarkerStyle(20);
        graph_mean_adc_sum_vs_energy->SetMarkerSize(1.0);
        graph_mean_adc_sum_vs_energy->SetLineWidth(2);
        graph_mean_adc_sum_vs_energy->SetLineColor(kBlue);

        graph_mean_adc_sum_vs_energy->Draw("AP");

        graph_mean_adc_sum_vs_energy->GetYaxis()->SetTitle("Mean ADC sum");
        graph_mean_adc_sum_vs_energy->GetYaxis()->SetTitleFont(kFont);
        graph_mean_adc_sum_vs_energy->GetYaxis()->SetLabelFont(kFont);
        graph_mean_adc_sum_vs_energy->GetYaxis()->SetTitleSize(titleSize);
        graph_mean_adc_sum_vs_energy->GetYaxis()->SetLabelSize(labelSize);
        graph_mean_adc_sum_vs_energy->GetYaxis()->SetTitleOffset(yTitleOffsetTop);

        graph_mean_adc_sum_vs_energy->GetXaxis()->SetLimits(xmin, xmax);
        graph_mean_adc_sum_vs_energy->GetXaxis()->SetRangeUser(xmin, xmax);
        graph_mean_adc_sum_vs_energy->GetXaxis()->SetLabelSize(0);   // 隐藏上图 x 刻度文字
        graph_mean_adc_sum_vs_energy->GetXaxis()->SetTitleSize(0);   // 隐藏上图 x 轴标题

        linear_fit->SetRange(xmin, xmax);
        linear_fit->SetLineColor(kRed);
        linear_fit->SetLineWidth(2);
        linear_fit->Draw("SAME");

        TLegend* leg = new TLegend(0.45, 0.12, 0.88, 0.32, "", "brNDC");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(kFont);
        leg->SetTextSize(0.040);
        leg->AddEntry(linear_fit,
                    fmt::format("Fit: offset = {:.1f}, slope = {:.1f}",
                                linear_fit->GetParameter(0), linear_fit->GetParameter(1)).c_str(),
                    "l");
        leg->Draw();

        
        
        pad_dev->cd();

        // draw an empty graph to set axis ranges
        TGraph* empty_graph = new TGraph();
        empty_graph->SetPoint(0, xmin, -5);
        empty_graph->SetPoint(1, xmax, 5);
        empty_graph->GetXaxis()->SetLimits(xmin, xmax);
        empty_graph->GetXaxis()->SetRangeUser(xmin, xmax);
        empty_graph->GetYaxis()->SetRangeUser(-5, 5);
        empty_graph->Draw("AP");
        empty_graph->GetXaxis()->SetTitle("Beam energy [GeV]");
        empty_graph->GetYaxis()->SetTitle("Deviation [%]");
        // show less y ticks
        empty_graph->GetYaxis()->SetNdivisions(505);
        empty_graph->GetXaxis()->SetTitleFont(kFont);
        empty_graph->GetXaxis()->SetLabelFont(kFont);
        empty_graph->GetXaxis()->SetTitleSize(titleSize/bottomZoomRatio);
        empty_graph->GetXaxis()->SetLabelSize(labelSize/bottomZoomRatio);
        empty_graph->GetXaxis()->SetTitleOffset(xTitleOffsetBot);
        empty_graph->GetYaxis()->SetTitleFont(kFont);
        empty_graph->GetYaxis()->SetLabelFont(kFont);
        empty_graph->GetYaxis()->SetTitleSize(titleSize/bottomZoomRatio);
        empty_graph->GetYaxis()->SetLabelSize(labelSize/bottomZoomRatio);
        // empty_graph->GetYaxis()->SetTitleOffset(yTitleOffsetBot);


        TGraph* nonlinear_deviation = new TGraph();
        nonlinear_deviation->SetName("nonlinear_deviation");
        nonlinear_deviation->SetTitle("");
        for (int i = 0; i < graph_mean_adc_sum_vs_energy->GetN(); i++) {
            double energy, mean_adc_sum;
            graph_mean_adc_sum_vs_energy->GetPoint(i, energy, mean_adc_sum);
            double mean_adc_sum_fit = linear_fit->Eval(energy);
            double deviation = (mean_adc_sum - mean_adc_sum_fit) / mean_adc_sum_fit * 100;
            nonlinear_deviation->SetPoint(i, energy, deviation);
        }
        nonlinear_deviation->SetMarkerStyle(21);
        nonlinear_deviation->SetMarkerSize(1);
        nonlinear_deviation->SetLineColor(kBlue);
        nonlinear_deviation->SetLineWidth(2);
        nonlinear_deviation->GetYaxis()->SetRangeUser(-5, 5);
        // set the x range if beam energies provided
        if (found_run_energies) {
            nonlinear_deviation->GetXaxis()->SetRangeUser(xmin, xmax);
            nonlinear_deviation->GetXaxis()->SetLimits(xmin, xmax);
        }
        nonlinear_deviation->SetMarkerStyle(21);
        nonlinear_deviation->SetMarkerSize(1.0);
        nonlinear_deviation->SetLineWidth(2);
        nonlinear_deviation->SetLineColor(kBlue);

        nonlinear_deviation->Draw("P");

        // nonlinear_deviation->GetXaxis()->SetTitle("Beam energy [GeV]");
        // nonlinear_deviation->GetYaxis()->SetTitle("Deviation from linear [%]");

        // nonlinear_deviation->GetXaxis()->SetLimits(xmin, xmax);
        // nonlinear_deviation->GetXaxis()->SetTitleFont(kFont);
        // nonlinear_deviation->GetXaxis()->SetLabelFont(kFont);
        // nonlinear_deviation->GetXaxis()->SetTitleSize(titleSize/bottomZoomRatio);
        // nonlinear_deviation->GetXaxis()->SetLabelSize(labelSize/bottomZoomRatio);
        // nonlinear_deviation->GetXaxis()->SetTitleOffset(xTitleOffsetBot);

        // nonlinear_deviation->GetYaxis()->SetTitleFont(kFont);
        // nonlinear_deviation->GetYaxis()->SetLabelFont(kFont);
        // nonlinear_deviation->GetYaxis()->SetTitleSize(titleSize/bottomZoomRatio);
        // nonlinear_deviation->GetYaxis()->SetLabelSize(labelSize/bottomZoomRatio);
        // nonlinear_deviation->GetYaxis()->SetTitleOffset(yTitleOffsetBot);
        // nonlinear_deviation->GetYaxis()->SetRangeUser(-5, 5);

        // draw a horizontal line at y=0
        TLine* line_y0 = new TLine(xmin, 0, xmax, 0);
        line_y0->SetLineStyle(2);
        line_y0->SetLineColor(kBlack);
        line_y0->Draw("SAME");
    }

    // * Draw resolution vs beam energy if gaussian fit enabled and beam energies provided
    TCanvas* canvas_resolution_vs_energy = new TCanvas("canvas_resolution_vs_energy", "Canvas Resolution vs Beam Energy", 800, 600);
    TLegend *legend_resolution_vs_energy = new TLegend(0.12, 0.12, 0.4, 0.3);
    legend_resolution_vs_energy->SetBorderSize(0);
    legend_resolution_vs_energy->SetTextSize(0.03);
    if (enable_gaussian_fit && found_run_energies) {
        TGraphErrors* graph_resolution_vs_energy = new TGraphErrors();
        graph_resolution_vs_energy->SetName("graph_resolution_vs_energy");
        graph_resolution_vs_energy->SetTitle("");
        for (int _file_index = 0; _file_index < run_adc_sum_hist1d_list.size(); _file_index++) {
            auto* hist = run_adc_sum_hist1d_list[_file_index];
            // use all the fitting results to calculate the average resolution
            double resolution = 0;
            double resolution_err = 0;
            std::vector<double> resolutions;

            int n_fits = run_adc_sum_hist1d_gauss_fit_list[_file_index].size();
            for (auto* gauss_fit : run_adc_sum_hist1d_gauss_fit_list[_file_index]) {
                double sigma = gauss_fit->GetParameter(2);
                double mean = gauss_fit->GetParameter(1);
                if (mean != 0) {
                    resolution += sigma / mean;
                }
                resolutions.push_back(sigma / mean);
            }
            if (n_fits > 0) {
                resolution /= n_fits;
                // calculate standard deviation of the resolutions
                double sum_sq_diff = 0;
                for (const auto& res : resolutions) {
                    sum_sq_diff += (res - resolution) * (res - resolution);
                }
                resolution_err = std::sqrt(sum_sq_diff / n_fits);
            } else {
                resolution = 0;
                resolution_err = 0;
            }
            // assume 5% error on beam energy
            graph_resolution_vs_energy->SetPoint(_file_index, run_energies[_file_index], resolution);
            graph_resolution_vs_energy->SetPointError(_file_index, run_energies[_file_index] * 0.05, resolution_err);
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
    float text_start_y = 0.85;
    float text_start_x = 0.12;
    latex_adc_sum.DrawLatex(text_start_x, text_start_y, annotation_canvas_title.c_str());
    latex_adc_sum.SetTextSize(0.035);
    latex_adc_sum.SetTextFont(42);
    latex_adc_sum.DrawLatex(text_start_x, text_start_y - 0.05, annotation_testbeam_title.c_str());
    // write run number, date time
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream date_stream;
    date_stream << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    latex_adc_sum.DrawLatex(text_start_x, text_start_y - 0.09, config_description.c_str());
    // hardon or electron beam
    latex_adc_sum.DrawLatex(text_start_x, text_start_y - 0.13, ("Beam: " + beam_type).c_str());
    latex_adc_sum.DrawLatex(text_start_x, text_start_y - 0.17, date_stream.str().c_str());
    canvas_adc_sum->Write();
    canvas_adc_sum->Close();

    if (enable_gaussian_fit && found_run_energies) {
        canvas_mean_adc_sum_vs_energy->cd();
        text_start_x = 0.13;
        text_start_y = 0.89;
        TLatex mean_adc_sum_vs_energy_latex;
        mean_adc_sum_vs_energy_latex.SetTextColor(kGray+2);
        mean_adc_sum_vs_energy_latex.SetNDC();
        mean_adc_sum_vs_energy_latex.SetTextSize(0.05);
        mean_adc_sum_vs_energy_latex.SetTextFont(62);
        mean_adc_sum_vs_energy_latex.DrawLatex(text_start_x, text_start_y, annotation_canvas_title.c_str());
        mean_adc_sum_vs_energy_latex.SetTextSize(0.035);
        mean_adc_sum_vs_energy_latex.SetTextFont(42);
        mean_adc_sum_vs_energy_latex.DrawLatex(text_start_x, text_start_y - 0.05, annotation_testbeam_title.c_str());
        mean_adc_sum_vs_energy_latex.DrawLatex(text_start_x, text_start_y - 0.09, config_description.c_str());
        mean_adc_sum_vs_energy_latex.DrawLatex(text_start_x, text_start_y - 0.13, ("Beam: " + beam_type).c_str());
        mean_adc_sum_vs_energy_latex.DrawLatex(text_start_x, text_start_y - 0.17, date_stream.str().c_str());

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
        resolution_vs_energy_latex.DrawLatex(0.12, 0.76, config_description.c_str());
        resolution_vs_energy_latex.DrawLatex(0.12, 0.72, ("Beam: " + beam_type).c_str());
        resolution_vs_energy_latex.DrawLatex(0.12, 0.68, date_stream.str().c_str());
        canvas_resolution_vs_energy->Write();
        canvas_resolution_vs_energy->Close();
    }

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}