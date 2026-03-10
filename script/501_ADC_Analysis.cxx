#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "TRandom3.h"
#include "CommonParams.hxx"

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

    const int example_channel = CommonParams::example_channel;
    const int CM_channel = 9; // this is for each half of the chip, not the global channel number

        cxxopts::Options options(script_name, "Generate heatmaps from machine gun data");

    options.add_options()
        ("f,file", "Input .root file", cxxopts::value<std::string>())
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

    // dump/102_EventMatch/beamtests/RunXXXX.root
    script_input_file  = parsed["file"].as<std::string>();
    std::string script_input_run_number = script_input_file.substr(script_input_file.find_last_of("Run") + 1).append(4, '0').substr(0, 4);
    script_output_file = parsed["output"].as<std::string>();
    script_n_events    = parsed["events"].as<int>();
    script_verbose     = parsed["verbose"].as<bool>();

    configure_logger(script_verbose);

    spdlog::info("Input file: {}", script_input_file);

    if (access(script_input_file.c_str(), F_OK) == -1) {
        spdlog::error("Input file {} does not exist!", script_input_file);
        return 1;
    }
    if (script_input_file.substr(script_input_file.find_last_of(".") + 1) != "root") {
        spdlog::error("Input file {} should end with .root!", script_input_file);
        return 1;
    }
    if (script_output_file.substr(script_output_file.find_last_of(".") + 1) != "root") {
        spdlog::error("Output file {} should end with .root!", script_output_file);
        return 1;
    }

    std::string mapping_json_file = "config/mapping_Feb2026.json";
    std::ifstream mapping_json_ifs(mapping_json_file);
    if (!mapping_json_ifs.is_open()) {
        spdlog::error("Failed to open mapping json file {}", mapping_json_file);
        return 1;
    }
    json mapping_json;
    mapping_json_ifs >> mapping_json;
    const auto& sipm_board      = mapping_json.at("SiPM_Board");
    const auto& board_loc       = mapping_json.at("Board_Loc");
    const auto& board_rotation  = mapping_json.at("Board_Rotation");
    const auto& board_flip      = mapping_json.at("Board_Flip");

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
    spdlog::info("Input file: {}", script_input_file);
    spdlog::info("Output file: {} in {}", script_output_file, script_output_folder);
    spdlog::info("Number of events: {}", script_n_events);

    // * --- Read input file ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    const int machine_gun_samples = 16;
    const int vldb_number = 2;

    TFile *input_root = new TFile(script_input_file.c_str(), "READ");
    if (input_root->IsZombie()) {
        spdlog::error("Failed to open input file {}", script_input_file);
        return 1;
    }
    TTree *input_tree = (TTree*) input_root->Get("data_tree");
    if (input_tree == nullptr) {
        spdlog::error("Failed to get data tree from input file {}", script_input_file);
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

    int entry_max = input_tree->GetEntries();

    if (script_n_events > 0 && script_n_events < entry_max) {
        entry_max = script_n_events;
    }

    spdlog::info("Processing {} events", entry_max);

    if (entry_max < 1) {
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

    const int adc_peak_min_index = CommonParams::adc_peak_min_index;
    const int adc_peak_max_index = CommonParams::adc_peak_max_index;

    std::vector<double> adc_peak_sums_per_event;
    adc_peak_sums_per_event.reserve(entry_max);

    std::vector<std::vector<std::vector<double>>> adc_samples_in_CM_channels_all_vldb(
        vldb_number, std::vector<std::vector<double>>(4, std::vector<double>(machine_gun_samples, 0.0)));
    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        double adc_peak_sum = 0.0;
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            std::vector<std::vector<double>>& adc_samples_in_CM_channels = adc_samples_in_CM_channels_all_vldb[vldb_id];
            for (int half = 0; half < 4; half++) {
                int channel = half * (FPGA_CHANNEL_NUMBER / 4) + CM_channel;
                std::vector<UInt_t> adc_pedestal_samples; // only take the first 3 samples
                std::vector<UInt_t> adc_samples; // all samples for common-mode calculation
                adc_pedestal_samples.reserve(3);
                adc_samples.reserve(machine_gun_samples);
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    UInt_t adc_value = val0_list_pools[vldb_id][0][idx];
                    if (sample < 3) {
                        adc_pedestal_samples.push_back(adc_value);
                    }
                    adc_samples.push_back(adc_value);
                } // end of sample loop

                double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    adc_samples_in_CM_channels[half][sample] = static_cast<double>(adc_samples[sample]) - adc_pedestal;
                }
            } // end of half loop for CM
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int valid_channel_number = get_valid_fpga_channel(channel);
                if (valid_channel_number < 0) {
                    continue; // skip invalid channels
                }
                int channel_half_index = (channel / (FPGA_CHANNEL_NUMBER / 4)) % 4;
                if ((channel_half_index == 2 || channel_half_index == 3) && vldb_id == 0) {
                    continue; // skip CM channels
                }
                std::vector<double> CM_samples = adc_samples_in_CM_channels[channel_half_index];
                std::vector<double> adc_pedestal_samples; // only take the first 3 samples
                std::vector<double> adc_samples; // all samples for common-mode calculation
                int adc_peak_index = 0;
                double adc_peak_value = 0.0;

                adc_pedestal_samples.reserve(3);
                adc_samples.reserve(machine_gun_samples);
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    double adc_value = static_cast<double>(val0_list_pools[vldb_id][0][idx]);
                    UInt_t tot_value = val1_list_pools[vldb_id][0][idx];
                    UInt_t toa_value = val2_list_pools[vldb_id][0][idx];

                    adc_value -= CM_samples[sample]; // common-mode subtraction
                    if (sample < 3) {
                        adc_pedestal_samples.push_back(adc_value);
                    }
                    adc_samples.push_back(adc_value);
                    if (adc_value > adc_peak_value && sample >= adc_peak_min_index && sample <= adc_peak_max_index) {
                        adc_peak_value = adc_value;
                        adc_peak_index = sample;
                    }
                } // end of sample loop
                // subtract the pedestal from the peak value
                double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
                double adc_peak_value_pede_sub = adc_peak_value - adc_pedestal;
                if (adc_peak_value_pede_sub < 0) {
                    adc_peak_value_pede_sub = 0.0; // set negative values to zero
                }
                adc_peak_sum += adc_peak_value_pede_sub;
            } // end of channel loop
        } // end of vldb loop
        adc_peak_sums_per_event.push_back(adc_peak_sum);
    } // end of event loop

    input_root->Close();

    const int NX = 16, NY = 12;
    const int board_cols = 8, board_rows = 4;

    auto chan2pad = build_chan2pad_LUT(
        vldb_number, FPGA_CHANNEL_NUMBER,
        NX, NY, board_cols, board_rows,
        sipm_board, board_loc, board_rotation, board_flip
    );

    const std::vector<double> fit_range_sigmas = {2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5}; // in sigma
    const std::vector<double> fit_range_offsets = {-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5}; // in sigma

    MosaicTopology topo_wave;
    topo_wave.NX = NX;
    topo_wave.NY = NY;
    topo_wave.vldb_number = vldb_number;
    topo_wave.channels_per_vldb = FPGA_CHANNEL_NUMBER;
    topo_wave.reverse_row = true;
    topo_wave.minimalist_axis = true;
    topo_wave.th2_logz = true;
    topo_wave.chan2pad = chan2pad;

    MosaicTopology topo_ped_median = topo_wave;
    topo_ped_median = topo_wave;
    topo_ped_median.th2_logz = true;

    auto output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    output_root->cd();

    // create the histogram for ADC peak sum per event
    const double adc_peak_sum_hist_min = 0.0;
    const double adc_peak_sum_hist_bin_size = 200.0; // adjust the bin size as needed
    // define the max as the 95% percentile of the adc_peak_sums_per_event to avoid outliers dominating the histogram
    std::vector<double> adc_peak_sums_sorted = adc_peak_sums_per_event;
    std::sort(adc_peak_sums_sorted.begin(), adc_peak_sums_sorted.end());
    double adc_peak_sum_hist_max = adc_peak_sums_sorted[static_cast<int>(adc_peak_sums_sorted.size() * 0.95)] * 1.5;
    auto canvas_adc_peak_sum = new TCanvas("canvas", "ADC Peak Sum per Event", 800, 600);
    auto h1d_adc_peak_sum = new TH1D("adc_peak_sum_hist", "ADC Peak Sum per Event;ADC Peak Sum;Number of Events", static_cast<int>((adc_peak_sum_hist_max - adc_peak_sum_hist_min) / adc_peak_sum_hist_bin_size), adc_peak_sum_hist_min, adc_peak_sum_hist_max);

    auto legend_adc_peak_sum = new TLegend(0.62, 0.4, 0.9, 0.9);
    legend_adc_peak_sum->SetBorderSize(0);
    legend_adc_peak_sum->SetFillStyle(0);
    legend_adc_peak_sum->SetTextSize(0.025);

    for (const auto& adc_peak_sum : adc_peak_sums_per_event) {
        if (adc_peak_sum <= 0.0) {
            continue; // skip non-positive values
        }
        h1d_adc_peak_sum->Fill(adc_peak_sum);
    }
    std::string canvas_info_sub = "ADC Peak Sum, Run " + script_input_run_number;
    format_1d_hist_canvas(canvas_adc_peak_sum, h1d_adc_peak_sum, kBlue, annotation_canvas_title, annotation_testbeam_title, canvas_info_sub, false);
    legend_adc_peak_sum->AddEntry(h1d_adc_peak_sum, "Data", "l");

    // do the fitting
    int gaus_fit_color = kCyan+2;
    double gaus_fit_mean_nom, gaus_fit_mean_err_stat, gaus_fit_mean_err_syst;
    double gaus_fit_sigma_nom, gaus_fit_sigma_err_stat, gaus_fit_sigma_err_syst;
    TF1* nominal_fit_func = nullptr;
    gaussian_fit_th1d_nom(*canvas_adc_peak_sum, *h1d_adc_peak_sum, fit_range_sigmas, fit_range_offsets, gaus_fit_color, gaus_fit_mean_nom, gaus_fit_sigma_nom, gaus_fit_mean_err_stat, gaus_fit_sigma_err_stat, gaus_fit_mean_err_syst, gaus_fit_sigma_err_syst, nominal_fit_func);
    legend_adc_peak_sum->AddEntry(nominal_fit_func, "Gaussian Fit", "l");
    // add fit results to the legend
    auto format_sig3 = [](double value) {
        std::ostringstream oss;
        oss << std::setprecision(3) << std::defaultfloat << value;
        return oss.str();
    };
    auto format_sig5 = [](double value) {
        std::ostringstream oss;
        oss << std::setprecision(5) << std::defaultfloat << value;
        return oss.str();
    };
    // std::string gaus_fit_info = "#mu = " + format_sig3(gaus_fit_mean_nom) + " #pm " + format_sig3(gaus_fit_mean_err_stat) + " (stat) #pm " + format_sig3(gaus_fit_mean_err_syst) + " (syst)";
    // std::string gaus_fit_info_sigma = "#sigma = " + format_sig3(gaus_fit_sigma_nom) + " #pm " + format_sig3(gaus_fit_sigma_err_stat) + " (stat) #pm " + format_sig3(gaus_fit_sigma_err_syst) + " (syst)";
    // legend_adc_peak_sum->AddEntry((TObject*)nullptr, gaus_fit_info.c_str(), "");
    // legend_adc_peak_sum->AddEntry((TObject*)nullptr, gaus_fit_info_sigma.c_str(), "");
    std::string gaus_fit_mean_value_info = "#mu = " + format_sig5(gaus_fit_mean_nom);
    std::string gaus_fit_sigma_value_info = "#sigma = " + format_sig5(gaus_fit_sigma_nom);
    std::string gaus_fit_mean_err_info = "#pm " + format_sig3(gaus_fit_mean_err_stat) + " (stat) #pm " + format_sig3(gaus_fit_mean_err_syst) + " (syst)";
    std::string gaus_fit_sigma_err_info = "#pm " + format_sig3(gaus_fit_sigma_err_stat) + " (stat) #pm " + format_sig3(gaus_fit_sigma_err_syst) + " (syst)";
    legend_adc_peak_sum->AddEntry((TObject*)nullptr, gaus_fit_mean_value_info.c_str(), "");
    legend_adc_peak_sum->AddEntry((TObject*)nullptr, gaus_fit_mean_err_info.c_str(), "");
    legend_adc_peak_sum->AddEntry((TObject*)nullptr, gaus_fit_sigma_value_info.c_str(), "");
    legend_adc_peak_sum->AddEntry((TObject*)nullptr, gaus_fit_sigma_err_info.c_str(), "");

    // do crystal ball fit
    int cb_fit_color = kPink+2;
    double cb_fit_mean_nom, cb_fit_sigma_nom;
    double cb_fit_mean_err_stat, cb_fit_sigma_err_stat;
    double cb_fit_mean_err_syst, cb_fit_sigma_err_syst;
    TF1* cb_nominal_fit_func = nullptr;
    // inline void crystalball_fit_th1d_nom(TCanvas& canvas, TH1D& hist, const std::vector<double>& fit_range_sigma, const std::vector<double>& fit_range_offset, int color, double& mean_nom, double& sigma_nom, double& mean_err_stat, double& sigma_err_stat, double& mean_err_sys, double& sigma_err_sys)
    crystalball_fit_th1d_nom(*canvas_adc_peak_sum, *h1d_adc_peak_sum, fit_range_sigmas, fit_range_offsets, cb_fit_color, cb_fit_mean_nom, cb_fit_sigma_nom, cb_fit_mean_err_stat, cb_fit_sigma_err_stat, cb_fit_mean_err_syst, cb_fit_sigma_err_syst, cb_nominal_fit_func);
    legend_adc_peak_sum->AddEntry(cb_nominal_fit_func, "Crystal Ball Fit", "l");
    // std::string cb_fit_info = "#mu = " + format_sig3(cb_fit_mean_nom) + " #pm " + format_sig3(cb_fit_mean_err_stat) + " (stat) #pm " + format_sig3(cb_fit_mean_err_syst) + " (syst)";
    // std::string cb_fit_info_sigma = "#sigma = " + format_sig3(cb_fit_sigma_nom) + " #pm " + format_sig3(cb_fit_sigma_err_stat) + " (stat) #pm " + format_sig3(cb_fit_sigma_err_syst) + " (syst)";
    // legend_adc_peak_sum->AddEntry((TObject*)nullptr, cb_fit_info.c_str(), "");
    // legend_adc_peak_sum->AddEntry((TObject*)nullptr, cb_fit_info_sigma.c_str(), "");
    std::string cb_fit_mean_value_info = "#mu = " + format_sig5(cb_fit_mean_nom);
    std::string cb_fit_sigma_value_info = "#sigma = " + format_sig5(cb_fit_sigma_nom);
    std::string cb_fit_mean_err_info = "#pm " + format_sig3(cb_fit_mean_err_stat) + " (stat) #pm " + format_sig3(cb_fit_mean_err_syst) + " (syst)";
    std::string cb_fit_sigma_err_info = "#pm " + format_sig3(cb_fit_sigma_err_stat) + " (stat) #pm " + format_sig3(cb_fit_sigma_err_syst) + " (syst)";
    legend_adc_peak_sum->AddEntry((TObject*)nullptr, cb_fit_mean_value_info.c_str(), "");
    legend_adc_peak_sum->AddEntry((TObject*)nullptr, cb_fit_mean_err_info.c_str(), "");
    legend_adc_peak_sum->AddEntry((TObject*)nullptr, cb_fit_sigma_value_info.c_str(), "");
    legend_adc_peak_sum->AddEntry((TObject*)nullptr, cb_fit_sigma_err_info.c_str(), "");

    double cb_gaus_mix_mean_nom, cb_gaus_mix_mean_err_stat, cb_gaus_mix_mean_err_syst;
    double cb_gaus_mix_sigma_nom, cb_gaus_mix_sigma_err_stat, cb_gaus_mix_sigma_err_syst;
    cb_gaus_combine(cb_fit_mean_nom, cb_fit_mean_err_stat, cb_fit_mean_err_syst, gaus_fit_mean_nom, gaus_fit_mean_err_stat, gaus_fit_mean_err_syst, cb_gaus_mix_mean_nom, cb_gaus_mix_mean_err_stat, cb_gaus_mix_mean_err_syst);
    cb_gaus_combine(cb_fit_sigma_nom, cb_fit_sigma_err_stat, cb_fit_sigma_err_syst, gaus_fit_sigma_nom, gaus_fit_sigma_err_stat, gaus_fit_sigma_err_syst, cb_gaus_mix_sigma_nom, cb_gaus_mix_sigma_err_stat, cb_gaus_mix_sigma_err_syst);
    double resolution, resolution_err_stat, resolution_err_syst;
    resolution_calculator(cb_gaus_mix_mean_nom, cb_gaus_mix_mean_err_stat, cb_gaus_mix_mean_err_syst, cb_gaus_mix_sigma_nom, cb_gaus_mix_sigma_err_stat, cb_gaus_mix_sigma_err_syst, resolution, resolution_err_stat, resolution_err_syst);
    // std::string resolution_info = "Resolution = " + format_sig3(resolution) + "% #pm " + format_sig3(resolution_err_stat) + "% (stat) #pm " + format_sig3(resolution_err_syst) + "% (syst)";
    // legend_adc_peak_sum->AddEntry((TObject*)nullptr, resolution_info.c_str(), "");

    legend_adc_peak_sum->Draw();
    canvas_adc_peak_sum->Modified();
    canvas_adc_peak_sum->Update();
    canvas_adc_peak_sum->Write();
    canvas_adc_peak_sum->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_adc_peak_sum_distribution.pdf").c_str());
    canvas_adc_peak_sum->Close();

    // write the six values to the output file
    // double cb_gaus_mix_mean_nom, cb_gaus_mix_mean_err_stat, cb_gaus_mix_mean_err_syst;
    // double cb_gaus_mix_sigma_nom, cb_gaus_mix_sigma_err_stat, cb_gaus_mix_sigma_err_syst;
    auto* cb_gaus_mix_mean_value_obj     = new TParameter<double>("cb_gaus_mix_mean_value",     cb_gaus_mix_mean_nom);
    auto* cb_gaus_mix_mean_err_stat_obj  = new TParameter<double>("cb_gaus_mix_mean_err_stat",  cb_gaus_mix_mean_err_stat);
    auto* cb_gaus_mix_mean_err_syst_obj  = new TParameter<double>("cb_gaus_mix_mean_err_syst",  cb_gaus_mix_mean_err_syst);
    auto* cb_gaus_mix_sigma_value_obj    = new TParameter<double>("cb_gaus_mix_sigma_value",    cb_gaus_mix_sigma_nom);
    auto* cb_gaus_mix_sigma_err_stat_obj = new TParameter<double>("cb_gaus_mix_sigma_err_stat", cb_gaus_mix_sigma_err_stat);
    auto* cb_gaus_mix_sigma_err_syst_obj = new TParameter<double>("cb_gaus_mix_sigma_err_syst", cb_gaus_mix_sigma_err_syst);
    cb_gaus_mix_mean_value_obj->Write();
    cb_gaus_mix_mean_err_stat_obj->Write();
    cb_gaus_mix_mean_err_syst_obj->Write();
    cb_gaus_mix_sigma_value_obj->Write();
    cb_gaus_mix_sigma_err_stat_obj->Write();
    cb_gaus_mix_sigma_err_syst_obj->Write();

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}