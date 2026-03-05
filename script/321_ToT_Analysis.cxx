#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "H2GCROC_ToT.hxx"
#include "TRandom3.h"
#include "CommonParams.hxx"

double reverse_Hill(double tot_decoded) {
    double A = 3640.0;
    double n = 0.62;
    double x_1_2 = 0.13;

    double x_0 = 5.8;
    double x_zoom = 9.95;

    double u = x_1_2 * pow((tot_decoded / (A - tot_decoded)), (1/n));
    double x = (u / x_zoom) + x_0;
    return x;
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

    std::string example_tot_lut_file = "config/ToTScan5.root_LUT_Channel_52.txt";
    TotToAdcLUT lut = LoadTotToAdcLUT(example_tot_lut_file);

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
    const int adc_peak_min_index = CommonParams::adc_peak_min_index;
    const int adc_peak_max_index = CommonParams::adc_peak_max_index;

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

    // * === Define output data structures ===================================================
    // * =====================================================================================
    const double tot_bin_max = 5000.0; // this is the maximum possible ToT value (12 bits)
    const double tot_bin_min = 0.0;
    const int tot_bin_number = 500; // we want to keep the full resolution

    const double tot_recovered_adc_max = 5e3; // this is an estimated maximum ADC value after ToT to ADC conversion, used for histogram range
    const double tot_recovered_adc_min = 0.0;
    const int tot_recovered_adc_bin_number = 500; // we want to keep the

    const double adc_bin_max = 1024.0;
    const double adc_bin_min = 0.0;
    const int adc_bin_number = 256;

    const double hill_bin_max = 10.0;
    const double hill_bin_min = 5.5;
    const int hill_bin_number = 200;

    std::vector<TH1D*> h1_adc_peak_channel_list(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h1_name = "h1_adc_peak_channel_" + std::to_string(vldb_id) + "_ch" + std::to_string(channel);
            std::string h1_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " ADC Peak Distribution;ADC Peak Value;Counts";
            h1_adc_peak_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH1D(
                h1_name.c_str(),
                h1_title.c_str(),
                adc_bin_number, adc_bin_min, adc_bin_max
            );
            h1_adc_peak_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<TH1D*> h1_tot_channel_list(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h1_name = "h1_tot_vldb" + std::to_string(vldb_id) + "_channel" + std::to_string(channel);
            std::string h1_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " ToT Distribution;ToT Value;Counts";
            h1_tot_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH1D(
                h1_name.c_str(),
                h1_title.c_str(),
                tot_bin_number, tot_bin_min, tot_bin_max
            );
            h1_tot_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<TH1D*> h1_tot_recovered_adc_channel_list(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h1_name = "h1_tot_recovered_adc_vldb" + std::to_string(vldb_id) + "_channel" + std::to_string(channel);
            std::string h1_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " Recovered ADC from ToT;ADC Value;Counts";
            h1_tot_recovered_adc_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH1D(
                h1_name.c_str(),
                h1_title.c_str(),
                tot_recovered_adc_bin_number, tot_recovered_adc_min, tot_recovered_adc_max
            );
            h1_tot_recovered_adc_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<TH1D*> h1_hill_recovered_channel_list(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h1_name = "h1_hill_recovered_adc_vldb" + std::to_string(vldb_id) + "_channel" + std::to_string(channel);
            std::string h1_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " Recovered ADC from ToT with Hill Function;ADC Value;Counts";
            h1_hill_recovered_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH1D(
                h1_name.c_str(),
                h1_title.c_str(),
                hill_bin_number, hill_bin_min, hill_bin_max
            );
            h1_hill_recovered_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<TH2D*> h2_hill_adc_average_correlation_list(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h2_name = "h2_hill_adc_average_correlation_vldb" + std::to_string(vldb_id) + "_channel" + std::to_string(channel);
            std::string h2_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " Recovered ADC from ToT with Hill Function vs ADC Average Correlation;Recovered ADC Value;ADC Average Value";
            h2_hill_adc_average_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH2D(
                h2_name.c_str(),
                h2_title.c_str(),
                hill_bin_number, hill_bin_min, hill_bin_max,
                adc_bin_number, adc_bin_min, adc_bin_max
            );
            h2_hill_adc_average_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<TH2D*> h2_tot_adc_correlation_list(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h2_name = "h2_tot_adc_correlation_vldb" + std::to_string(vldb_id) + "_channel" + std::to_string(channel);
            std::string h2_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " ToT vs ADC Correlation;ToT Value;ADC Value";
            h2_tot_adc_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH2D(
                h2_name.c_str(),
                h2_title.c_str(),
                tot_bin_number, tot_bin_min, tot_bin_max,
                adc_bin_number, adc_bin_min, adc_bin_max
            );
            h2_tot_adc_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<TH2D*> h2_recovered_adc_adc_correlation_list(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h2_name = "h2_recovered_adc_adc_correlation_vldb" + std::to_string(vldb_id) + "_channel" + std::to_string(channel);
            std::string h2_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " Recovered ADC vs ADC Correlation;Recovered ADC Value;ADC Value";
            h2_recovered_adc_adc_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH2D(
                h2_name.c_str(),
                h2_title.c_str(),
                tot_recovered_adc_bin_number, tot_recovered_adc_min, tot_recovered_adc_max,
                adc_bin_number, adc_bin_min, adc_bin_max
            );
            h2_recovered_adc_adc_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<TH2D*> h2_adc_peak_adc_average_correlation_list(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h2_name = "h2_adc_peak_adc_average_correlation_vldb" + std::to_string(vldb_id) + "_channel" + std::to_string(channel);
            std::string h2_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " ADC Peak vs ADC Average Correlation;ADC Peak Value;ADC Average Value";
            h2_adc_peak_adc_average_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH2D(
                h2_name.c_str(),
                h2_title.c_str(),
                adc_bin_number, adc_bin_min, adc_bin_max,
                adc_bin_number, adc_bin_min, adc_bin_max
            );
            h2_adc_peak_adc_average_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }
    
    std::vector<std::vector<double>> adc_peak_value_in_channel_list(vldb_number*FPGA_CHANNEL_NUMBER);
    std::vector<std::vector<double>> adc_peak_raw_value_in_channel_list(vldb_number*FPGA_CHANNEL_NUMBER); // for saturation detection
    std::vector<std::vector<double>> tot_value_in_channel_list(vldb_number*FPGA_CHANNEL_NUMBER);
    std::vector<std::vector<double>> hill_in_channel_list(vldb_number*FPGA_CHANNEL_NUMBER);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            adc_peak_value_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].reserve(entry_max);
            adc_peak_raw_value_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].reserve(entry_max);
            tot_value_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].reserve(entry_max);
            hill_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].reserve(entry_max);
        }
    }

    std::vector<std::vector<std::vector<double>>> adc_samples_in_CM_channels_all_vldb(
        vldb_number, std::vector<std::vector<double>>(4, std::vector<double>(machine_gun_samples, 0.0)));

    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
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
            } // end of half loop for CM subtraction
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int valid_channel_number = get_valid_fpga_channel(channel);
                int channel_half_index = (channel / (FPGA_CHANNEL_NUMBER / 4)) % 4;
                std::vector<double> CM_samples = adc_samples_in_CM_channels[channel_half_index];

                double channel_tot_value = 0.0;
                double channel_toa_value = 0.0;
                int channel_tot_index = -1;
                int channel_toa_index = -1;

                std::vector<UInt_t> adc_pedestal_samples; // only take the first 3 samples
                std::vector<UInt_t> adc_samples; // all samples for common-mode calculation
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    UInt_t adc_value = val0_list_pools[vldb_id][0][idx];
                    UInt_t tot_value = val1_list_pools[vldb_id][0][idx];
                    UInt_t toa_value = val2_list_pools[vldb_id][0][idx];

                    adc_samples.push_back(adc_value);
                     if (sample < 3) {
                        adc_pedestal_samples.push_back(adc_value);
                    }
                    if (channel_tot_value == 0.0 && tot_value > 0) {
                        channel_tot_value = static_cast<double>(tot_value);
                        if (channel_tot_value >= 512.0) {
                            channel_tot_value = (channel_tot_value - 512.0) * 8.0;
                        }
                        channel_tot_index = sample;
                    }
                    if (channel_toa_value == 0.0 && toa_value > 0) {
                        channel_toa_value = static_cast<double>(toa_value);
                        channel_toa_index = sample;
                    }
                } // end of sample loop
                double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
                double adc_peak = 0.0;
                for (int sample = adc_peak_min_index; sample <= adc_peak_max_index; sample++) {
                    double adc_value = static_cast<double>(adc_samples[sample]) - adc_pedestal - CM_samples[sample];
                    if (adc_value > adc_peak) {
                        adc_peak = adc_value;
                    }
                }
                if (adc_peak < 0.0) {
                    adc_peak = 0.0;
                }
                double adc_sum = 0.0;
                double adc_count = 0.0;
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    double adc_value = static_cast<double>(adc_samples[sample]) - adc_pedestal - CM_samples[sample];
                    adc_sum += adc_value;
                    adc_count += 1.0;
                }
                double adc_average = adc_sum / adc_count;

                if (adc_peak < 880.0) {
                    h2_adc_peak_adc_average_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(adc_peak, adc_average);
                }

                adc_peak_value_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].push_back(adc_peak);
                adc_peak_raw_value_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].push_back(adc_samples[adc_peak_max_index]);

                h1_adc_peak_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(adc_peak);
                if (channel_tot_value > 0.0 && channel_tot_value < 4080.0) {
                    h1_tot_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(channel_tot_value);
                    double tot_recovered_adc = lut.Eval(channel_tot_value);
                    h1_tot_recovered_adc_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(tot_recovered_adc);
                    
                    h2_tot_adc_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(channel_tot_value, adc_average);
                    h2_recovered_adc_adc_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(tot_recovered_adc, adc_average);

                    auto hill_estimated_adc = reverse_Hill(channel_tot_value);
                    h1_hill_recovered_channel_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(hill_estimated_adc);
                    h2_hill_adc_average_correlation_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(hill_estimated_adc, adc_average);

                    tot_value_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].push_back(channel_tot_value);
                    hill_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].push_back(hill_estimated_adc);
                } else {
                    tot_value_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].push_back(0.0);
                    hill_in_channel_list[vldb_id*FPGA_CHANNEL_NUMBER + channel].push_back(0.0);
                }
            } // end of channel loop
        } // end of VLDB loop
    } // end of event loop
    input_root->Close();
    spdlog::info("Finished processing {} events", entry_max);

    const int NX = 16, NY = 12;
    const int board_cols = 8, board_rows = 4;

    auto chan2pad = build_chan2pad_LUT(
        vldb_number, FPGA_CHANNEL_NUMBER,
        NX, NY, board_cols, board_rows,
        sipm_board, board_loc, board_rotation, board_flip
    );

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
    const std::string out_pdf = script_output_file.substr(0, script_output_file.find_last_of(".")) + ".pdf";

    auto canvas_start_dummy = new TCanvas("canvas_start_dummy", "Start Dummy Canvas", 800, 600);
    canvas_start_dummy->Print((out_pdf + "[").c_str()); // begin of pdf
    canvas_start_dummy->Close();

    auto canvas_adc_peak_channel = new TCanvas("canvas_adc_peak_channel", "ADC Peak Distribution per Channel", 1600, 1200);
    draw_mosaic_fixed(*canvas_adc_peak_channel, h1_adc_peak_channel_list, topo_ped_median);
    canvas_adc_peak_channel->Modified();
    canvas_adc_peak_channel->Update();
    canvas_adc_peak_channel->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_adc_peak_dist.pdf").c_str());
    canvas_adc_peak_channel->Write();
    canvas_adc_peak_channel->Close();
    delete canvas_adc_peak_channel;

    auto canvas_map_tot_channel = new TCanvas("canvas_map_tot_channel", "ToT Distribution per Channel", 1600, 1200);
    draw_mosaic_fixed(*canvas_map_tot_channel, h1_tot_channel_list, topo_ped_median);
    canvas_map_tot_channel->Modified();
    canvas_map_tot_channel->Update();
    canvas_map_tot_channel->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_tot_dist.pdf").c_str());
    canvas_map_tot_channel->Write();
    canvas_map_tot_channel->Close();
    delete canvas_map_tot_channel;

    auto canvas_map_tot_recovered_adc_channel = new TCanvas("canvas_map_tot_recovered_adc_channel", "Recovered ADC from ToT Distribution per Channel", 1600, 1200);
    draw_mosaic_fixed(*canvas_map_tot_recovered_adc_channel, h1_tot_recovered_adc_channel_list, topo_ped_median);
    canvas_map_tot_recovered_adc_channel->Modified();
    canvas_map_tot_recovered_adc_channel->Update();
    canvas_map_tot_recovered_adc_channel->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_tot_recovered_adc_dist.pdf").c_str());
    canvas_map_tot_recovered_adc_channel->Write();
    canvas_map_tot_recovered_adc_channel->Close();
    delete canvas_map_tot_recovered_adc_channel;

    auto canvas_map_tot_adc_correlation = new TCanvas("canvas_map_tot_adc_correlation", "ToT vs ADC Correlation per Channel", 1600, 1200);
    draw_mosaic_fixed(*canvas_map_tot_adc_correlation, h2_tot_adc_correlation_list, topo_wave);
    canvas_map_tot_adc_correlation->Modified();
    canvas_map_tot_adc_correlation->Update();
    canvas_map_tot_adc_correlation->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_tot_adc_correlation.pdf").c_str());
    canvas_map_tot_adc_correlation->Write();
    canvas_map_tot_adc_correlation->Close();
    delete canvas_map_tot_adc_correlation;

    // take the correlation plot of channel 72 and do a linear fit
    auto canvas_example_tot_adc_correlation = new TCanvas("canvas_example_tot_adc_correlation", "Example ToT vs ADC Correlation", 800, 600);
    int example_channel_global = 33 + FPGA_CHANNEL_NUMBER; // this is the global channel number, which corresponds to channel 9 in the second VLDB
    int example_vldb_id = example_channel_global / FPGA_CHANNEL_NUMBER;
    int example_channel_local = example_channel_global % FPGA_CHANNEL_NUMBER;
    TH2D* example_hist = h2_tot_adc_correlation_list[example_vldb_id * FPGA_CHANNEL_NUMBER + example_channel_local];
    example_hist->GetXaxis()->SetTitle("ToT Value");
    example_hist->GetYaxis()->SetTitle("ADC Value");
    std::string example_canvas_info = "Channel " + std::to_string(example_channel_global) + " ToT vs ADC Correlation";
    format_2d_hist_canvas(canvas_example_tot_adc_correlation, example_hist, kBlue+2, annotation_canvas_title, annotation_testbeam_title, example_canvas_info);
    example_hist->Fit("pol1", "Q");
    TF1* fit_func = example_hist->GetFunction("pol1");
    if (fit_func != nullptr) {
        double fit_slope = fit_func->GetParameter(1);
        double fit_intercept = fit_func->GetParameter(0);
        TLatex fit_info;
        fit_info.SetNDC();
        fit_info.SetTextSize(0.03);
        fit_info.SetTextColor(kRed+2);
        std::string fit_info_text = "Fit: ADC = " + std::to_string(fit_slope) + " * ToT + " + std::to_string(fit_intercept);
        fit_info.DrawLatex(0.40, 0.13, fit_info_text.c_str());
    }
    canvas_example_tot_adc_correlation->Modified();
    canvas_example_tot_adc_correlation->Update();
    canvas_example_tot_adc_correlation->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_example_tot_adc_correlation.pdf").c_str());
    canvas_example_tot_adc_correlation->Write();
    canvas_example_tot_adc_correlation->Close();
    delete canvas_example_tot_adc_correlation; 

    auto canvas_map_adc_peak_adc_average_correlation = new TCanvas("canvas_map_adc_peak_adc_average_correlation", "ADC Peak vs ADC Average Correlation per Channel", 1600, 1200);
    draw_mosaic_fixed(*canvas_map_adc_peak_adc_average_correlation, h2_adc_peak_adc_average_correlation_list, topo_wave);
    canvas_map_adc_peak_adc_average_correlation->Modified();
    canvas_map_adc_peak_adc_average_correlation->Update();
    canvas_map_adc_peak_adc_average_correlation->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_adc_peak_adc_average_correlation.pdf").c_str());
    canvas_map_adc_peak_adc_average_correlation->Write();
    canvas_map_adc_peak_adc_average_correlation->Close();
    delete canvas_map_adc_peak_adc_average_correlation;

    auto canvas_map_hill_recovered_adc_channel = new TCanvas("canvas_map_hill_recovered_adc_channel", "Recovered ADC from ToT with Hill Function per Channel", 1600, 1200);
    draw_mosaic_fixed(*canvas_map_hill_recovered_adc_channel, h1_hill_recovered_channel_list, topo_ped_median);
    canvas_map_hill_recovered_adc_channel->Modified();
    canvas_map_hill_recovered_adc_channel->Update();
    canvas_map_hill_recovered_adc_channel->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_hill_recovered_adc_dist.pdf").c_str());
    canvas_map_hill_recovered_adc_channel->Write();
    canvas_map_hill_recovered_adc_channel->Close();
    delete canvas_map_hill_recovered_adc_channel;

    auto canvas_map_hill_adc_average_correlation = new TCanvas("canvas_map_hill_adc_average_correlation", "Recovered ADC from ToT with Hill Function vs ADC Average Correlation per Channel", 1600, 1200);
    draw_mosaic_fixed(*canvas_map_hill_adc_average_correlation, h2_hill_adc_average_correlation_list, topo_wave);
    canvas_map_hill_adc_average_correlation->Modified();
    canvas_map_hill_adc_average_correlation->Update();
    canvas_map_hill_adc_average_correlation->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_hill_adc_average_correlation.pdf").c_str());
    canvas_map_hill_adc_average_correlation->Write();
    canvas_map_hill_adc_average_correlation->Close();
    delete canvas_map_hill_adc_average_correlation;

    // fit for the example channel
    auto canvas_example_adc_peak_adc_average_correlation = new TCanvas("canvas_example_adc_peak_adc_average_correlation", "Example ADC Peak vs ADC Average Correlation", 800, 600);
    TH2D* example_hist_adc_peak_average = h2_adc_peak_adc_average_correlation_list[example_vldb_id * FPGA_CHANNEL_NUMBER + example_channel_local];
    example_hist_adc_peak_average->GetXaxis()->SetTitle("ADC Peak Value");
    example_hist_adc_peak_average->GetYaxis()->SetTitle("ADC Average Value");
    std::string example_canvas_info_adc_peak_average = "Channel " + std::to_string(example_channel_global) + " ADC Peak vs ADC Average Correlation";
    format_2d_hist_canvas(canvas_example_adc_peak_adc_average_correlation, example_hist_adc_peak_average, kGreen+2, annotation_canvas_title, annotation_testbeam_title, example_canvas_info_adc_peak_average);
    example_hist_adc_peak_average->Fit("pol1", "Q");
    TF1* fit_func_adc_peak_average = example_hist_adc_peak_average->GetFunction("pol1");
    if (fit_func_adc_peak_average != nullptr) {
        double fit_slope = fit_func_adc_peak_average->GetParameter(1);
        double fit_intercept = fit_func_adc_peak_average->GetParameter(0);
        TLatex fit_info;
        fit_info.SetNDC();
        fit_info.SetTextSize(0.03);
        fit_info.SetTextColor(kRed+2);
        std::string fit_info_text = "Fit: ADC Average = " + std::to_string(fit_slope) + " * ADC Peak + " + std::to_string(fit_intercept);
        fit_info.DrawLatex(0.40, 0.13, fit_info_text.c_str());
    }
    canvas_example_adc_peak_adc_average_correlation->Modified();
    canvas_example_adc_peak_adc_average_correlation->Update();
    canvas_example_adc_peak_adc_average_correlation->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_example_adc_peak_adc_average_correlation.pdf").c_str());
    canvas_example_adc_peak_adc_average_correlation->Write();
    canvas_example_adc_peak_adc_average_correlation->Close();
    delete canvas_example_adc_peak_adc_average_correlation;

    auto canvas_map_recovered_adc_adc_correlation = new TCanvas("canvas_map_recovered_adc_adc_correlation", "Recovered ADC vs ADC Correlation per Channel", 1600, 1200);
    draw_mosaic_fixed(*canvas_map_recovered_adc_adc_correlation, h2_recovered_adc_adc_correlation_list, topo_wave);
    canvas_map_recovered_adc_adc_correlation->Modified();
    canvas_map_recovered_adc_adc_correlation->Update();
    canvas_map_recovered_adc_adc_correlation->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_recovered_adc_adc_correlation.pdf").c_str());
    canvas_map_recovered_adc_adc_correlation->Write();
    canvas_map_recovered_adc_adc_correlation->Close();
    delete canvas_map_recovered_adc_adc_correlation;

    // ! do the scan over different hill - adc slope and intercept to find the best parameters
    const int scan_channel_global = 33 + FPGA_CHANNEL_NUMBER;
    double scan_slope_min       = 100;
    double scan_slope_max       = 8000;
    double scan_intercept_min   = 2.0;
    double scan_intercept_max   = 6.0;
    int scan_slope_steps        = 200;
    int scan_intercept_steps    = 100;
    double adc_bin_size         = 8.0; 

    double example_slope = 7100.0;
    double example_intercept = 5.72;

    int example_slope_step = static_cast<int>((example_slope - scan_slope_min) / (scan_slope_max - scan_slope_min) * (scan_slope_steps - 1));
    int example_intercept_step = static_cast<int>((example_intercept - scan_intercept_min) / (scan_intercept_max - scan_intercept_min) * (scan_intercept_steps - 1));

    auto& scan_adc_list      = adc_peak_value_in_channel_list[scan_channel_global];
    auto& scan_adc_raw_list  = adc_peak_raw_value_in_channel_list[scan_channel_global];
    auto& scan_tot_list      = tot_value_in_channel_list[scan_channel_global];
    auto& scan_hill_list     = hill_in_channel_list[scan_channel_global];

    // create a 2D histogram to evaluate the combination of slope and intercept
    TH2D* scan_result_hist = new TH2D("scan_result_hist", "Scan Result for Hill-ADC Slope and Intercept;Slope;Intercept", scan_slope_steps, scan_slope_min, scan_slope_max, scan_intercept_steps, scan_intercept_min, scan_intercept_max);
    TH2D* scan_prob_hist = new TH2D("scan_prob_hist", "Scan Result for Hill-ADC Slope and Intercept;Slope;Intercept", scan_slope_steps, scan_slope_min, scan_slope_max, scan_intercept_steps, scan_intercept_min, scan_intercept_max);

    for (int slope_step = 0; slope_step < scan_slope_steps; slope_step++) {
        double slope = scan_slope_min + (scan_slope_max - scan_slope_min) * slope_step / (scan_slope_steps - 1);
        for (int intercept_step = 0; intercept_step < scan_intercept_steps; intercept_step++) {
            double intercept = scan_intercept_min + (scan_intercept_max - scan_intercept_min) * intercept_step / (scan_intercept_steps - 1);
            // estimate the maximum adc from the combination
            double max_hill_allowed = 8.0;
            double max_adc_from_hill = (max_hill_allowed-intercept) * slope;
            int adc_bin_number = static_cast<int>(std::ceil(max_adc_from_hill / adc_bin_size));
            TH1D* combined_adc_hist = new TH1D("combined_adc_hist", "Combined ADC Histogram for Scan;ADC Value;Counts", adc_bin_number, 0, max_adc_from_hill);
            int total_entries = scan_adc_list.size();
            int tot_used_entries = 0;
            for (size_t i = 0; i < scan_adc_list.size(); i++) {
                auto& adc_peak  = scan_adc_list[i];
                auto& adc_raw   = scan_adc_raw_list[i];
                auto& tot   = scan_tot_list[i];
                auto& hill  = scan_hill_list[i];
                if (adc_raw < 1023.0) { // only consider non-saturated events for the scan
                    combined_adc_hist->Fill(adc_peak);
                } else if (tot > 0.0 && tot < 4080.0) { // for saturated events, use the hill estimation if available
                    double hill_estimated_adc = (hill - intercept) * slope;
                    if (hill_estimated_adc > 0.0 && hill_estimated_adc < max_adc_from_hill) {
                        combined_adc_hist->Fill(hill_estimated_adc);
                        tot_used_entries++;
                    } else {
                        combined_adc_hist->Fill(adc_peak);
                    }
                } else {
                    combined_adc_hist->Fill(adc_peak);
                }
            } // end of event loop

            LOG(INFO) << "Slope: " << slope << ", Intercept: " << intercept << ", Total Entries: " << total_entries << ", Tot Used Entries: " << tot_used_entries;
            double mean_adc = combined_adc_hist->GetMean();
            double rms_adc  = combined_adc_hist->GetRMS();

            // guard against empty / pathological hist
            double integral = combined_adc_hist->Integral();
            if (integral <= 0 || !std::isfinite(mean_adc) || !std::isfinite(rms_adc) || rms_adc <= 0) {
                scan_result_hist->SetBinContent(slope_step + 1, intercept_step + 1, 9999.0);
                scan_prob_hist->SetBinContent  (slope_step + 1, intercept_step + 1, 0.0);
                delete combined_adc_hist;
                continue;
            }

            double smoothness = SmoothnessRoughness2(combined_adc_hist);
            scan_result_hist->SetBinContent(slope_step + 1, intercept_step + 1, smoothness);
            // perform a Kolmogorov test against the example histogram with


            if (slope_step == example_slope_step && intercept_step == example_intercept_step) {
                TCanvas* example_canvas = new TCanvas("example_canvas", "Example Combined ADC Distribution with Hill-ADC Parameters", 800, 600);
                combined_adc_hist->GetXaxis()->SetTitle("ADC Value");
                combined_adc_hist->GetYaxis()->SetTitle("Counts");
                std::string example_canvas_info = "Example Combined ADC Distribution with Hill-ADC Parameters: Slope = " + std::to_string(slope) + ", Intercept = " + std::to_string(intercept);
                auto combined_adc_hist_clone = (TH1D*)combined_adc_hist->Clone("combined_adc_hist_clone");
                format_1d_hist_canvas(example_canvas, combined_adc_hist_clone, kMagenta+2, annotation_canvas_title, annotation_testbeam_title, example_canvas_info);
                example_canvas->Modified();
                example_canvas->Update();
                example_canvas->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_example_combined_adc.pdf").c_str());
                example_canvas->Write();
                example_canvas->Close();
            }

            delete combined_adc_hist;
        } // end of intercept loop
    } // end of slope loop
    // draw the scan result
    auto canvas_scan_result = new TCanvas("canvas_scan_result", "Scan Result for Hill-ADC Slope and Intercept", 800, 600);
    canvas_scan_result->SetLogz();
    scan_result_hist->GetXaxis()->SetTitle("Hill-ADC Slope");
    scan_result_hist->GetYaxis()->SetTitle("Hill-ADC Intercept");
    scan_result_hist->GetZaxis()->SetTitle("Mean ADC Value");
    scan_result_hist->SetStats(false);
    scan_result_hist->Draw("COLZ");
    canvas_scan_result->Modified();
    canvas_scan_result->Update();
    canvas_scan_result->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_hill_adc_scan_result.pdf").c_str());
    canvas_scan_result->Write();
    canvas_scan_result->Close();
    delete canvas_scan_result;

    auto canvas_scan_prob = new TCanvas("canvas_scan_prob", "Scan Probability for Hill-ADC Slope and Intercept", 800, 600);
    canvas_scan_prob->SetLogz();
    scan_prob_hist->GetXaxis()->SetTitle("Hill-ADC Slope");
    scan_prob_hist->GetYaxis()->SetTitle("Hill-ADC Intercept");
    scan_prob_hist->GetZaxis()->SetTitle("Kolmogorov Test Probability");
    scan_prob_hist->SetStats(false);
    scan_prob_hist->Draw("COLZ");
    canvas_scan_prob->Modified();
    canvas_scan_prob->Update();
    canvas_scan_prob->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_hill_adc_scan_prob.pdf").c_str());
    canvas_scan_prob->Write();
    canvas_scan_prob->Close();
    delete canvas_scan_prob;

    auto canvas_dummy = new TCanvas("canvas_dummy", "Dummy Canvas", 800, 600);
    canvas_dummy->Print((out_pdf + "]").c_str()); // end of pdf
    canvas_dummy->Close();
    delete canvas_dummy;

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);
    return 0;
}
