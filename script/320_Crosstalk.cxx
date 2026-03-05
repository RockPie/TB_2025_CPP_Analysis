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

    // * === Define output data structures ===================================================
    // * =====================================================================================
    // * --- Accumulate waveform for each channel ---
    std::vector<TH2D*> h2_adc_waveforms(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h2_name = "h2_adc_waveform_vldb" + std::to_string(vldb_id) + "_ch" + std::to_string(channel);
            std::string h2_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " ADC Waveform;Sample Index;ADC Value";
            h2_adc_waveforms[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH2D(
                h2_name.c_str(),
                h2_title.c_str(),
                machine_gun_samples, -0.5, machine_gun_samples - 0.5,
                512, 0, 1024
            );
            h2_adc_waveforms[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    // * --- ADC Peak Position Histogram ---
    TH1D* h1_adc_peak_position = new TH1D("h1_adc_peak_position", "ADC Peak Position;Sample Index;Counts", machine_gun_samples, -0.5, machine_gun_samples - 0.5);
    h1_adc_peak_position->SetDirectory(nullptr);

    // * --- ADC Peak Value - Channel Correlation Histogram ---
    TH2D *h2_adc_peak_channel_correlation = new TH2D(
        "h2_adc_peak_channel_correlation",
        "ADC Peak Value vs Channel Correlation;Channel;ADC Peak Value",
        FPGA_CHANNEL_NUMBER, -0.5, FPGA_CHANNEL_NUMBER - 0.5,
        512, 0, 1024
    );
    h2_adc_peak_channel_correlation->SetDirectory(nullptr);

    // * --- ADC Sum distribution Histogram ---
    std::vector<double> adc_sum_list;
    adc_sum_list.reserve(entry_max);

    std::vector<std::vector<double>> CM_example_channel_samples(entry_max, std::vector<double>(machine_gun_samples, 0.0));
    std::vector<std::vector<double>> example_channel_samples(entry_max, std::vector<double>(machine_gun_samples, 0.0));
    // get which ASIC half the example channel belongs to
    int example_channel_vldb = example_channel / FPGA_CHANNEL_NUMBER;
    int example_channel_local = example_channel % FPGA_CHANNEL_NUMBER;
    int half_size = FPGA_CHANNEL_NUMBER / 4; // 4 halves per VLDB
    int example_channel_half = example_channel_local / half_size;
    int CM_channel_global = example_channel_vldb * FPGA_CHANNEL_NUMBER + example_channel_half * half_size + CM_channel;

    // const int adc_peak_min_index = CommonParams::adc_peak_min_index;
    // const int adc_peak_max_index = CommonParams::adc_peak_max_index;
    const int adc_peak_min_index = 5;
    const int adc_peak_max_index = 7;

    // print the adc peak index range
    spdlog::info("Using ADC peak index range: {} to {}", adc_peak_min_index, adc_peak_max_index);

    // first loop to find the dynamic range of each channel and whether it is a saturated channel
    std::vector<std::vector<UInt_t>> adc_peak_values_per_channel;
    adc_peak_values_per_channel.resize(vldb_number * FPGA_CHANNEL_NUMBER);
    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            // channel loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int valid_channel_number = get_valid_fpga_channel(channel);
                std::vector<UInt_t> adc_pedestal_samples; // only take the first 3 samples
                int adc_peak_index = -1;
                UInt_t adc_peak_value = 0;
                int adc_peak_ranged_index = -1;
                UInt_t adc_peak_ranged_value = 0;
                UInt_t adc_min_value = 1024;
                adc_pedestal_samples.reserve(3);
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    UInt_t adc_value = val0_list_pools[vldb_id][0][idx];
                    UInt_t tot_value = val1_list_pools[vldb_id][0][idx];
                    UInt_t toa_value = val2_list_pools[vldb_id][0][idx];

                    if (sample < 3) {
                        adc_pedestal_samples.push_back(adc_value);
                    }
                    if (adc_value > adc_peak_value) {
                        adc_peak_value = adc_value;
                        adc_peak_index = sample;
                    }
                    if (sample < 8) {
                        if (adc_value < adc_min_value) {
                            adc_min_value = adc_value;
                        }
                    }
                    if (sample >= 0 && sample <= 15) {  
                        // relaxed range because I only care about saturation, not the peak position
                        if (adc_value > adc_peak_ranged_value) {
                            adc_peak_ranged_value = adc_value;
                            adc_peak_ranged_index = sample;
                        }
                    }
                } // end of sample loop
                // double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
                double adc_pedestal = static_cast<double>(adc_min_value);
                double adc_peak_value_pede_sub = static_cast<double>(adc_peak_ranged_value) - adc_pedestal;
                // double adc_peak_value_CM_sub = adc_peak_value_pede_sub - static_cast<double>(CM_samples[adc_peak_ranged_index]);
                adc_peak_values_per_channel[vldb_id * FPGA_CHANNEL_NUMBER + channel].push_back(adc_peak_value_pede_sub);
                // adc_peak_values_per_channel[vldb_id * FPGA_CHANNEL_NUMBER + channel].push_back(adc_peak_value_CM_sub);
            }
        }
    }
    UInt_t global_min_dynamic_range = 1024;
    for (int channel_global = 0; channel_global < vldb_number * FPGA_CHANNEL_NUMBER; channel_global++) {
        auto& values = adc_peak_values_per_channel[channel_global];
        std::sort(values.begin(), values.end());
        // if more than 1% of the values are at the maximum ADC value, consider this channel as saturated
        int n_saturated = 0;
        UInt_t max_adc_value = values.back();
        for (const auto& val : values) {
            if (val >= max_adc_value) {
                n_saturated++;
            }
        }
        double saturated_fraction = static_cast<double>(n_saturated) / static_cast<double>(values.size());
        if (saturated_fraction > 0.0002 && max_adc_value > 800) {
            spdlog::warn("Channel {} is considered saturated ({}% of events at max ADC value {})", channel_global, format_decimal(saturated_fraction * 100.0, 2), max_adc_value);
            if (max_adc_value < global_min_dynamic_range) {
                global_min_dynamic_range = max_adc_value;
            }
        }
    }
    spdlog::info("Global minimum dynamic range across all channels is {}", global_min_dynamic_range);
    // Initialize common-mode samples storage for all VLDBs
    std::vector<std::vector<std::vector<double>>> adc_samples_in_CM_channels_all_vldb(
        vldb_number, std::vector<std::vector<double>>(4, std::vector<double>(machine_gun_samples, 0.0)));
    // start event loop
    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        double adc_sum = 0.0;
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

                // if it is example CM channel, store the samples for later plotting
                if (channel + vldb_id * FPGA_CHANNEL_NUMBER == CM_channel_global) {
                    for (int sample = 0; sample < machine_gun_samples; sample++) {
                        CM_example_channel_samples[entry][sample] = adc_samples_in_CM_channels[half][sample];
                    }
                }
            }
            // channel loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int valid_channel_number = get_valid_fpga_channel(channel);
                int channel_half_index = (channel / (FPGA_CHANNEL_NUMBER / 4)) % 4;
                // skip
                std::vector<double> CM_samples = adc_samples_in_CM_channels[channel_half_index];
                std::vector<double> adc_pedestal_samples; // only take the first 3 samples
                std::vector<double> adc_raw_samples; // all samples for later use
                int adc_peak_index = -1;
                double adc_peak_value = 0.0;
                int adc_peak_ranged_index = -1;
                double adc_peak_ranged_value = 0.0;
                double adc_min_value = 1024.0;
                adc_pedestal_samples.reserve(3);
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    double adc_value = static_cast<double>(val0_list_pools[vldb_id][0][idx]);
                    UInt_t tot_value = val1_list_pools[vldb_id][0][idx];
                    UInt_t toa_value = val2_list_pools[vldb_id][0][idx];
                    adc_raw_samples.push_back(adc_value);
                    // spdlog::info("Event {} VLDB {} Channel {} Sample {}: ADC {}, ToT {}, ToA {}", entry, vldb_id, channel, sample, adc_value, tot_value, toa_value);
                    adc_value -= CM_samples[sample]; // common-mode subtraction
                    if (sample < 3) {
                        adc_pedestal_samples.push_back(adc_value);
                    }
                    if (adc_value > adc_peak_value) {
                        adc_peak_value = adc_value;
                        adc_peak_index = sample;
                    }
                    if (sample < adc_peak_max_index) {
                        if (adc_value < adc_min_value) {
                            adc_min_value = adc_value;
                        }
                    }
                    // if (sample >= 5 && sample <= 8) {
                    if (sample >= adc_peak_min_index && sample <= adc_peak_max_index) {
                        if (adc_value > adc_peak_ranged_value) {
                            adc_peak_ranged_value = adc_value;
                            adc_peak_ranged_index = sample;
                        }
                    }
                    // fill ADC waveform histogram
                    h2_adc_waveforms[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(sample, adc_value);
                    h1_adc_peak_position->Fill(adc_peak_index);
                } // end of sample loop
                // calculate the pedestal
                // double adc_pedestal = pedestal_median_of_first3(adc_pedestal_samples);
                //double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
                double adc_pedestal = static_cast<double>(adc_min_value);
                double adc_peak_value_pede_sub = static_cast<double>(adc_peak_ranged_value) - adc_pedestal;
                // double adc_peak_value_CM_sub = adc_peak_value_pede_sub - static_cast<double>(CM_samples[adc_peak_ranged_index]);
                // if (adc_peak_value_pede_sub > global_min_dynamic_range) {
                //     adc_peak_value_pede_sub = static_cast<double>(global_min_dynamic_range);
                // }
                if (adc_peak_value_pede_sub > 0) {
                    bool do_sum_flag = true;
                    if (valid_channel_number == -1) {
                        do_sum_flag = false; // skip channels that are not in the valid channel list
                    }
                    if (valid_channel_number % 8 == 0) {
                        do_sum_flag = false; // skip the first channel of each half, which is used for common-mode calculation and usually has lower dynamic range
                    }
                    if (vldb_id == 0 && channel >= 76) {
                        do_sum_flag = false; // skip the last 4 channels of VLDB 0, which are known to be noisy and have lower dynamic range
                    }
                    if (do_sum_flag) {
                        h2_adc_peak_channel_correlation->Fill(channel, adc_peak_value_pede_sub);
                        adc_sum += adc_peak_value_pede_sub;
                    }
                }

                // if the channel is the example channel, store the samples for later plotting
                if (channel  + vldb_id * FPGA_CHANNEL_NUMBER == example_channel) {
                    for (int sample = 0; sample < machine_gun_samples; sample++) {
                        example_channel_samples[entry][sample] = adc_raw_samples[sample];
                    }
                }
                // if (adc_peak_value_CM_sub > 0) {
                //     h2_adc_peak_channel_correlation->Fill(channel, adc_peak_value_CM_sub);
                //     adc_sum += adc_peak_value_CM_sub;
                // }
            } // end of channel loop
        } // end of vldb loop
        adc_sum_list.push_back(adc_sum); // because we don't know the max adc sum value beforehand
    } // end of event loop
    input_root->Close();

    spdlog::info("Finished processing {} events", entry_max);

    // * === Plot to output file =============================================================
    // * =====================================================================================
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
    
    // draw the accumulated ADC waveforms for example channel
    auto canvas_example_channel = new TCanvas("canvas_example_channel", "Example Channel ADC Waveforms", 800, 600);
    auto th2d_example_channel_raw = new TH2D("", "Example Channel Raw ADC Waveforms;Sample Index;ADC Value", machine_gun_samples, -0.5, machine_gun_samples - 0.5, 512, 0, 1024);
    for (int entry = 0; entry < entry_max; entry++) {
        for (int sample = 0; sample < machine_gun_samples; sample++) {
            th2d_example_channel_raw->Fill(sample, example_channel_samples[entry][sample]);
        }
    }
    std::string example_canvas_info = "Channel " + std::to_string(example_channel);
    format_2d_hist_canvas(canvas_example_channel, th2d_example_channel_raw, kBlue+2, annotation_canvas_title, annotation_testbeam_title, example_canvas_info, true, false, false);
    canvas_example_channel->Print((out_pdf + "(").c_str()); // open pdf
    canvas_example_channel->Write();
    // write to a seprate pdf file
    canvas_example_channel->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_example_channel.pdf").c_str());
    canvas_example_channel->Close();

    // draw the example CM channel
    auto canvas_example_CM_channel = new TCanvas("canvas_example_CM_channel", "Example CM Channel ADC Waveforms", 800, 600);
    auto th2d_example_CM_channel = new TH2D("", "Example CM Channel ADC Waveforms;Sample Index;ADC Value", machine_gun_samples, -0.5, machine_gun_samples - 0.5,   576, -128, 1024); 
    for (int entry = 0; entry < entry_max; entry++) {
        for (int sample = 0; sample < machine_gun_samples; sample++) {
            th2d_example_CM_channel->Fill(sample, CM_example_channel_samples[entry][sample]);
        }
    }
    std::string example_CM_canvas_info = "CM Channel " + std::to_string(CM_channel_global) + ", pedestal subtracted";
    format_2d_hist_canvas(canvas_example_CM_channel, th2d_example_CM_channel, kRed+2, annotation_canvas_title, annotation_testbeam_title, example_CM_canvas_info, true, false, false);
    canvas_example_CM_channel->Print((out_pdf).c_str()); // append to pdf
    canvas_example_CM_channel->Write();
    // write to a seprate pdf file
    canvas_example_CM_channel->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_example_CM_channel.pdf").c_str());
    canvas_example_CM_channel->Close();

    // draw the example channel minus CM channel
    auto canvas_example_channel_CM_sub = new TCanvas("canvas_example_channel_CM_sub", "Example Channel minus CM Channel ADC Waveforms", 800, 600);
    auto th2d_example_channel_CM_sub = new TH2D("", "Example Channel minus CM Channel ADC Waveforms;Sample Index;ADC Value", machine_gun_samples, -0.5, machine_gun_samples - 0.5, 512, 0, 1024); 
    for (int entry = 0; entry < entry_max; entry++) {
        for (int sample = 0; sample < machine_gun_samples; sample++) {
            double CM_value = CM_example_channel_samples[entry][sample];
            double channel_value = example_channel_samples[entry][sample];
            double channel_CM_sub_value = channel_value - CM_value;
            th2d_example_channel_CM_sub->Fill(sample, channel_CM_sub_value);
        }
    }
    std::string example_channel_CM_sub_canvas_info = "Channel " + std::to_string(example_channel) + " minus CM Channel " + std::to_string(CM_channel_global);
    format_2d_hist_canvas(canvas_example_channel_CM_sub, th2d_example_channel_CM_sub, kMagenta+2, annotation_canvas_title, annotation_testbeam_title, example_channel_CM_sub_canvas_info, true, false, false);
    canvas_example_channel_CM_sub->Print((out_pdf).c_str()); // append to pdf
    canvas_example_channel_CM_sub->Write();
    // write to a seprate pdf file
    canvas_example_channel_CM_sub->SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_example_channel_CM_sub.pdf").c_str());
    canvas_example_channel_CM_sub->Close();

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}
