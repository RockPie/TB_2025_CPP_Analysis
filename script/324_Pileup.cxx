#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "TRandom3.h"
#include "CommonParams.hxx"

double pedestal_rms_of_first3(std::vector<double>& pedestal_samples) {
    if (pedestal_samples.size() != 3) {
        spdlog::error("Pedestal samples size is {}, expected 3", pedestal_samples.size());
        return 0.0;
    }
    double average = 0.0;
    for (const auto& sample : pedestal_samples) {
        average += static_cast<double>(sample);
    }
    average /= static_cast<double>(pedestal_samples.size());
    double sq_err_sum = 0.0;
    for (const auto& sample : pedestal_samples) {
        double diff = static_cast<double>(sample) - average;
        sq_err_sum += diff * diff;
    }
    double rms = std::sqrt(sq_err_sum / static_cast<double>(pedestal_samples.size()));
    return rms;
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

    const int adc_peak_min_index = 5;
    const int adc_peak_max_index = 7;

    const int adc_value_bin_size = 4;
    const double adc_value_max   = 1024.0;
    const double adc_value_min   = 0.0;
    const int adc_value_n_bins   = static_cast<int>((adc_value_max - adc_value_min) / adc_value_bin_size);

    const double pedestal_bin_size = 2.0;
    const double pedestal_max = 256.0;
    const double pedestal_min = 0.0;
    const int pedestal_n_bins = static_cast<int>((pedestal_max - pedestal_min) / pedestal_bin_size);

    std::vector<TH1D*> h1d_pedestal_distributions_list(vldb_number * FPGA_CHANNEL_NUMBER);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int valid_channel_number = get_valid_fpga_channel(channel);
            if (valid_channel_number == -1) {
                continue;
            }
            h1d_pedestal_distributions_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH1D(Form("h1d_pedestal_distribution_vldb%d_channel%d", vldb_id, channel), Form("Pedestal Distribution VLDB %d Channel %d;Pedestal (ADC);Entries", vldb_id, channel), pedestal_n_bins, pedestal_min, pedestal_max);
            h1d_pedestal_distributions_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<TH1D*> h1d_pedestal_rms_list(vldb_number * FPGA_CHANNEL_NUMBER);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int valid_channel_number = get_valid_fpga_channel(channel);
            if (valid_channel_number == -1) {
                continue;
            }
            h1d_pedestal_rms_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH1D(Form("h1d_pedestal_rms_vldb%d_channel%d", vldb_id, channel), Form("Pedestal RMS VLDB %d Channel %d;Pedestal RMS (ADC);Entries", vldb_id, channel), 50, 0.0, 25.0);
            h1d_pedestal_rms_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<TH2D*> h2d_adc_waveform_list(vldb_number * FPGA_CHANNEL_NUMBER);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int valid_channel_number = get_valid_fpga_channel(channel);
            if (valid_channel_number == -1) {
                continue;
            }
            h2d_adc_waveform_list[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH2D(Form("h2d_adc_waveform_vldb%d_channel%d", vldb_id, channel), Form("ADC Waveform VLDB %d Channel %d;Sample Index;ADC Value", vldb_id, channel), machine_gun_samples, -0.5, machine_gun_samples - 0.5, adc_value_n_bins, adc_value_min, adc_value_max);
            h2d_adc_waveform_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    std::vector<std::vector<std::vector<double>>> adc_samples_in_CM_channels_all_vldb(vldb_number, std::vector<std::vector<double>>(4, std::vector<double>(machine_gun_samples, 0.0)));
    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        if (entry % 5000 == 0) {
            spdlog::info("Processing entry {}/{}", entry, entry_max);
        }
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
            } // end of crosstalk correction loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int valid_channel_number = get_valid_fpga_channel(channel);
                int channel_half_index = (channel / (FPGA_CHANNEL_NUMBER / 4)) % 4;
                if (entry < 10) {
                    // spdlog::info("Entry {} VLDB {} Channel {}: Valid channel number {}, Channel half index {}", entry, vldb_id, channel, valid_channel_number, channel_half_index);
                }
                if (valid_channel_number == -1) {
                    continue;
                } // skip invalid channels

                std::vector<double> CM_samples = adc_samples_in_CM_channels[channel_half_index];
                std::vector<double> adc_pedestal_samples; // only take the first 3 samples
                adc_pedestal_samples.reserve(3);
                double  adc_peak_value = 0.0;
                int     adc_peak_index = -1;
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    UInt_t adc_raw_value = val0_list_pools[vldb_id][0][idx];
                    UInt_t tot_raw_value = val1_list_pools[vldb_id][0][idx];
                    UInt_t toa_raw_value = val2_list_pools[vldb_id][0][idx];

                    double adc_value = static_cast<double>(adc_raw_value);
                    double adc_value_CM_subtracted = adc_value - CM_samples[sample];

                    if (sample < 3) adc_pedestal_samples.push_back(adc_value_CM_subtracted);
                    if (adc_value_CM_subtracted > adc_peak_value && sample >= adc_peak_min_index && sample <= adc_peak_max_index) {
                        adc_peak_value = adc_value_CM_subtracted;
                        adc_peak_index = sample;
                    }
                } // end of sample loop
                double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
                double adc_pedestal_rms = pedestal_rms_of_first3(adc_pedestal_samples);

                double adc_peak_value_pede_sub = adc_peak_value - adc_pedestal;

                h1d_pedestal_distributions_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(adc_pedestal);
                h1d_pedestal_rms_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(adc_pedestal_rms);
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    UInt_t adc_raw_value = val0_list_pools[vldb_id][0][idx];
                    double adc_value = static_cast<double>(adc_raw_value);
                    double adc_value_CM_subtracted = adc_value - CM_samples[sample];
                    h2d_adc_waveform_list[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(sample, adc_value_CM_subtracted);
                }
            } // end of channel loop
        } // end of VLDB loop
    } // end of event loop

    input_root->Close();

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

    // draw the raw pedestal distribution for each channel
    TCanvas canvas_pedestal_distribution("canvas_pedestal_distribution", "Pedestal Distribution Canvas", 1200, 800);
    canvas_pedestal_distribution.SetLogy();
    draw_mosaic_fixed(canvas_pedestal_distribution, h1d_pedestal_distributions_list, topo_ped_median);
    canvas_pedestal_distribution.Modified();
    canvas_pedestal_distribution.Update();
    canvas_pedestal_distribution.Write();
    canvas_pedestal_distribution.Close();

    // draw the pedestal RMS distribution for each channel
    TCanvas canvas_pedestal_rms("canvas_pedestal_rms", "Pedestal RMS Canvas", 1200, 800);
    canvas_pedestal_rms.SetLogy();
    draw_mosaic_fixed(canvas_pedestal_rms, h1d_pedestal_rms_list, topo_ped_median);
    canvas_pedestal_rms.Modified();
    canvas_pedestal_rms.Update();
    canvas_pedestal_rms.Write();
    canvas_pedestal_rms.Close();

    // draw the ADC waveform for each channel
    TCanvas canvas_adc_waveform("canvas_adc_waveform", "ADC Waveform Canvas", 1200, 800);
    draw_mosaic_fixed(canvas_adc_waveform, h2d_adc_waveform_list, topo_wave);
    canvas_adc_waveform.Modified();
    canvas_adc_waveform.Update();
    canvas_adc_waveform.Write();
    canvas_adc_waveform.Close();

    // // draw the channel map of all adc peak distributions for quick visualization
    // TCanvas canvas_adc_peak_distribution("canvas_adc_peak_distribution", "ADC Peak Distribution Canvas", 1200, 800);
    // draw_mosaic_fixed(canvas_adc_peak_distribution, h1_adc_peak_distribution_channels, topo_ped_median);
    // canvas_adc_peak_distribution.Modified();
    // canvas_adc_peak_distribution.Update();
    // canvas_adc_peak_distribution.Write();
    // canvas_adc_peak_distribution.Close();


    output_root->Close();

    return 0;
}