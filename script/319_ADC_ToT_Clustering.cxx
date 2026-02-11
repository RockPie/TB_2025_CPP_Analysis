#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "H2GCROC_Clustering.hxx"
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
    const double tot_adc_slope      = 0.667;
    const double tot_adc_intercept  = 0;
    const double tot_baseline       = 1500.0;
    const double adc_saturation_threshold = 1023.0;

    double slope_scan_range_min     = 0.5;
    double slope_scan_range_max     = 8.0;
    double slope_scan_range_step    = 0.1;

    // const double cluster_cell_pitch = 1.0; // cm
    // const double cluster_seed_threshold = 1000.0; // ADC counts
    // const double cluster_neighbor_threshold = 40.0; // ADC counts

    std::vector<double> slope_scan_values;
    for (double slope = slope_scan_range_min; slope <= slope_scan_range_max; slope += slope_scan_range_step) {
        slope_scan_values.push_back(slope);
    }

    double intercept_scan_range_min = -1500.0;
    double intercept_scan_range_max = 1000.0;
    double intercept_scan_range_step = 50.0;
    std::vector<double> intercept_scan_values;
    for (double intercept = intercept_scan_range_min; intercept <= intercept_scan_range_max; intercept += intercept_scan_range_step) {
        intercept_scan_values.push_back(intercept);
    }
    
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
    TH1D* h1_tot_first_position = new TH1D("h1_tot_first_position", "First ToT Position;Sample Index;Counts", machine_gun_samples, -0.5, machine_gun_samples - 0.5);
    h1_tot_first_position->SetDirectory(nullptr);

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
    std::vector<double> tot_first_decoded_list;     // first tot value 
    tot_first_decoded_list.reserve(entry_max);
    std::vector<double> adc_in_tot_channels_list;   // sum of adc in the channels that have tot > 0
    adc_in_tot_channels_list.reserve(entry_max);
    std::vector<int> tot_valid_channel_count_list;
    tot_valid_channel_count_list.reserve(entry_max);

    const int adc_peak_min_index = CommonParams::adc_peak_min_index;
    const int adc_peak_max_index = CommonParams::adc_peak_max_index;

    // start event loop
    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        double adc_sum = 0.0;
        double tot_first_decoded = 0.0;
        double adc_in_tot_channels = 0.0;
        int tot_valid_channel_count = 0;

        std::vector<double> adc_peak_values_current_event;
        std::vector<int> adc_peak_values_chn_number;
        adc_peak_values_current_event.reserve(vldb_number * FPGA_CHANNEL_NUMBER);
        adc_peak_values_chn_number.reserve(vldb_number * FPGA_CHANNEL_NUMBER);
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            // channel loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int channel_global = vldb_id * FPGA_CHANNEL_NUMBER + channel;
                std::vector<UInt_t> adc_pedestal_samples; // only take the first 3 samples
                int adc_peak_index = -1;
                UInt_t adc_peak_value = 0;
                int adc_peak_ranged_index = -1;
                UInt_t adc_peak_ranged_value = 0;
                adc_pedestal_samples.reserve(3);
                UInt_t first_tot_value = 0;
                int tot_first_index = -1;
                int tot_count = 0;
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    UInt_t adc_value = val0_list_pools[vldb_id][0][idx];
                    UInt_t tot_value = val1_list_pools[vldb_id][0][idx];
                    UInt_t toa_value = val2_list_pools[vldb_id][0][idx];
                    // spdlog::info("Event {} VLDB {} Channel {} Sample {}: ADC {}, ToT {}, ToA {}", entry, vldb_id, channel, sample, adc_value, tot_value, toa_value);
                    if (sample < 3) {
                        adc_pedestal_samples.push_back(adc_value);
                    }
                    if (adc_value > adc_peak_value) {
                        adc_peak_value = adc_value;
                        adc_peak_index = sample;
                    }
                    if (sample >= adc_peak_min_index && sample <= adc_peak_max_index) {
                        if (adc_value > adc_peak_ranged_value) {
                            adc_peak_ranged_value = adc_value;
                            adc_peak_ranged_index = sample;
                        }
                    }
                    if (first_tot_value == 0 && tot_value > 0) {
                        first_tot_value = tot_value;
                        tot_first_index = sample;
                    }
                    if (tot_value > 0) {
                        tot_count++;
                    }
                    // fill ADC waveform histogram
                    h2_adc_waveforms[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(sample, adc_value);
                    h1_adc_peak_position->Fill(adc_peak_index);
                    h1_tot_first_position->Fill(tot_first_index);
                } // end of sample loop
                // calculate the pedestal
                // double adc_pedestal = pedestal_median_of_first3(adc_pedestal_samples);
                double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
                double adc_peak_value_pede_sub = static_cast<double>(adc_peak_ranged_value) - adc_pedestal;

                if (adc_peak_value_pede_sub > 0) {
                    h2_adc_peak_channel_correlation->Fill(channel, adc_peak_value_pede_sub);
                    adc_sum += adc_peak_value_pede_sub;
                }
                // * when the tot is valid and adc is saturated
                if (first_tot_value > 0 && adc_peak_ranged_value >= adc_saturation_threshold) {
                    UInt_t tot_decoded = decode_tot_value(first_tot_value);
                    if (tot_decoded > tot_baseline) {
                        tot_first_decoded += static_cast<double>(tot_decoded);
                        double adc_in_channel = adc_peak_value_pede_sub;
                        adc_in_tot_channels += adc_in_channel;
                        tot_valid_channel_count++;

                        double tot_adc_equivalent = static_cast<double>(tot_decoded) * tot_adc_slope + tot_adc_intercept;
                        adc_peak_values_current_event.push_back(tot_adc_equivalent);
                        adc_peak_values_chn_number.push_back(vldb_id * FPGA_CHANNEL_NUMBER + channel);
                    } 
                    else {
                        adc_peak_values_current_event.push_back(adc_peak_value_pede_sub);
                        adc_peak_values_chn_number.push_back(vldb_id * FPGA_CHANNEL_NUMBER + channel);
                    }
                }
                else {
                    adc_peak_values_current_event.push_back(adc_peak_value_pede_sub);
                    adc_peak_values_chn_number.push_back(vldb_id * FPGA_CHANNEL_NUMBER + channel);
                }
            } // end of channel loop
        } // end of vldb loop
        // adc_sum_list.push_back(adc_sum);
        // ! -- Do the clustering here if needed -- !
        std::vector<int> channel_x_positions;
        std::vector<int> channel_y_positions;
        std::vector<double> channel_value_positions;
        mapping_chn_to_xy(
            adc_peak_values_current_event,
            adc_peak_values_chn_number,
            topo_wave,
            channel_x_positions,
            channel_y_positions,
            &channel_value_positions
        );

        double cluster_energy_sum = 0.0;
        if (true) {
            TH2D* h2_event_display = event_map(
                channel_value_positions,
                channel_x_positions,
                channel_y_positions,
                NX,
                NY
            );
            if (!h2_event_display) {
                spdlog::warn("Failed to build event map for entry {} due to inconsistent pixel data", entry);
                continue;
            }
            const std::string canvas_name = "event_display_canvas_" + std::to_string(entry);
            TCanvas event_display_canvas(canvas_name.c_str(), "Event Display Canvas", 800, 600);
            h2_event_display->GetXaxis()->SetTitle("X Position");
            h2_event_display->GetYaxis()->SetTitle("Y Position");
            std::string canvas_info = "Event " + std::to_string(entry) + " ADC Peak Map";
            format_2d_hist_canvas(&event_display_canvas, h2_event_display, kBlue+2, CANVAS_TITLE, TESTBEAM_TITLE, canvas_info, true, false, false);

            std::vector<int> cluster_x_list, cluster_y_list, cluster_id_list;
            std::vector<double> cluster_com_x, cluster_com_y;

            auto clusters = cluster_uniform_grid_topcells(
                channel_x_positions,
                channel_y_positions,
                channel_value_positions,
                NX,
                NY,
                1.0,    // pitch in cm
                500.0,  // threshold in ADC units
                20.0,   // neighbor threshold in ADC units
                false,  // use diagonal neighbors
                200.0,  // min cluster energy
                256,    // max cluster size
                1       // max cluster number
            );

            if (clusters.empty()) {
                h2_event_display->Delete();
                event_display_canvas.Close();
                continue;
            }
            for (const auto& cluster : clusters) {
                auto cluster_xs = cluster.xs;
                auto cluster_ys = cluster.ys;
                cluster_x_list.insert(cluster_x_list.end(), cluster_xs.begin(), cluster_xs.end());
                cluster_y_list.insert(cluster_y_list.end(), cluster_ys.begin(), cluster_ys.end());
                cluster_com_x.push_back(cluster.x_cm_cell);
                cluster_com_y.push_back(cluster.y_cm_cell);
                cluster_energy_sum += cluster.sumE;
                break; // only consider the first cluster
            }
    
            label_pixels(&event_display_canvas, cluster_x_list, cluster_y_list);
            label_markers(&event_display_canvas, cluster_com_x, cluster_com_y);

            if (entry < 10) {
                output_root->cd();
                event_display_canvas.Write();
            }
            h2_event_display->Delete();
            event_display_canvas.Close();

        } // end of first 10 events
        if (cluster_energy_sum > 0.0) {
            adc_sum_list.push_back(cluster_energy_sum); // because we don't know the max adc sum value beforehand
        }

        tot_first_decoded_list.push_back(tot_first_decoded);
        adc_in_tot_channels_list.push_back(adc_in_tot_channels);
        tot_valid_channel_count_list.push_back(tot_valid_channel_count);
    } // end of event loop
    input_root->Close();

    spdlog::info("Finished processing {} events", entry_max);

    // * === Plot to output file =============================================================
    // * =====================================================================================
    // const int NX = 16, NY = 12;
    // const int board_cols = 8, board_rows = 4;

    // auto chan2pad = build_chan2pad_LUT(
    //     vldb_number, FPGA_CHANNEL_NUMBER,
    //     NX, NY, board_cols, board_rows,
    //     sipm_board, board_loc, board_rotation, board_flip
    // );

    // MosaicTopology topo_wave;
    // topo_wave.NX = NX;
    // topo_wave.NY = NY;
    // topo_wave.vldb_number = vldb_number;
    // topo_wave.channels_per_vldb = FPGA_CHANNEL_NUMBER;
    // topo_wave.reverse_row = true;
    // topo_wave.minimalist_axis = true;
    // topo_wave.th2_logz = true;
    // topo_wave.chan2pad = chan2pad;

    // MosaicTopology topo_ped_median = topo_wave;
    // topo_ped_median = topo_wave;
    // topo_ped_median.th2_logz = true;

    // auto output_root = new TFile(script_output_file.c_str(), "RECREATE");
    // if (output_root->IsZombie()) {
    //     spdlog::error("Failed to create output file {}", script_output_file);
    //     return 1;
    // }

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    output_root->cd();
    const std::string out_pdf = script_output_file.substr(0, script_output_file.find_last_of(".")) + ".pdf";

    // --- Write to output file ------------------------------------------------------------
    TCanvas adc_waveform_canvas("adc_waveform_canvas", "ADC Waveform Canvas", 1600, 1200);
    draw_mosaic_fixed(adc_waveform_canvas, h2_adc_waveforms, topo_wave);
    output_root->cd();
    adc_waveform_canvas.Modified();
    adc_waveform_canvas.Update();
    adc_waveform_canvas.Print((out_pdf + "[").c_str()); // begin of pdf
    adc_waveform_canvas.Write();
    adc_waveform_canvas.Close();

    auto example_waveform_hist = h2_adc_waveforms[example_channel];
    TCanvas example_waveform_canvas("example_waveform_canvas", "Example Waveform Canvas", 800, 600);
    example_waveform_hist->GetXaxis()->SetTitle("Sample Index");
    example_waveform_hist->GetYaxis()->SetTitle("ADC Value");

    std::string example_canvas_info = "Channel " + std::to_string(example_channel) + " Raw Waveform";
    format_2d_hist_canvas(&example_waveform_canvas, example_waveform_hist, kBlue+2, annotation_canvas_title, annotation_testbeam_title, example_canvas_info, true, false, false);
    example_waveform_canvas.Print(out_pdf.c_str());
    // save to a separate pdf file
    example_waveform_canvas.SaveAs((script_output_file.substr(0, script_output_file.find_last_of(".")) + "_example_waveform.pdf").c_str());
    example_waveform_canvas.Write();
    example_waveform_canvas.Close();

    TCanvas adc_peak_position_canvas("adc_peak_position_canvas", "ADC Peak Position Canvas", 800, 600);
    std::string canvas_info = "ADC Peak Position Distribution";
    format_1d_hist_canvas(&adc_peak_position_canvas, h1_adc_peak_position, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
    adc_peak_position_canvas.Print(out_pdf.c_str());
    adc_peak_position_canvas.Write();
    adc_peak_position_canvas.Close();

    TCanvas tot_first_position_canvas("tot_first_position_canvas", "ToT First Position Canvas", 800, 600);
    canvas_info = "ToT First Position Distribution";
    format_1d_hist_canvas(&tot_first_position_canvas, h1_tot_first_position, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
    tot_first_position_canvas.Print(out_pdf.c_str());
    tot_first_position_canvas.Write();
    tot_first_position_canvas.Close();

    TCanvas adc_peak_channel_correlation_canvas("adc_peak_channel_correlation_canvas", "ADC Peak Channel Correlation Canvas", 800, 600);
    canvas_info = "ADC Peak Value vs Channel Correlation";
    format_2d_hist_canvas(&adc_peak_channel_correlation_canvas, h2_adc_peak_channel_correlation, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
    adc_peak_channel_correlation_canvas.Print(out_pdf.c_str());
    adc_peak_channel_correlation_canvas.Write();
    adc_peak_channel_correlation_canvas.Close();

    double mean_avg = 0.0;
    double mean_err_sys = 0.0;
    double mean_err_stat = 0.0;
    double sigma_avg = 0.0;
    double sigma_err_sys = 0.0;
    double sigma_err_stat = 0.0;
    double sigma_err = 0.0;
    double resolution = 0.0;
    double resolution_err = 0.0;

    double cb_mean_avg = 0.0;
    double cb_mean_err_sys = 0.0;
    double cb_mean_err_stat = 0.0;
    double cb_sigma_avg = 0.0;
    double cb_sigma_err_sys = 0.0;
    double cb_sigma_err_stat = 0.0;
    double cb_effective_sigma_avg = 0.0;
    double cb_effective_sigma_err_sys = 0.0;
    double cb_effective_sigma_err_stat = 0.0;
    double cb_resolution_avg = 0.0;
    double cb_resolution_err = 0.0;
    double cb_effective_resolution = 0.0;
    double cb_effective_resolution_err = 0.0;
    double cb_sigma_err = 0.0;

    double cb_tot_mean_avg = 0.0;
    double cb_tot_mean_err_sys = 0.0;
    double cb_tot_mean_err_stat = 0.0;
    double cb_tot_sigma_avg = 0.0;
    double cb_tot_sigma_err_sys = 0.0;
    double cb_tot_sigma_err_stat = 0.0;
    double cb_tot_effective_sigma_avg = 0.0;
    double cb_tot_effective_sigma_err_sys = 0.0;
    double cb_tot_effective_sigma_err_stat = 0.0;
    double cb_tot_resolution_avg = 0.0;
    double cb_tot_resolution_err = 0.0;
    double cb_tot_effective_resolution = 0.0;

    auto hist_color = kMagenta+2;

    // ! -- Raw ADC Sum distribution Histogram ---
    // ! =====================================================================================
    if (adc_sum_list.empty()) {
        spdlog::warn("ADC sum list is empty, skipping ADC/ToT summary plots");
    } else {
        TCanvas adc_sum_distribution_canvas("adc_sum_distribution_canvas", "ADC Sum Distribution Canvas", 800, 600);
        canvas_info = "ADC Sum Distribution, Run " + script_input_run_number;
    // determine 90% percentile max value for better visualization
        auto sorted_adc_sum_list = adc_sum_list;
        std::sort(sorted_adc_sum_list.begin(), sorted_adc_sum_list.end());
        double adc_sum_90pct_max = sorted_adc_sum_list[static_cast<size_t>(0.9 * sorted_adc_sum_list.size()) - 1];
    // create histogram and fill
    double x_hist_max = adc_sum_90pct_max * 3.6;
    double x_hist_min = 0.0;
    const double bin_width = 300.0;
    int n_bins = static_cast<int>((x_hist_max - x_hist_min) / bin_width);
    if (n_bins < 100) n_bins = 100; // at least 100 bins
    TH1D* h1_adc_sum_distribution = new TH1D("h1_adc_sum_distribution", "ADC Sum Distribution;ADC Sum;Counts", n_bins, x_hist_min, x_hist_max);
    for (const auto& adc_sum_value : adc_sum_list) {
        if (adc_sum_value <= 0.0) continue; // skip invalid values
        h1_adc_sum_distribution->Fill(adc_sum_value);
    }
    format_1d_hist_canvas(&adc_sum_distribution_canvas, h1_adc_sum_distribution, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
    // add additional statistics
    double adc_sum_mean = h1_adc_sum_distribution->GetMean();
    double adc_sum_entry = h1_adc_sum_distribution->GetEntries();
    double adc_sum_rms = h1_adc_sum_distribution->GetRMS();
    TPaveText *pave_stats = new TPaveText(0.6,0.7,0.88,0.88,"NDC");
    pave_stats->SetFillColorAlpha(kWhite,0.0);
    pave_stats->SetBorderSize(0);
    // right align text
    pave_stats->SetTextAlign(31);
    pave_stats->SetTextFont(42);
    pave_stats->SetTextSize(0.03);
    pave_stats->SetTextColor(kBlue+2);
    pave_stats->AddText(("Entries: " + std::to_string(static_cast<int>(adc_sum_entry))).c_str());
    pave_stats->AddText(("Mean: " + format_decimal(adc_sum_mean)).c_str());
    pave_stats->AddText(("RMS: " + format_decimal(adc_sum_rms)).c_str());
    pave_stats->Draw();

    adc_sum_distribution_canvas.Print(out_pdf.c_str());
    adc_sum_distribution_canvas.Write();
    adc_sum_distribution_canvas.Close();

    // ! -- Fitted ADC Sum distribution Histogram ---
    TCanvas adc_sum_distribution_fitted_canvas("adc_sum_distribution_fitted_canvas", "ADC Sum Distribution Fitted Canvas", 800, 600);
    canvas_info = "Fitted ADC Sum Distribution, Run " + script_input_run_number;
    TH1D* h1_adc_sum_distribution_fitted = (TH1D*) h1_adc_sum_distribution->Clone("h1_adc_sum_distribution_fitted");
    format_1d_hist_canvas(&adc_sum_distribution_fitted_canvas, h1_adc_sum_distribution_fitted, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
    // pre-fit for initial values
    double fit_min = adc_sum_mean - 2 * adc_sum_rms;
    double fit_max = adc_sum_mean + 3 * adc_sum_rms;
    if (fit_min < 0) fit_min = 0;
    double fit_amp_init = h1_adc_sum_distribution_fitted->GetMaximum();
    TF1 *gaus_fit_pre = new TF1("gaus_fit_pre", "gaus", fit_min, fit_max);
    // set initial values
    gaus_fit_pre->SetParameter(0, fit_amp_init);
    gaus_fit_pre->SetParameter(1, adc_sum_mean);
    gaus_fit_pre->SetParameter(2, adc_sum_rms);
    h1_adc_sum_distribution_fitted->Fit(gaus_fit_pre, "RQ");
    double pre_fit_amp = gaus_fit_pre->GetParameter(0);
    double pre_fit_mean = gaus_fit_pre->GetParameter(1);
    double pre_fit_sigma = gaus_fit_pre->GetParameter(2);
    double pre_fit_chi2 = gaus_fit_pre->GetChisquare();
    double pre_fit_ndf = gaus_fit_pre->GetNDF();
    // draw fitted function
    gaus_fit_pre->SetLineColor(kRed);
    gaus_fit_pre->SetLineWidth(2);
    gaus_fit_pre->Draw("same");

    // do the second round of fit with initial values from pre-fit
    fit_min = pre_fit_mean - 1.5 * pre_fit_sigma;
    fit_max = pre_fit_mean + 2.5 * pre_fit_sigma;
    if (fit_min < 0) fit_min = 0;
    TF1 *gaus_fit_final = new TF1("gaus_fit_final", "gaus", fit_min, fit_max);
    // set initial values
    gaus_fit_final->SetParameter(0, pre_fit_amp);
    gaus_fit_final->SetParameter(1, pre_fit_mean);
    gaus_fit_final->SetParameter(2, pre_fit_sigma);
    h1_adc_sum_distribution_fitted->Fit(gaus_fit_final, "RQ+");

    double pre_fit_2_amp = gaus_fit_final->GetParameter(0);
    double pre_fit_2_mean = gaus_fit_final->GetParameter(1);
    double pre_fit_2_sigma = gaus_fit_final->GetParameter(2);
    double pre_fit_2_chi2 = gaus_fit_final->GetChisquare();
    double pre_fit_2_ndf = gaus_fit_final->GetNDF();
    // draw final fitted function
    gaus_fit_final->SetLineColor(kGreen+2);
    gaus_fit_final->SetLineWidth(2);
    gaus_fit_final->Draw("same");

    // add additional statistics
    TPaveText *pave_stats_pre = new TPaveText(0.55,0.73,0.90,0.88,"NDC");
    pave_stats_pre->SetFillColorAlpha(kWhite,0.0);
    pave_stats_pre->SetBorderSize(0);
    pave_stats_pre->SetTextAlign(31);
    pave_stats_pre->SetTextFont(42);
    pave_stats_pre->SetTextSize(0.03);
    pave_stats_pre->SetTextColor(kRed);
    pave_stats_pre->AddText("Pre-Fit:");
    pave_stats_pre->AddText(("  Mean: " + format_decimal(pre_fit_mean)).c_str());
    pave_stats_pre->AddText(("  Sigma: " + format_decimal(pre_fit_sigma)).c_str());
    pave_stats_pre->AddText(("  Chi2/NDF: " + format_decimal(pre_fit_chi2) + "/" + std::to_string(static_cast<int>(pre_fit_ndf))).c_str());
    pave_stats_pre->Draw();

    TPaveText *pave_stats_pre_2 = new TPaveText(0.55,0.58,0.90,0.73,"NDC");
    pave_stats_pre_2->SetFillColorAlpha(kWhite,0.0);
    pave_stats_pre_2->SetBorderSize(0);
    pave_stats_pre_2->SetTextAlign(31);
    pave_stats_pre_2->SetTextFont(42);
    pave_stats_pre_2->SetTextSize(0.03);
    pave_stats_pre_2->SetTextColor(kGreen+2);
    pave_stats_pre_2->AddText("Pre-Fit 2:");
    pave_stats_pre_2->AddText(("  Mean: " + format_decimal(pre_fit_2_mean)).c_str());
    pave_stats_pre_2->AddText(("  Sigma: " + format_decimal(pre_fit_2_sigma)).c_str());
    pave_stats_pre_2->AddText(("  Chi2/NDF: " + format_decimal(pre_fit_2_chi2) + "/" + std::to_string(static_cast<int>(pre_fit_2_ndf))).c_str());
    pave_stats_pre_2->Draw();

    adc_sum_distribution_fitted_canvas.Modified();
    adc_sum_distribution_fitted_canvas.Update();
    adc_sum_distribution_fitted_canvas.Print(out_pdf.c_str());
    adc_sum_distribution_fitted_canvas.Write();
    adc_sum_distribution_fitted_canvas.Close();

    // ! --- Do the real gaussian fit ---
    TCanvas adc_sum_distribution_gausfit_canvas("adc_sum_distribution_gausfit_canvas", "ADC Sum Distribution Gaus Fit Canvas", 800, 600);
    canvas_info = "Gaussian Fitted ADC, Run " + script_input_run_number;
    TH1D* h1_adc_sum_distribution_gausfit = (TH1D*) h1_adc_sum_distribution->Clone("h1_adc_sum_distribution_gausfit");
    format_1d_hist_canvas(&adc_sum_distribution_gausfit_canvas, h1_adc_sum_distribution_gausfit, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);

    // define the fitting range
    std::vector<double> fit_range_sigmas = {2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5}; // in sigma
    std::vector<double> fit_range_offsets = {-0.4, -0.2, 0.0, 0.2, 0.4}; // in sigma
    std::vector<double> fit_results_means;
    std::vector<double> fit_results_mean_errs;
    std::vector<double> fit_results_sigmas;
    std::vector<double> fit_results_sigma_errs;
    std::vector<double> fit_results_chi2s;
    std::vector<double> fit_results_ndfs;

    for (const auto& range_sigma : fit_range_sigmas) {
        for (const auto& range_offset : fit_range_offsets) {
            double fit_min = pre_fit_2_mean - range_sigma * pre_fit_2_sigma + range_offset * pre_fit_2_sigma;
            double fit_max = pre_fit_2_mean + range_sigma * pre_fit_2_sigma + range_offset * pre_fit_2_sigma;
            if (fit_min < 0) fit_min = 0;
            TF1 *gaus_fit_real = new TF1("gaus_fit_real", "gaus", fit_min, fit_max);
            // set initial values
            gaus_fit_real->SetParameter(0, pre_fit_2_amp);
            gaus_fit_real->SetParameter(1, pre_fit_2_mean);
            gaus_fit_real->SetParameter(2, pre_fit_2_sigma);
            h1_adc_sum_distribution_gausfit->Fit(gaus_fit_real, "RNQ+");

            double fit_mean = gaus_fit_real->GetParameter(1);
            double fit_mean_err = gaus_fit_real->GetParError(1);
            double fit_sigma = gaus_fit_real->GetParameter(2);
            double fit_sigma_err = gaus_fit_real->GetParError(2);

            double fit_chi2 = gaus_fit_real->GetChisquare();
            double fit_ndf = gaus_fit_real->GetNDF();

            fit_results_means.push_back(fit_mean);
            fit_results_mean_errs.push_back(fit_mean_err);
            fit_results_sigmas.push_back(fit_sigma);
            fit_results_sigma_errs.push_back(fit_sigma_err);
            fit_results_chi2s.push_back(fit_chi2);
            fit_results_ndfs.push_back(fit_ndf);

            gaus_fit_real->SetLineColorAlpha(kCyan+2, 0.2);
            gaus_fit_real->SetLineWidth(2);
            gaus_fit_real->Draw("same");

            spdlog::info("Gaus Fit Range Sigma: {}, Offset: {} => Mean: {}, Sigma: {}, Chi2/NDF: {}/{}", 
                range_sigma, range_offset, fit_mean, fit_sigma, fit_chi2, fit_ndf
            );
        }
    }

    spdlog::info("Calculating mean and sigma from multiple fit results...");
    mean_sigma_list_calculator(
        fit_results_means,
        fit_results_mean_errs,
        fit_results_sigmas,
        fit_results_sigma_errs,
        mean_avg,
        mean_err_sys,
        mean_err_stat,
        sigma_avg,
        sigma_err_sys,
        sigma_err_stat
    );
    spdlog::info("Mean: {} +/- {} (stat) +/- {} (sys)", mean_avg, mean_err_stat, mean_err_sys);
    spdlog::info("Sigma: {} +/- {} (stat) +/- {} (sys)", sigma_avg, sigma_err_stat, sigma_err_sys);

    resolution = (sigma_avg / mean_avg) * 100.0;
    double mean_err = std::sqrt(mean_err_stat * mean_err_stat + mean_err_sys * mean_err_sys);
    sigma_err = std::sqrt(sigma_err_stat * sigma_err_stat + sigma_err_sys * sigma_err_sys);
    resolution_err = resolution * std::sqrt( (sigma_err / sigma_avg) * (sigma_err / sigma_avg) + (mean_err / mean_avg) * (mean_err / mean_avg) );
    spdlog::info("Calculated Resolution: {} % +/- {} %", resolution, resolution_err);

    // add additional statistics
    TPaveText *pave_stats_gausfit = new TPaveText(0.55,0.6,0.90,0.88,"NDC");
    pave_stats_gausfit->SetFillColorAlpha(kWhite,0.0);
    pave_stats_gausfit->SetBorderSize(0);
    pave_stats_gausfit->SetTextAlign(31);
    pave_stats_gausfit->SetTextFont(42);
    pave_stats_gausfit->SetTextSize(0.03);
    pave_stats_gausfit->SetTextColor(kCyan+2);
    pave_stats_gausfit->AddText("Gaussian Fit Results:");
    
    std::string mean_str = format_decimal(mean_avg);
    std::string mean_err_stat_str = format_decimal(mean_err_stat);
    std::string mean_err_sys_str = format_decimal(mean_err_sys);
    std::string sigma_str = format_decimal(sigma_avg);
    std::string sigma_err_stat_str = format_decimal(sigma_err_stat);
    std::string sigma_err_sys_str = format_decimal(sigma_err_sys);
    
    std::string mean_text = "  Mean: " + mean_str + " #pm " + mean_err_stat_str + " (stat)";
    pave_stats_gausfit->AddText(mean_text.c_str());
    std::string mean_sys_text = "        #pm " + mean_err_sys_str + " (sys)";
    pave_stats_gausfit->AddText(mean_sys_text.c_str());

    std::string sigma_text = "  Sigma: " + sigma_str + " #pm " + sigma_err_stat_str + " (stat)";
    pave_stats_gausfit->AddText(sigma_text.c_str());
    std::string sigma_sys_text = "         #pm " + sigma_err_sys_str + " (sys)";
    pave_stats_gausfit->AddText(sigma_sys_text.c_str());

    std::string resolution_str = format_decimal(resolution);
    std::string resolution_err_str = format_decimal(resolution_err);
    std::string resolution_text = "  Resolution: " + resolution_str + " % #pm " + resolution_err_str + " %";
    pave_stats_gausfit->AddText(resolution_text.c_str());
    pave_stats_gausfit->Draw();

    adc_sum_distribution_gausfit_canvas.Modified();
    adc_sum_distribution_gausfit_canvas.Update();
    adc_sum_distribution_gausfit_canvas.Print(out_pdf.c_str());
    adc_sum_distribution_gausfit_canvas.Write();
    adc_sum_distribution_gausfit_canvas.Close();
    
     // ! --- Do the real crystalball fit ---
    TCanvas adc_sum_distribution_cbfit_canvas("adc_sum_distribution_cbfit_canvas", "ADC Sum Distribution CB Fit Canvas", 800, 600);
    canvas_info = "Crystal Ball Fitted ADC, Run " + script_input_run_number;
    TH1D* h1_adc_sum_distribution_cbfit = (TH1D*) h1_adc_sum_distribution->Clone("h1_adc_sum_distribution_cbfit");
    format_1d_hist_canvas(&adc_sum_distribution_cbfit_canvas, h1_adc_sum_distribution_cbfit, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);

    crystalball_fit_th1d(adc_sum_distribution_cbfit_canvas, *h1_adc_sum_distribution_cbfit, fit_range_sigmas, fit_range_offsets, hist_color, cb_mean_avg, cb_mean_err_sys, cb_mean_err_stat, cb_sigma_avg, cb_sigma_err_sys, cb_sigma_err_stat, cb_resolution_avg, cb_resolution_err);

    cb_effective_resolution = (cb_effective_sigma_avg / cb_mean_avg) * 100.0;
    cb_effective_resolution_err = cb_effective_resolution * std::sqrt( (cb_effective_sigma_err_stat / cb_effective_sigma_avg) * (cb_effective_sigma_err_stat / cb_effective_sigma_avg) + (cb_mean_err_stat / cb_mean_avg) * (cb_mean_err_stat / cb_mean_avg) );
    cb_sigma_err = std::sqrt(cb_sigma_err_stat * cb_sigma_err_stat + cb_sigma_err_sys * cb_sigma_err_sys);

    // add additional statistics
    TPaveText *pave_stats_cbfit = new TPaveText(0.55,0.4,0.90,0.88,"NDC");
    pave_stats_cbfit->SetFillColorAlpha(kWhite,0.0);
    pave_stats_cbfit->SetBorderSize(0);
    pave_stats_cbfit->SetTextAlign(31);
    pave_stats_cbfit->SetTextFont(42);
    pave_stats_cbfit->SetTextSize(0.03);
    pave_stats_cbfit->SetTextColor(kMagenta+2);
    pave_stats_cbfit->AddText("Crystal Ball Fit Results:");
    std::string cb_mean_str = format_decimal(cb_mean_avg);
    std::string cb_mean_err_stat_str = format_decimal(cb_mean_err_stat);
    std::string cb_mean_err_sys_str = format_decimal(cb_mean_err_sys);
    std::string cb_sigma_str = format_decimal(cb_sigma_avg);
    std::string cb_sigma_err_stat_str = format_decimal(cb_sigma_err_stat);
    std::string cb_sigma_err_sys_str = format_decimal(cb_sigma_err_sys);
    std::string cb_mean_text = "  Mean: " + cb_mean_str + " #pm " + cb_mean_err_stat_str + " (stat)";
    pave_stats_cbfit->AddText(cb_mean_text.c_str());
    std::string cb_mean_sys_text = "        #pm " + cb_mean_err_sys_str + " (sys)";
    pave_stats_cbfit->AddText(cb_mean_sys_text.c_str());
    std::string cb_sigma_text = "  Sigma: " + cb_sigma_str + " #pm " + cb_sigma_err_stat_str + " (stat)";
    pave_stats_cbfit->AddText(cb_sigma_text.c_str());
    std::string cb_sigma_sys_text = "         #pm " + cb_sigma_err_sys_str + " (sys)";
    pave_stats_cbfit->AddText(cb_sigma_sys_text.c_str());
    std::string cb_resolution_str = format_decimal(cb_resolution_avg);
    std::string cb_resolution_err_str = format_decimal(cb_resolution_err);
    std::string cb_resolution_text = "  Core Resolution: " + cb_resolution_str + " % #pm " + cb_resolution_err_str + " %";
    pave_stats_cbfit->AddText(cb_resolution_text.c_str());
    std::string cb_effective_sigma_str = format_decimal(cb_effective_sigma_avg);
    std::string cb_effective_sigma_err_stat_str = format_decimal(cb_effective_sigma_err_stat);
    std::string cb_effective_sigma_err_sys_str = format_decimal(cb_effective_sigma_err_sys);
    std::string cb_effective_sigma_text = "  Eff Sigma: " + cb_effective_sigma_str;
    pave_stats_cbfit->AddText(cb_effective_sigma_text.c_str());
    std::string cb_effective_resolution_str = format_decimal(cb_effective_resolution);
    std::string cb_effective_resolution_text = "  Eff Resolution: " + cb_effective_resolution_str + " %";
    pave_stats_cbfit->AddText(cb_effective_resolution_text.c_str());
    pave_stats_cbfit->Draw();

    adc_sum_distribution_cbfit_canvas.Modified();
    adc_sum_distribution_cbfit_canvas.Update();
    adc_sum_distribution_cbfit_canvas.Print(out_pdf.c_str());
    std::string cbfit_pic_file_name = script_output_file.substr(0, script_output_file.find_last_of(".")) + "_cbfit.pdf";
    adc_sum_distribution_cbfit_canvas.SaveAs(cbfit_pic_file_name.c_str());
    adc_sum_distribution_cbfit_canvas.Write();
    adc_sum_distribution_cbfit_canvas.Close();
    // ! =====================================================================================

    // ! do the same crystal ball fit but for tot and adc combined
    TCanvas adc_in_tot_channels_distribution_cbfit_canvas("adc_in_tot_channels_distribution_cbfit_canvas", "ADC in ToT Channels Distribution CB Fit Canvas", 800, 600);
    canvas_info = "Crystal Ball Fitted ADC and ToT Combined, Run " + script_input_run_number;
    // create histogram and fill
    double x_hist_max_adc_in_tot = adc_sum_90pct_max * 3.6;
    double x_hist_min_adc_in_tot = 0.0;
    TH1D* h1_adc_tot_combined = new TH1D("h1_adc_tot_combined", "ADC in ToT Channels Distribution;ADC Sum in ToT Channels;Counts", n_bins, x_hist_min_adc_in_tot, x_hist_max_adc_in_tot);
    TH1D* h1_adc_sums_compared = new TH1D("h1_adc_sums_compared", "ADC Sum Comparison;ADC Sum;Counts", n_bins, x_hist_min, x_hist_max);
    for (int i = 0; i < adc_sum_list.size(); i++) {
        auto adc_in_tot_value = adc_in_tot_channels_list[i];
        auto adc_sum_value = adc_sum_list[i];
        auto tot_value = tot_first_decoded_list[i];
        auto tot_valid_channel_count = tot_valid_channel_count_list[i];
        double adc_tot_combined = double(adc_sum_value) - double(adc_in_tot_value) + double(tot_value) * tot_adc_slope + tot_adc_intercept * tot_valid_channel_count;
        h1_adc_tot_combined->Fill(adc_tot_combined);
        h1_adc_sums_compared->Fill(adc_sum_value);
    }

    // rescale the y axis
    double y_max = h1_adc_tot_combined->GetMaximum();
    if (h1_adc_sums_compared->GetMaximum() > y_max) {
        y_max = h1_adc_sums_compared->GetMaximum();
    }
    h1_adc_sums_compared->SetMaximum(y_max * 1.3);
    h1_adc_tot_combined->SetMaximum(y_max * 1.3);
    
    TH1D* h1_adc_tot_combined_cbfit = (TH1D*) h1_adc_tot_combined->Clone("h1_adc_tot_combined_cbfit");
    format_1d_hist_canvas(&adc_in_tot_channels_distribution_cbfit_canvas, h1_adc_tot_combined_cbfit, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
        crystalball_fit_th1d(adc_in_tot_channels_distribution_cbfit_canvas, *h1_adc_tot_combined_cbfit, fit_range_sigmas, fit_range_offsets, hist_color, cb_tot_mean_avg, cb_tot_mean_err_sys, cb_tot_mean_err_stat, cb_tot_sigma_avg, cb_tot_sigma_err_sys, cb_tot_sigma_err_stat, cb_tot_resolution_avg, cb_tot_resolution_err);
        if (cb_tot_mean_avg != 0.0) {
            cb_tot_effective_resolution = (cb_tot_effective_sigma_avg / cb_tot_mean_avg) * 100.0;
        } else {
            cb_tot_effective_resolution = 0.0;
        }


    // draw adc only for comparison
    TH1D* h1_adc_sums_compared_cbfit = (TH1D*) h1_adc_sums_compared->Clone("h1_adc_sums_compared_cbfit");
    h1_adc_sums_compared_cbfit->SetLineColor(kGreen+2);;
    h1_adc_sums_compared_cbfit->SetLineWidth(2);
    h1_adc_sums_compared_cbfit->Draw("same");

    // add additional statistics
    TPaveText *pave_stats_cbfit_tot = new TPaveText(0.55,0.4,0.90,0.88,"NDC");
    pave_stats_cbfit_tot->SetFillColorAlpha(kWhite,0.0);
    pave_stats_cbfit_tot->SetBorderSize(0);
    pave_stats_cbfit_tot->SetTextAlign(31);
    pave_stats_cbfit_tot->SetTextFont(42);
    pave_stats_cbfit_tot->SetTextSize(0.03);    
    pave_stats_cbfit_tot->SetTextColor(kMagenta+2);
    pave_stats_cbfit_tot->AddText("Crystal Ball Fit Results (ADC + ToT):");
    std::string cb_tot_mean_str = format_decimal(cb_tot_mean_avg);
    std::string cb_tot_mean_err_stat_str = format_decimal(cb_tot_mean_err_stat);
    std::string cb_tot_mean_err_sys_str = format_decimal(cb_tot_mean_err_sys);
    std::string cb_tot_sigma_str = format_decimal(cb_tot_sigma_avg);
    std::string cb_tot_sigma_err_stat_str = format_decimal(cb_tot_sigma_err_stat);
    std::string cb_tot_sigma_err_sys_str = format_decimal(cb_tot_sigma_err_sys);
    std::string cb_tot_mean_text = "  Mean: " + cb_tot_mean_str + " #pm " + cb_tot_mean_err_stat_str + " (stat)";
    pave_stats_cbfit_tot->AddText(cb_tot_mean_text.c_str());
    std::string cb_tot_mean_sys_text = "        #pm " + cb_tot_mean_err_sys_str + " (sys)"; 
    pave_stats_cbfit_tot->AddText(cb_tot_mean_sys_text.c_str());
    std::string cb_tot_sigma_text = "  Sigma: " + cb_tot_sigma_str + " #pm " + cb_tot_sigma_err_stat_str + " (stat)";
    pave_stats_cbfit_tot->AddText(cb_tot_sigma_text.c_str());
    std::string cb_tot_sigma_sys_text = "         #pm " + cb_tot_sigma_err_sys_str + " (sys)";
    pave_stats_cbfit_tot->AddText(cb_tot_sigma_sys_text.c_str());
    std::string cb_tot_resolution_str = format_decimal(cb_tot_resolution_avg);
    std::string cb_tot_resolution_err_str = format_decimal(cb_tot_resolution_err);
    std::string cb_tot_resolution_text = "  Core Resolution: " + cb_tot_resolution_str + " % #pm " + cb_tot_resolution_err_str + " %";
    pave_stats_cbfit_tot->AddText(cb_tot_resolution_text.c_str());
    pave_stats_cbfit_tot->Draw();

        adc_in_tot_channels_distribution_cbfit_canvas.Modified();
        adc_in_tot_channels_distribution_cbfit_canvas.Update();
        adc_in_tot_channels_distribution_cbfit_canvas.Print(out_pdf.c_str());
        std::string cbfit_tot_pic_file_name = script_output_file.substr(0, script_input_file.find_last_of(".")) + "_cbfit_adc_tot.pdf";
        adc_in_tot_channels_distribution_cbfit_canvas.SaveAs(cbfit_tot_pic_file_name.c_str());
        adc_in_tot_channels_distribution_cbfit_canvas.Write();
        adc_in_tot_channels_distribution_cbfit_canvas.Close();

        // save h1_adc_tot_combined to the output root file
        output_root->cd();
        h1_adc_tot_combined->Write();
    }
    // scan over different ToT ADC slope and intercept values

    // end of pdf, create dummy canvas
    TCanvas end_canvas("end_canvas", "End Canvas", 800, 600);
    end_canvas.Print((out_pdf + "]").c_str());
    end_canvas.Close();

    // write important values for comparison
    TParameter<double> param_gaus_fit_mean("gaus_fit_mean", mean_avg); 
    TParameter<double> param_gaus_fit_sigma("gaus_fit_sigma", sigma_avg);
    TParameter<double> param_gaus_fit_resolution("gaus_fit_resolution", resolution);
    TParameter<double> param_gaus_fit_mean_err_stat("gaus_fit_mean_err_stat", mean_err_stat);
    TParameter<double> param_gaus_fit_mean_err_sys("gaus_fit_mean_err_sys", mean_err_sys);
    TParameter<double> param_gaus_fit_sigma_err_stat("gaus_fit_sigma_err_stat", sigma_err);
    TParameter<double> param_gaus_fit_sigma_err_sys("gaus_fit_sigma_err_sys", sigma_err);
    TParameter<double> param_gaus_fit_resolution_err("gaus_fit_resolution_err", resolution_err);

    TParameter<double> param_cb_fit_mean("cb_fit_mean", cb_mean_avg);
    TParameter<double> param_cb_fit_sigma("cb_fit_sigma", cb_sigma_avg);
    TParameter<double> param_cb_fit_effective_sigma("cb_fit_effective_sigma", cb_effective_sigma_avg);
    TParameter<double> param_cb_fit_resolution("cb_fit_resolution", cb_resolution_avg);
    TParameter<double> param_cb_fit_effective_resolution("cb_fit_effective_resolution", cb_effective_resolution);
    TParameter<double> param_cb_fit_mean_err_stat("cb_fit_mean_err_stat", cb_mean_err_stat);
    TParameter<double> param_cb_fit_mean_err_sys("cb_fit_mean_err_sys", cb_mean_err_sys);
    TParameter<double> param_cb_fit_sigma_err_stat("cb_fit_sigma_err_stat", cb_sigma_err);
    TParameter<double> param_cb_fit_sigma_err_sys("cb_fit_sigma_err_sys", cb_sigma_err);
    TParameter<double> param_cb_fit_effective_sigma_err_stat("cb_fit_effective_sigma_err_stat", cb_effective_sigma_err_stat);
    TParameter<double> param_cb_fit_effective_sigma_err_sys("cb_fit_effective_sigma_err_sys", cb_effective_sigma_err_sys);
    TParameter<double> param_cb_fit_resolution_err("cb_fit_resolution_err", cb_resolution_err);
    TParameter<double> param_cb_fit_effective_resolution_err("cb_fit_effective_resolution_err", cb_effective_resolution_err);

    // double cb_tot_mean_avg=0.0;
    // double cb_tot_mean_err_sys=0.0;
    // double cb_tot_mean_err_stat=0.0;
    // double cb_tot_sigma_avg=0.0;
    // double cb_tot_sigma_err_sys=0.0;
    // double cb_tot_sigma_err_stat=0.0;
    // double cb_tot_effective_sigma_avg=0.0;
    // double cb_tot_effective_sigma_err_sys=0.0;
    // double cb_tot_effective_sigma_err_stat=0.0;
    // double cb_tot_resolution_avg=0.0;
    // double cb_tot_resolution_err=0.0;
    TParameter<double> param_tot_cb_fit_mean("tot_cb_fit_mean", cb_tot_mean_avg);
    TParameter<double> param_tot_cb_fit_sigma("tot_cb_fit_sigma", cb_tot_sigma_avg);
    TParameter<double> param_tot_cb_fit_effective_sigma("tot_cb_fit_effective_sigma", cb_tot_effective_sigma_avg);
    TParameter<double> param_tot_cb_fit_resolution("tot_cb_fit_resolution", cb_tot_resolution_avg);
    TParameter<double> param_tot_cb_fit_effective_resolution("tot_cb_fit_effective_resolution", cb_tot_effective_resolution);
    TParameter<double> param_tot_cb_fit_mean_err_stat("tot_cb_fit_mean_err_stat", cb_tot_mean_err_stat);
    TParameter<double> param_tot_cb_fit_mean_err_sys("tot_cb_fit_mean_err_sys", cb_tot_mean_err_sys);
    TParameter<double> param_tot_cb_fit_sigma_err_stat("tot_cb_fit_sigma_err_stat", cb_tot_sigma_err_stat);
    TParameter<double> param_tot_cb_fit_sigma_err_sys("tot_cb_fit_sigma_err_sys", cb_tot_sigma_err_sys);
    TParameter<double> param_tot_cb_fit_effective_sigma_err_stat("tot_cb_fit_effective_sigma_err_stat", cb_tot_effective_sigma_err_stat);
    TParameter<double> param_tot_cb_fit_effective_sigma_err_sys("tot_cb_fit_effective_sigma_err_sys", cb_tot_effective_sigma_err_sys);
    TParameter<double> param_tot_cb_fit_resolution_err("tot_cb_fit_resolution_err", cb_tot_resolution_err);

    output_root->cd();
    param_gaus_fit_mean.Write();
    param_gaus_fit_sigma.Write();
    param_gaus_fit_resolution.Write();
    param_gaus_fit_mean_err_stat.Write();
    param_gaus_fit_mean_err_sys.Write();
    param_gaus_fit_sigma_err_stat.Write();
    param_gaus_fit_sigma_err_sys.Write();
    param_gaus_fit_resolution_err.Write();

    param_cb_fit_mean.Write();
    param_cb_fit_sigma.Write();
    param_cb_fit_effective_sigma.Write();
    param_cb_fit_resolution.Write();
    param_cb_fit_effective_resolution.Write();
    param_cb_fit_mean_err_stat.Write();
    param_cb_fit_mean_err_sys.Write();
    param_cb_fit_sigma_err_stat.Write();
    param_cb_fit_sigma_err_sys.Write();
    param_cb_fit_effective_sigma_err_stat.Write();
    param_cb_fit_effective_sigma_err_sys.Write();
    param_cb_fit_resolution_err.Write();
    param_cb_fit_effective_resolution_err.Write();

    param_tot_cb_fit_mean.Write();
    param_tot_cb_fit_sigma.Write();
    param_tot_cb_fit_effective_sigma.Write();
    param_tot_cb_fit_resolution.Write();
    param_tot_cb_fit_effective_resolution.Write();
    param_tot_cb_fit_mean_err_stat.Write();
    param_tot_cb_fit_mean_err_sys.Write();
    param_tot_cb_fit_sigma_err_stat.Write();
    param_tot_cb_fit_sigma_err_sys.Write();
    param_tot_cb_fit_effective_sigma_err_stat.Write();
    param_tot_cb_fit_effective_sigma_err_sys.Write();
    param_tot_cb_fit_resolution_err.Write();

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}

