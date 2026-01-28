#include <fstream>
#include <limits>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
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

    std::string mapping_json_file = "config/mapping_Oct2025.json";
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

    spdlog::info("Script name: {}", script_name);
    spdlog::info("Input file: {}", script_input_file);
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
    int machine_gun_samples = 16;
    int vldb_number = 2;
    int chn_example = CommonParams::example_channel;

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

    size_t num_single_ToA = 0;
    size_t num_multi_ToA = 0;

    TH1I *h1i_toa_frequency = new TH1I("h1i_toa_frequency", "h1i_toa_frequency;ToA (25 ns bins);Counts", machine_gun_samples, 0, machine_gun_samples);
    h1i_toa_frequency->SetDirectory(nullptr);

    // * --- Raw ToA ns histograms ----------------------------------------------------------
    std::vector<TH1D*> h1d_raw_toa_ns_list;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h1d_raw_toa_ns_ch" + std::to_string(i);
        TH1D *h1d_raw_toa_ns = new TH1D(hist_name.c_str(), (hist_name + ";ToA (ns);Counts").c_str(), 256, 0, 25.0*static_cast<double>(machine_gun_samples));
        h1d_raw_toa_ns->SetDirectory(nullptr);
        h1d_raw_toa_ns_list.push_back(h1d_raw_toa_ns);
    }

    // * --- ToA ns v.s. ADC max correction histograms --------------------------------------
    std::vector<TH2D*> h2d_toa_adc_max_corr_list;;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h2d_toa_adc_max_corr_ch" + std::to_string(i);
        TH2D *h2d_toa_adc_max_corr = new TH2D(hist_name.c_str(), (hist_name + ";Corrected ToA (ns);ADC Peak Value (pedestal subtracted)").c_str(),
                                             256, 0, 1024,
                                             256, 0, 25.0*static_cast<double>(machine_gun_samples));
        h2d_toa_adc_max_corr->SetDirectory(nullptr);
        h2d_toa_adc_max_corr_list.push_back(h2d_toa_adc_max_corr);
    }

    // * --- ToA ns v.s. ToA code correction histograms --------------------------------------
    std::vector<TH2D*> h2d_toa_code_corr_list;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h2d_toa_code_corr_ch" + std::to_string(i);
        TH2D *h2d_toa_code_corr = new TH2D(hist_name.c_str(), (hist_name + ";Corrected ToA (ns);ToA Code").c_str(),
                                           256, 0, 1024,
                                           256, 0, 25.0*static_cast<double>(machine_gun_samples));
        h2d_toa_code_corr->SetDirectory(nullptr);
        h2d_toa_code_corr_list.push_back(h2d_toa_code_corr);
    }

    std::vector<std::vector<UInt_t>> toa_code_matrix(FPGA_CHANNEL_NUMBER * vldb_number, std::vector<UInt_t>());
    for (auto & vec : toa_code_matrix) {
        vec.reserve(entry_max);
    }
    std::vector<std::vector<double>> toa_ns_matrix(FPGA_CHANNEL_NUMBER * vldb_number, std::vector<double>());
    for (auto & vec : toa_ns_matrix) {
        vec.reserve(entry_max);
    }

    const int adc_peak_min_index = CommonParams::adc_peak_min_index;
    const int adc_peak_max_index = CommonParams::adc_peak_max_index;

    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        double adc_sum = 0.0;
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            // channel loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int channel_index_in_h2g = channel % 76;
                int channel_index_in_h2g_half = channel_index_in_h2g % 38;
                int asic_id = channel / 76;
                int half_id = channel_index_in_h2g / 38;

                std::vector<UInt_t> adc_pedestal_samples; // only take the first 3 samples
                int adc_peak_index = -1;
                UInt_t adc_peak_value = 0;
                int adc_peak_ranged_index = -1;
                UInt_t adc_peak_ranged_value = 0;
                adc_pedestal_samples.reserve(3);

                int toa_showup_times = 0;
                int tot_showup_times = 0;

                UInt_t toa_first_value = 0;
                int toa_first_index = -1;

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
                    if (sample >= adc_peak_min_index && sample <= adc_peak_max_index) {
                        if (adc_value > adc_peak_ranged_value) {
                            adc_peak_ranged_value = adc_value;
                            adc_peak_ranged_index = sample;
                        }
                    }
                    if (toa_value > 0) {
                        toa_showup_times++;
                        if (toa_first_index == -1) {
                            toa_first_index = sample;
                            toa_first_value = toa_value;
                        }
                    }
                    if (tot_value > 0) {
                        tot_showup_times++;
                    }
                } // end of sample loop
                // double adc_pedestal = pedestal_median_of_first3(adc_pedestal_samples);
                double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
                double adc_peak_value_pede_sub = static_cast<double>(adc_peak_ranged_value) - adc_pedestal;
                adc_sum += adc_peak_value_pede_sub;

                h1i_toa_frequency->Fill(tot_showup_times);
                if (toa_showup_times >= 1) {
                    num_single_ToA++;
                    // only one ToA hit
                    // apply correction here if needed
                    double toa_ns = decode_toa_value_ns(toa_first_value);
                    toa_ns += 25.0 * static_cast<double>(toa_first_index);
                    h1d_raw_toa_ns_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_ns);
                    h2d_toa_adc_max_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(adc_peak_value_pede_sub, toa_ns);
                    h2d_toa_code_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_first_value, toa_ns);
                    toa_code_matrix[channel + vldb_id * FPGA_CHANNEL_NUMBER].push_back(toa_first_value);
                    toa_ns_matrix[channel + vldb_id * FPGA_CHANNEL_NUMBER].push_back(toa_ns);
                }
                // } else if (toa_showup_times > 1) {
                //     num_multi_ToA++;
                    // multiple ToA hits
                    // currently, just take the first one
                    // double toa_ns = decode_toa_value_ns(toa_first_value);
                    // toa_ns += 25.0 * static_cast<double>(toa_first_index);
                    // h1d_raw_toa_ns_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_ns);
                // }
            } // end of channel loop
        } // end of vldb loop
    } // end of event loop

    input_root->Close();

    spdlog::info("Number of single ToA hits: {}", num_single_ToA);
    spdlog::info("Number of multiple ToA hits: {}", num_multi_ToA);
    spdlog::info("Percentage of single ToA hits: {:.2f}%", 100.0 * static_cast<double>(num_single_ToA) / static_cast<double>(num_single_ToA + num_multi_ToA));  

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    const std::string out_pdf = script_output_file.substr(0, script_output_file.find_last_of(".")) + ".pdf";
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
    
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream date_stream;
    date_stream << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");

    output_root->cd();

    TCanvas *canvas_toa_frequency = new TCanvas("canvas_toa_frequency", "ToA Frequency", 800, 600);
    canvas_toa_frequency->cd();
    format_1i_hist_canvas(canvas_toa_frequency, h1i_toa_frequency, kBlue, annotation_canvas_title, annotation_testbeam_title, "Generated on " + date_stream.str());
    canvas_toa_frequency->Print((out_pdf + "[").c_str());
    canvas_toa_frequency->Write();
    canvas_toa_frequency->Close();

        // ! do the bx slip searching for each channel
    UInt_t max_bx_slip_threshold = 1024;
    UInt_t min_bx_slip_threshold = 0;
    UInt_t bx_slip_step = 1;
    std::vector<double> bx_slip_values;
    std::vector<std::vector<double>> smooth_level_matrix(FPGA_CHANNEL_NUMBER * vldb_number, std::vector<double>());
    for (UInt_t slip = min_bx_slip_threshold; slip <= max_bx_slip_threshold; slip += bx_slip_step) {
        bx_slip_values.push_back(static_cast<double>(slip));
        for (int ch = 0; ch < FPGA_CHANNEL_NUMBER * vldb_number; ch++) {
            std::vector<UInt_t> & toa_codes = toa_code_matrix[ch];
            std::vector<double> & toa_ns_values = toa_ns_matrix[ch];
            // if the sample size is less than 100, skip
            if (toa_codes.size() < 100) {
                smooth_level_matrix[ch].push_back(std::numeric_limits<double>::quiet_NaN());
                continue;
            }
            // make the th1d with the current slip value
            TH1D h1d_toa_code_corr_slip("h1d_toa_code_corr_slip", "h1d_toa_code_corr_slip;Corrected ToA (ns);Counts", 256, 0, 25.0*static_cast<double>(machine_gun_samples));
            for (size_t i = 0; i < toa_codes.size(); i++) {
                UInt_t toa_code = toa_codes[i];
                double toa_ns = toa_ns_values[i];
                // apply slip
                if (toa_code >= slip) {
                    toa_ns -= 25.0;
                }
                h1d_toa_code_corr_slip.Fill(toa_ns);
            }
            auto smooth_level = FastSmoothMetrics(&h1d_toa_code_corr_slip);
            smooth_level_matrix[ch].push_back(smooth_level.chi2_adj);
        }
    }
    // draw the TGraphs
    TCanvas *canvas_bx_slip = new TCanvas("canvas_bx_slip", "Bx Slip Scan", 1200, 800);
    canvas_bx_slip->cd();
    std::vector<TGraph*> tg_bx_slip_list;
    std::vector<TH1D*> h1d_half_optimal_bx_slip_list;
    std::vector<std::pair<int, double>> optimal_bx_slip_values;
    optimal_bx_slip_values.reserve(FPGA_CHANNEL_NUMBER * vldb_number);
    for (int half = 0; half < 8; half++) {
        TH1D *h1d_half_optimal_bx_slip = new TH1D(("h1d_half_optimal_bx_slip_half" + std::to_string(half)).c_str(),
                                                  ("h1d_half_optimal_bx_slip_half" + std::to_string(half) + ";Optimal Bx Slip Value;Counts").c_str(),
                                                  256, 0, 1024);
        h1d_half_optimal_bx_slip->SetDirectory(nullptr);
        h1d_half_optimal_bx_slip_list.push_back(h1d_half_optimal_bx_slip);
    }
    std::vector<Color_t> half_colors = {kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange+7, kViolet, kPink+9};
    for (int ch = 0; ch < FPGA_CHANNEL_NUMBER * vldb_number; ch++) {
        int channel_in_vldb = ch % FPGA_CHANNEL_NUMBER;
        int vldb_id = ch / FPGA_CHANNEL_NUMBER;
        int asic_id = channel_in_vldb / 76;
        int half_id = (channel_in_vldb % 76) / 38;
        int global_half_index = vldb_id * 4 + asic_id * 2 + half_id;
        TGraph *tg_bx_slip = new TGraph();
        tg_bx_slip->SetName(("tg_bx_slip_ch" + std::to_string(ch)).c_str());
        tg_bx_slip->SetTitle(("bx slip scan ch" + std::to_string(ch) + ";Bx Slip Threshold;Smooth Level (TV Norm)").c_str());
        tg_bx_slip->SetMarkerStyle(20);
        tg_bx_slip->SetMarkerSize(0.5);
        tg_bx_slip->SetMarkerColor(half_colors[global_half_index % half_colors.size()]);
        // if the toa_codes size is less than 100, make dummy graph
        if (toa_code_matrix[ch].size() < 100) {
            for (size_t i = 0; i < bx_slip_values.size(); i++) {
                double slip_value = bx_slip_values[i];
                tg_bx_slip->SetPoint(static_cast<int>(i), slip_value, std::numeric_limits<double>::quiet_NaN());
            }
            tg_bx_slip_list.push_back(tg_bx_slip);
            continue;
        }
        for (size_t i = 0; i < bx_slip_values.size(); i++) {
            double slip_value = bx_slip_values[i];
            double smooth_level = smooth_level_matrix[ch][i];
            tg_bx_slip->SetPoint(static_cast<int>(i), slip_value, smooth_level);
        }
        // find the minimal point
        double min_smooth_level = std::numeric_limits<double>::max();
        double optimal_bx_slip = 0.0;
        for (size_t i = 0; i < bx_slip_values.size(); i++) {
            double smooth_level = smooth_level_matrix[ch][i];
            if (smooth_level < min_smooth_level) {
                min_smooth_level = smooth_level;
                optimal_bx_slip = bx_slip_values[i];
            }
        }
        
        spdlog::info("Channel {} (VLDB {}, ASIC {}, Half {}) optimal bx slip: {}, smooth level: {}", ch, vldb_id, asic_id, half_id, optimal_bx_slip, min_smooth_level);
        if (global_half_index >= 0 && global_half_index < static_cast<int>(h1d_half_optimal_bx_slip_list.size())) {
            h1d_half_optimal_bx_slip_list[global_half_index]->Fill(optimal_bx_slip);
        } else {
            spdlog::warn("Computed global_half_index {} out of range for channel {}", global_half_index, ch);
        }
        tg_bx_slip_list.push_back(tg_bx_slip);
        optimal_bx_slip_values.emplace_back(ch, optimal_bx_slip);

        // save channel 299 as a pdf file
        if (ch == chn_example) {
            TCanvas* canvas_bx_slip_example = new TCanvas("canvas_bx_slip_example", "BX Slip vs Optimal BX Slip", 800, 600);
            canvas_bx_slip_example->SetLeftMargin(0.15);
            tg_bx_slip->SetTitle("");
            // {\chi^2}_\Delta [ToA]
            tg_bx_slip->GetYaxis()->SetTitle("#chi^{2}_{#Delta} [ToA]");
            tg_bx_slip->GetYaxis()->SetTitleOffset(1.5);
            tg_bx_slip->GetXaxis()->SetTitle("BX Slip Threshold [ToA]");
            tg_bx_slip->Draw("AP");
            // add TLatex annotation on the bottom right
            TLatex *latex_bxslip = new TLatex();
            latex_bxslip->SetNDC();
            latex_bxslip->SetTextAlign(31);
            latex_bxslip->SetTextSize(0.05);
            latex_bxslip->SetTextFont(62);
            latex_bxslip->DrawLatex(0.88, 0.3, annotation_canvas_title.c_str());
            latex_bxslip->SetTextSize(0.035);
            latex_bxslip->SetTextFont(42);
            latex_bxslip->DrawLatex(0.88, 0.25, annotation_testbeam_title.c_str());
            latex_bxslip->DrawLatex(0.88, 0.2, ("Channel " + std::to_string(ch)).c_str());
            auto now = std::chrono::system_clock::now();
            std::time_t now_c = std::chrono::system_clock::to_time_t(now);
            char date_str[100];
            std::strftime(date_str, sizeof(date_str), "%d-%m-%Y", std::localtime(&now_c));
            latex_bxslip->DrawLatex(0.88, 0.15, date_str);
            
            
            canvas_bx_slip_example->SaveAs((out_pdf + "_chn299.pdf").c_str());
            canvas_bx_slip_example->Close();
        }
    }

    // if the channel is not in the mapping, set the optimal bx slip to average of the half
    // create a directory to save the TParameters
    output_root->mkdir("Optimal_Bx_Slip_Parameters")->cd();
    std::vector<double> half_optimal_averages;
    for (size_t half = 0; half < 2*vldb_number*2; half++) {
        half_optimal_averages.push_back(h1d_half_optimal_bx_slip_list[half]->GetMean());
    }
    for (int ch = 0; ch < FPGA_CHANNEL_NUMBER * vldb_number; ch++) {
        int channel_in_vldb = ch % FPGA_CHANNEL_NUMBER;
        int vldb_id = ch / FPGA_CHANNEL_NUMBER;
        int asic_id = channel_in_vldb / 76;
        int half_id = (channel_in_vldb % 76) / 38;
        int global_half_index = vldb_id * 4 + asic_id * 2 + half_id;
        // check if the channel is in optimal_bx_slip_values
        auto it = std::find_if(optimal_bx_slip_values.begin(), optimal_bx_slip_values.end(),
                               [ch](const std::pair<int, double>& p) { return p.first == ch; });
        if (it == optimal_bx_slip_values.end()) {
            double average_bx_slip = half_optimal_averages[global_half_index];
            spdlog::info("Channel {} (VLDB {}, ASIC {}, Half {}) not in mapping, set optimal bx slip to half average: {}", ch, vldb_id, asic_id, half_id, average_bx_slip);
            optimal_bx_slip_values.emplace_back(ch, average_bx_slip);
        }
        // write to output file as TParameter
        double optimal_bx_slip = 0.0;
        for (const auto& p : optimal_bx_slip_values) {
            if (p.first == ch) {
                optimal_bx_slip = p.second;
                break;
            }
        }
        TParameter<double> param_optimal_bx_slip(("optimal_bx_slip_ch" + std::to_string(ch)).c_str(), optimal_bx_slip);
        param_optimal_bx_slip.Write();
    }
    output_root->cd();
    draw_mosaic_fixed(*canvas_bx_slip, tg_bx_slip_list, topo_ped_median);
    canvas_bx_slip->Modified();
    canvas_bx_slip->Update();
    canvas_bx_slip->Print(out_pdf.c_str());
    canvas_bx_slip->Write();
    canvas_bx_slip->Close();

    TCanvas *canvas_half_optimal_bx_slip = new TCanvas("canvas_half_optimal_bx_slip", "Half Optimal Bx Slip", 1200, 800);
    // draw the half optimal bx slip histograms

    TLegend *legend_half_optimal_bx_slip = new TLegend(0.7, 0.7, 0.89, 0.89);
    legend_half_optimal_bx_slip->SetFillStyle(0);
    legend_half_optimal_bx_slip->SetBorderSize(0);
    for (int half = 0; half < 8; half++) {
        TH1D *h1d_half_optimal_bx_slip = h1d_half_optimal_bx_slip_list[half];
        h1d_half_optimal_bx_slip->SetTitle("");
        h1d_half_optimal_bx_slip->SetStats(0);
        h1d_half_optimal_bx_slip->SetLineColor(half_colors[half]);
        h1d_half_optimal_bx_slip->SetLineWidth(2);
        h1d_half_optimal_bx_slip->Draw("HIST SAME");
        legend_half_optimal_bx_slip->AddEntry(h1d_half_optimal_bx_slip, ("Half " + std::to_string(half)).c_str(), "l");
        // find the peak position
        int max_bin = h1d_half_optimal_bx_slip->GetMaximumBin();
        double peak_bx_slip = h1d_half_optimal_bx_slip->GetBinCenter(max_bin);
        double peak_counts = h1d_half_optimal_bx_slip->GetBinContent(max_bin);
        spdlog::info("Half {} optimal bx slip peak: {}, counts: {}", half, peak_bx_slip, peak_counts);
        TParameter<double> param_peak_bx_slip(("peak_bx_slip_half" + std::to_string(half)).c_str(), peak_bx_slip);
        param_peak_bx_slip.Write();
    }
    legend_half_optimal_bx_slip->Draw();
    canvas_half_optimal_bx_slip->Modified();
    canvas_half_optimal_bx_slip->Update();
    canvas_half_optimal_bx_slip->Print(out_pdf.c_str());
    canvas_half_optimal_bx_slip->Write();
    canvas_half_optimal_bx_slip->Close();
    

    TCanvas *canvas_raw_toa_ns = new TCanvas("canvas_raw_toa_ns", "Raw ToA in ns", 1200, 800);
    draw_mosaic_fixed(*canvas_raw_toa_ns, h1d_raw_toa_ns_list, topo_ped_median);
    canvas_raw_toa_ns->Modified();
    canvas_raw_toa_ns->Update();
    canvas_raw_toa_ns->Print(out_pdf.c_str());
    canvas_raw_toa_ns->Write();
    canvas_raw_toa_ns->Close();

    TCanvas *canvas_toa_adc_max_corr = new TCanvas("canvas_toa_adc_max_corr", "ToA ADC Max Correction", 1200, 800);
    draw_mosaic_fixed(*canvas_toa_adc_max_corr, h2d_toa_adc_max_corr_list, topo_wave);
    canvas_toa_adc_max_corr->Modified();
    canvas_toa_adc_max_corr->Update();
    canvas_toa_adc_max_corr->Print(out_pdf.c_str());
    canvas_toa_adc_max_corr->Write();
    canvas_toa_adc_max_corr->Close();

    TCanvas *canvas_toa_adc_max_corr_example = new TCanvas("canvas_toa_adc_max_corr_example", "ToA ADC Max Correction Example Channel", 800, 600);
    auto h2d_example = h2d_toa_adc_max_corr_list[chn_example];
    h2d_example->GetXaxis()->SetTitle("ADC Peak (pedestal subtracted)");
    h2d_example->GetYaxis()->SetTitle("Raw ToA [ns]");
    format_2d_hist_canvas(canvas_toa_adc_max_corr_example, h2d_example, kBlue+2, annotation_canvas_title, annotation_testbeam_title, "Channel_" + std::to_string(chn_example));
    canvas_toa_adc_max_corr_example->Print(out_pdf.c_str());
    canvas_toa_adc_max_corr_example->Write();
    canvas_toa_adc_max_corr_example->Close();

    TCanvas *canvas_toa_code_corr = new TCanvas("canvas_toa_code_corr", "ToA Code Correction", 1200, 800);
    draw_mosaic_fixed(*canvas_toa_code_corr, h2d_toa_code_corr_list, topo_wave);
    canvas_toa_code_corr->Modified();
    canvas_toa_code_corr->Update();
    canvas_toa_code_corr->Print(out_pdf.c_str());
    canvas_toa_code_corr->Write();
    canvas_toa_code_corr->Close();

    // draw dummy page to close the pdf
    TCanvas *canvas_dummy = new TCanvas("canvas_dummy", "Dummy Canvas", 800, 600);
    canvas_dummy->Print((out_pdf + "]").c_str());

    output_root->Close();

    delete canvas_toa_frequency;
    delete canvas_raw_toa_ns;
    delete canvas_toa_adc_max_corr;
    delete canvas_toa_code_corr;
    delete canvas_toa_adc_max_corr_example;
    delete canvas_bx_slip;
    delete canvas_half_optimal_bx_slip;
    delete canvas_dummy;

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}