#include <fstream>
#include "H2GCROC_Common.hxx"



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

    script_input_file  = parsed["file"].as<std::string>();
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

    spdlog::info("Script name: {}", script_name);
    spdlog::info("Input file: {}", script_input_file);
    spdlog::info("Output file: {} in {}", script_output_file, script_output_folder);
    spdlog::info("Number of events: {}", script_n_events);

    std::string mapping_json_file = "config/mapping_Oct2025.json";
    std::ifstream mapping_json_ifs(mapping_json_file);
    if (!mapping_json_ifs.is_open()) {
        spdlog::error("Failed to open mapping json file {}", mapping_json_file);
        return 1;
    }
    json mapping_json;
    mapping_json_ifs >> mapping_json;
    const auto& sipm_board = mapping_json.at("SiPM_Board");
    const auto& board_loc = mapping_json.at("Board_Loc");
    const auto& board_rotation = mapping_json.at("Board_Rotation");
    const auto& board_flip = mapping_json.at("Board_Flip");

    // * --- Read input file ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    int machine_gun_samples = 16;
    int vldb_number = 2;

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

    // create hist2d for heatmap
    const int board_cols = 8;
    const int board_rows = 4;
    int heat_map_x_bins = 16;
    int heat_map_y_bins = 12;
    TH2D* heat_map = new TH2D("", "Heat Map;X;Y;Counts", heat_map_x_bins, 0, heat_map_x_bins, heat_map_y_bins, 0, heat_map_y_bins);
    // no title for the histogram
    heat_map->SetTitle("");

    // start event loop
    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            // channel loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int channel_index_in_h2g = channel % 76;
                int channel_index_in_h2g_half = channel_index_in_h2g % 38;
                int asic_id = channel / 76;
                int half_id = channel_index_in_h2g / 38;
                
                if (channel_index_in_h2g_half == 19)
                    continue;
                if (channel_index_in_h2g_half == 0)
                    continue;

                int channel_index_no_CM_Calib = channel_index_in_h2g_half;
                
                if (channel_index_in_h2g_half > 19)
                    channel_index_no_CM_Calib -= 2;
                else if (channel_index_in_h2g_half > 0)
                    channel_index_no_CM_Calib -= 1;

                // if (channel_index_no_CM_Calib % 8 == 0)
                //     continue;

                std::string chn_key = std::to_string(channel_index_no_CM_Calib);
                int col = -1, row = -1;
                if (sipm_board.contains(chn_key)) {
                    col = sipm_board.at(chn_key).at("col").get<int>();
                    row = sipm_board.at(chn_key).at("row").get<int>();
                }

                std::string board_key = std::to_string(half_id + asic_id * 2 + vldb_id * 4);
                int board_col = -1, board_row = -1;
                // print board_key
                // spdlog::info("Board Key: {} for VLDB {} Channel {}", board_key, vldb_id, channel);
                int board_rotated = 0;
                int board_flipped = 0;
                if (board_loc.contains(board_key)) {
                    board_col = board_loc.at(board_key).at("col").get<int>();
                    board_row = board_loc.at(board_key).at("row").get<int>();
                }
                if (board_rotation.contains(board_key)) {
                    board_rotated = board_rotation.at(board_key).get<int>();
                }
                if (board_flip.contains(board_key)) {
                    board_flipped = board_flip.at(board_key).get<int>();
                }

                int uni_col = -1, uni_row = -1;
                if (board_rotated == 0) {
                    uni_col = board_col * board_cols + col;
                    uni_row = board_row * board_rows + row;
                } else {
                    // rotate 180 degrees
                    uni_col = board_col * board_cols + (board_cols - 1 - col);
                    uni_row = board_row * board_rows + (board_rows - 1 - row);
                }

                UInt_t adc_pedestal = 1024;
                UInt_t adc_peak = 0;
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    auto &adc = val0_list_pools[vldb_id][0][channel + sample * FPGA_CHANNEL_NUMBER];
                    if (sample == 0) {
                        adc_pedestal = adc;
                    }
                    if (adc > adc_peak) {
                        adc_peak = adc;
                    }
                }
                auto adc_subtracted = int(adc_peak) - int(adc_pedestal);

                // fill the heatmap with the adc_subtracted value
                if (uni_col >= 0 && uni_col < heat_map_x_bins && uni_row >= 0 && uni_row < heat_map_y_bins) {
                    heat_map->Fill(uni_col + 0.5, uni_row + 0.5, adc_subtracted);
                }
                // spdlog::info("VLDB {} Channel {} ADC Pedestal {} Peak {} Subtracted {}", vldb_id, channel, adc_pedestal, adc_peak, adc_subtracted);
            }
        }
    }

    auto output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }
    output_root->cd();
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    heat_map->SetStats(0);
    heat_map->Draw("COLZ");
    
    // remove stat box
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderSize(0);
    c1->SetBorderSize(0);
    c1->SetTickx(0);
    c1->SetTicky(0);
    c1->Update();

    // write annotation
    TLatex latex;
    latex.SetTextColor(kGray+2);
    latex.SetNDC();
    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;

    latex.SetTextSize(0.05);
    latex.SetTextFont(62);
    latex.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());

    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
    // write run number, date time
    
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream date_stream;
    date_stream << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    latex.DrawLatex(0.12, 0.76, ("Run " + std::to_string(0)).c_str());
    latex.DrawLatex(0.12, 0.72, date_stream.str().c_str());

    c1->Write();

    output_root->Close();

    input_root->Close();


    spdlog::info("Output file {} has been saved.", script_output_file);


    //     std::vector <ULong64_t*> branch_timestamps_list;
    //     std::vector <UInt_t*> branch_daqh_list_list;
    //     std::vector <Bool_t*> branch_tc_list_list;
    //     std::vector <Bool_t*> branch_tp_list_list;
    //     std::vector <UInt_t*> branch_val0_list_list;
    //     std::vector <UInt_t*> branch_val1_list_list;
    //     std::vector <UInt_t*> branch_val2_list_list;
    //     std::vector <UInt_t*> branch_crc32_list_list;
    //     std::vector <UInt_t*> branch_last_heartbeat_list;

    //     for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
    //         auto _fpga_id = _legal_fpga_id_list[_fpga_index];
    //         auto *branch_timestamps = new ULong64_t[machine_gun_samples];  // 64 bits
    //         auto *branch_daqh_list  = new UInt_t[4 * machine_gun_samples];    // 32 bits
    //         auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    //         auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    //         auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    //         auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    //         auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    //         auto *branch_crc32_list = new UInt_t[4 * machine_gun_samples];
    //         auto *branch_last_heartbeat = new UInt_t[machine_gun_samples];
    
    //         _input_tree->SetBranchAddress(("timestamps_"    + std::to_string(_fpga_id)).c_str(), branch_timestamps);
    //         _input_tree->SetBranchAddress(("daqh_list_"     + std::to_string(_fpga_id)).c_str(), branch_daqh_list);
    //         _input_tree->SetBranchAddress(("tc_list_"       + std::to_string(_fpga_id)).c_str(), branch_tc_list);
    //         _input_tree->SetBranchAddress(("tp_list_"       + std::to_string(_fpga_id)).c_str(), branch_tp_list);
    //         _input_tree->SetBranchAddress(("val0_list_"     + std::to_string(_fpga_id)).c_str(), branch_val0_list);
    //         _input_tree->SetBranchAddress(("val1_list_"     + std::to_string(_fpga_id)).c_str(), branch_val1_list);
    //         _input_tree->SetBranchAddress(("val2_list_"     + std::to_string(_fpga_id)).c_str(), branch_val2_list);
    //         _input_tree->SetBranchAddress(("crc32_list_"    + std::to_string(_fpga_id)).c_str(), branch_crc32_list);
    //         _input_tree->SetBranchAddress(("last_heartbeat_"+ std::to_string(_fpga_id)).c_str(), branch_last_heartbeat);
    
    //         branch_timestamps_list.push_back(branch_timestamps);
    //         branch_daqh_list_list.push_back(branch_daqh_list);
    //         branch_tc_list_list.push_back(branch_tc_list);
    //         branch_tp_list_list.push_back(branch_tp_list);
    //         branch_val0_list_list.push_back(branch_val0_list);
    //         branch_val1_list_list.push_back(branch_val1_list);
    //         branch_val2_list_list.push_back(branch_val2_list);
    //         branch_crc32_list_list.push_back(branch_crc32_list);
    //         branch_last_heartbeat_list.push_back(branch_last_heartbeat);
    //     }

    //     int _input_double_tot_counter = 0;
    //     int _input_double_toa_counter = 0;

    //     int _input_valid_tot_counter = 0;
    //     int _input_valid_toa_counter = 0;

    //     auto _run_adc_sum_hist1d = run_adc_sum_hist1d[_run_index];

    //     for (int _entry = 0; _entry < entry_max; _entry++) {
    //         _input_tree->GetEntry(_entry);

    //         double _run_adc_sum = 0;

    //         for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
    //             auto _fpga_id    = _legal_fpga_id_list[_fpga_index];
    //             auto _timestamp  = branch_timestamps_list[_fpga_index][0];
    //             auto _daqh_list  = branch_daqh_list_list[_fpga_index];
    //             auto _tc_list    = branch_tc_list_list[_fpga_index];
    //             auto _tp_list    = branch_tp_list_list[_fpga_index];
    //             auto _val0_list  = branch_val0_list_list[_fpga_index];
    //             auto _val1_list  = branch_val1_list_list[_fpga_index];
    //             auto _val2_list  = branch_val2_list_list[_fpga_index];

    //             std::vector <bool> _hamming_code_pass_list;
    //             for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
    //                 bool _hamming_code_pass = true;
    //                 for (int _daqh_index = 0; _daqh_index < 4; _daqh_index++){
    //                     auto _daqh = _daqh_list[_sample_index * 4 + _daqh_index];
    //                     auto _h1h2h3 = (_daqh >> 4) & 0x7;
    //                     if (_h1h2h3 != 0x00){
    //                         _hamming_code_pass = false;
    //                     }
    //                 }
    //                 _hamming_code_pass_list.push_back(_hamming_code_pass);
    //             }

    //             std::vector <int> _channel_val1_max_list;
    //             std::vector <int> _channel_val2_max_list;
    //             std::vector <int> _channel_val1_max_index_list;
    //             std::vector <int> _channel_val2_max_index_list;

    //             for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
    //                 auto _channel_valid = get_valid_fpga_channel(_channel_index);
    //                 if (_channel_valid == -1){
    //                     continue;
    //                 }
    //                 auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
    //                 auto _hist_index = channel_wise_hists_unified_channel_map[_unified_valid_channel_number];
    //                 int _val1_max = -1;
    //                 int _val1_max_index = -1;
    //                 int _val2_max = -1;
    //                 int _val2_max_index = -1;

    //                 for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
    //                     auto _val1 = _val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
    //                     auto _val2 = _val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
    //                     if (_val1 > 0){
    //                         if (_val1_max == -1) {
    //                             _val1_max = _val1;
    //                             _val1_max_index = _sample_index;
    //                             _input_valid_tot_counter++;
    //                         } else {
    //                             _input_double_tot_counter++;
    //                         }
    //                     }
    //                     if (_val2 > 0) {
    //                         if (_val2_max == -1) {
    //                             _val2_max = _val2;
    //                             _val2_max_index = _sample_index;
    //                             _input_valid_toa_counter++;
    //                         } else {
    //                             _input_double_toa_counter++;
    //                         }
    //                     }
    //                 }
    //                 _channel_val1_max_list.push_back(_val1_max);
    //                 _channel_val2_max_list.push_back(_val2_max);
    //                 _channel_val1_max_index_list.push_back(_val1_max_index);
    //                 _channel_val2_max_index_list.push_back(_val2_max_index);
    //             }

    //             for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
    //                 auto _channel_valid = get_valid_fpga_channel(_channel_index);
    //                 auto _val1_max = _channel_val1_max_list[_channel_index];
    //                 auto _val2_max = _channel_val2_max_list[_channel_index];
    //                 auto _val1_max_index = _channel_val1_max_index_list[_channel_index];
    //                 auto _val2_max_index = _channel_val2_max_index_list[_channel_index];

    //                 if (_channel_valid == -1){
    //                     continue;
    //                 }
    //                 auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
    //                 auto _hist_index = channel_wise_hists_unified_channel_map[_unified_valid_channel_number];
    //                 auto _channel_sample_hist2d = run_chn_hist2d[_run_index][_hist_index];
    //                 auto _channel_toatime_toacode_hist2d = run_chn_toatime_toacode_hist2d[_run_index][_hist_index];
    //                 auto _channel_toatime_adcmax_hist2d = run_chn_toatime_adcmax_hist2d[_run_index][_hist_index];

    //                 int _val0_max = -1;
    //                 int _pedestal_event = -1;

    //                 double _channel_adc_value = 0;

    //                 for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
    //                     auto _val0 = _val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER];
    //                     if (_sample_index == 0) {
    //                         _pedestal_event = int(_val0);
    //                     }
    //                     if (int(_val0) > _val0_max) {
    //                         _val0_max = int(_val0);
    //                     }

    //                     if (_hamming_code_pass_list[_sample_index] && _val2_max > 0) {
    //                         _channel_sample_hist2d->Fill(_sample_index * sample_time, _val0);
    //                     }
    //                 }

    //                 // * --- Calculate pedestal subtracted ADC value -------------------------
    //                 // * --------------------------------------------------------------------
    //                 double _pedestal_mean1       = channel_pede_mean1[_unified_valid_channel_number];
    //                 double _pedestal_error1      = channel_pede_error1[_unified_valid_channel_number];
    //                 double _pedestal_mean2       = channel_pede_mean2[_unified_valid_channel_number]; 
    //                 double _pedestal_error2      = channel_pede_error2[_unified_valid_channel_number];
    //                 int _pedestal_peak_counts    = int(channel_pede_peak_counts[_unified_valid_channel_number]);

    //                 if (_val2_max > 0) {
    //                     double _toa_time = decode_toa_value_ns(_val2_max) + _val2_max_index * sample_time;
    //                     _channel_toatime_toacode_hist2d->Fill(_toa_time, _val2_max);
    //                     _channel_toatime_adcmax_hist2d->Fill(_val0_max, _toa_time);
    //                     // LOG(DEBUG) << "ADC Max: " << _val0_max << " TOA Time: " << _toa_time << " TOA Code: " << _val2_max;
    //                 }

    //                 if (_pedestal_peak_counts == 1){
    //                     auto _pedestal = _pedestal_mean1;
    //                     _channel_adc_value = _val0_max - _pedestal;
    //                     _run_adc_sum += _channel_adc_value;
    //                 } else if (_pedestal_peak_counts == 2){
    //                     auto _pedestal1_dist = std::abs(_pedestal_event - _pedestal_mean1);
    //                     auto _pedestal2_dist = std::abs(_pedestal_event - _pedestal_mean2);
    //                     auto _pedestal = 0;
    //                     if (_pedestal1_dist < _pedestal2_dist){
    //                         _pedestal = _pedestal_mean1;
    //                     } else {
    //                         _pedestal = _pedestal_mean2;
    //                     }
    //                     _channel_adc_value = _val0_max - _pedestal;
    //                     _run_adc_sum += _channel_adc_value;
    //                 }
    //             }
    //         }
    //         _run_adc_sum_hist1d->Fill(_run_adc_sum);
    //     }
    //     LOG(INFO) << "Run " << _run_number << " with " << _input_valid_tot_counter << " valid TOTs (" << (100.0 * _input_valid_tot_counter / (entry_max * fpga_count * FPGA_CHANNEL_NUMBER)) << "%) and " << _input_double_tot_counter << " double TOTs (" << (100.0 * _input_double_tot_counter / (entry_max * fpga_count * FPGA_CHANNEL_NUMBER)) << "%)";
    //     LOG(INFO) << "Run " << _run_number << " with " << _input_valid_toa_counter << " valid TOAs (" << (100.0 * _input_valid_toa_counter / (entry_max * fpga_count * FPGA_CHANNEL_NUMBER)) << "%) and " << _input_double_toa_counter << " double TOAs (" << (100.0 * _input_double_toa_counter / (entry_max * fpga_count * FPGA_CHANNEL_NUMBER)) << "%)";
    //     _input_root->Close();
    // }

    // TColor *softBlue   = new TColor(4001, 0.35, 0.55, 0.75);  // Steel Blue
    // TColor *softRed    = new TColor(4002, 0.75, 0.35, 0.35);  // Soft Coral
    // TColor *softGreen  = new TColor(4003, 0.45, 0.65, 0.45);  // Sage Green
    // TColor *softPurple = new TColor(4004, 0.55, 0.45, 0.65);  // Lavender Gray
    // TColor *softTeal   = new TColor(4005, 0.35, 0.65, 0.65);  // Light Teal
    // TColor *softOrange = new TColor(4006, 0.85, 0.55, 0.35);  // Warm Orange
    // TColor *softBrown  = new TColor(4007, 0.65, 0.55, 0.45);  // Sand Brown
    // TColor *softGray   = new TColor(4008, 0.55, 0.55, 0.55);  // Mid Gray
    // TColor *softOlive  = new TColor(4009, 0.65, 0.65, 0.45);  // Olive Green
    // std::vector <Color_t> color_list = {kBlack, 4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008, 4009};
    // std::string pdf_file_name = opts.output_file.substr(0, opts.output_file.find_last_of('.')) + ".pdf";

    // // * --- Save the histograms --------------------------------------------------------
    // // * --------------------------------------------------------------------------------
    // bool is_first_pdf_page_written = false;
    // LOG(INFO) << "Saving channel sample histograms ...";
    // for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
    //     auto _channel_hist2d = run_chn_hist2d[_run_index];
    //     run_chn_sample_hist2d_folder->cd();
    //     for (int _channel_index = 0; _channel_index < _channel_hist2d.size(); _channel_index++) {
    //         auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_channel_index];
    //         auto _hist2d = _channel_hist2d[_channel_index];
    //         _hist2d->SetStats(0);
    //         _hist2d->Write();
    //     }
    //     if (_run_index == run_index_highest_energy || _run_index == run_index_lowest_energy){
    //         global_painter->draw_global_channel_hists2D(_channel_hist2d, channel_wise_hists_unified_channel_map, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnSampleHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
    //         auto _canvas = global_painter->get_canvas();
    //         if (!is_first_pdf_page_written){
    //             _canvas->SaveAs((pdf_file_name + "(").c_str());
    //             is_first_pdf_page_written = true;
    //         } else {
    //             _canvas->SaveAs(pdf_file_name.c_str());
    //         }
    //         output_root->cd();
    //         _canvas->Write();
    //     }
    // }

    // LOG(INFO) << "Saving toa time vs toa code histograms ...";
    // for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
    //     auto _channel_toatime_toacode_hist2d = run_chn_toatime_toacode_hist2d[_run_index];
    //     run_chn_toatime_toacode_hist2d_folder->cd();
    //     for (int _channel_index = 0; _channel_index < _channel_toatime_toacode_hist2d.size(); _channel_index++) {
    //         auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_channel_index];
    //         auto _hist2d = _channel_toatime_toacode_hist2d[_channel_index];
    //         _hist2d->GetXaxis()->SetTitle("TOA Time [ns]");
    //         _hist2d->GetYaxis()->SetTitle("TOA Code");
    //         _hist2d->SetStats(0);
    //         _hist2d->Write();
    //     }
    //     if (_run_index == run_index_lowest_energy || _run_index == run_index_highest_energy){
    //         global_painter->draw_global_channel_hists2D(_channel_toatime_toacode_hist2d, channel_wise_hists_unified_channel_map, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnToaTimeToaCodeHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
    //         auto _canvas = global_painter->get_canvas();
    //         _canvas->SaveAs(pdf_file_name.c_str());
    //         output_root->cd();
    //         _canvas->Write();
    //     }
    // }

    // LOG(INFO) << "Saving toa time vs adc max histograms ...";
    // for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
    //     auto _channel_toatime_adcmax_hist2d = run_chn_toatime_adcmax_hist2d[_run_index];
    //     run_chn_toatime_adcmax_hist2d_folder->cd();
    //     for (int _channel_index = 0; _channel_index < _channel_toatime_adcmax_hist2d.size(); _channel_index++) {
    //         auto _unified_valid_channel_number = channel_wise_hists_unified_channel_map[_channel_index];
    //         auto _hist2d = _channel_toatime_adcmax_hist2d[_channel_index];
    //         _hist2d->GetXaxis()->SetTitle("ADC Max");
    //         _hist2d->GetYaxis()->SetTitle("TOA Time [ns]");
    //         _hist2d->SetStats(0);
    //         _hist2d->Write();
    //     }
    //     if (_run_index == run_index_lowest_energy || _run_index == run_index_highest_energy){
    //         global_painter->draw_global_channel_hists2D(_channel_toatime_adcmax_hist2d, channel_wise_hists_unified_channel_map, ("Run"+std::to_string(config_run_numbers[_run_index])+"ChnToaTimeAdcMaxHist2D").c_str(), ("Run " + std::to_string(config_run_numbers[_run_index])).c_str());
    //         auto _canvas = global_painter->get_canvas();
    //         _canvas->Print(pdf_file_name.c_str());
    //         output_root->cd();
    //         _canvas->Write();
    //     }
    // }

    // int run_adc_sum_hist_max = 0;
    // for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
    //     auto _adc_sum_hist1d = run_adc_sum_hist1d[_run_index];
    //     if (_adc_sum_hist1d->GetMaximum() > run_adc_sum_hist_max) {
    //         run_adc_sum_hist_max = _adc_sum_hist1d->GetMaximum();
    //     }
    //     run_adc_sum_hist1d_folder->cd();
    //     _adc_sum_hist1d->SetStats(0);
    //     _adc_sum_hist1d->Write();
    // }
    // auto adc_sum_hist1d_canvas = new TCanvas("adc_sum_hist1d_canvas", "adc_sum_hist1d_canvas", 800, 600);
    // auto adc_sum_legend = new TLegend(0.42, 0.5, 0.89, 0.89);
    // adc_sum_legend->SetBorderSize(0);
    // adc_sum_legend->SetFillStyle(0);
    // adc_sum_legend->SetTextSize(0.02);
    // adc_sum_legend->SetTextFont(102);

    // std::vector <double> adc_sum_fit_mu_list;
    // std::vector <double> adc_sum_fit_sigma_list;
    // std::vector <double> adc_sum_fit_mu_error_list;
    // std::vector <double> adc_sum_fit_sigma_error_list;
    // std::vector <double> adc_sum_beam_energy_list;
    // std::vector <double> adc_sum_beam_energy_error_list;

    // for (int _run_index = 0; _run_index < config_run_numbers.size(); _run_index++) {
    //     auto _adc_sum_hist1d = run_adc_sum_hist1d[_run_index];
    //     _adc_sum_hist1d->SetMaximum(run_adc_sum_hist_max * 1.3);
    //     _adc_sum_hist1d->SetLineColor(color_list[_run_index % color_list.size()]);
    //     if (_run_index == 0) {
    //         _adc_sum_hist1d->SetTitle("ADC Sum for Different Beam Energies");
    //         _adc_sum_hist1d->GetXaxis()->SetTitle("ADC Sum");
    //         _adc_sum_hist1d->GetYaxis()->SetTitle("Counts");
    //         _adc_sum_hist1d->Draw();
    //     } else {
    //         _adc_sum_hist1d->Draw("SAME");
    //     }
    //     // * do the gaussian fit
    //     std::vector <double> _adc_sum_fit_range = {1.6, 1.8, 2.0, 2.2, 2.4, 2.6}; // unit: sigma
    //     std::vector <double> _adc_sum_fit_offsets = {-1.0, -0.5, 0.0, 0.5, 1.0};
    //     double _adc_sum_fit_range_min = 0;
    //     double _adc_sum_fit_range_max = 0;
    //     double _adc_sum_fit_offset_min = 0;
    //     double _adc_sum_fit_offset_max = 0;
    //     for (auto _range: _adc_sum_fit_range) {
    //         if (_range < _adc_sum_fit_range_min) {
    //             _adc_sum_fit_range_min = _range;
    //         }
    //         if (_range > _adc_sum_fit_range_max) {
    //             _adc_sum_fit_range_max = _range;
    //         }
    //     }
    //     for (auto _offset: _adc_sum_fit_offsets) {
    //         if (_offset < _adc_sum_fit_offset_min) {
    //             _adc_sum_fit_offset_min = _offset;
    //         }
    //         if (_offset > _adc_sum_fit_offset_max) {
    //             _adc_sum_fit_offset_max = _offset;
    //         }
    //     }

    //     // * find initial parameters
    //     double _adc_sum_fit_mean = _adc_sum_hist1d->GetBinCenter(_adc_sum_hist1d->GetMaximumBin());
    //     double _adc_sum_fit_rms = _adc_sum_hist1d->GetRMS();
    //     double _adc_sum_fit_initial_fit_min = _adc_sum_fit_mean - _adc_sum_fit_rms * 2;
    //     double _adc_sum_fit_initial_fit_max = _adc_sum_fit_mean + _adc_sum_fit_rms * 2;
    //     if (_adc_sum_fit_initial_fit_min < 0) {
    //         _adc_sum_fit_initial_fit_min = 0;
    //     }
    //     if (_adc_sum_fit_initial_fit_max > adc_sum_hist_max) {
    //         _adc_sum_fit_initial_fit_max = adc_sum_hist_max;
    //     }
    //     // * do the first round of fitting to get the mu and sigma
    //     TF1 *_adc_sum_fit_function = new TF1("adc_sum_fit_function", "gaus", _adc_sum_fit_initial_fit_min, _adc_sum_fit_initial_fit_max);
    //     _adc_sum_hist1d->Fit(_adc_sum_fit_function, "RQN");
    //     double _adc_sum_fit_mu = _adc_sum_fit_function->GetParameter(1);
    //     double _adc_sum_fit_sigma = _adc_sum_fit_function->GetParameter(2);

    //     // * do the fitting for pedestal peak
    //     TF1 *_adc_sum_pede_fit_function = new TF1("adc_sum_pede_fit_function", "gaus", adc_sum_hist_min, adc_sum_hist_min + (adc_sum_hist_max - adc_sum_hist_min) * 0.06);
    //     // set parameter range
    //     _adc_sum_pede_fit_function->SetParLimits(1, adc_sum_hist_min, adc_sum_hist_min + (adc_sum_hist_max - adc_sum_hist_min) * 0.1);
    //     _adc_sum_hist1d->Fit(_adc_sum_pede_fit_function, "RQN");
    //     double _adc_sum_pede_mu = _adc_sum_pede_fit_function->GetParameter(1);
    //     double _adc_sum_pede_sigma = _adc_sum_pede_fit_function->GetParameter(2);
    //     // draw the pede mu
    //     TLine *_adc_sum_pede_mu_line = new TLine(_adc_sum_pede_mu, 0, _adc_sum_pede_mu, run_adc_sum_hist_max * 1.3);
    //     _adc_sum_pede_mu_line->SetLineColor(color_list[_run_index % color_list.size()]);
    //     _adc_sum_pede_mu_line->SetLineStyle(2);
    //     _adc_sum_pede_mu_line->Draw("SAME");

    //     std::vector <double> _adc_sum_fit_res_mean;
    //     std::vector <double> _adc_sum_fit_res_sigma;
    //     std::vector <double> _adc_sum_fit_res_mean_error;
    //     std::vector <double> _adc_sum_fit_res_sigma_error;

    //     for (auto _fit_range: _adc_sum_fit_range) {
    //         for (auto _fit_offset: _adc_sum_fit_offsets) {
    //             double _adc_sum_fit_min = _adc_sum_fit_mu - _fit_range * _adc_sum_fit_sigma + _fit_offset * _adc_sum_fit_sigma;
    //             double _adc_sum_fit_max = _adc_sum_fit_mu + _fit_range * _adc_sum_fit_sigma + _fit_offset * _adc_sum_fit_sigma;
    //             if (_adc_sum_fit_min < 0) {
    //                 _adc_sum_fit_min = 0;
    //             }
    //             if (_adc_sum_fit_max > adc_sum_hist_max) {
    //                 _adc_sum_fit_max = adc_sum_hist_max;
    //             }
    //             auto _adc_sum_fit_function = new TF1(("adc_sum_fit_function_"+std::to_string(_run_index)+"_"+std::to_string(_fit_range)+"_"+std::to_string(_fit_offset)).c_str(), "gaus", _adc_sum_fit_min, _adc_sum_fit_max);
    //             _adc_sum_fit_function->SetParameters(_adc_sum_hist1d->GetMaximum(), _adc_sum_fit_mu, _adc_sum_fit_sigma);
    //             _adc_sum_hist1d->Fit(_adc_sum_fit_function, "RQN+");
    //             _adc_sum_fit_function->SetLineColorAlpha(color_list[_run_index % color_list.size()], 0.2);
    //             _adc_sum_fit_function->Draw("SAME");

    //             _adc_sum_fit_res_mean.push_back(_adc_sum_fit_function->GetParameter(1));
    //             _adc_sum_fit_res_sigma.push_back(_adc_sum_fit_function->GetParameter(2));
    //             _adc_sum_fit_res_mean_error.push_back(_adc_sum_fit_function->GetParError(1));
    //             _adc_sum_fit_res_sigma_error.push_back(_adc_sum_fit_function->GetParError(2));
    //         }
    //     }

    //     // * calculate the mean and sigma
    //     double _adc_sum_fit_res_mean_weighted  = 0;
    //     double _adc_sum_fit_res_sigma_weighted = 0;
    //     double _adc_sum_fit_res_mean_err_sys   = 0;
    //     double _adc_sum_fit_res_sigma_err_sys  = 0;
    //     double _adc_sum_fit_res_mean_err_stat  = 0;
    //     double _adc_sum_fit_res_sigma_err_stat = 0;

    //     // * calculate the weighted mean and sigma
    //     double _adc_sum_fit_res_mean_weighted_sum = 0;
    //     double _adc_sum_fit_res_sigma_weighted_sum = 0;
    //     double _adc_sum_fit_res_mean_weighted_sum_weight = 0;
    //     double _adc_sum_fit_res_sigma_weighted_sum_weight = 0;

    //     for (int _fit_index = 0; _fit_index < _adc_sum_fit_res_mean.size(); _fit_index++) {
    //         double _mean = _adc_sum_fit_res_mean[_fit_index];
    //         double _sigma = _adc_sum_fit_res_sigma[_fit_index];
    //         double _mean_error = _adc_sum_fit_res_mean_error[_fit_index];
    //         double _sigma_error = _adc_sum_fit_res_sigma_error[_fit_index];
    //         double _weight_mean = 1.0 / (_mean_error * _mean_error);
    //         double _weight_sigma = 1.0 / (_sigma_error * _sigma_error);
    //         _adc_sum_fit_res_mean_weighted_sum += _mean * _weight_mean;
    //         _adc_sum_fit_res_sigma_weighted_sum += _sigma * _weight_sigma;
    //         _adc_sum_fit_res_mean_weighted_sum_weight += _weight_mean;
    //         _adc_sum_fit_res_sigma_weighted_sum_weight += _weight_sigma;
    //     }

    //     _adc_sum_fit_res_mean_weighted = _adc_sum_fit_res_mean_weighted_sum / _adc_sum_fit_res_mean_weighted_sum_weight;
    //     _adc_sum_fit_res_sigma_weighted = _adc_sum_fit_res_sigma_weighted_sum / _adc_sum_fit_res_sigma_weighted_sum_weight;

    //     // * calculate the systematic error
    //     double _adc_sum_fit_res_mean_err_sys_sum = 0;
    //     double _adc_sum_fit_res_sigma_err_sys_sum = 0;
    //     for (int _fit_index = 0; _fit_index < _adc_sum_fit_res_mean.size(); _fit_index++) {
    //         double _mean = _adc_sum_fit_res_mean[_fit_index];
    //         double _sigma = _adc_sum_fit_res_sigma[_fit_index];
    //         double _mean_error = _adc_sum_fit_res_mean_error[_fit_index];
    //         double _sigma_error = _adc_sum_fit_res_sigma_error[_fit_index];
    //         _adc_sum_fit_res_mean_err_sys_sum += (_mean - _adc_sum_fit_res_mean_weighted) * (_mean - _adc_sum_fit_res_mean_weighted);
    //         _adc_sum_fit_res_sigma_err_sys_sum += (_sigma - _adc_sum_fit_res_sigma_weighted) * (_sigma - _adc_sum_fit_res_sigma_weighted);
    //     }
    //     _adc_sum_fit_res_mean_err_sys = std::sqrt(_adc_sum_fit_res_mean_err_sys_sum / (_adc_sum_fit_res_mean.size() - 1));
    //     _adc_sum_fit_res_sigma_err_sys = std::sqrt(_adc_sum_fit_res_sigma_err_sys_sum / (_adc_sum_fit_res_sigma.size() - 1));

    //     // subtract the pedestal
    //     // _adc_sum_fit_res_mean_weighted -= _adc_sum_pede_mu;

    //     // * calculate the statistical error
    //     double _adc_sum_fit_res_mean_err_stat_sum = 0;
    //     double _adc_sum_fit_res_sigma_err_stat_sum = 0;
    //     for (int _fit_index = 0; _fit_index < _adc_sum_fit_res_mean.size(); _fit_index++) {
    //         double _mean = _adc_sum_fit_res_mean[_fit_index];
    //         double _sigma = _adc_sum_fit_res_sigma[_fit_index];
    //         double _mean_error = _adc_sum_fit_res_mean_error[_fit_index];
    //         double _sigma_error = _adc_sum_fit_res_sigma_error[_fit_index];
    //         _adc_sum_fit_res_mean_err_stat_sum += 1.0 / (_mean_error * _mean_error);
    //         _adc_sum_fit_res_sigma_err_stat_sum += 1.0 / (_sigma_error * _sigma_error);
    //     }
    //     _adc_sum_fit_res_mean_err_stat = 1.0 / std::sqrt(_adc_sum_fit_res_mean_err_stat_sum);
    //     _adc_sum_fit_res_sigma_err_stat = 1.0 / std::sqrt(_adc_sum_fit_res_sigma_err_stat_sum);

    //     // fixed length string
    //     std::string _adc_sum_mean_str = std::to_string(_adc_sum_fit_res_mean_weighted).substr(0, 5);
    //     std::string _adc_sum_sigma_str = std::to_string(_adc_sum_fit_res_sigma_weighted).substr(0, 4);
    //     std::string _adc_sum_mean_err_stat_str = std::to_string(_adc_sum_fit_res_mean_err_stat).substr(0, 4);
    //     std::string _adc_sum_sigma_err_stat_str = std::to_string(_adc_sum_fit_res_sigma_err_stat).substr(0, 4);
    //     std::string _adc_sum_mean_err_sys_str = std::to_string(_adc_sum_fit_res_mean_err_sys).substr(0, 5);
    //     std::string _adc_sum_sigma_err_sys_str = std::to_string(_adc_sum_fit_res_sigma_err_sys).substr(0, 5);
    //     std::string _adc_sum_beam_energy_str;
    //     if (int(config_beam_energies[_run_index]) < 100) {
    //         _adc_sum_beam_energy_str = " " + std::to_string(int(config_beam_energies[_run_index]));
    //     } else {
    //         _adc_sum_beam_energy_str = std::to_string(int(config_beam_energies[_run_index]));
    //     }
    //     std::string _adc_sum_beam_energy_dummy_str = "     ";
    //     for (int _beam_energy_index = 0; _beam_energy_index < 3 - _adc_sum_beam_energy_str.size(); _beam_energy_index++) {
    //         _adc_sum_beam_energy_dummy_str += " ";
    //     }

    //     adc_sum_fit_mu_list.push_back(_adc_sum_fit_res_mean_weighted);
    //     adc_sum_fit_sigma_list.push_back(_adc_sum_fit_res_sigma_weighted);
    //     adc_sum_fit_mu_error_list.push_back(sqrt(_adc_sum_fit_res_mean_err_stat * _adc_sum_fit_res_mean_err_stat + _adc_sum_fit_res_mean_err_sys * _adc_sum_fit_res_mean_err_sys));
    //     adc_sum_fit_sigma_error_list.push_back(sqrt(_adc_sum_fit_res_sigma_err_stat * _adc_sum_fit_res_sigma_err_stat + _adc_sum_fit_res_sigma_err_sys * _adc_sum_fit_res_sigma_err_sys));
    //     adc_sum_beam_energy_list.push_back(config_beam_energies[_run_index]);
    //     adc_sum_beam_energy_error_list.push_back(config_beam_energies[_run_index] * beam_energy_relative_error);

    //     // * write the results to legend
    //     std::string _adc_sum_fit_res_mean_str = (_adc_sum_beam_energy_str + "GeV mu: " + _adc_sum_mean_str + " #pm " + _adc_sum_mean_err_stat_str + "(stat) #pm " + _adc_sum_mean_err_sys_str + "(sys)").c_str();
    //     std::string _adc_sum_fit_res_sigma_str = (_adc_sum_beam_energy_dummy_str + "sigma: " + _adc_sum_sigma_str + " #pm " + _adc_sum_sigma_err_stat_str + "(stat) #pm " + _adc_sum_sigma_err_sys_str + "(sys)").c_str();
    //     auto dummy_hist = new TH1D("dummy_hist", "dummy_hist", 1, 0, 1);
    //     // no line for the dummy hist
    //     dummy_hist->SetLineColor(kWhite);
    //     adc_sum_legend->AddEntry(_adc_sum_hist1d, _adc_sum_fit_res_mean_str.c_str(), "l");
    //     adc_sum_legend->AddEntry(dummy_hist, _adc_sum_fit_res_sigma_str.c_str(), "l");
        
    //     // adc_sum_legend->AddEntry(_adc_sum_hist1d, (std::to_string(int(config_beam_energies[_run_index])) + " GeV " + config_beam_particles[_run_index]).c_str(), "l");
    // }
    // adc_sum_legend->Draw();
    // auto adc_sum_latex = new TLatex();
    // const double _text_line_height = 0.04;
    // const double _text_line_start = 0.85;
    // const double _text_line_left = 0.13;
    // adc_sum_latex->SetNDC();
    // adc_sum_latex->SetTextSize(0.04);
    // adc_sum_latex->SetTextFont(62);
    // adc_sum_latex->DrawLatex(_text_line_left, _text_line_start, (config_plot_info[0].c_str()));
    // adc_sum_latex->SetTextSize(0.03);
    // adc_sum_latex->SetTextFont(42);
    // for (int _info_index = 1; _info_index < config_plot_info.size(); _info_index++) {
    //     adc_sum_latex->DrawLatex(_text_line_left, _text_line_start - _info_index * _text_line_height, (config_plot_info[_info_index].c_str()));
    // }
    // if (enable_working_in_progress){
    //     adc_sum_latex->SetTextFont(52);
    //     adc_sum_latex->SetTextColor(kGray+3);
    //     adc_sum_latex->DrawLatex(_text_line_left, _text_line_start - (config_plot_info.size()) * _text_line_height, "Work in progress");
    // }
    
    // output_root->cd();
    // adc_sum_hist1d_canvas->Print(pdf_file_name.c_str());
    // adc_sum_hist1d_canvas->Write();

    // // * --- Draw linearity plot --------------------------------------------------------
    // // * --------------------------------------------------------------------------------
    // auto linearity_canvas = new TCanvas("linearity_canvas", "linearity_canvas", 800, 600);
    // linearity_canvas->SetMargin(0.15, 0.1, 0.15, 0.1);
    // auto linearity_dummy_graph = new TGraphErrors(1); // only for setting the axis
    // linearity_dummy_graph->SetTitle("Linearity of ADC Sum");
    // linearity_dummy_graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    // linearity_dummy_graph->GetYaxis()->SetTitle("ADC Sum Mean");
    // linearity_dummy_graph->GetXaxis()->SetRangeUser(0, 400);
    // linearity_dummy_graph->GetXaxis()->SetLimits(0, 400);
    // linearity_dummy_graph->GetYaxis()->SetRangeUser(0, adc_sum_hist_max);
    // linearity_dummy_graph->GetYaxis()->SetLimits(0, adc_sum_hist_max);  

    // linearity_dummy_graph->Draw("AP");

    // auto linearity_graph = new TGraphErrors(adc_sum_beam_energy_list.size(), &adc_sum_beam_energy_list[0], &adc_sum_fit_mu_list[0], &adc_sum_beam_energy_error_list[0], &adc_sum_fit_mu_error_list[0]);
    // linearity_graph->SetMarkerStyle(20);
    // linearity_graph->SetMarkerSize(0.5);
    // linearity_graph->SetLineWidth(2);

    // auto linearity_fit_function = new TF1("linearity_fit_function", "[0]*x + [1]", 0, 400);
    // linearity_fit_function->SetParameter(0, 400.0/adc_sum_hist_max);
    // linearity_fit_function->SetParameter(1, 0);
    // linearity_graph->Fit(linearity_fit_function, "RQN");
    // linearity_graph->Draw("PE SAME");

    // linearity_fit_function->SetLineColor(kCyan+3);
    // linearity_fit_function->SetLineStyle(2);
    // linearity_fit_function->Draw("SAME");

    // auto linearity_fit_confidence_band = new TH1F("linearity_fit_confidence_band", "linearity_fit_confidence_band", 100, 0, 400);
    // TVirtualFitter::GetFitter()->GetConfidenceIntervals(linearity_fit_confidence_band, 0.68);
    // linearity_fit_confidence_band->SetFillColorAlpha(kCyan+3, 0.5);
    // linearity_fit_confidence_band->SetMarkerColorAlpha(kCyan+3, 0.0);
    // linearity_fit_confidence_band->Draw("E3 SAME");

    // auto linearity_legend = new TLegend(0.5, 0.7, 0.89, 0.89);
    // linearity_legend->SetBorderSize(0);
    // linearity_legend->SetFillStyle(0);
    // linearity_legend->SetTextSize(0.02);
    // linearity_legend->SetTextFont(102);

    // double _linearity_fit_slope = linearity_fit_function->GetParameter(0);
    // double _linearity_fit_slope_error = linearity_fit_function->GetParError(0);
    // double _linearity_fit_intercept = linearity_fit_function->GetParameter(1);
    // double _linearity_fit_intercept_error = linearity_fit_function->GetParError(1);
    // double _linearity_fit_chi2 = linearity_fit_function->GetChisquare();
    // double _linearity_fit_ndf = linearity_fit_function->GetNDF();

    // double nonlinearity = 0;
    // for (int _beam_energy_index = 0; _beam_energy_index < adc_sum_beam_energy_list.size(); _beam_energy_index++) {
    //     double _beam_energy = adc_sum_beam_energy_list[_beam_energy_index];
    //     double _adc_sum_mu = adc_sum_fit_mu_list[_beam_energy_index];
    //     double _reconstructed_adc_sum_mu = _linearity_fit_slope * _beam_energy + _linearity_fit_intercept;
    //     double _nonlinearity = (_adc_sum_mu - _reconstructed_adc_sum_mu) / _reconstructed_adc_sum_mu;
    //     if (_nonlinearity > nonlinearity) {
    //         nonlinearity = _nonlinearity;
    //     }
    // }

    // nonlinearity = nonlinearity * 100;


    // auto linearity_dummy_hist = new TH1D("linearity_dummy_hist", "linearity_dummy_hist", 1, 0, 1);
    // linearity_dummy_hist->SetLineColor(kWhite);
    // linearity_legend->AddEntry(linearity_graph, "Gaussian Fit Mean", "ep");
    // linearity_legend->AddEntry(linearity_fit_function, ("Fit: slope     = " + std::to_string(_linearity_fit_slope).substr(0,6) + " #pm " + std::to_string(_linearity_fit_slope_error).substr(0,3)).c_str(), "l");
    // linearity_legend->AddEntry(linearity_dummy_hist, ("     intercept = " + std::to_string(_linearity_fit_intercept).substr(0,6) + " #pm " + std::to_string(_linearity_fit_intercept_error).substr(0,3)).c_str(), "l");
    // linearity_legend->AddEntry(linearity_dummy_hist, ("     #chi^{2}/NDF     = " + std::to_string(_linearity_fit_chi2).substr(0,4) + "/" + std::to_string(int(_linearity_fit_ndf))).c_str(), "l");
    // linearity_legend->AddEntry(linearity_dummy_hist, ("     Nonlinearity = " + std::to_string(nonlinearity).substr(0,4) + "%").c_str(), "l");
    
    // linearity_legend->Draw();

    // auto linearity_latex = new TLatex();
    // linearity_latex->SetNDC();
    // linearity_latex->SetTextSize(0.04);
    // linearity_latex->SetTextFont(62);
    // linearity_latex->DrawLatex(_text_line_left+0.05, _text_line_start, (config_plot_info[0].c_str()));
    // linearity_latex->SetTextSize(0.03);
    // linearity_latex->SetTextFont(42);
    // for (int _info_index = 1; _info_index < config_plot_info.size(); _info_index++) {
    //     linearity_latex->DrawLatex(_text_line_left+0.05, _text_line_start - _info_index * _text_line_height, (config_plot_info[_info_index].c_str()));
    // }
    // if (enable_working_in_progress){
    //     linearity_latex->SetTextFont(52);
    //     linearity_latex->SetTextColor(kGray+3);
    //     linearity_latex->DrawLatex(_text_line_left+0.05, _text_line_start - (config_plot_info.size()) * _text_line_height, "Work in progress");
    // }

    // linearity_canvas->Print(pdf_file_name.c_str());
    // linearity_canvas->Write();

    // // * --- Draw resolution plot --------------------------------------------------------
    // // * --------------------------------------------------------------------------------
    // auto resolution_canvas = new TCanvas("resolution_canvas", "resolution_canvas", 800, 600);
    // auto resolution_dummy_graph = new TGraphErrors(1); // only for setting the axis
    // resolution_dummy_graph->SetTitle("Resolution of ADC Sum");
    // resolution_dummy_graph->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    // resolution_dummy_graph->GetYaxis()->SetTitle("#sigma_{E}/E");
    // resolution_dummy_graph->GetXaxis()->SetRangeUser(0, 400);
    // resolution_dummy_graph->GetYaxis()->SetRangeUser(0.1, 0.4);
    // resolution_dummy_graph->GetYaxis()->SetLimits(0.1, 0.4);
    // resolution_dummy_graph->GetXaxis()->SetLimits(0, 400);
    // resolution_dummy_graph->Draw("AP");

    // std::vector <double> adc_sum_fit_sigmaE_E_list;
    // std::vector <double> adc_sum_fit_sigmaE_E_error_list;
    // std::vector <double> reconstructed_beam_energy_list;
    // std::vector <double> reconstructed_beam_energy_list_error;
    // for (int _beam_energy_index = 0; _beam_energy_index < adc_sum_beam_energy_list.size(); _beam_energy_index++) {
    //     double _beam_energy = adc_sum_beam_energy_list[_beam_energy_index];
    //     double _adc_sum_sigma = adc_sum_fit_sigma_list[_beam_energy_index];
    //     double _adc_sum_sigma_error = adc_sum_fit_sigma_error_list[_beam_energy_index];
    //     double _adc_sum_mu = adc_sum_fit_mu_list[_beam_energy_index];
    //     double _adc_sum_mu_error = adc_sum_fit_mu_error_list[_beam_energy_index];
    //     double _adc_sum_sigmaE_E = _adc_sum_sigma / _adc_sum_mu;
    //     double _adc_sum_sigmaE_E_error = _adc_sum_sigmaE_E * std::sqrt((_adc_sum_sigma_error / _adc_sum_sigma) * (_adc_sum_sigma_error / _adc_sum_sigma) + (_adc_sum_mu_error / _adc_sum_mu) * (_adc_sum_mu_error / _adc_sum_mu));
    //     adc_sum_fit_sigmaE_E_list.push_back(_adc_sum_sigmaE_E);
    //     adc_sum_fit_sigmaE_E_error_list.push_back(_adc_sum_sigmaE_E_error);
    //     double _reconstructed_beam_energy = (_adc_sum_mu - _linearity_fit_intercept)/_linearity_fit_slope;
    //     double _reconstructed_beam_energy_error = _reconstructed_beam_energy * beam_energy_relative_error;
    //     reconstructed_beam_energy_list.push_back(_reconstructed_beam_energy);
    //     reconstructed_beam_energy_list_error.push_back(_reconstructed_beam_energy_error);
    // }
    // auto resolution_legend = new TLegend(0.5, 0.65, 0.89, 0.89);
    // resolution_legend->SetBorderSize(0);
    // resolution_legend->SetFillStyle(0);
    // resolution_legend->SetTextSize(0.02);
    // resolution_legend->SetTextFont(102);

    // // * draw reference resolution plots
    // std::vector<double> _resolution_mc_resolution = {0.164829,0.157141,0.147176,0.132711,0.121076,0.113879,0.108920,0.105994};
    // std::vector<double> _resolution_mc_energy = {60, 80, 100, 150, 200, 250, 300, 350};
    // std::vector<double> _resolution_mc_energy_error = {0,0,0,0,0,0,0,0};
    // std::vector<double> _resolution_mc_resolution_error = {0,0,0,0,0,0,0,0};
    // auto resolution_mc_graph = new TGraphErrors(_resolution_mc_energy.size(), &_resolution_mc_energy[0], &_resolution_mc_resolution[0], &_resolution_mc_energy_error[0], &_resolution_mc_resolution_error[0]);
    // // resolution_mc_graph->SetTitle("Resolution of ADC Sum");
    // // resolution_mc_graph->GetXaxis()->SetTitle("Beam Energy [GeV]");
    // // resolution_mc_graph->GetXaxis()->SetRangeUser(0, 400);
    // // resolution_mc_graph->GetYaxis()->SetTitle("#sigma_{E}/E");
    // // resolution_mc_graph->GetYaxis()->SetRangeUser(0.1, 0.4);
    // resolution_mc_graph->SetMarkerStyle(20);
    // resolution_mc_graph->SetMarkerSize(0.5);
    // resolution_mc_graph->SetMarkerColor(kGreen+3);
    // resolution_mc_graph->SetLineWidth(2);
    // resolution_mc_graph->SetLineColor(kGreen+3);
    // resolution_mc_graph->Draw("PE same");

    // resolution_legend->AddEntry(resolution_mc_graph, "Geant4 Simulation", "ep");

    // auto resolution_mc_fit_function = new TF1("resolution_mc_fit_function", "[0]/sqrt(x) + [1]", 0, 400);
    // resolution_mc_fit_function->SetParameter(0, 0.1);
    // resolution_mc_fit_function->SetParameter(1, 0.01);
    // resolution_mc_graph->Fit(resolution_mc_fit_function, "RQN");
    // resolution_mc_fit_function->SetLineColor(kGreen+3);
    // resolution_mc_fit_function->SetLineStyle(2);
    // resolution_mc_fit_function->Draw("SAME");

    // auto resolution_mc_fit_confidence_band = new TH1F("resolution_mc_fit_confidence_band", "resolution_mc_fit_confidence_band", 400, 0, 400);
    // TVirtualFitter *fitter_mc = TVirtualFitter::GetFitter();
    // if (fitter_mc) {
    //     fitter_mc->GetConfidenceIntervals(resolution_mc_fit_confidence_band, 0.68);
    // } else {
    //     LOG(WARNING) << "Fitter is not initialized. Confidence intervals will not be drawn.";
    // }
    // resolution_mc_fit_confidence_band->SetFillColorAlpha(kGreen+3, 0.3);
    // resolution_mc_fit_confidence_band->SetLineColor(kGreen+3);
    // resolution_mc_fit_confidence_band->SetMarkerColorAlpha(kGreen+3, 0.0);
    // resolution_mc_fit_confidence_band->SetLineStyle(2);
    // resolution_mc_fit_confidence_band->Draw("e3 SAME");

    // resolution_legend->AddEntry(resolution_mc_fit_function, ("#sigma_{E}/E   = #frac{" + std::to_string(resolution_mc_fit_function->GetParameter(0)).substr(0,4) + " #pm " + std::to_string(resolution_mc_fit_function->GetParError(0)).substr(0,3) + "}{#sqrt{E}} #oplus (" + std::to_string(resolution_mc_fit_function->GetParameter(1)).substr(0,4) + " #pm " + std::to_string(resolution_mc_fit_function->GetParError(1)).substr(0,4) + ")").c_str(), "l");

    // std::vector<double> _resolution_caen_resolution = {0.202903,0.193606,0.180118,0.157828,0.143362,0.135984,0.128683,0.124429};
    // std::vector<double> _resolution_caen_energy = {60, 80, 100, 150, 200, 250, 300, 350};
    // std::vector<double> _resolution_caen_energy_error = {1.2, 1.6, 2, 3, 4, 5, 6, 7};
    // std::vector<double> _resolution_caen_resolution_error = {0.0223685347754385,
    //     0.0276764080400618,
    //     0.0256984028686609,
    //     0.0235984254347615,
    //     0.016989855708628,
    //     0.0151618517338747,
    //     0.0109004483394033,
    //     0.0102402795860269
    // };

    // auto resolution_caen_graph = new TGraphErrors(_resolution_caen_energy.size(), &_resolution_caen_energy[0], &_resolution_caen_resolution[0], &_resolution_caen_energy_error[0], &_resolution_caen_resolution_error[0]);
    // resolution_caen_graph->SetMarkerStyle(20);
    // resolution_caen_graph->SetMarkerSize(0.5);
    // resolution_caen_graph->SetMarkerColor(kOrange - 7);
    // resolution_caen_graph->SetLineWidth(2);
    // resolution_caen_graph->SetLineColor(kOrange - 7);
    // resolution_caen_graph->Draw("PE same");

    // resolution_legend->AddEntry(resolution_caen_graph, "CAEN Data", "ep");

    // auto resolution_caen_fit_function = new TF1("resolution_caen_fit_function", "[0]/sqrt(x) + [1]", 0, 400);
    // resolution_caen_fit_function->SetParameter(0, 0.1);
    // resolution_caen_fit_function->SetParameter(1, 0.01);
    // resolution_caen_graph->Fit(resolution_caen_fit_function, "RQN");
    // resolution_caen_graph->Draw("PE SAME");

    // resolution_caen_fit_function->SetLineColor(kOrange - 7);
    // resolution_caen_fit_function->SetLineStyle(2);
    // resolution_caen_fit_function->Draw("SAME");

    // resolution_legend->AddEntry(resolution_caen_graph, ("#sigma_{E}/E   = #frac{" + std::to_string(resolution_caen_fit_function->GetParameter(0)).substr(0,4) + " #pm " + std::to_string(resolution_caen_fit_function->GetParError(0)).substr(0,3) + "}{#sqrt{E}} #oplus (" + std::to_string(resolution_caen_fit_function->GetParameter(1)).substr(0,4) + " #pm " + std::to_string(resolution_caen_fit_function->GetParError(1)).substr(0,4) + ")").c_str(), "l");

    // auto resolution_caen_fit_confidence_band = new TH1F("resolution_caen_fit_confidence_band", "resolution_caen_fit_confidence_band", 400, 0, 400);
    // TVirtualFitter *fitter_caen = TVirtualFitter::GetFitter();
    // if (fitter_caen) {
    //     fitter_caen->GetConfidenceIntervals(resolution_caen_fit_confidence_band, 0.68);
    // } else {
    //     LOG(WARNING) << "Fitter is not initialized. Confidence intervals will not be drawn.";
    // }
    // resolution_caen_fit_confidence_band->SetFillColorAlpha(kOrange - 7, 0.3);
    // resolution_caen_fit_confidence_band->SetMarkerColorAlpha(kOrange - 7, 0.0);
    // resolution_caen_fit_confidence_band->SetLineColor(kOrange - 7);
    // resolution_caen_fit_confidence_band->SetLineStyle(2);
    // resolution_caen_fit_confidence_band->Draw("e3 SAME");

    // auto resolution_graph = new TGraphErrors(reconstructed_beam_energy_list.size(), &reconstructed_beam_energy_list[0], &adc_sum_fit_sigmaE_E_list[0], &reconstructed_beam_energy_list_error[0], &adc_sum_fit_sigmaE_E_error_list[0]);
    // resolution_graph->SetMarkerStyle(20);
    // resolution_graph->SetMarkerSize(0.5);
    // resolution_graph->SetMarkerColor(kMagenta - 3);
    // resolution_graph->SetLineWidth(2);
    // resolution_graph->SetLineColor(kMagenta - 3);
    
    // auto resolution_fit_function = new TF1("resolution_fit_function", "[0]/sqrt(x) + [1]", 0, 400);
    // resolution_fit_function->SetParameter(0, 0.1);
    // resolution_fit_function->SetParameter(1, 0.01);
    // resolution_graph->Fit(resolution_fit_function, "RQN");
    // resolution_graph->Draw("PE SAME");

    // resolution_fit_function->SetLineColor(kMagenta - 3);
    // resolution_fit_function->SetLineStyle(2);
    // resolution_fit_function->Draw("SAME");

    // resolution_legend->AddEntry(resolution_graph, "H2GCROC3A Data", "ep");

    // auto resolution_fit_confidence_band = new TH1F("resolution_fit_confidence_band", "resolution_fit_confidence_band", 400, 0, 400);
    // TVirtualFitter *fitter = TVirtualFitter::GetFitter();
    // if (fitter) {
    //     fitter->GetConfidenceIntervals(resolution_fit_confidence_band, 0.68);
    // } else {
    //     LOG(WARNING) << "Fitter is not initialized. Confidence intervals will not be drawn.";
    // }
    // resolution_fit_confidence_band->SetFillColorAlpha(kMagenta - 3, 0.3);
    // resolution_fit_confidence_band->SetLineColor(kMagenta - 3);
    // resolution_fit_confidence_band->SetMarkerColorAlpha(kMagenta - 3, 0.0);
    // resolution_fit_confidence_band->SetLineStyle(2);
    // resolution_fit_confidence_band->Draw("e3 SAME");

    // double _resolution_fit_a = resolution_fit_function->GetParameter(0);
    // double _resolution_fit_a_error = resolution_fit_function->GetParError(0);
    // double _resolution_fit_b = resolution_fit_function->GetParameter(1);
    // double _resolution_fit_b_error = resolution_fit_function->GetParError(1);
    // double _resolution_fit_chi2 = resolution_fit_function->GetChisquare();
    // double _resolution_fit_ndf = resolution_fit_function->GetNDF();

    // auto resolution_dummy_hist = new TH1D("resolution_dummy_hist", "resolution_dummy_hist", 1, 0, 1);
    // resolution_dummy_hist->SetLineColor(kWhite);

    // resolution_legend->AddEntry(resolution_fit_function, ("#sigma_{E}/E   = #frac{" + std::to_string(_resolution_fit_a).substr(0,4) + " #pm " + std::to_string(_resolution_fit_a_error).substr(0,3) + "}{#sqrt{E}} #oplus (" + std::to_string(_resolution_fit_b).substr(0,4) + " #pm " + std::to_string(_resolution_fit_b_error).substr(0,4) + ")").c_str(), "l");
    // // resolution_legend->AddEntry(resolution_dummy_hist, ("#chi^{2}/NDF = " + std::to_string(_resolution_fit_chi2).substr(0,4) + "/" + std::to_string(int(_resolution_fit_ndf))).c_str(), "l");
    
    // resolution_legend->Draw();

    // auto resolution_latex = new TLatex();
    // resolution_latex->SetNDC();
    // resolution_latex->SetTextSize(0.04);
    // resolution_latex->SetTextFont(62);
    // resolution_latex->DrawLatex(_text_line_left, _text_line_start, (config_plot_info[0].c_str()));
    // resolution_latex->SetTextSize(0.03);
    // resolution_latex->SetTextFont(42);
    // for (int _info_index = 1; _info_index < config_plot_info.size(); _info_index++) {
    //     resolution_latex->DrawLatex(_text_line_left, _text_line_start - _info_index * _text_line_height, (config_plot_info[_info_index].c_str()));
    // }

    // if (enable_working_in_progress){
    //     resolution_latex->SetTextFont(52);
    //     resolution_latex->SetTextColor(kGray+3);
    //     resolution_latex->DrawLatex(_text_line_left, _text_line_start - (config_plot_info.size()) * _text_line_height, "Work in progress");
    // }

    // resolution_canvas->Print(pdf_file_name.c_str());
    // resolution_canvas->Write();

    // // * --- Close the files ------------------------------------------------------------
    // // * --------------------------------------------------------------------------------
    // TCanvas* dummy_canvas = new TCanvas("dummy_canvas", "dummy_canvas", 800, 600);
    // dummy_canvas->SaveAs((pdf_file_name + ")").c_str());
    // dummy_canvas->Close();

    // output_root->Close();
    // LOG(INFO) << "Output file " << opts.output_file << " has been saved.";
    // LOG(INFO) << "PDF file " << pdf_file_name << " has been saved.";

    return 0;
}