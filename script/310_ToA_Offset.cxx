#include <fstream>
#include <limits>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"

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
        ("c,correction", "Input ToA correction .root file", cxxopts::value<std::string>())
        ("d,dnl", "Input DNL .root file", cxxopts::value<std::string>())
        ("t,timewalk", "Input Timewalk correction .root file", cxxopts::value<std::string>()->default_value("dump/ToA_Timewalk_Correction.root"))
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
    spdlog::info("Input bx correction file: {}", parsed["correction"].as<std::string>());
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
    int chn_example = 299; // print this channel seprately

    TFile *input_correction_root = new TFile(parsed["correction"].as<std::string>().c_str(), "READ");
    if (input_correction_root->IsZombie()) {
        spdlog::error("Failed to open input correction file {}", parsed["correction"].as<std::string>());
        return 1;
    }
    // read all channel correction optimal values
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
    // go to "Optimal_Bx_Slip_Parameters" directory
    input_correction_root->cd("Optimal_Bx_Slip_Parameters");
    TDirectory* param_dir = gDirectory;
    std::map<int, double> toa_correction_map;
    for (int chn=0; chn<FPGA_CHANNEL_NUMBER * vldb_number; chn++) {
        std::string param_name = "optimal_bx_slip_ch" + std::to_string(chn);
        double optimal_value = getParam(param_dir, param_name.c_str());
        toa_correction_map[chn] = optimal_value;
        spdlog::info("Channel {}: optimal_bx_slip = {}", chn, optimal_value);
    }
    input_correction_root->Close();

    TFile *input_timewalk_root = new TFile(parsed["timewalk"].as<std::string>().c_str(), "READ");
    if (input_timewalk_root->IsZombie()) {
        spdlog::error("Failed to open input timewalk file {}", parsed["timewalk"].as<std::string>());
        return 1;
    }
    input_timewalk_root->cd("Optimal_ToA_ADC_Max_Correction_Parameters");
    param_dir = gDirectory;
    std::map<int, double> timewalk_param0_map;
    std::map<int, double> timewalk_param1_map;
    std::map<int, double> timewalk_param2_map;
    std::map<int, double> timewalk_param3_map;
    for (int chn=0; chn<FPGA_CHANNEL_NUMBER * vldb_number; chn++) {
        std::string param_name0 = "toa_adc_max_corr_parm0_ch" + std::to_string(chn);
        std::string param_name1 = "toa_adc_max_corr_parm1_ch" + std::to_string(chn);
        std::string param_name2 = "toa_adc_max_corr_parm2_ch" + std::to_string(chn);
        std::string param_name3 = "toa_adc_max_corr_parm3_ch" + std::to_string(chn);
        double parm0 = getParam(param_dir, param_name0.c_str());
        double parm1 = getParam(param_dir, param_name1.c_str());
        double parm2 = getParam(param_dir, param_name2.c_str());
        double parm3 = getParam(param_dir, param_name3.c_str());
        timewalk_param0_map[chn] = parm0;
        timewalk_param1_map[chn] = parm1;
        timewalk_param2_map[chn] = parm2;
        timewalk_param3_map[chn] = parm3;
        spdlog::info("Channel {}: ToA ADC Max Correction Parameters = {}, {}, {}, {}",
                     chn, parm0, parm1, parm2, parm3);
    }
    input_timewalk_root->Close();

    TFile *input_dnl_root = new TFile(parsed["dnl"].as<std::string>().c_str(), "READ");
    if (input_dnl_root->IsZombie()) {
        spdlog::error("Failed to open input DNL file {}", parsed["dnl"].as<std::string>());
        return 1;
    }

    std::vector<std::vector<double>> dnl_coarse_tdc_list(vldb_number* FPGA_CHANNEL_NUMBER);
    std::vector<std::vector<double>> dnl_fine_tdc_list(vldb_number* FPGA_CHANNEL_NUMBER);
    input_dnl_root->cd("DNL_Calibration");
    for (int chn=0; chn<FPGA_CHANNEL_NUMBER * vldb_number; chn++) {
        std::string dnl_coarse_name = "dnl_coarse_tdc_ch" + std::to_string(chn);
        TVectorD* dnl_coarse_tdc = nullptr;
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
        dnl_coarse_tdc = gDirectory->Get<TVectorD>(dnl_coarse_name.c_str());
    #else
        gDirectory->GetObject(dnl_coarse_name.c_str(), dnl_coarse_tdc);
    #endif
        if (!dnl_coarse_tdc) {
            spdlog::error("Missing TVectorD key '{}' in directory '{}'",
                        dnl_coarse_name, gDirectory->GetPath());
            return 1;
        }
        std::vector<double> coarse_dnl(dnl_coarse_tdc->GetNoElements());
        for (int idx = 0; idx < dnl_coarse_tdc->GetNoElements(); ++idx) {
            coarse_dnl[idx] = (*dnl_coarse_tdc)[idx];
        }
        dnl_coarse_tdc_list[chn] = std::move(coarse_dnl);
        std::string dnl_fine_name = "dnl_fine_tdc_ch" + std::to_string(chn);
        TVectorD* dnl_fine_tdc = nullptr;
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
        dnl_fine_tdc = gDirectory->Get<TVectorD>(dnl_fine_name.c_str());
    #else
        gDirectory->GetObject(dnl_fine_name.c_str(), dnl_fine_tdc);
    #endif
        if (!dnl_fine_tdc) {
            spdlog::error("Missing TVectorD key '{}' in directory '{}'",
                        dnl_fine_name, gDirectory->GetPath());
            return 1;
        }
        std::vector<double> fine_dnl(dnl_fine_tdc->GetNoElements());
        for (int idx = 0; idx < dnl_fine_tdc->GetNoElements(); ++idx) {
            fine_dnl[idx] = (*dnl_fine_tdc)[idx];
        }
        dnl_fine_tdc_list[chn] = std::move(fine_dnl);
        // spdlog::info("Channel {}: DNL coarse TDC size = {}, fine TDC size = {}",
        //              chn,
        //              dnl_coarse_tdc->GetNoElements(),
        //              dnl_fine_tdc->GetNoElements());
        // spdlog::info("Coarse DNL: ");
        // for (const auto& val : dnl_coarse_tdc_list[chn]) {
        //     spdlog::info("{:.6f}, ", val);
        // }
        // spdlog::info("Fine DNL: ");
        // for (const auto& val : dnl_fine_tdc_list[chn]) {
        //     spdlog::info("{:.6f}, ", val);
        // }

    }
    input_dnl_root->Close();

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

    // std::vector<std::vector<UInt_t>> toa_code_matrix(FPGA_CHANNEL_NUMBER * vldb_number, std::vector<UInt_t>());
    // for (auto & vec : toa_code_matrix) {
    //     vec.reserve(entry_max);
    // }
    // std::vector<std::vector<double>> toa_ns_matrix(FPGA_CHANNEL_NUMBER * vldb_number, std::vector<double>());
    // for (auto & vec : toa_ns_matrix) {
    //     vec.reserve(entry_max);
    // }

    const int event_display_number = 20; // show the first 10 events for example

    const int NX = 16, NY = 12;
    const int board_cols = 8, board_rows = 4;
    const int NPIX = NX * NY;
    const int TOTAL_CH = FPGA_CHANNEL_NUMBER * vldb_number;

    std::vector<int> ch2pid(TOTAL_CH, -1);
    std::vector<int> pid2ch(NPIX, -1);

    std::vector<TH2D*> event_display_adc_hitmap;
    std::vector<TH2D*> event_display_toa_hitmap;
    for (int ev=0; ev<event_display_number; ev++) {
        std::string adc_hist_name = "event_" + std::to_string(ev) + "_adc_hitmap";
        TH2D *adc_hitmap = new TH2D(adc_hist_name.c_str(), (adc_hist_name + ";Channel X;Channel Y;ADC Peak Value (pedestal subtracted)").c_str(),
                                    NX, 0, NX,
                                    NY, 0, NY);
        adc_hitmap->SetDirectory(nullptr);
        event_display_adc_hitmap.push_back(adc_hitmap);

        std::string toa_hist_name = "event_" + std::to_string(ev) + "_toa_hitmap";
        TH2D *toa_hitmap = new TH2D(toa_hist_name.c_str(), (toa_hist_name + ";Channel X;Channel Y;ToA (ns)").c_str(),
                                    NX, 0, NX,
                                    NY, 0, NY);
        toa_hitmap->SetDirectory(nullptr);
        event_display_toa_hitmap.push_back(toa_hitmap);
    }

    TH1D *h1d_event_average_toa = new TH1D("h1d_event_average_toa", "h1d_event_average_toa;Average ToA per Event (ns);Counts", 256, 0, 25.0*static_cast<double>(machine_gun_samples));
    h1d_event_average_toa->SetDirectory(nullptr);

    const int adc_peak_min_index = 4;
    const int adc_peak_max_index = 8;

    std::vector<std::vector<double>> event_chn_adc_matrix;
    std::vector<std::vector<double>> event_chn_toa_matrix;
    event_chn_adc_matrix.resize(entry_max, std::vector<double>(FPGA_CHANNEL_NUMBER * vldb_number, NAN)); // event -> channel
    event_chn_toa_matrix.resize(entry_max, std::vector<double>(FPGA_CHANNEL_NUMBER * vldb_number, NAN)); // event -> channel

    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
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

                if (uni_col < 0 || uni_col >= NX || uni_row < 0 || uni_row >= NY)
                    continue;

                const int pid = uni_row*NX + uni_col;
                const int gch = channel + vldb_id*FPGA_CHANNEL_NUMBER;

                ch2pid[gch] = pid;
                // 若一个 pid 对应多个通道（少见），按你的硬件定义选择其一；这里只保留第一次
                if (pid2ch[pid] == -1) pid2ch[pid] = gch;
        }
    }

    const int DX[4] = {+1,-1,0,0};
    const int DY[4] = {0,0,+1,-1};
    std::vector<std::array<int,4>> neigh(NPIX);
    std::vector<int> deg(NPIX, 0);

    for (int r=0; r<NY; ++r){
        for (int c=0; c<NX; ++c){
            int p = r*NX + c;
            int k = 0;
            for (int m=0; m<4; ++m){
                int cc = c + DX[m], rr = r + DY[m];
                int q = (rr>=0 && rr<NY && cc>=0 && cc<NX) ? rr*NX + cc : -1;
                neigh[p][m] = q;
                if (q!=-1) ++k;
            }
            deg[p] = k;
        }
    }

    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        double adc_sum = 0.0;
        double toa_weighted_sum = 0.0;
        double toa_weight_total = 0.0;
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            // channel loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                
                int gch = channel + vldb_id * FPGA_CHANNEL_NUMBER;
                int pid = ch2pid[gch];
                if (pid == -1)
                    continue;
                int uni_col = pid % NX;
                int uni_row = pid / NX;

                if (uni_col < 0 || uni_col >= NX || uni_row < 0 || uni_row >= NY)
                    continue;

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
                double adc_pedestal = pedestal_median_of_first3(adc_pedestal_samples);
                double adc_peak_value_pede_sub = static_cast<double>(adc_peak_ranged_value) - adc_pedestal;
                adc_sum += adc_peak_value_pede_sub;

                h1i_toa_frequency->Fill(tot_showup_times);
                if (toa_showup_times >= 1) {
                    num_single_ToA++;
                    // only one ToA hit
                    // apply correction here if needed
                    auto& coarse_dnl_lookup = dnl_coarse_tdc_list[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    auto& fine_dnl_lookup = dnl_fine_tdc_list[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    // double toa_ns = decode_toa_value_ns(toa_first_value);
                    double toa_ns = decode_toa_value_ns_with_dnl(
                        toa_first_value,
                        coarse_dnl_lookup,
                        fine_dnl_lookup
                    );
                    toa_ns += 25.0 * static_cast<double>(toa_first_index);
                    // do the bx slip correction
                    double channel_optimal_bx_slip = toa_correction_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    if (toa_first_value >= static_cast<UInt_t>(channel_optimal_bx_slip)) {
                        toa_ns -= 25.0;
                    }
                    // do the timewalk correction
                    double param0 = timewalk_param0_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    double param1 = timewalk_param1_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    double param2 = timewalk_param2_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    double param3 = timewalk_param3_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    // "[0] + [1]/pow(x-[3],[2])"
                    double timewalk_correction_ns = param1 / std::pow(adc_peak_value_pede_sub - param2, param3);
                    toa_ns -= timewalk_correction_ns;
                    h1d_raw_toa_ns_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_ns);
                    h2d_toa_adc_max_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(adc_peak_value_pede_sub, toa_ns);
                    h2d_toa_code_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_first_value, toa_ns);
                    if (toa_ns > 0.0) {
                        toa_weighted_sum += toa_ns * adc_peak_value_pede_sub;
                        toa_weight_total += adc_peak_value_pede_sub;
                        event_chn_toa_matrix[entry][channel + vldb_id * FPGA_CHANNEL_NUMBER] = toa_ns;
                        event_chn_adc_matrix[entry][channel + vldb_id * FPGA_CHANNEL_NUMBER] = adc_peak_value_pede_sub;
                    }
                    if (entry < event_display_number) {
                        event_display_adc_hitmap[entry]->Fill(uni_col + 0.5, uni_row + 0.5, adc_peak_value_pede_sub);
                        if (toa_ns > 0.0)
                            event_display_toa_hitmap[entry]->Fill(uni_col + 0.5, uni_row + 0.5, toa_ns);
                    }
                }
            } // end of channel loop
        } // end of vldb loop
        if (toa_weight_total > 0.0) {
            double event_average_toa = toa_weighted_sum / toa_weight_total;
            h1d_event_average_toa->Fill(event_average_toa);
        }
    } // end of event loop

    input_root->Close();

    spdlog::info("Number of single ToA hits: {}", num_single_ToA);
    spdlog::info("Number of multiple ToA hits: {}", num_multi_ToA);
    spdlog::info("Percentage of single ToA hits: {:.2f}%", 100.0 * static_cast<double>(num_single_ToA) / static_cast<double>(num_single_ToA + num_multi_ToA)); 
    
    
    auto output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }

    std::vector<double> toa_offsets;
    // diff_toa_offset_calculator_mapped(
    //     event_chn_toa_matrix,
    //     event_chn_adc_matrix,
    //     toa_offsets,
    //     neigh,
    //     deg,
    //     pid2ch
    // );
    diff_toa_offset_calculator_2(
        event_chn_toa_matrix,
        event_chn_adc_matrix,
        toa_offsets,
        neigh,
        deg
    );
    // create directory for ToA Offsets
    TDirectory* toa_offset_dir = output_root->mkdir("ToA_Offsets");
    toa_offset_dir->cd();
    spdlog::info("Calculated ToA offsets for each channel:");
    for (size_t chn = 0; chn < toa_offsets.size(); chn++) {
        spdlog::info("Channel {}: ToA Offset = {:.3f} ns", chn, toa_offsets[chn]);
        TParameter<double> *p_offset = new TParameter<double>(("toa_offset_ch" + std::to_string(chn)).c_str(), toa_offsets[chn]);
        p_offset->Write();
    }
    output_root->cd();

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    const std::string out_pdf = script_output_file.substr(0, script_output_file.find_last_of(".")) + ".pdf";

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

    // print the example channel
    TCanvas *canvas_toa_adc_max_corr_example = new TCanvas("canvas_toa_adc_max_corr_example", "ToA ADC Max Correction Example Channel", 800, 600);
    auto h2d_example = h2d_toa_adc_max_corr_list[chn_example];
    h2d_example->GetXaxis()->SetTitle("ADC Peak (pedestal subtracted)");
    h2d_example->GetYaxis()->SetTitle("Corrected ToA [ns]");
    format_2d_hist_canvas(canvas_toa_adc_max_corr_example, h2d_example, kBlue+2, annotation_canvas_title, annotation_testbeam_title, "Channel_" + std::to_string(chn_example));
    canvas_toa_adc_max_corr_example->Print(out_pdf.c_str());
    canvas_toa_adc_max_corr_example->Write();
    canvas_toa_adc_max_corr_example->Close();

    TCanvas *canvas_event_toa_average = new TCanvas("canvas_event_toa_average", "Event Average ToA", 800, 600);
    canvas_event_toa_average->cd();
    format_1d_hist_canvas(canvas_event_toa_average, h1d_event_average_toa, kBlue, annotation_canvas_title, annotation_testbeam_title, "Event Average ToA");
    canvas_event_toa_average->Print(out_pdf.c_str());
    canvas_event_toa_average->Write();
    canvas_event_toa_average->Close();

    TCanvas *canvas_toa_code_corr = new TCanvas("canvas_toa_code_corr", "ToA Code Correction", 1200, 800);
    draw_mosaic_fixed(*canvas_toa_code_corr, h2d_toa_code_corr_list, topo_wave);
    canvas_toa_code_corr->Modified();
    canvas_toa_code_corr->Update();
    canvas_toa_code_corr->Print(out_pdf.c_str());
    canvas_toa_code_corr->Write();
    canvas_toa_code_corr->Close();

    for (int ev=0; ev<event_display_number; ev++) {
        TCanvas *canvas_event_adc = new TCanvas(("canvas_event_" + std::to_string(ev) + "_adc").c_str(),
                                                ("Event " + std::to_string(ev) + " ADC Hitmap").c_str(),
                                                800, 600);
        auto h2d_adc = event_display_adc_hitmap[ev];
        
        h2d_adc->GetXaxis()->SetTitle("X");
        h2d_adc->GetYaxis()->SetTitle("Y");
        format_2d_hist_canvas(canvas_event_adc, h2d_adc, kBlue+2, annotation_canvas_title, annotation_testbeam_title, "Event_" + std::to_string(ev), false);
        canvas_event_adc->Print(out_pdf.c_str());
        canvas_event_adc->Write();
        canvas_event_adc->Close();

        TCanvas *canvas_event_toa = new TCanvas(("canvas_event_" + std::to_string(ev) + "_toa").c_str(),
                                                ("Event " + std::to_string(ev) + " ToA Hitmap").c_str(),
                                                800, 600);
        auto h2d_toa = event_display_toa_hitmap[ev];
        // fill the bins with the offset
        for (int ix = 1; ix <= h2d_toa->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= h2d_toa->GetNbinsY(); iy++) {
                int pid = (iy - 1) * NX + (ix - 1);
                int chn = pid2ch[pid];
                if (chn >= 0 && chn < static_cast<int>(toa_offsets.size())) {
                    double offset = toa_offsets[chn];
                    double original_val = h2d_toa->GetBinContent(ix, iy);
                    if (original_val > 0.0) {
                        h2d_toa->SetBinContent(ix, iy, original_val - offset);
                    }
                }
            }
        }
        h2d_toa->GetXaxis()->SetTitle("X");
        h2d_toa->GetYaxis()->SetTitle("Y");
        format_2d_hist_canvas(canvas_event_toa, h2d_toa, kBlue+2, annotation_canvas_title, annotation_testbeam_title, "Event_" + std::to_string(ev), false);
        canvas_event_toa->Print(out_pdf.c_str());
        canvas_event_toa->Write();
        canvas_event_toa->Close();

        delete canvas_event_adc;
        delete canvas_event_toa;
    }

    // draw dummy page to close the pdf
    TCanvas *canvas_dummy = new TCanvas("canvas_dummy", "Dummy Canvas", 800, 600);
    canvas_dummy->Print((out_pdf + "]").c_str());
    
    output_root->Close();

    delete canvas_toa_frequency;
    delete canvas_raw_toa_ns;
    delete canvas_toa_adc_max_corr;
    delete canvas_toa_code_corr;
    delete canvas_toa_adc_max_corr_example;
    delete canvas_dummy;

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}