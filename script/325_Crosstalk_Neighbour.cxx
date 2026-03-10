#include <fstream>
#include <limits>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "H2GCROC_TimewalkLUT.hxx"
#include "H2GCROC_ToA_Align.hxx"
#include "CommonParams.hxx"

std::vector<int> find_seed_channel(std::vector<double> &adc_vec, std::vector<double> &toa_vec, std::vector<double> &tot_vec, int number_of_seed_channels = 1, int candidate_channel_number = 0) {
    int size = adc_vec.size();
    if (size == 0 || toa_vec.size() != size || tot_vec.size() != size) {
        spdlog::error("Input vectors must have the same non-zero size");
        return {};
    }
    if (candidate_channel_number <= 0) {
        candidate_channel_number = number_of_seed_channels + 2;
        if (candidate_channel_number > size) {
            candidate_channel_number = size;
        }
    }
    
    // Start with ALL channel indices
    std::vector<int> candidate_channels(size);
    std::iota(candidate_channels.begin(), candidate_channels.end(), 0);
    
    // Step 1: Find top candidate_channel_number channels with highest ADC values
    std::partial_sort(candidate_channels.begin(), 
                      candidate_channels.begin() + candidate_channel_number, 
                      candidate_channels.end(),
                      [&adc_vec](int a, int b) {
                          return adc_vec[a] > adc_vec[b];
                      });

    // Step 2: From those candidates, find number_of_seed_channels with minimal ToA values
    std::partial_sort(candidate_channels.begin(), 
                      candidate_channels.begin() + number_of_seed_channels, 
                      candidate_channels.begin() + candidate_channel_number,
                      [&toa_vec](int a, int b) {
                          return toa_vec[a] < toa_vec[b];
                      });
    
    // Return the seed channels
    std::vector<int> seed_channels(number_of_seed_channels);
    std::copy(candidate_channels.begin(), 
              candidate_channels.begin() + number_of_seed_channels, 
              seed_channels.begin());
    
    return seed_channels;
}

// Find seed channel by maximum ToT, fallback to maximum ADC if no ToT
std::vector<int> find_seed_channel_tot(std::vector<double> &adc_vec, std::vector<double> &tot_vec, int number_of_seed_channels = 1, int candidate_channel_number = 0) {
    int size = adc_vec.size();
    if (size == 0 || tot_vec.size() != size) {
        spdlog::error("Input vectors must have the same non-zero size");
        return {};
    }
    
    if (candidate_channel_number <= 0) {
        candidate_channel_number = number_of_seed_channels + 2;
        if (candidate_channel_number > size) {
            candidate_channel_number = size;
        }
    }
    
    // Check if there are any non-zero ToT values
    bool has_tot = false;
    for (const auto& tot : tot_vec) {
        if (tot > 0) {
            has_tot = true;
            break;
        }
    }
    
    // Initialize candidate channels with all indices
    std::vector<int> candidate_channels(size);
    std::iota(candidate_channels.begin(), candidate_channels.end(), 0);
    
    if (has_tot) {
        // Step 1: Find top candidate_channel_number channels with highest ToT values
        std::partial_sort(candidate_channels.begin(), 
                          candidate_channels.begin() + candidate_channel_number, 
                          candidate_channels.end(),
                          [&tot_vec](int a, int b) {
                              return tot_vec[a] > tot_vec[b];
                          });
    } else {
        // No ToT available, use ADC instead
        spdlog::warn("No ToT values detected; using ADC for seed channel selection");
        std::partial_sort(candidate_channels.begin(), 
                          candidate_channels.begin() + candidate_channel_number, 
                          candidate_channels.end(),
                          [&adc_vec](int a, int b) {
                              return adc_vec[a] > adc_vec[b];
                          });
    }
    
    // Return top number_of_seed_channels from candidates
    std::vector<int> seed_channels(std::min(number_of_seed_channels, candidate_channel_number));
    std::copy(candidate_channels.begin(), 
              candidate_channels.begin() + seed_channels.size(), 
              seed_channels.begin());
    
    return seed_channels;
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
    
    cxxopts::Options options(script_name, "Generate heatmaps from machine gun data");

    options.add_options()
        ("f,file", "Input .root file", cxxopts::value<std::string>())
        ("c,correction", "Input ToA correction .root file", cxxopts::value<std::string>())
        ("d,dnl", "Input DNL .root file", cxxopts::value<std::string>())
        ("t,timewalk", "Input Timewalk correction .root file", cxxopts::value<std::string>())
        ("O,offset", "Input ToA offset .root file", cxxopts::value<std::string>())
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
    std::string script_input_run_number = script_input_file.substr(script_input_file.find_last_of("Run") + 1).append(4, '0').substr(0, 4);
    script_output_file = parsed["output"].as<std::string>();
    script_n_events    = parsed["events"].as<int>();
    script_verbose     = parsed["verbose"].as<bool>();

    configure_logger(script_verbose);

    spdlog::info("Input file: {}", script_input_file);

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

    const int CM_channel = 9;
    int machine_gun_samples = 16;
    int vldb_number = 2;
    int chn_example = CommonParams::example_channel;

    const int adc_peak_min_index = CommonParams::adc_peak_min_index;
    const int adc_peak_max_index = CommonParams::adc_peak_max_index;
    const int undershoot_min_index = 3;
    const int undershoot_max_index = 7;
    
    const int event_display_number = 10;
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

    auto chan2pad = build_chan2pad_LUT(
        vldb_number,
        FPGA_CHANNEL_NUMBER,
        NX,
        NY,
        board_cols,
        board_rows,
        sipm_board,
        board_loc,
        board_rotation,
        board_flip
    );

    if (static_cast<int>(chan2pad.size()) != TOTAL_CH) {
        spdlog::error("Channel-to-pad LUT size mismatch: expected {} entries, got {}", TOTAL_CH, chan2pad.size());
        return 1;
    }

    for (int ch = 0; ch < TOTAL_CH; ++ch) {
        int pad_linear = chan2pad[ch];
        ch2pid[ch] = pad_linear;
        if (pad_linear < 0 || pad_linear >= NPIX) {
            continue;
        }
        if (pid2ch[pad_linear] != -1 && pid2ch[pad_linear] != ch) {
            spdlog::warn("Pixel {} already assigned to channel {} when mapping channel {}", pad_linear, pid2ch[pad_linear], ch);
            continue;
        }
        pid2ch[pad_linear] = ch;
    }


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
    spdlog::info("Input bx correction file: {}", parsed["correction"].as<std::string>());
    spdlog::info("Output file: {} in {}", script_output_file, script_output_folder);
    spdlog::info("Number of events: {}", script_n_events);



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

    // Try to load timewalk LUTs, but handle gracefully if directory/keys don't exist
    std::map<int, TGraph*> timewalk_lut_map;
    TDirectory* lut_dir = input_timewalk_root->GetDirectory("Timewalk_LUTs");
    
    if (lut_dir != nullptr) {
        lut_dir->cd();
        for (int chn=0; chn<FPGA_CHANNEL_NUMBER * vldb_number; chn++) {
            std::string lut_name = "timewalk_lut_ch" + std::to_string(chn);
            TGraph* lut = nullptr;
        #if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
            lut = gDirectory->Get<TGraph>(lut_name.c_str());
        #else
            gDirectory->GetObject(lut_name.c_str(), lut);
        #endif
            if (!lut) {
                spdlog::warn("Missing TGraph key '{}' in directory '{}'; skipping this channel",
                            lut_name, gDirectory->GetPath());
                continue;
            }
            timewalk_lut_map[chn] = lut;
            spdlog::info("Channel {}: Timewalk LUT with {} points loaded",
                         chn, lut->GetN());
        }
        if (timewalk_lut_map.empty()) {
            spdlog::warn("No timewalk LUTs loaded; all channels were missing");
        }
    } else {
        spdlog::warn("Timewalk_LUTs directory not found in input timewalk file");
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
    }
    input_dnl_root->Close();

    TFile *input_offset_root = new TFile(parsed["offset"].as<std::string>().c_str(), "READ");
    if (input_offset_root->IsZombie()) {
        spdlog::error("Failed to open input offset file {}", parsed["offset"].as<std::string>());
        return 1;
    }

    input_offset_root->cd("ToA_Offsets");
    std::map<int, double> toa_offset_map;
    for (int chn=0; chn<FPGA_CHANNEL_NUMBER * vldb_number; chn++) {
        auto valid_aligner_channel_index = from_total_channel_to_valid_channel(chn);
        if (valid_aligner_channel_index < 0)
            continue;
        std::string param_name = "toa_offset_ch" + std::to_string(valid_aligner_channel_index);
        double offset_value = getParam(gDirectory, param_name.c_str());
        toa_offset_map[valid_aligner_channel_index] = offset_value;
        spdlog::info("Channel {}: ToA offset = {}", chn, offset_value);
    }
    input_offset_root->Close();

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

    const double adc_average_zoom_factor = 1.5;

    double adc_hist_min = 0.0;
    double adc_hist_max = 1024.0;
    double adc_hist_bin_width = 4.0;
    int adc_hist_bins = static_cast<int>((adc_hist_max - adc_hist_min) / adc_hist_bin_width);

    double adc_channel_average_min = 0.0;
    double adc_channel_average_max = 200.0;
    double adc_channel_average_bin_width = 0.5;
    int adc_channel_average_bins = static_cast<int>((adc_channel_average_max - adc_channel_average_min) / adc_channel_average_bin_width);

    double adc_channel_undershoot_min = 0.0;
    double adc_channel_undershoot_max = 20.0;
    double adc_channel_undershoot_bin_width = 0.05;
    int adc_channel_undershoot_bins = static_cast<int>((adc_channel_undershoot_max - adc_channel_undershoot_min) / adc_channel_undershoot_bin_width);

    std::vector<std::vector<std::vector<double>>> adc_samples_in_CM_channels_all_vldb(
        vldb_number, std::vector<std::vector<double>>(4, std::vector<double>(machine_gun_samples, 0.0)));
    
    auto th1d_seed_channel_distribution = new TH1D("th1d_seed_channel_distribution", "Seed Channel Distribution;Channel Number;Entries", FPGA_CHANNEL_NUMBER*vldb_number, 0, FPGA_CHANNEL_NUMBER*vldb_number);
    th1d_seed_channel_distribution->SetDirectory(nullptr);

    auto th2d_chn221_chn217_adc_correlation = new TH2D("th2d_chn221_chn217_adc_correlation", "ADC Correlation between Chn 221 and Chn 217;Chn 221 ADC;Chn 217 ADC", adc_hist_bins, adc_hist_min, adc_hist_max, adc_hist_bins, adc_hist_min, adc_hist_max);
    th2d_chn221_chn217_adc_correlation->SetDirectory(nullptr);

    auto th2d_chn221_chn217_undershoot_correlation = new TH2D("th2d_chn221_chn217_undershoot_correlation", "Undershoot Correlation between Chn 221 and Chn 217;Chn 221 ADC;Chn 217 Undershoot", adc_hist_bins, adc_hist_min, adc_hist_max, 128, 0, 128);
    th2d_chn221_chn217_undershoot_correlation->SetDirectory(nullptr);

    auto th2d_adc_average_vs_adc_correlation = new TH2D("th2d_adc_average_vs_adc_correlation", "Correlation between ADC Average and ADC Peak;ADC Average;ADC Peak", adc_channel_average_bins, adc_channel_average_min, adc_channel_average_max, adc_channel_average_bins, adc_channel_average_min, adc_channel_average_max);
    th2d_adc_average_vs_adc_correlation->SetDirectory(nullptr);

    auto th2d_undershoot_vs_adc_correlation = new TH2D("th2d_undershoot_vs_adc_correlation", "Correlation between Undershoot and ADC Peak;Undershoot;ADC Peak", adc_channel_undershoot_bins, adc_channel_undershoot_min, adc_channel_undershoot_max, adc_channel_average_bins, adc_channel_average_min, adc_channel_average_max);
    th2d_undershoot_vs_adc_correlation->SetDirectory(nullptr);


    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        std::vector<double> toa_list_this_event;
        std::vector<double> adc_list_this_event;
        std::vector<double> tot_list_this_event;
        std::vector<double> undershoot_adc_list_this_event;
        std::vector<double> adc_average_list_this_event;
        toa_list_this_event.reserve(FPGA_CHANNEL_NUMBER * vldb_number);
        adc_list_this_event.reserve(FPGA_CHANNEL_NUMBER * vldb_number);
        tot_list_this_event.reserve(FPGA_CHANNEL_NUMBER * vldb_number);
        undershoot_adc_list_this_event.reserve(FPGA_CHANNEL_NUMBER * vldb_number);
        adc_average_list_this_event.reserve(FPGA_CHANNEL_NUMBER * vldb_number);
        for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
            toa_list_this_event.push_back(0.0);
            adc_list_this_event.push_back(0.0);
            tot_list_this_event.push_back(0.0);
            undershoot_adc_list_this_event.push_back(0.0);
            adc_average_list_this_event.push_back(0.0);
        }

        int seed_channel_number = -1; // -1 for not found
        double seed_channel_adc = 0.0;
        double seed_channel_toa = 0.0;
        double seed_channel_tot = 0.0;
        double seed_channel_adc_average = 0.0;

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
            } // end of half CM loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int valid_channel_number = get_valid_fpga_channel(channel);
                if (valid_channel_number < 0) {
                    continue; // skip invalid channels
                }
                int channel_half_index = (channel / (FPGA_CHANNEL_NUMBER / 4)) % 4;
                if (channel_half_index == 2 || channel_half_index == 3) {
                    continue; // skip CM channels
                }
                std::vector<double> CM_samples = adc_samples_in_CM_channels[channel_half_index];
                std::vector<double> adc_samples; // all samples
                std::vector<double> adc_pedestal_samples; // only take the first 3 samples
                int adc_peak_index = 0;
                double adc_peak_value = 0.0;
                int undershoot_peak_index = 0;
                double undershoot_peak_value = 1024.0; // start with a large value since we are looking for minimum
                int toa_first_index = -1;
                double toa_peak_value = 0.0;
                int tot_first_index = -1;
                double tot_peak_value = 0.0;
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
                    if (adc_value < undershoot_peak_value && sample >= undershoot_min_index && sample <= undershoot_max_index) {
                        undershoot_peak_value = adc_value;
                        undershoot_peak_index = sample;
                    }
                    if (toa_first_index < 0 && tot_value > 0) {
                        toa_first_index = sample;
                        toa_peak_value = static_cast<double>(toa_value);
                        tot_first_index = sample;
                        tot_peak_value = static_cast<double>(tot_value);
                    }
                } // end of sample loop
                double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
                double adc_peak_pede_subtracted = adc_peak_value - adc_pedestal;
                if (adc_peak_pede_subtracted < 0) {
                    adc_peak_pede_subtracted = 0.0; // set negative values to zero
                }
                double undershoot_adc_pede_subtracted = adc_pedestal - undershoot_peak_value; // note the subtraction order for undershoot
                if (undershoot_adc_pede_subtracted < 0) {
                    undershoot_adc_pede_subtracted = 0.0; // set negative values to zero
                }
                adc_list_this_event[channel + vldb_id*FPGA_CHANNEL_NUMBER] = adc_peak_pede_subtracted;
                undershoot_adc_list_this_event[channel + vldb_id*FPGA_CHANNEL_NUMBER] = undershoot_adc_pede_subtracted;
                // calculate average ADC of all the event
                double adc_sum = 0.0;
                for (double sample_adc : adc_samples) {
                    adc_sum += sample_adc;
                }
                double adc_average = adc_sum / static_cast<double>(adc_samples.size());
                adc_average_list_this_event[channel + vldb_id*FPGA_CHANNEL_NUMBER] = adc_average;
                if (toa_first_index >= 0) {
                    // * step 1: dnl correction
                    auto& coarse_dnl_lookup = dnl_coarse_tdc_list[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    auto& fine_dnl_lookup = dnl_fine_tdc_list[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    double toa_ns = decode_toa_value_ns_with_dnl(
                        toa_peak_value,
                        coarse_dnl_lookup,
                        fine_dnl_lookup
                    );
                    // * step 2: include bx time
                    toa_ns += 25.0 * static_cast<double>(toa_first_index);
                    // * step 3: optimal bx slip correction
                    double channel_optimal_bx_slip = toa_correction_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    if (toa_peak_value >= static_cast<UInt_t>(channel_optimal_bx_slip)) {
                        toa_ns -= 25.0;
                    }
                    // * step 4: timewalk correction
                    int ch_index = channel + vldb_id * FPGA_CHANNEL_NUMBER;
                    auto it = timewalk_lut_map.find(ch_index);
                    if (it != timewalk_lut_map.end() && it->second != nullptr) {
                        auto* timewalk_lut = it->second;
                        double timewalk_correction_ns = timewalk_lut->Eval(adc_peak_pede_subtracted);
                        toa_ns -= timewalk_correction_ns;
                    }
                    // If no timewalk LUT available, skip correction (already warned during loading)
                    toa_list_this_event[channel + vldb_id*FPGA_CHANNEL_NUMBER] = toa_ns;
                } else {
                    toa_list_this_event[channel + vldb_id*FPGA_CHANNEL_NUMBER] = 0.0;
                }
                if (tot_first_index >= 0) {
                    double tot_decoded = tot_peak_value;
                    if (tot_decoded > 512.0) {
                        tot_decoded = 8.0 * (tot_decoded - 512.0);
                    }
                    tot_list_this_event[channel + vldb_id*FPGA_CHANNEL_NUMBER] = tot_decoded;
                } else {
                    tot_list_this_event[channel + vldb_id*FPGA_CHANNEL_NUMBER] = 0.0;
                }
            } // end of channel loop
        } // end of VLDB loop
        // find the seed channel
        // std::vector<int> seed_channels = find_seed_channel(adc_list_this_event, toa_list_this_event, tot_list_this_event, 1);
        std::vector<int> seed_channels_no_toa = find_seed_channel_tot(adc_list_this_event, tot_list_this_event, 1);
        if (!seed_channels_no_toa.empty()) {
            seed_channel_number = seed_channels_no_toa[0];
            seed_channel_adc = adc_list_this_event[seed_channel_number];
            seed_channel_toa = toa_list_this_event[seed_channel_number];
            seed_channel_tot = tot_list_this_event[seed_channel_number];
            seed_channel_adc_average = adc_average_list_this_event[seed_channel_number] * adc_average_zoom_factor;  
            th1d_seed_channel_distribution->Fill(seed_channel_number);
            // if (seed_channel_number == 221) {
            if (true) {
                int chn217_index = 192;
                double chn217_adc = adc_list_this_event[chn217_index];
                double chn217_adc_average = adc_average_list_this_event[chn217_index] * adc_average_zoom_factor;
                th2d_chn221_chn217_adc_correlation->Fill(seed_channel_adc_average, chn217_adc_average);
                double chn217_undershoot_adc = undershoot_adc_list_this_event[chn217_index];
                th2d_chn221_chn217_undershoot_correlation->Fill(seed_channel_adc_average, chn217_undershoot_adc);
            }
        }

        // calculate event ADC peak sum and event undershoot sum
        double event_adc_peak_sum = 0.0;
        double event_undershoot_sum = 0.0;
        double event_adc_average_sum = 0.0;
        double event_nonzero_adc_channel_count = 0.0;
        for (int ch = 0; ch < FPGA_CHANNEL_NUMBER * vldb_number; ch++) {
            if (ch < 152 || ch >= 152+38) {
                continue; // focus on the central region channels for average calculation
            }
            event_adc_peak_sum += adc_list_this_event[ch];
            event_undershoot_sum += undershoot_adc_list_this_event[ch];
            event_adc_average_sum += adc_average_list_this_event[ch];
            if (adc_list_this_event[ch] > 0.0) {
                event_nonzero_adc_channel_count += 1.0;
            }
        }
        double event_adc_peak_average = event_adc_peak_sum / (event_nonzero_adc_channel_count > 0 ? event_nonzero_adc_channel_count : 1.0);
        double event_undershoot_average = event_undershoot_sum / (event_nonzero_adc_channel_count > 0 ? event_nonzero_adc_channel_count : 1.0);
        double event_adc_average = event_adc_average_sum / (event_nonzero_adc_channel_count > 0 ? event_nonzero_adc_channel_count : 1.0);

        th2d_adc_average_vs_adc_correlation->Fill(event_adc_average, event_adc_peak_average);
        th2d_undershoot_vs_adc_correlation->Fill(event_undershoot_average, event_adc_peak_average);
    } // end of event loop

    input_root->Close();

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    const std::string out_pdf = script_output_file.substr(0, script_output_file.find_last_of(".")) + ".pdf";

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

    auto canvas_seed_channel_distribution = new TCanvas("canvas_seed_channel_distribution", "Seed Channel Distribution", 800, 600);
    format_1d_hist_canvas(canvas_seed_channel_distribution, th1d_seed_channel_distribution, kBlue, annotation_canvas_title, annotation_testbeam_title, "Seed Channel Number");
    canvas_seed_channel_distribution->Modified();
    canvas_seed_channel_distribution->Update();
    canvas_seed_channel_distribution->Write();
    canvas_seed_channel_distribution->Close();
    
    auto canvas_chn221_chn217_adc_correlation = new TCanvas("canvas_chn221_chn217_adc_correlation", "ADC Correlation between Chn 221 and Chn 217", 800, 600);
    format_2d_hist_canvas(canvas_chn221_chn217_adc_correlation, th2d_chn221_chn217_adc_correlation, kRed, annotation_canvas_title, annotation_testbeam_title, "Channel 221 v.s. 217 ADC Correlation");
    canvas_chn221_chn217_adc_correlation->Modified();
    canvas_chn221_chn217_adc_correlation->Update();
    canvas_chn221_chn217_adc_correlation->Write();
    canvas_chn221_chn217_adc_correlation->Close();

    auto canvas_chn221_chn217_undershoot_correlation = new TCanvas("canvas_chn221_chn217_undershoot_correlation", "Undershoot Correlation between Chn 221 and Chn 217", 800, 600);
    format_2d_hist_canvas(canvas_chn221_chn217_undershoot_correlation, th2d_chn221_chn217_undershoot_correlation, kRed, annotation_canvas_title, annotation_testbeam_title, "Channel 221 ADC v.s. Channel 217 Undershoot Correlation");
    canvas_chn221_chn217_undershoot_correlation->Modified();
    canvas_chn221_chn217_undershoot_correlation->Update();
    canvas_chn221_chn217_undershoot_correlation->Write();
    canvas_chn221_chn217_undershoot_correlation->Close();

    auto canvas_adc_average_vs_adc_correlation = new TCanvas("canvas_adc_average_vs_adc_correlation", "Correlation between ADC Average and ADC Peak", 800, 600);
    format_2d_hist_canvas(canvas_adc_average_vs_adc_correlation, th2d_adc_average_vs_adc_correlation, kRed, annotation_canvas_title, annotation_testbeam_title, "Event ADC Average vs. ADC Peak Correlation");
    canvas_adc_average_vs_adc_correlation->Modified();
    canvas_adc_average_vs_adc_correlation->Update();
    canvas_adc_average_vs_adc_correlation->Write();
    canvas_adc_average_vs_adc_correlation->Close();

    auto canvas_undershoot_vs_adc_correlation = new TCanvas("canvas_undershoot_vs_adc_correlation", "Correlation between Undershoot and ADC Peak", 800, 600);
    format_2d_hist_canvas(canvas_undershoot_vs_adc_correlation, th2d_undershoot_vs_adc_correlation, kRed, annotation_canvas_title, annotation_testbeam_title, "Event Undershoot vs. ADC Peak Correlation");
    canvas_undershoot_vs_adc_correlation->Modified();
    canvas_undershoot_vs_adc_correlation->Update();
    canvas_undershoot_vs_adc_correlation->Write();
    canvas_undershoot_vs_adc_correlation->Close();

    output_root->Close();

    return 0;
}