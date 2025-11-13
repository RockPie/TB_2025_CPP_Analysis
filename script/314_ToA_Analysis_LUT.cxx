#include <fstream>
#include <limits>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "H2GCROC_TimewalkLUT.hxx"
#include "H2GCROC_ToA_Align.hxx"

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

    // * --- Read input file ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    int machine_gun_samples = 16;
    int vldb_number = 2;
    int chn_example = 299; // print this channel seprately
    const double channel_toa_min_ns = 0.0;
    const double channel_toa_max_ns = 300.0;
    const double event_toa_min_ns = 0.0;
    const double event_toa_max_ns = 300.0;
    // const double channel_toa_min_ns = -200.0;
    // const double channel_toa_max_ns = 200.0;
    // const double event_toa_min_ns = -120.0;
    // const double event_toa_max_ns = 120.0;

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

    // TFile *input_timewalk_root = new TFile(parsed["timewalk"].as<std::string>().c_str(), "READ");
    // if (input_timewalk_root->IsZombie()) {
    //     spdlog::error("Failed to open input timewalk file {}", parsed["timewalk"].as<std::string>());
    //     return 1;
    // }

    // input_timewalk_root->cd("Optimal_ToA_ADC_Max_Correction_Parameters");
    // param_dir = gDirectory;
    // std::map<int, double> timewalk_param0_map;
    // std::map<int, double> timewalk_param1_map;
    // std::map<int, double> timewalk_param2_map;
    // std::map<int, double> timewalk_param3_map;
    // for (int chn=0; chn<FPGA_CHANNEL_NUMBER * vldb_number; chn++) {
    //     std::string param_name0 = "toa_adc_max_corr_parm0_ch" + std::to_string(chn);
    //     std::string param_name1 = "toa_adc_max_corr_parm1_ch" + std::to_string(chn);
    //     std::string param_name2 = "toa_adc_max_corr_parm2_ch" + std::to_string(chn);
    //     std::string param_name3 = "toa_adc_max_corr_parm3_ch" + std::to_string(chn);
    //     double parm0 = getParam(param_dir, param_name0.c_str());
    //     double parm1 = getParam(param_dir, param_name1.c_str());
    //     double parm2 = getParam(param_dir, param_name2.c_str());
    //     double parm3 = getParam(param_dir, param_name3.c_str());
    //     timewalk_param0_map[chn] = parm0;
    //     timewalk_param1_map[chn] = parm1;
    //     timewalk_param2_map[chn] = parm2;
    //     timewalk_param3_map[chn] = parm3;
    //     spdlog::info("Channel {}: ToA ADC Max Correction Parameters = {}, {}, {}, {}",
    //                  chn, parm0, parm1, parm2, parm3);
    // }
    // input_timewalk_root->Close();

    TFile *input_timewalk_root = new TFile(parsed["timewalk"].as<std::string>().c_str(), "READ");
    if (input_timewalk_root->IsZombie()) {
        spdlog::error("Failed to open input timewalk file {}", parsed["timewalk"].as<std::string>());
        return 1;
    }

    input_timewalk_root->cd("Timewalk_LUTs");
    std::map<int, TGraph*> timewalk_lut_map;
    for (int chn=0; chn<FPGA_CHANNEL_NUMBER * vldb_number; chn++) {
        std::string lut_name = "timewalk_lut_ch" + std::to_string(chn);
        TGraph* lut = nullptr;
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
        lut = gDirectory->Get<TGraph>(lut_name.c_str());
    #else
        gDirectory->GetObject(lut_name.c_str(), lut);
    #endif
        if (!lut) {
            spdlog::error("Missing TGraph key '{}' in directory '{}'",
                        lut_name, gDirectory->GetPath());
            return 1;
        }
        timewalk_lut_map[chn] = lut;
        spdlog::info("Channel {}: Timewalk LUT with {} points loaded",
                     chn, lut->GetN());
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

    size_t num_single_ToA = 0;
    size_t num_multi_ToA = 0;

    TH1I *h1i_toa_frequency = new TH1I("h1i_toa_frequency", "h1i_toa_frequency;ToA (25 ns bins);Counts", machine_gun_samples, 0, machine_gun_samples);
    h1i_toa_frequency->SetDirectory(nullptr);

    // * --- Raw ToA ns histograms ----------------------------------------------------------
    std::vector<TH1D*> h1d_raw_toa_ns_list;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h1d_raw_toa_ns_ch" + std::to_string(i);
        TH1D *h1d_raw_toa_ns = new TH1D(hist_name.c_str(), (hist_name + ";ToA (ns);Counts").c_str(), 256, channel_toa_min_ns, channel_toa_max_ns);
        h1d_raw_toa_ns->SetDirectory(nullptr);
        h1d_raw_toa_ns_list.push_back(h1d_raw_toa_ns);
    }

    // * --- ToA ns v.s. ADC max correction histograms --------------------------------------
    std::vector<TH2D*> h2d_toa_adc_max_corr_list;;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h2d_toa_adc_max_corr_ch" + std::to_string(i);
        TH2D *h2d_toa_adc_max_corr = new TH2D(hist_name.c_str(), (hist_name + ";Corrected ToA (ns);ADC Peak Value (pedestal subtracted)").c_str(),
                                             256, 0, 1024,
                                             256, channel_toa_min_ns, channel_toa_max_ns);
        h2d_toa_adc_max_corr->SetDirectory(nullptr);
        h2d_toa_adc_max_corr_list.push_back(h2d_toa_adc_max_corr);
    }

    // * --- ToA ns v.s. ToA code correction histograms --------------------------------------
    std::vector<TH2D*> h2d_toa_code_corr_list;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h2d_toa_code_corr_ch" + std::to_string(i);
        TH2D *h2d_toa_code_corr = new TH2D(hist_name.c_str(), (hist_name + ";Corrected ToA (ns);ToA Code").c_str(),
                                           256, 0, 1024,
                                           256, channel_toa_min_ns, channel_toa_max_ns);
        h2d_toa_code_corr->SetDirectory(nullptr);
        h2d_toa_code_corr_list.push_back(h2d_toa_code_corr);
    }

    // * --- Channel-wise waveform -------------------------------------------------
    std::vector<TH2D*> h2d_channel_waveform_list;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h2d_channel_waveform_ch" + std::to_string(i);
        TH2D *h2d_channel_waveform = new TH2D(hist_name.c_str(), (hist_name + ";Time Sample;ADC Value").c_str(),
                                              machine_gun_samples*30, -5*machine_gun_samples, machine_gun_samples*25,
                                              256, 0, 1024);
        h2d_channel_waveform->SetDirectory(nullptr);
        h2d_channel_waveform_list.push_back(h2d_channel_waveform);
    }

    // * --- Channel-seed channel ToA correlation histogram ------------------------------
    std::vector<TH2D*> h2d_channel_seed_corr_list;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h2d_channel_seed_corr_ch" + std::to_string(i);
        TH2D *h2d_channel_seed_corr = new TH2D(hist_name.c_str(), (hist_name + ";Channel ToA (ns);Seed Channel ToA (ns)").c_str(),
                                               256, channel_toa_min_ns, channel_toa_max_ns,
                                               256, channel_toa_min_ns, channel_toa_max_ns);
        h2d_channel_seed_corr->SetDirectory(nullptr);
        // fill one dummy entry to avoid empty histogram issue
        h2d_channel_seed_corr->Fill(channel_toa_min_ns, channel_toa_min_ns);
        h2d_channel_seed_corr_list.push_back(h2d_channel_seed_corr);
    }

    const int event_display_number = 10; // show the first 10 events for example

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

    // TH1D *h1d_event_average_toa = new TH1D("h1d_event_average_toa", "h1d_event_average_toa;Average ToA per Event (ns);Counts", 256, 0, 25.0*static_cast<double>(machine_gun_samples));
    // h1d_event_average_toa->SetDirectory(nullptr);

    const int adc_peak_min_index = 6;
    const int adc_peak_max_index = 7;

    std::vector<double> event_adc_list;
    std::vector<double> event_adc_toa_on_list;
    std::vector<double> event_leading_toa_list;
    std::vector<double> event_adc_weighted_toa_list;
    event_adc_list.reserve(entry_max);
    event_adc_toa_on_list.reserve(entry_max);
    event_leading_toa_list.reserve(entry_max);
    event_adc_weighted_toa_list.reserve(entry_max);

    auto dup_list = map_channels(
        vldb_number, FPGA_CHANNEL_NUMBER,
        NX, NY, board_cols, board_rows,
        sipm_board, board_loc, board_rotation, board_flip,
        ch2pid, pid2ch
    );

    const int seed_channel = 219; // use seed channel to do two-channel toa correlation
    int seed_channel_vldb_id = seed_channel / FPGA_CHANNEL_NUMBER;
    int seed_channel_local = seed_channel % FPGA_CHANNEL_NUMBER;
    int seed_channel_pid = ch2pid[seed_channel];
    if (seed_channel_pid == -1) {
        spdlog::error("Seed channel {} has no mapped pixel ID!", seed_channel);
        return 1;
    }
    int seed_channel_col = seed_channel_pid % NX;
    int seed_channel_row = seed_channel_pid / NX;
    spdlog::info("Seed channel {} mapped to pixel ID {} at (col={}, row={})", seed_channel, seed_channel_pid, seed_channel_col, seed_channel_row);

    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        double adc_sum = 0.0;
        double adc_toa_on_sum = 0.0;
        double adc_weighted_toa_sum = 0.0;
        double adc_weight_sum = 0.0;
        std::vector<double> toa_list_this_event;
        std::vector<double> adc_list_this_event;
        toa_list_this_event.reserve(FPGA_CHANNEL_NUMBER * vldb_number);
        adc_list_this_event.reserve(FPGA_CHANNEL_NUMBER * vldb_number);

        // get the toa of the seed channel
        double seed_channel_toa_ns = -1.0;
        UInt_t seed_channel_toa_code = 0;
        UInt_t seed_channel_adc_peak_value = 0;
        std::vector<UInt_t> seed_channel_pedestal_samples; // only take the first 3 samples
        seed_channel_pedestal_samples.reserve(3);
        int seed_channel_toa_sample_index = -1;
        for (int sample = 0; sample < machine_gun_samples; sample++) {
            int idx = sample*FPGA_CHANNEL_NUMBER + seed_channel_local;
            UInt_t toa_value = val2_list_pools[seed_channel_vldb_id][0][idx];
            UInt_t adc_value = val0_list_pools[seed_channel_vldb_id][0][idx];
            if (sample < 3) {
                seed_channel_pedestal_samples.push_back(adc_value);
            }
            if (sample >= adc_peak_min_index && sample <= adc_peak_max_index) {
                if (adc_value > seed_channel_adc_peak_value) {
                    seed_channel_adc_peak_value = adc_value;
                }
            }
            if (toa_value > 0 && seed_channel_toa_sample_index == -1) {
                seed_channel_toa_sample_index = sample;
                seed_channel_toa_code = toa_value;
            }
        }
        double seed_channel_pedestal = pedestal_median_of_first3(seed_channel_pedestal_samples);
        double seed_channel_adc_peak_value_pede_sub = static_cast<double>(seed_channel_adc_peak_value) - seed_channel_pedestal;
        if (seed_channel_toa_sample_index != -1) {
            // * --- DNL correction ---
            auto& coarse_dnl_lookup = dnl_coarse_tdc_list[seed_channel];
            auto& fine_dnl_lookup = dnl_fine_tdc_list[seed_channel];
            double toa_ns = decode_toa_value_ns_with_dnl(
                seed_channel_toa_code,
                coarse_dnl_lookup,
                fine_dnl_lookup
            );
            toa_ns += 25.0 * static_cast<double>(seed_channel_toa_sample_index);
            // * --- bx slip correction ---
            double channel_optimal_bx_slip = toa_correction_map[seed_channel];
            if (seed_channel_toa_code >= static_cast<UInt_t>(channel_optimal_bx_slip)) {
                toa_ns -= 25.0;
            }
            // * --- timewalk correction ---
            auto* timewalk_lut = timewalk_lut_map[seed_channel];
            double timewalk_correction_ns = timewalk_lut->Eval(seed_channel_adc_peak_value_pede_sub);
            double toa_before_timewalk = toa_ns;
            toa_ns -= timewalk_correction_ns;
            seed_channel_toa_ns = toa_ns;
        } else {
            // continue;
            // spdlog::warn("Event {}: Seed channel {} has no ToA!", entry, seed_channel);
        }
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
                std::vector<UInt_t> adc_samples;
                adc_samples.reserve(machine_gun_samples);
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
                    adc_samples.push_back(adc_value);
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
                    // * --- DNL correction ---
                    auto& coarse_dnl_lookup = dnl_coarse_tdc_list[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    auto& fine_dnl_lookup = dnl_fine_tdc_list[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    double toa_ns = decode_toa_value_ns_with_dnl(
                        toa_first_value,
                        coarse_dnl_lookup,
                        fine_dnl_lookup
                    );
                    toa_ns += 25.0 * static_cast<double>(toa_first_index);
                    // * --- bx slip correction ---
                    double channel_optimal_bx_slip = toa_correction_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    if (toa_first_value >= static_cast<UInt_t>(channel_optimal_bx_slip)) {
                        toa_ns -= 25.0;
                    }
                    // * --- timewalk correction ---
                    // double param0 = timewalk_param0_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    // double param1 = timewalk_param1_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    // double param2 = timewalk_param2_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    // double param3 = timewalk_param3_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    // double timewalk_correction_ns = param1 / std::pow(adc_peak_value_pede_sub - param2, param3);
                    // toa_ns -= timewalk_correction_ns;
                    auto* timewalk_lut = timewalk_lut_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    double timewalk_correction_ns = timewalk_lut->Eval(adc_peak_value_pede_sub);
                    double timewalk_correction_ns_max = timewalk_lut->Eval(1023.0);
                    double toa_before_timewalk = toa_ns;
                    toa_ns -= timewalk_correction_ns;
                    // * --- channel offset correction ---
                    auto valid_aligner_channel_index = from_total_channel_to_valid_channel(channel + vldb_id * FPGA_CHANNEL_NUMBER);
                    if (valid_aligner_channel_index >= 0){
                        double channel_toa_offset = toa_offset_map[valid_aligner_channel_index];
                        toa_ns -= channel_toa_offset;
                    }

                    h1d_raw_toa_ns_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_ns);
                    h2d_toa_adc_max_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(adc_peak_value_pede_sub, toa_ns);
                    h2d_toa_code_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_first_value, toa_ns);

                    if (toa_first_value > 0) {
                        toa_list_this_event.push_back(toa_ns);
                        adc_list_this_event.push_back(adc_peak_value_pede_sub);
                        adc_weighted_toa_sum += toa_ns * adc_peak_value_pede_sub;
                        adc_weight_sum += adc_peak_value_pede_sub;
                        adc_toa_on_sum += adc_peak_value_pede_sub;

                        for (int _sample_machinegun = 0; _sample_machinegun < machine_gun_samples; _sample_machinegun++) {
                            double _sample_time_ns = static_cast<double>(_sample_machinegun) * 25.0 - toa_ns;
                            double _adc_value = static_cast<double>(adc_samples[_sample_machinegun]);
                            h2d_channel_waveform_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(_sample_time_ns, _adc_value);
                        }

                        h2d_channel_seed_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(
                            seed_channel_toa_ns, toa_ns
                        );

                    }
                    if (entry < event_display_number) {
                        event_display_adc_hitmap[entry]->Fill(uni_col + 0.5, uni_row + 0.5, adc_peak_value_pede_sub);
                        if (toa_first_value > 0 && toa_ns < 300.0)
                            event_display_toa_hitmap[entry]->Fill(uni_col + 0.5, uni_row + 0.5, toa_ns);
                    }
                }
            } // end of channel loop
        } // end of vldb loop
        double leading_toa = find_leading_toa(toa_list_this_event, adc_list_this_event, 4, 1);
        double adc_weighted_event_toa = adc_weight_sum > 0.0 ? (adc_weighted_toa_sum / adc_weight_sum) : NAN;
        event_adc_weighted_toa_list.push_back(adc_weighted_event_toa);
        event_adc_toa_on_list.push_back(adc_toa_on_sum);
        event_adc_list.push_back(adc_sum);
        event_leading_toa_list.push_back(leading_toa);
        // if (adc_weight_sum > 0.0) {
        //     adc_weighted_event_toa = adc_weighted_toa_sum / adc_weight_sum;
        //     // h1d_event_average_toa->Fill(adc_weighted_event_toa);
        //     event_adc_weighted_toa_list.push_back(adc_weighted_event_toa);
        // }
        // if (leading_toa != NAN) {
        //     // h2d_leading_toa_vs_adc->Fill(adc_sum, leading_toa);
        //     event_adc_list.push_back(adc_sum);
        //     event_leading_toa_list.push_back(leading_toa);
        // }
    } // end of event loop

    input_root->Close();

    spdlog::info("Number of single ToA hits: {}", num_single_ToA);
    spdlog::info("Number of multiple ToA hits: {}", num_multi_ToA);
    spdlog::info("Percentage of single ToA hits: {:.2f}%", 100.0 * static_cast<double>(num_single_ToA) / static_cast<double>(num_single_ToA + num_multi_ToA));  

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

    // skip if single ToA histogram is empty
    if (h1i_toa_frequency->GetEntries() < 1) {
        spdlog::warn("No ToA hits found, skipping output histograms.");
        output_root->Close();
        return 0;
    }

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

    TCanvas *canvas_waveform = new TCanvas("canvas_waveform", "Channel Waveform", 1200, 800);
    draw_mosaic_fixed(*canvas_waveform, h2d_channel_waveform_list, topo_wave);
    canvas_waveform->Modified();
    canvas_waveform->Update();
    canvas_waveform->Print(out_pdf.c_str());
    canvas_waveform->Write();
    canvas_waveform->Close();

    // only draw if the total entry > 100 to avoid empty histograms
    if (num_single_ToA > 100) {
        TCanvas *canvas_channel_seed_corr = new TCanvas("canvas_channel_seed_corr", "Channel-Seed ToA Correlation", 1200, 800);
        draw_mosaic_fixed(*canvas_channel_seed_corr, h2d_channel_seed_corr_list, topo_wave);
        canvas_channel_seed_corr->Modified();
        canvas_channel_seed_corr->Update();
        canvas_channel_seed_corr->Print(out_pdf.c_str());
        canvas_channel_seed_corr->Write();
        canvas_channel_seed_corr->Close();
    }

    TCanvas *canvas_leading_toa_vs_adc = new TCanvas("canvas_leading_toa_vs_adc", "Leading ToA vs ADC Sum", 800, 600);
    canvas_leading_toa_vs_adc->cd();
    auto sorted_event_adc_list = event_adc_list;
    std::sort(sorted_event_adc_list.begin(), sorted_event_adc_list.end());
    double adc_90_percentile = sorted_event_adc_list[static_cast<size_t>(0.9 * static_cast<double>(sorted_event_adc_list.size()))];
    TH2D *h2d_leading_toa_vs_adc = new TH2D("h2d_leading_toa_vs_adc", "h2d_leading_toa_vs_adc;ADC Peak Value;Leading ToA [ns]", 256, 0, adc_90_percentile*1.3, 250, event_toa_min_ns, event_toa_max_ns);
    for (size_t i = 0; i < event_adc_list.size(); i++) {
        h2d_leading_toa_vs_adc->Fill(event_adc_list[i], event_leading_toa_list[i]);
    }
    format_2d_hist_canvas(canvas_leading_toa_vs_adc, h2d_leading_toa_vs_adc, kBlue+2, annotation_canvas_title, annotation_testbeam_title, "Leading ToA vs ADC Sum");
    canvas_leading_toa_vs_adc->Print(out_pdf.c_str());
    canvas_leading_toa_vs_adc->Write();
    canvas_leading_toa_vs_adc->Close();

    TCanvas *canvas_event_toa = new TCanvas("canvas_event_toa", "Event ToA vs ADC Weighted ToA", 800, 600);
    canvas_event_toa->cd();
    TH1D *h1d_event_adc_weighted_toa = new TH1D("h1d_event_adc_weighted_toa", "h1d_event_adc_weighted_toa;ADC Weighted ToA per Event (ns);Counts", 250, event_toa_min_ns, event_toa_max_ns);
    TH1D *h1d_event_leading_toa = new TH1D("h1d_event_leading_toa", "h1d_event_leading_toa;Leading ToA per Event (ns);Counts", 250, event_toa_min_ns, event_toa_max_ns);
    for (size_t i = 0; i < event_adc_weighted_toa_list.size(); i++) {
        h1d_event_adc_weighted_toa->Fill(event_adc_weighted_toa_list[i]);
    }
    for (size_t i = 0; i < event_leading_toa_list.size(); i++) {
        h1d_event_leading_toa->Fill(event_leading_toa_list[i]);
    }
    double max_y_adc_weighted = h1d_event_adc_weighted_toa->GetMaximum();
    double max_y_leading = h1d_event_leading_toa->GetMaximum();
    double max_y = std::max(max_y_adc_weighted, max_y_leading);
    h1d_event_adc_weighted_toa->SetMaximum(max_y * 1.3);
    format_1d_hist_canvas(canvas_event_toa, h1d_event_adc_weighted_toa, kBlue, annotation_canvas_title, annotation_testbeam_title, "Event ADC Weighted ToA");
    double adc_weighted_pre_fit_min = h1d_event_adc_weighted_toa->GetMean() - 2.0 * h1d_event_adc_weighted_toa->GetRMS();
    double adc_weighted_pre_fit_max = h1d_event_adc_weighted_toa->GetMean() + 2.0 * h1d_event_adc_weighted_toa->GetRMS();
    TF1 *pre_fit_adc_weighted = new TF1("pre_fit_adc_weighted", "gaus", adc_weighted_pre_fit_min, adc_weighted_pre_fit_max);
    h1d_event_adc_weighted_toa->Fit(pre_fit_adc_weighted, "RQ");
    double adc_weighted_fit_mean = pre_fit_adc_weighted->GetParameter(1);
    double adc_weighted_fit_sigma = pre_fit_adc_weighted->GetParameter(2);
    TF1 *fit_adc_weighted = new TF1("fit_adc_weighted", "gaus", adc_weighted_fit_mean - 2.0 * adc_weighted_fit_sigma, adc_weighted_fit_mean + 2.0 * adc_weighted_fit_sigma);
    h1d_event_adc_weighted_toa->Fit(fit_adc_weighted, "RQ");
    delete pre_fit_adc_weighted;
    fit_adc_weighted->SetLineColor(kBlue);
    fit_adc_weighted->Draw("SAME");

    h1d_event_leading_toa->SetLineColor(kRed);
    h1d_event_leading_toa->Draw("SAME");
    double leading_pre_fit_min = h1d_event_leading_toa->GetMean() - 2.0 * h1d_event_leading_toa->GetRMS();
    double leading_pre_fit_max = h1d_event_leading_toa->GetMean() + 2.0 * h1d_event_leading_toa->GetRMS();
    TF1 *pre_fit_leading = new TF1("pre_fit_leading", "gaus", leading_pre_fit_min, leading_pre_fit_max);
    h1d_event_leading_toa->Fit(pre_fit_leading, "RQ");
    double leading_fit_mean = pre_fit_leading->GetParameter(1);
    double leading_fit_sigma = pre_fit_leading->GetParameter(2);
    TF1 *fit_leading = new TF1("fit_leading", "gaus", leading_fit_mean - 2.0 * leading_fit_sigma, leading_fit_mean + 2.0 * leading_fit_sigma);
    h1d_event_leading_toa->Fit(fit_leading, "RQ");
    delete pre_fit_leading;
    fit_leading->SetLineColor(kRed);
    fit_leading->Draw("SAME");
    TLegend *legend_event_toa = new TLegend(0.6, 0.6, 0.88, 0.88);
    legend_event_toa->SetBorderSize(0);
    legend_event_toa->SetFillStyle(0);
    legend_event_toa->AddEntry(h1d_event_adc_weighted_toa, "ADC Weighted ToA", "l");
    legend_event_toa->AddEntry(fit_adc_weighted, Form("Fit: #mu=%.2f ns, #sigma=%.2f ns", adc_weighted_fit_mean, adc_weighted_fit_sigma), "l");
    legend_event_toa->AddEntry(h1d_event_leading_toa, "Leading ToA", "l");
    legend_event_toa->AddEntry(fit_leading, Form("Fit: #mu=%.2f ns, #sigma=%.2f ns", leading_fit_mean, leading_fit_sigma), "l");
    legend_event_toa->Draw();
    canvas_event_toa->Print(out_pdf.c_str());
    canvas_event_toa->Write();
    canvas_event_toa->Close();

    // ! Compare different adc sum
    std::vector<double> fit_range_sigmas = {2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5}; // in sigma
    std::vector<double> fit_range_offsets = {-0.4, -0.2, 0.0, 0.2, 0.4}; // in sigma

    TCanvas *canvas_adc_sum_comparison = new TCanvas("canvas_adc_sum_comparison", "ADC Sum Comparison", 800, 600);
    canvas_adc_sum_comparison->cd();
    auto sorted_event_adc_toa_on_list = event_adc_toa_on_list;
    std::sort(sorted_event_adc_toa_on_list.begin(), sorted_event_adc_toa_on_list.end());
    double adc_toa_on_90_percentile = sorted_event_adc_toa_on_list[static_cast<size_t>(0.9 * static_cast<double>(sorted_event_adc_toa_on_list.size()))];
    auto sorted_event_raw_adc_list = event_adc_list;
    std::sort(sorted_event_raw_adc_list.begin(), sorted_event_raw_adc_list.end());
    double raw_adc_90_percentile = sorted_event_raw_adc_list[static_cast<size_t>(0.9 * static_cast<double>(sorted_event_raw_adc_list.size()))];
    double x_max = std::max(adc_toa_on_90_percentile, raw_adc_90_percentile) * 1.8;
    TH1D *h1d_event_raw_adc = new TH1D("h1d_event_raw_adc", "h1d_event_raw_adc;ADC Sum;Counts", 256, 0, x_max);
    TH1D *h1d_event_adc_toa_on = new TH1D("h1d_event_adc_toa_on", "h1d_event_adc_toa_on;ADC Sum (ToA on);Counts", 256, 0, x_max);
    TH1D *h1d_event_adc_25ns_central_window = new TH1D("h1d_event_adc_25ns_central_window", "h1d_event_adc_25ns_central_window;ADC Sum (25 ns central window);Counts", 256, 0, x_max);
    TH1D *h1d_event_adc_10ns_central_window = new TH1D("h1d_event_adc_10ns_central_window", "h1d_event_adc_10ns_central_window;ADC Sum (10 ns central window);Counts", 256, 0, x_max);
    for (size_t i = 0; i < event_adc_list.size(); i++) {
        h1d_event_raw_adc->Fill(event_adc_list[i]);
        if (event_adc_toa_on_list[i] > 0.0)
            h1d_event_adc_toa_on->Fill(event_adc_toa_on_list[i]);
        if (event_leading_toa_list[i] >= leading_fit_mean - 12.5 && event_leading_toa_list[i] <= leading_fit_mean + 12.5)
            h1d_event_adc_25ns_central_window->Fill(event_adc_list[i]);
        if (event_leading_toa_list[i] >= leading_fit_mean - 5.0 && event_leading_toa_list[i] <= leading_fit_mean + 5.0)
            h1d_event_adc_10ns_central_window->Fill(event_adc_list[i]);
    }
    // normalize to number of events
    h1d_event_raw_adc->Scale(1.0 / static_cast<double>(h1d_event_raw_adc->GetEntries()));
    h1d_event_adc_toa_on->Scale(1.0 / static_cast<double>(h1d_event_adc_toa_on->GetEntries()));
    h1d_event_adc_25ns_central_window->Scale(1.0 / static_cast<double>(h1d_event_adc_25ns_central_window->GetEntries()));
    h1d_event_adc_10ns_central_window->Scale(1.0 / static_cast<double>(h1d_event_adc_10ns_central_window->GetEntries()));
    double max_y_raw_adc = h1d_event_raw_adc->GetMaximum();
    double max_y_adc_toa_on = h1d_event_adc_toa_on->GetMaximum();
    double max_y_adc = std::max(max_y_raw_adc, max_y_adc_toa_on);
    h1d_event_raw_adc->SetMaximum(max_y_adc * 1.3);
    format_1d_hist_canvas(canvas_adc_sum_comparison, h1d_event_raw_adc, kBlue, annotation_canvas_title, annotation_testbeam_title, "Event Raw ADC Sum");
    h1d_event_adc_toa_on->SetLineColor(kRed);
    h1d_event_adc_toa_on->Draw("SAME");
    h1d_event_adc_25ns_central_window->SetLineColor(kGreen);
    h1d_event_adc_25ns_central_window->Draw("SAME");
    h1d_event_adc_10ns_central_window->SetLineColor(kMagenta);
    h1d_event_adc_10ns_central_window->Draw("SAME");

    double raw_adc_fit_mean = 0.0;
    double raw_adc_fit_mean_err_stat = 0.0;
    double raw_adc_fit_mean_err_syst = 0.0;
    double raw_adc_fit_sigma = 0.0;
    double raw_adc_fit_sigma_err_stat = 0.0;
    double raw_adc_fit_sigma_err_syst = 0.0;
    double raw_adc_fit_resolution = 0.0;
    double raw_adc_fit_resolution_err = 0.0;
    crystalball_fit_th1d(*canvas_adc_sum_comparison, *h1d_event_raw_adc, fit_range_sigmas, fit_range_offsets, kBlue, raw_adc_fit_mean, raw_adc_fit_mean_err_stat, raw_adc_fit_mean_err_syst, raw_adc_fit_sigma, raw_adc_fit_sigma_err_stat, raw_adc_fit_sigma_err_syst, raw_adc_fit_resolution, raw_adc_fit_resolution_err);

    double adc_toa_on_fit_mean = 0.0;
    double adc_toa_on_fit_mean_err_stat = 0.0;
    double adc_toa_on_fit_mean_err_syst = 0.0;
    double adc_toa_on_fit_sigma = 0.0;
    double adc_toa_on_fit_sigma_err_stat = 0.0;
    double adc_toa_on_fit_sigma_err_syst = 0.0;
    double adc_toa_on_fit_resolution = 0.0;
    double adc_toa_on_fit_resolution_err = 0.0;
    crystalball_fit_th1d(*canvas_adc_sum_comparison, *h1d_event_adc_toa_on, fit_range_sigmas, fit_range_offsets, kRed, adc_toa_on_fit_mean, adc_toa_on_fit_mean_err_stat, adc_toa_on_fit_mean_err_syst, adc_toa_on_fit_sigma, adc_toa_on_fit_sigma_err_stat, adc_toa_on_fit_sigma_err_syst, adc_toa_on_fit_resolution, adc_toa_on_fit_resolution_err);

    double adc_25ns_fit_mean = 0.0;
    double adc_25ns_fit_mean_err_stat = 0.0;
    double adc_25ns_fit_mean_err_syst = 0.0;
    double adc_25ns_fit_sigma = 0.0;
    double adc_25ns_fit_sigma_err_stat = 0.0;
    double adc_25ns_fit_sigma_err_syst = 0.0;
    double adc_25ns_fit_resolution = 0.0;
    double adc_25ns_fit_resolution_err = 0.0;
    crystalball_fit_th1d(*canvas_adc_sum_comparison, *h1d_event_adc_25ns_central_window, fit_range_sigmas, fit_range_offsets, kGreen, adc_25ns_fit_mean, adc_25ns_fit_mean_err_stat, adc_25ns_fit_mean_err_syst, adc_25ns_fit_sigma, adc_25ns_fit_sigma_err_stat, adc_25ns_fit_sigma_err_syst, adc_25ns_fit_resolution, adc_25ns_fit_resolution_err);

    double adc_10ns_fit_mean = 0.0;
    double adc_10ns_fit_mean_err_stat = 0.0;
    double adc_10ns_fit_mean_err_syst = 0.0;
    double adc_10ns_fit_sigma = 0.0;
    double adc_10ns_fit_sigma_err_stat = 0.0;
    double adc_10ns_fit_sigma_err_syst = 0.0;
    double adc_10ns_fit_resolution = 0.0;
    double adc_10ns_fit_resolution_err = 0.0;
    crystalball_fit_th1d(*canvas_adc_sum_comparison, *h1d_event_adc_10ns_central_window, fit_range_sigmas, fit_range_offsets, kMagenta, adc_10ns_fit_mean, adc_10ns_fit_mean_err_stat, adc_10ns_fit_mean_err_syst, adc_10ns_fit_sigma, adc_10ns_fit_sigma_err_stat, adc_10ns_fit_sigma_err_syst, adc_10ns_fit_resolution, adc_10ns_fit_resolution_err);

    TLegend *legend_adc_sum = new TLegend(0.6, 0.3, 0.89, 0.89);
    legend_adc_sum->SetBorderSize(0);
    legend_adc_sum->SetFillStyle(0);
    legend_adc_sum->AddEntry(h1d_event_raw_adc, Form("ADC Sum (Raw) (N=%.0f)", h1d_event_raw_adc->GetEntries()), "l");
    legend_adc_sum->AddEntry(h1d_event_raw_adc, Form("Resolution: %.2f%% #pm %.2f%%", raw_adc_fit_resolution, raw_adc_fit_resolution_err), "l");
    legend_adc_sum->AddEntry(h1d_event_raw_adc, Form("Mean: %.2f #pm %.2f (stat) #pm %.2f (syst)", raw_adc_fit_mean, raw_adc_fit_mean_err_stat, raw_adc_fit_mean_err_syst), "l");
    legend_adc_sum->AddEntry(h1d_event_raw_adc, Form("Sigma: %.2f #pm %.2f (stat) #pm %.2f (syst)", raw_adc_fit_sigma, raw_adc_fit_sigma_err_stat, raw_adc_fit_sigma_err_syst), "l");

    legend_adc_sum->AddEntry(h1d_event_adc_toa_on, Form("ADC Sum (ToA on) (N=%.0f)", h1d_event_adc_toa_on->GetEntries()), "l");
    legend_adc_sum->AddEntry(h1d_event_adc_toa_on, Form("Resolution: %.2f%% #pm %.2f%%", adc_toa_on_fit_resolution, adc_toa_on_fit_resolution_err), "l");
    legend_adc_sum->AddEntry(h1d_event_adc_toa_on, Form("Mean: %.2f #pm %.2f (stat) #pm %.2f (syst)", adc_toa_on_fit_mean, adc_toa_on_fit_mean_err_stat, adc_toa_on_fit_mean_err_syst), "l");
    legend_adc_sum->AddEntry(h1d_event_adc_toa_on, Form("Sigma: %.2f #pm %.2f (stat) #pm %.2f (syst)", adc_toa_on_fit_sigma, adc_toa_on_fit_sigma_err_stat, adc_toa_on_fit_sigma_err_syst), "l");


    legend_adc_sum->AddEntry(h1d_event_adc_25ns_central_window, Form("ADC Sum (25 ns window) (N=%.0f)", h1d_event_adc_25ns_central_window->GetEntries()), "l");
    legend_adc_sum->AddEntry(h1d_event_adc_25ns_central_window, Form("Resolution: %.2f%% #pm %.2f%%", adc_25ns_fit_resolution, adc_25ns_fit_resolution_err), "l");
    legend_adc_sum->AddEntry(h1d_event_adc_25ns_central_window, Form("Mean: %.2f #pm %.2f (stat) #pm %.2f (syst)", adc_25ns_fit_mean, adc_25ns_fit_mean_err_stat, adc_25ns_fit_mean_err_syst), "l");
    legend_adc_sum->AddEntry(h1d_event_adc_25ns_central_window, Form("Sigma: %.2f #pm %.2f (stat) #pm %.2f (syst)", adc_25ns_fit_sigma, adc_25ns_fit_sigma_err_stat, adc_25ns_fit_sigma_err_syst), "l");

    legend_adc_sum->AddEntry(h1d_event_adc_10ns_central_window, Form("ADC Sum (10 ns window) (N=%.0f)", h1d_event_adc_10ns_central_window->GetEntries()), "l");
    legend_adc_sum->AddEntry(h1d_event_adc_10ns_central_window, Form("Resolution: %.2f%% #pm %.2f%%", adc_10ns_fit_resolution, adc_10ns_fit_resolution_err), "l");
    legend_adc_sum->AddEntry(h1d_event_adc_10ns_central_window, Form("Mean: %.2f #pm %.2f (stat) #pm %.2f (syst)", adc_10ns_fit_mean, adc_10ns_fit_mean_err_stat, adc_10ns_fit_mean_err_syst), "l");
    legend_adc_sum->AddEntry(h1d_event_adc_10ns_central_window, Form("Sigma: %.2f #pm %.2f (stat) #pm %.2f (syst)", adc_10ns_fit_sigma, adc_10ns_fit_sigma_err_stat, adc_10ns_fit_sigma_err_syst), "l");

    legend_adc_sum->Draw();
    canvas_adc_sum_comparison->Print(out_pdf.c_str());
    canvas_adc_sum_comparison->Write();
    canvas_adc_sum_comparison->Close();

    TCanvas *canvas_example_chn_waveform = new TCanvas("canvas_example_chn_waveform", "Example Channel Waveform", 800, 600);
    auto h2d_example_waveform = h2d_channel_waveform_list[chn_example];
    format_2d_hist_canvas(canvas_example_chn_waveform, h2d_example_waveform, kBlue+2, annotation_canvas_title, annotation_testbeam_title, "Channel_" + std::to_string(chn_example));
    canvas_example_chn_waveform->Print(out_pdf.c_str());
    canvas_example_chn_waveform->Write();
    canvas_example_chn_waveform->Close();

    // print the example channel
    TCanvas *canvas_toa_adc_max_corr_example = new TCanvas("canvas_toa_adc_max_corr_example", "ToA ADC Max Correction Example Channel", 800, 600);
    auto h2d_example = h2d_toa_adc_max_corr_list[chn_example];
    h2d_example->GetXaxis()->SetTitle("ADC Peak (pedestal subtracted)");
    h2d_example->GetYaxis()->SetTitle("Corrected ToA [ns]");
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