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

    AlignConfig cfg;
    cfg.NCH         = 76*3;
    cfg.DT_MIN      = -50.0;
    cfg.DT_MAX      = +50.0;
    cfg.DT_NBINS    = 800;
    cfg.MIN_ENTRIES = 80;
    cfg.Tbc_ns      = 0.0; 

    cfg.out_root    = nullptr;

    ToAAligner aligner(cfg);


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
        spdlog::info("Channel {}: Timewalk LUT with {} points loaded", chn, lut->GetN());
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

    if (!dup_list.empty()) {
        spdlog::error("Found {} duplicate locations.", dup_list.size());
        for (auto [pid,gch] : dup_list) {
            int r = pid / NX, c = pid % NX;
            spdlog::error("  pid={} → uni({},{}) claimed by gch={} and gch={}",
                        pid, c, r, pid2ch[pid], gch);
        }
    }

    const int seed_channel = 74; // use seed channel to do two-channel toa correlation
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
        std::vector<std::pair<int,double>> hits; // for toa alignment
        hits.reserve(FPGA_CHANNEL_NUMBER * vldb_number);
        std::vector<double> toa_list_this_event;
        std::vector<double> adc_list_this_event;
        toa_list_this_event.reserve(FPGA_CHANNEL_NUMBER * vldb_number);
        adc_list_this_event.reserve(FPGA_CHANNEL_NUMBER * vldb_number);

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
                // double adc_pedestal = pedestal_median_of_first3(adc_pedestal_samples);
                double adc_pedestal = pedestal_average_of_first3(adc_pedestal_samples);
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
                    auto* timewalk_lut = timewalk_lut_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    double timewalk_correction_ns = timewalk_lut->Eval(adc_peak_value_pede_sub);
                    double timewalk_correction_ns_max = timewalk_lut->Eval(1023.0);
                    double toa_before_timewalk = toa_ns;
                    toa_ns -= timewalk_correction_ns;
                    // * --- channel offset correction ---
                    // double channel_toa_offset = toa_offset_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    // toa_ns -= channel_toa_offset;

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

                        // h2d_channel_seed_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(
                        //     seed_channel_toa_ns, toa_ns
                        // );
                        auto valid_aligner_channel_index = from_total_channel_to_valid_channel(channel + vldb_id * FPGA_CHANNEL_NUMBER);
                        if (valid_aligner_channel_index != -1)
                            hits.emplace_back(valid_aligner_channel_index, toa_ns);

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
        aligner.addEvent(hits);
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

    auto offsets = aligner.solveAll(false);

    double prev_mean_abs = 1e99;
    for (int iter=0; iter<5; ++iter) {
        aligner.irlsHuberUpdate(offsets);
        auto new_offsets = aligner.solveAll(false);

        // 计算平均绝对残差
        double mean_abs = 0; int cnt=0;
        for (const auto& p : aligner.pairs) if (p.w>0) {
            double r = (new_offsets[p.i]-new_offsets[p.j]) - p.mu;
            mean_abs += std::abs(r); ++cnt;
        }
        mean_abs /= std::max(1, cnt);
        std::cout << "IRLS iter " << iter+1 << " mean|res|=" << mean_abs << " ns\n";

        if (std::abs(prev_mean_abs - mean_abs) < 0.01) { // 0.01ns 阈值
            offsets = std::move(new_offsets);
            break;
        }
        prev_mean_abs = mean_abs;
        offsets = std::move(new_offsets);
    }

    std::cout << "First 16 recovered offsets (ns):\n";
    for (int i=0;i<16;++i) {
        std::cout << i << ": " << std::fixed << std::setprecision(3) << offsets[i] << "\n";
    }



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

        // 保存
    aligner.cfg.out_root = output_root;
    aligner.saveResults(offsets);

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

    // draw dummy page to close the pdf
    TCanvas *canvas_dummy = new TCanvas("canvas_dummy", "Dummy Canvas", 800, 600);
    canvas_dummy->Print((out_pdf + "]").c_str());
    
    output_root->Close();

    delete canvas_toa_frequency;
    delete canvas_raw_toa_ns;
    delete canvas_toa_adc_max_corr;
    delete canvas_dummy;

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}