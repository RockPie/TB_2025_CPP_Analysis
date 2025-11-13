#include <fstream>
#include <limits>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "H2GCROC_TimewalkLUT.hxx"

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
        spdlog::info("Channel {}: DNL coarse TDC size = {}, fine TDC size = {}",
                     chn,
                     dnl_coarse_tdc->GetNoElements(),
                     dnl_fine_tdc->GetNoElements());
        spdlog::info("Coarse DNL: ");
        for (const auto& val : dnl_coarse_tdc_list[chn]) {
            spdlog::info("{:.6f}, ", val);
        }
        spdlog::info("Fine DNL: ");
        for (const auto& val : dnl_fine_tdc_list[chn]) {
            spdlog::info("{:.6f}, ", val);
        }
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

    std::vector<std::vector<UInt_t>> toa_code_matrix(FPGA_CHANNEL_NUMBER * vldb_number, std::vector<UInt_t>());
    for (auto & vec : toa_code_matrix) {
        vec.reserve(entry_max);
    }

    std::vector<std::vector<double>> toa_ns_matrix(FPGA_CHANNEL_NUMBER * vldb_number, std::vector<double>());
    for (auto & vec : toa_ns_matrix) {
        vec.reserve(entry_max);
    }

    const int adc_peak_min_index = 4;
    const int adc_peak_max_index = 8;

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
                    double channel_optimal_bx_slip = toa_correction_map[channel + vldb_id * FPGA_CHANNEL_NUMBER];
                    if (toa_first_value >= static_cast<UInt_t>(channel_optimal_bx_slip)) {
                        toa_ns -= 25.0;
                    }
                    h1d_raw_toa_ns_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_ns);
                    h2d_toa_adc_max_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(adc_peak_value_pede_sub, toa_ns);
                    h2d_toa_code_corr_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_first_value, toa_ns);
                    toa_code_matrix[channel + vldb_id * FPGA_CHANNEL_NUMBER].push_back(toa_first_value);
                    toa_ns_matrix[channel + vldb_id * FPGA_CHANNEL_NUMBER].push_back(toa_ns);
                }
                    
                // } else if (toa_showup_times > 1) {
                //     num_multi_ToA++;
                //     // multiple ToA hits
                //     // currently, just take the first one
                //     // double toa_ns = decode_toa_value_ns(toa_first_value);
                //     // toa_ns += 25.0 * static_cast<double>(toa_first_index);
                //     // h1d_raw_toa_ns_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_ns);
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

    const int Nmin = 300;
    const int max_expand = 6;
    const bool enforce_monotonic = true;

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

    // do the timewalk LUT generation here
    // create a new directory
    output_root->cd();
    output_root->mkdir("Timewalk_LUTs");
    output_root->cd("Timewalk_LUTs");
    std::vector<TGraph*> lut_graphs;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * vldb_number; _chn++) {
        auto h2d_chn = h2d_toa_adc_max_corr_list[_chn];
        // rebin
        h2d_chn->RebinX(
            (h2d_chn->GetNbinsX() / 96) ?
             h2d_chn->GetNbinsX() / 96 : 1);
        auto lut = BuildFitFreeTimewalkLUT(
            h2d_chn,
            Nmin,
            max_expand,
            enforce_monotonic
        );
        auto g_lut = lut.AsGraph();
        g_lut->SetName(("timewalk_lut_ch" + std::to_string(_chn)).c_str());
        g_lut->SetTitle("Fit-free time-walk LUT;ADC;#Delta t_{median} [ps]");
        g_lut->Write();
        lut_graphs.push_back(g_lut);
    }
    output_root->cd();

    TCanvas *canvas_toa_adc_max_fit_corr = new TCanvas("canvas_toa_adc_max_fit_corr", "ToA ADC Max Fit Correction", 1200, 800);
    draw_mosaic_fixed(*canvas_toa_adc_max_fit_corr, h2d_toa_adc_max_corr_list, lut_graphs, topo_wave);
    canvas_toa_adc_max_fit_corr->Modified();
    canvas_toa_adc_max_fit_corr->Update();
    canvas_toa_adc_max_fit_corr->Print(out_pdf.c_str());
    canvas_toa_adc_max_fit_corr->Write();
    canvas_toa_adc_max_fit_corr->Close();

    TCanvas *canvas_toa_code_corr = new TCanvas("canvas_toa_code_corr", "ToA Code Correction", 1200, 800);
    draw_mosaic_fixed(*canvas_toa_code_corr, h2d_toa_code_corr_list, topo_wave);
    canvas_toa_code_corr->Modified();
    canvas_toa_code_corr->Update();
    canvas_toa_code_corr->Print(out_pdf.c_str());
    canvas_toa_code_corr->Write();
    canvas_toa_code_corr->Close();

    TCanvas *canvas_toa_adc_max_corr_example = new TCanvas("canvas_toa_adc_max_corr_example", "ToA ADC Max Correction Example Channel", 800, 600);
    auto h2d_example = h2d_toa_adc_max_corr_list[chn_example];
    h2d_example->GetXaxis()->SetTitle("ADC Peak (pedestal subtracted)");
    h2d_example->GetYaxis()->SetTitle("Corrected ToA [ns]");
    h2d_example->RebinX(
        (h2d_example->GetNbinsX() / 96) ?
         h2d_example->GetNbinsX() / 96 : 1);
    format_2d_hist_canvas(canvas_toa_adc_max_corr_example, h2d_example, kBlue+2, annotation_canvas_title, annotation_testbeam_title, "Channel_" + std::to_string(chn_example));
    TimewalkLUT lut = BuildFitFreeTimewalkLUT(h2d_example, Nmin, max_expand, enforce_monotonic);
    TLegend *legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    std::unique_ptr<TGraph> g_lut_example(lut.AsGraph());
    g_lut_example->SetLineColor(kRed+2);
    g_lut_example->SetLineWidth(2);
    g_lut_example->Draw("L SAME");
    legend->AddEntry(g_lut_example.get(), "Timewalk LUT", "l");
    legend->Draw();
    canvas_toa_adc_max_corr_example->Modified();
    canvas_toa_adc_max_corr_example->Update();

    canvas_toa_adc_max_corr_example->Print(out_pdf.c_str());
    canvas_toa_adc_max_corr_example->Write();
    canvas_toa_adc_max_corr_example->Close();

    auto g_lut = lut.AsGraph();
    g_lut->SetName(("timewalk_lut_ch" + std::to_string(chn_example)).c_str());
    g_lut->SetTitle("Fit-free time-walk LUT;ADC;#Delta t_{median} [ps]");
    g_lut->Write();

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