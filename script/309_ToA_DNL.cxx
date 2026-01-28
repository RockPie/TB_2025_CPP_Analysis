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

    // * --- Fine TDC Code distribution ----------------------------------------------------------
    // only 3 bits
    std::vector<TH1I*> h1d_fine_tdc_code_list;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h1d_fine_tdc_code_ch" + std::to_string(i);
        TH1I *h1d_fine_tdc_code = new TH1I(hist_name.c_str(), (hist_name + ";Fine TDC Code;Counts").c_str(), 8, 0, 8);
        h1d_fine_tdc_code->SetDirectory(nullptr);
        h1d_fine_tdc_code_list.push_back(h1d_fine_tdc_code);
    }

    // * --- Coarse TDC Code distribution ----------------------------------------------------------
    // 5 bits
    std::vector<TH1I*> h1d_coarse_tdc_code_list;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        std::string hist_name = "h1d_coarse_tdc_code_ch" + std::to_string(i);
        TH1I *h1d_coarse_tdc = new TH1I(hist_name.c_str(), (hist_name + ";Coarse TDC Code;Counts").c_str(), 32, 0, 32);
        h1d_coarse_tdc->SetDirectory(nullptr);
        h1d_coarse_tdc_code_list.push_back(h1d_coarse_tdc);
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
                    Int_t fine_tdc_code = toa_first_value & 0x07;
                    Int_t coarse_tdc_code = (toa_first_value >> 3) & 0x1F;
                    h1d_fine_tdc_code_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(fine_tdc_code);
                    h1d_coarse_tdc_code_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(coarse_tdc_code);
                    toa_ns += 25.0 * static_cast<double>(toa_first_index);
                    h1d_raw_toa_ns_list[channel + vldb_id * FPGA_CHANNEL_NUMBER]->Fill(toa_ns);
                    toa_code_matrix[channel + vldb_id * FPGA_CHANNEL_NUMBER].push_back(toa_first_value);
                    toa_ns_matrix[channel + vldb_id * FPGA_CHANNEL_NUMBER].push_back(toa_ns);
                }
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
    
    TCanvas *canvas_raw_toa_ns = new TCanvas("canvas_raw_toa_ns", "Raw ToA in ns", 1200, 800);
    draw_mosaic_fixed(*canvas_raw_toa_ns, h1d_raw_toa_ns_list, topo_ped_median);
    canvas_raw_toa_ns->Modified();
    canvas_raw_toa_ns->Update();
    canvas_raw_toa_ns->Print(out_pdf.c_str());
    canvas_raw_toa_ns->Write();
    canvas_raw_toa_ns->Close();

    TCanvas *canvas_fine_tdc_code = new TCanvas("canvas_fine_tdc_code", "Fine TDC Code", 1200, 800);
    draw_mosaic_fixed(*canvas_fine_tdc_code, h1d_fine_tdc_code_list, topo_ped_median);
    canvas_fine_tdc_code->Modified();
    canvas_fine_tdc_code->Update();
    canvas_fine_tdc_code->Print(out_pdf.c_str());
    canvas_fine_tdc_code->Write();
    canvas_fine_tdc_code->Close();

    TCanvas *canvas_coarse_tdc_code = new TCanvas("canvas_coarse_tdc_code", "Coarse TDC Code", 1200, 800);
    draw_mosaic_fixed(*canvas_coarse_tdc_code, h1d_coarse_tdc_code_list, topo_ped_median);
    canvas_coarse_tdc_code->Modified();
    canvas_coarse_tdc_code->Update();
    canvas_coarse_tdc_code->Print(out_pdf.c_str());
    canvas_coarse_tdc_code->Write();
    canvas_coarse_tdc_code->Close();

    output_root->cd();
    const int example_vldb_id = 1;
    const int example_asic_id = 0;
    const int example_half_id = 1;
    std::vector<int> color_list = {kBlue, kRed, kGreen+2, kMagenta+2, kCyan+2, kOrange+7, kViolet+7, kSpring+7, kTeal+7, kAzure+7, kPink+7};
    TCanvas *example_fine_tdc_codes_canvas = new TCanvas("example_fine_tdc_codes_canvas", "Example Fine TDC Codes", 800, 600);
    example_fine_tdc_codes_canvas->cd();

    // First pass: find global maximum
    double global_max = 0;
    for (int chn = 0; chn < FPGA_CHANNEL_NUMBER*vldb_number; chn++) {
        int vldb_id = chn / FPGA_CHANNEL_NUMBER;
        int chn_in_vldb = chn % FPGA_CHANNEL_NUMBER;
        int chn_in_asic = chn_in_vldb % 76;
        int asic_id = chn_in_vldb / 76;
        int half_id = chn_in_asic / 38;
        if (asic_id == example_asic_id && half_id == example_half_id && vldb_id == example_vldb_id) {
            double hist_max = h1d_fine_tdc_code_list[chn]->GetMaximum();
            if (hist_max > global_max) global_max = hist_max;
        }
    }

    // Second pass: draw with correct range
    bool first_hist_drawn = false;
    TLegend *legend_example_fine_tdc = new TLegend(0.8,0.5,0.89,0.89);
    legend_example_fine_tdc->SetBorderSize(0);
    legend_example_fine_tdc->SetFillStyle(0);
    for (int chn = 0; chn < FPGA_CHANNEL_NUMBER*vldb_number; chn++) {
        int vldb_id = chn / FPGA_CHANNEL_NUMBER;
        int chn_in_vldb = chn % FPGA_CHANNEL_NUMBER;
        int chn_in_asic = chn_in_vldb % 76;
        int asic_id = chn_in_vldb / 76;
        int half_id = chn_in_asic / 38;
        if (asic_id == example_asic_id && half_id == example_half_id && vldb_id == example_vldb_id) {
            TH1I* h1d_fine_tdc_copy = new TH1I(*h1d_fine_tdc_code_list[chn]);
            // normalize
            // h1d_fine_tdc_copy->Scale(1.0 / h1d_fine_tdc_copy->GetEntries());
            if (!first_hist_drawn) {
                // remove title and stat box for the first histogram
                // h1d_fine_tdc_copy->SetTitle("");
                // h1d_fine_tdc_copy->SetStats(0);
                // h1d_fine_tdc_copy->GetXaxis()->SetTitle("Fine TDC Code");
                // h1d_fine_tdc_copy->GetYaxis()->SetTitle("Counts");
                // h1d_fine_tdc_copy->GetYaxis()->SetRangeUser(0, global_max*1.5);  // Use global max
                // h1d_fine_tdc_copy->SetLineColor(kBlue);
                // h1d_fine_tdc_copy->Draw("HIST");

                format_1i_hist_canvas(example_fine_tdc_codes_canvas, h1d_fine_tdc_copy, kBlue, annotation_canvas_title, annotation_testbeam_title, "Fine TDC Code Distribution", static_cast<int>(global_max*1.5));
                first_hist_drawn = true;
            } else {
                h1d_fine_tdc_copy->SetLineColor(color_list[(chn % color_list.size())]);
                h1d_fine_tdc_copy->Draw("SAME");
            }
            legend_example_fine_tdc->AddEntry(h1d_fine_tdc_copy, ("Ch" + std::to_string(chn)).c_str(), "l");
        }
    }
    legend_example_fine_tdc->Draw();
    example_fine_tdc_codes_canvas->Modified();
    example_fine_tdc_codes_canvas->Update();
    example_fine_tdc_codes_canvas->Print(out_pdf.c_str());
    example_fine_tdc_codes_canvas->Write();
    example_fine_tdc_codes_canvas->Close();

    // Draw the example coarse TDC codes
    TCanvas *example_coarse_tdc_codes_canvas = new TCanvas("example_coarse_tdc_codes_canvas", "Example Coarse TDC Codes", 800, 600);
    example_coarse_tdc_codes_canvas->cd();
    // First pass: find global maximum
    global_max = 0;
    for (int chn = 0; chn < FPGA_CHANNEL_NUMBER*vldb_number; chn++) {
        int vldb_id = chn / FPGA_CHANNEL_NUMBER;
        int chn_in_vldb = chn % FPGA_CHANNEL_NUMBER;
        int chn_in_asic = chn_in_vldb % 76;
        int asic_id = chn_in_vldb / 76;
        int half_id = chn_in_asic / 38;
        if (asic_id == example_asic_id && half_id == example_half_id && vldb_id == example_vldb_id) {
            double hist_max = h1d_coarse_tdc_code_list[chn]->GetMaximum();
            if (hist_max > global_max) global_max = hist_max;
        }
    }

    // Second pass: draw with correct range
    first_hist_drawn = false;
    TLegend *legend_example_coarse_tdc = new TLegend(0.8,0.5,0.89,0.89);
    legend_example_coarse_tdc->SetBorderSize(0);
    legend_example_coarse_tdc->SetFillStyle(0);
    for (int chn = 0; chn < FPGA_CHANNEL_NUMBER*vldb_number; chn++) {
        int vldb_id = chn / FPGA_CHANNEL_NUMBER;
        int chn_in_vldb = chn % FPGA_CHANNEL_NUMBER;
        int chn_in_asic = chn_in_vldb % 76;
        int asic_id = chn_in_vldb / 76;
        int half_id = chn_in_asic / 38;
        if (asic_id == example_asic_id && half_id == example_half_id && vldb_id == example_vldb_id) {
            TH1I* h1d_coarse_tdc_copy = new TH1I(*h1d_coarse_tdc_code_list[chn]);
            // normalize
            // h1d_coarse_tdc_copy->Scale(1.0 / h1d_coarse_tdc_copy->GetEntries());
            if (!first_hist_drawn) {
                format_1i_hist_canvas(example_coarse_tdc_codes_canvas, h1d_coarse_tdc_copy, kBlue, annotation_canvas_title, annotation_testbeam_title, "Coarse TDC Code Distribution", static_cast<int>(global_max*1.5));
                first_hist_drawn = true;
            } else {
                h1d_coarse_tdc_copy->SetLineColor(color_list[(chn % color_list.size())]);
                h1d_coarse_tdc_copy->Draw("SAME");
            }
            legend_example_coarse_tdc->AddEntry(h1d_coarse_tdc_copy, ("Ch" + std::to_string(chn)).c_str(), "l");
        }
    }
    legend_example_coarse_tdc->Draw();
    example_coarse_tdc_codes_canvas->Modified();
    example_coarse_tdc_codes_canvas->Update();
    example_coarse_tdc_codes_canvas->Print(out_pdf.c_str());
    example_coarse_tdc_codes_canvas->Write();
    example_coarse_tdc_codes_canvas->Close();


    // calculate DNL from coarse TDC code histograms
    // create new directory in output root file
    output_root->cd();
    TDirectory* dnl_directory = output_root->mkdir("DNL_Calibration");
    dnl_directory->cd();
    for (int i = 0; i < FPGA_CHANNEL_NUMBER * vldb_number; i++) {
        TH1I* h1d_coarse_tdc = h1d_coarse_tdc_code_list[i];
        std::vector<double> code_centers;
        std::vector<double> normalized_counts;
        int n_bins = h1d_coarse_tdc->GetNbinsX();
        int total_counts = h1d_coarse_tdc->GetEntries();
        for (int bin = 1; bin <= n_bins; bin++) {
            double center = h1d_coarse_tdc->GetBinCenter(bin);
            double count = static_cast<double>(h1d_coarse_tdc->GetBinContent(bin));
            double normalized_count = (total_counts > 0) ? (count / static_cast<double>(total_counts)) : 0.0;
            normalized_counts.push_back(normalized_count);
        }
        double current_center = 0.0;
        for (int bin = 0; bin < n_bins; bin++) {
            current_center += normalized_counts[bin] * 0.5;
            code_centers.push_back(current_center);
            current_center += normalized_counts[bin] * 0.5;
        }
        // write to TVectorD
        TVectorD dnl_coarse_tdc(code_centers.size(), code_centers.data());
        const std::string dnl_coarse_name = "dnl_coarse_tdc_ch" + std::to_string(i);
        dnl_coarse_tdc.Write(dnl_coarse_name.c_str());
        
        TH1I* h1d_fine_tdc = h1d_fine_tdc_code_list[i];
        std::vector<double> fine_code_centers;
        std::vector<double> fine_normalized_counts;
        int n_fine_bins = h1d_fine_tdc->GetNbinsX();
        int total_fine_counts = h1d_fine_tdc->GetEntries();
        for (int bin = 1; bin <= n_fine_bins; bin++) {
            double center = h1d_fine_tdc->GetBinCenter(bin);
            double count = static_cast<double>(h1d_fine_tdc->GetBinContent(bin));
            double normalized_count = (total_fine_counts > 0) ? (count / static_cast<double>(total_fine_counts)) : 0.0;
            fine_normalized_counts.push_back(normalized_count);
        }
        double current_fine_center = 0.0;
        for (int bin = 0; bin < n_fine_bins; bin++) {
            current_fine_center += fine_normalized_counts[bin] * 0.5;
            fine_code_centers.push_back(current_fine_center);
            current_fine_center += fine_normalized_counts[bin] * 0.5;
        }

        TVectorD dnl_fine_tdc(fine_code_centers.size(), fine_code_centers.data());
        const std::string dnl_fine_name = "dnl_fine_tdc_ch" + std::to_string(i);
        dnl_fine_tdc.Write(dnl_fine_name.c_str());
    }
    output_root->cd();



    // draw dummy page to close the pdf
    TCanvas *canvas_dummy = new TCanvas("canvas_dummy", "Dummy Canvas", 800, 600);
    canvas_dummy->Print((out_pdf + "]").c_str());

    output_root->Close();

    delete canvas_toa_frequency;
    delete canvas_raw_toa_ns;
    delete canvas_dummy;

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}