#include <fstream>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <TStyle.h>
#include <TParameter.h>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"

int main(int argc, char **argv) {

    // * --- Basic setup ------------------------------------------------------
    // * ----------------------------------------------------------------------
    gROOT->SetBatch(kTRUE);
    std::string script_input_file, script_output_file;

    int script_n_events = -1;
    bool script_verbose = false;

    std::string script_name = __FILE__;
    script_name = script_name.substr(script_name.find_last_of("/\\") + 1).substr(0, script_name.find_last_of("."));
    const std::string script_version = "0.1";
    std::string script_output_folder;

    configure_logger(false);

    // * --- Handle input arguments -------------------------------------------
    // * ----------------------------------------------------------------------
    cxxopts::Options options(script_name, "Waveform level analysis");

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
    std::string script_input_run_number = script_input_file.substr(script_input_file.find_last_of("Run") + 1).append(4, '0').substr(0, 4);
    script_output_file = parsed["output"].as<std::string>();
    script_n_events    = parsed["events"].as<int>();
    script_verbose     = parsed["verbose"].as<bool>();

    configure_logger(script_verbose);

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
    // open the input root file
    TFile *input_root = TFile::Open(script_input_file.c_str(), "READ");
    if (input_root->IsZombie() || !input_root->IsOpen()) {
        spdlog::error("Failed to open input file {}", script_input_file);
        return 1;
    }

    // create the output root file
    TFile *output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
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
    const auto& sipm_board      = mapping_json.at("SiPM_Board");
    const auto& board_loc       = mapping_json.at("Board_Loc");
    const auto& board_rotation  = mapping_json.at("Board_Rotation");
    const auto& board_flip      = mapping_json.at("Board_Flip");

    // load the example timewalk calibration file
    std::string timewalk_json_file = "config/0416_example_adcmax_toa.json";
    std::ifstream timewalk_json_ifs(timewalk_json_file);
    if (!timewalk_json_ifs.is_open()) {
        spdlog::error("Failed to open timewalk json file {}", timewalk_json_file);
        return 1;
    }
    json timewalk_json;
    timewalk_json_ifs >> timewalk_json;
    
    // * --- Main code --------------------------------------------------------
    // * ----------------------------------------------------------------------
    const int machine_gun_samples = 16;
    const int vldb_number = 2;

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

    const double toa_ns_min = -0.5 * (double(machine_gun_samples) * 25.0);
    const double toa_ns_max = 0.5 * (double(machine_gun_samples) * 25.0);

    // create the 2D histograms for each SiPM channel
    std::vector<TH2D*> h2_waveforms;
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int hist_x_bins = machine_gun_samples;
            int hist_y_bins = 512;
            TH2D* h2_waveform = new TH2D(
                Form("h2_waveform_vldb%d_channel%d", vldb_id, channel),
                Form("Waveform (VLDB %d, Channel %d);Time (ns);ADC Counts", vldb_id, channel),
                hist_x_bins, 0, hist_x_bins,
                hist_y_bins, 0, 1024
            );
            h2_waveforms.push_back(h2_waveform);
        }
    }

    // create the 2D histograms for pedestal - adc max correlation
    std::vector<TH2D*> h2_pedestal_median_adcmax;
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int hist_x_bins = 512;
            int hist_y_bins = 256;
            TH2D* h2_ped_adcmax = new TH2D(
                Form("h2_pedestal_median_adcmax_vldb%d_channel%d", vldb_id, channel),
                Form("Pedestal vs ADC Max (VLDB %d, Channel %d);Pedestal;ADC Max", vldb_id, channel),
                hist_x_bins, 0, 1024,
                hist_y_bins, 0, 256
            );
            h2_pedestal_median_adcmax.push_back(h2_ped_adcmax);
        }
    }

    std::vector<TH2D*> h2_pedestal_median_pedestal_sigma;
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int hist_x_bins = 512;
            int hist_y_bins = 64;
            TH2D* h2_pede_sigma = new TH2D(
                Form("h2_pedestal_median_pedestal_sigma_vldb%d_channel%d", vldb_id, channel),
                Form("Pedestal vs Pedestal Sigma (VLDB %d, Channel %d);Pedestal;Pedestal Sigma", vldb_id, channel),
                hist_x_bins, 0, 512,
                hist_y_bins, 0, hist_y_bins
            );
            h2_pedestal_median_pedestal_sigma.push_back(h2_pede_sigma);
        }
    }

    std::vector<TH2D*> h2_adcmax_toa_correlation;
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int hist_x_bins = 256;
            int hist_y_bins = 25 * machine_gun_samples;
            TH2D* h2_adcmax_toa = new TH2D(
                Form("h2_adcmax_toa_correlation_vldb%d_channel%d", vldb_id, channel),
                Form("ADC Max vs ToA (VLDB %d, Channel %d);ToA;ADC Max", vldb_id, channel),
                hist_x_bins, 0, 1024,
                hist_y_bins, toa_ns_min, toa_ns_max
            );
            h2_adcmax_toa_correlation.push_back(h2_adcmax_toa);
        }
    }

    // 1d toa value in ns distribution
    std::vector<TH1D*> h1_toa_distributions;
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int hist_bins = 100;
            TH1D* h1_toa = new TH1D(
                Form("h1_toa_distribution_vldb%d_channel%d", vldb_id, channel),
                Form("ToA Distribution (VLDB %d, Channel %d);ToA (ns);Counts", vldb_id, channel),
                hist_bins, toa_ns_min, toa_ns_max
            );
            h1_toa_distributions.push_back(h1_toa);
        }
    }

    // toa code - toa ns correlation
    std::vector<TH2D*> h2_toa_correlation;
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int hist_x_bins = 256;
            int hist_y_bins = 25 * machine_gun_samples;
            TH2D* h2_toa_corr = new TH2D(
                Form("h2_toa_code_ns_correlation_vldb%d_channel%d", vldb_id, channel),
                Form("ToA Code vs ToA ns (VLDB %d, Channel %d);ToA Code;ToA (ns)", vldb_id, channel),
                hist_x_bins, 0, 1024,
                hist_y_bins, toa_ns_min, toa_ns_max
            );
            h2_toa_correlation.push_back(h2_toa_corr);
        }
    }

    std::vector<TH2D*> h2_waveform_toa_offset;
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            int hist_x_bins = machine_gun_samples * 25; // 1 bin per ns
            int hist_y_bins = 256;
            TH2D* h2_waveform_toa = new TH2D(
                Form("h2_waveform_toa_offset_vldb%d_channel%d", vldb_id, channel),
                Form("Waveform ToA Offset (VLDB %d, Channel %d);Sample;ToA Offset (ns)", vldb_id, channel),
                hist_x_bins, 0, 25.0 * machine_gun_samples,
                hist_y_bins, 0, 1024
            );
            h2_waveform_toa_offset.push_back(h2_waveform_toa);
        }
    }

    // monitor hamming code errors
    size_t total_hamming_code_errors = 0;
    TH1I* h1_hamming_code_errors = new TH1I("h1_hamming_code_errors", "Hamming Code Errors per Event;Hamming Code Errors;Counts", 8, 0, 8);
    size_t total_daqh_header_footer_errors = 0;

    size_t total_hamming_pass_events = 0;
    size_t total_non_zero_toa_counts = 0;
    size_t total_valid_toa_counts = 0;

    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        Int_t hamming_code_errors_in_event = 0;
        // check hamming code for each VLDB
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            for (int asic_id = 0; asic_id < 2; asic_id++) {
                for (int half_id = 0; half_id < 2; half_id++) {
                    auto daqh_value = daqh_list_pools[vldb_id][0][asic_id * 2 + half_id];
                    auto hamming_code = (daqh_value & 0x00000070) >> 4;
                    if (hamming_code != 0) {
                        total_hamming_code_errors++;
                        hamming_code_errors_in_event++;
                    }
                    auto daqh_header = (daqh_value & 0xF0000000) >> 28;
                    auto daqh_footer = (daqh_value & 0x0000000F);
                    if (daqh_header != 0xF || (daqh_footer != 0x5 && daqh_footer != 0x2)) {
                        total_daqh_header_footer_errors++;
                    }
                } // half_id
            } // asic_id
        } // vldb_id
        h1_hamming_code_errors->Fill(hamming_code_errors_in_event);
        if (hamming_code_errors_in_event > 0) {
            continue;
        }
        total_hamming_pass_events++;
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                const int hist_index = vldb_id * FPGA_CHANNEL_NUMBER + channel;
                std::vector<UInt_t> adc_pedestal_array;
                std::vector<UInt_t> adc_samples;
                UInt_t adc_max = 0;
                UInt_t toa_max = 0;
                int toa_max_index = -1;
                int non_zero_toa_counts = 0;
                int asic_id = channel / 76;
                int half_id = (channel % 76) / 38;
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int hist_index = vldb_id * FPGA_CHANNEL_NUMBER + channel;
                    UInt_t adc_value = val0_list_pools[vldb_id][0][channel + sample * FPGA_CHANNEL_NUMBER];
                    UInt_t toa_value = val2_list_pools[vldb_id][0][channel + sample * FPGA_CHANNEL_NUMBER];
                    h2_waveforms[hist_index]->Fill(sample, adc_value);

                    if (sample < 3) {
                        adc_pedestal_array.push_back(adc_value);
                    }
                    adc_samples.push_back(adc_value);

                    if (adc_value > adc_max) {
                        adc_max = adc_value;
                    }

                    if (toa_value > toa_max) {
                        toa_max = toa_value;
                        toa_max_index = sample;
                    }

                    if (toa_value != 0) {
                        non_zero_toa_counts++;
                        total_non_zero_toa_counts++;
                    }
                }
                auto pedestal_median = pedestal_median_of_first3(adc_pedestal_array);
                // use max-min as pedestal sigma
                auto pedestal_sigma = *std::max_element(adc_pedestal_array.begin(), adc_pedestal_array.end()) -
                                      *std::min_element(adc_pedestal_array.begin(), adc_pedestal_array.end());
                double adc_max_double = static_cast<double>(adc_max);
                double adc_max_minus_pedestal = adc_max_double - static_cast<double>(pedestal_median);

                h2_pedestal_median_adcmax[hist_index]->Fill(adc_max_double, pedestal_median);
                h2_pedestal_median_pedestal_sigma[hist_index]->Fill(pedestal_median, pedestal_sigma);

                // handle ToA value
                if (non_zero_toa_counts == 1 && toa_max != 0) {
                    total_valid_toa_counts++;
                    auto toa_value_ns = decode_toa_value_ns(toa_max);
                    toa_value_ns += static_cast<double>(toa_max_index) * 25.0;  // add the sample offset
                    toa_manual_correction(vldb_id, asic_id, half_id, toa_max, toa_value_ns);
                    auto timewalk_correction = find_closest_toa(timewalk_json, adc_max_minus_pedestal);
                    // spdlog::info("Event {}: VLDB {} Channel {}: ToA (corrected) = {}", entry, vldb_id, channel, toa_value_ns);
                    // debug: check the toa value
                    // spdlog::info("Event {}: VLDB {} Channel {}: ToA = {} (raw: {})", entry, vldb_id, channel, toa_value_ns, toa_max);
                    toa_value_ns -= timewalk_correction;
                    h2_adcmax_toa_correlation[hist_index]->Fill(adc_max_minus_pedestal, toa_value_ns);
                    h1_toa_distributions[hist_index]->Fill(toa_value_ns);
                    h2_toa_correlation[hist_index]->Fill(toa_max, toa_value_ns);

                    // fill waveform toa offset histogram
                    for (int sample = 0; sample < machine_gun_samples; sample++) {
                        double sample_time_ns = static_cast<double>(sample) * 25.0 - toa_value_ns;
                        double sample_amplitude = static_cast<double>(adc_samples[sample]) - static_cast<double>(pedestal_median) + 50.0;
                        h2_waveform_toa_offset[hist_index]->Fill(sample_time_ns, sample_amplitude);
                    }
                }
            }
        }
    }
    input_root->Close();

    // // parameterize the curve of vldb1, channel 106
    // auto h2_example = h2_adcmax_toa_correlation[1 * FPGA_CHANNEL_NUMBER + 106];
    // // create TGraph from TH2D
    // std::vector<double> x_vals;
    // std::vector<double> y_vals;
    // std::vector<double> y_errs;
    // for (int x_bin = 1; x_bin <= h2_example->GetNbinsX(); x_bin++) {
    //     double x_center = h2_example->GetXaxis()->GetBinCenter(x_bin);
    //     double y_sum = 0;
    //     double y_count = 0;
    //     for (int y_bin = 1; y_bin <= h2_example->GetNbinsY(); y_bin++) {
    //         double bin_content = h2_example->GetBinContent(x_bin, y_bin);
    //         if (bin_content > 0) {
    //             double y_center = h2_example->GetYaxis()->GetBinCenter(y_bin);
    //             y_sum += y_center * bin_content;
    //             y_count += bin_content;
    //         }
    //     }
    //     if (y_count > 0) {
    //         double y_mean = y_sum / y_count;
    //         x_vals.push_back(x_center);
    //         y_vals.push_back(y_mean);
    //         // estimate error as standard deviation
    //         double y_var_sum = 0;
    //         for (int y_bin = 1; y_bin <= h2_example->GetNbinsY(); y_bin++) {
    //             double bin_content = h2_example->GetBinContent(x_bin, y_bin);
    //             if (bin_content > 0) {
    //                 double y_center = h2_example->GetYaxis()->GetBinCenter(y_bin);
    //                 y_var_sum += bin_content * (y_center - y_mean) * (y_center - y_mean);
    //             }
    //         }
    //         double y_stddev = std::sqrt(y_var_sum / y_count);
    //         y_errs.push_back(y_stddev/std::sqrt(y_count));
    //     }
    // }
    // TGraphErrors* graph_example = new TGraphErrors(x_vals.size(), x_vals.data(), y_vals.data(), nullptr, y_errs.data());
    // graph_example->SetName("graph_example_vldb1_channel106");
    // graph_example->SetTitle("Example ADC Max vs ToA (VLDB 1, Channel 106);ADC Max - Pedestal;ToA (ns)");

    // // save the data as json file
    // std::string example_json_file = script_output_folder + "/" + script_input_run_number + "_example_adcmax_toa.json";
    // std::ofstream example_json_ofs(example_json_file);
    // if (!example_json_ofs.is_open()) {
    //     spdlog::error("Failed to open example json file {}", example_json_file);
    //     return 1;
    // }
    // json example_json;
    // example_json["adc_max_minus_pedestal"] = x_vals;
    // example_json["toa_ns"] = y_vals;
    // example_json_ofs << example_json.dump(4);
    // example_json_ofs.close();

    // output_root->cd();
    // graph_example->Write();

    // calculate QA values
    double hamming_code_error_rate = static_cast<double>(total_hamming_code_errors) / static_cast<double>(entry_max) / 4.0;
    double daqh_header_footer_error_rate = static_cast<double>(total_daqh_header_footer_errors) / static_cast<double>(entry_max) / 4.0;
    double toa_valid_fraction = static_cast<double>(total_valid_toa_counts) / static_cast<double>(vldb_number * FPGA_CHANNEL_NUMBER * total_hamming_pass_events);

    // print out total non-zero ToA counts and valid ToA counts
    spdlog::info("Total Hamming code errors rate: {} / {} = {:.6f}%", total_hamming_code_errors, entry_max, static_cast<double>(total_hamming_code_errors) / static_cast<double>(entry_max) * 25.0);
    spdlog::info("Total non-zero ToA counts: {}", total_non_zero_toa_counts);
    spdlog::info("Total valid ToA counts (only one non-zero ToA per channel): {}", total_valid_toa_counts);
    spdlog::info("Valid ToA fraction: {:.2f}%", static_cast<double>(total_valid_toa_counts) / static_cast<double>(total_non_zero_toa_counts) * 100.0);

    // * --- Handle output ----------------------------------------------------
    // * ----------------------------------------------------------------------
    const int NX = 16, NY = 12;
    const int board_cols = 8, board_rows = 4;

    output_root->cd();

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

    TCanvas canvas_waveform("canvas_waveform", "Waveforms", 1200, 800);
    draw_mosaic_fixed(canvas_waveform, h2_waveforms, topo_wave);
    output_root->cd();
    canvas_waveform.Modified();
    canvas_waveform.Update();
    canvas_waveform.Write();

    TCanvas canvas_pedestal_adcmax("canvas_pedestal_adcmax", "Pedestal vs ADC Max", 1200, 800);
    draw_mosaic_fixed(canvas_pedestal_adcmax, h2_pedestal_median_adcmax, topo_ped_median);
    output_root->cd();
    canvas_pedestal_adcmax.Modified();
    canvas_pedestal_adcmax.Update();
    canvas_pedestal_adcmax.Write();
    // save as png file too
    //canvas_pedestal_adcmax.SaveAs((script_output_folder + "/" + script_input_run_number + "_pedestal_median_adcmax.png").c_str());

    TCanvas canvas_pedestal_median_pedestal_sigma("canvas_pedestal_median_pedestal_sigma", "Pedestal vs Pedestal Sigma", 1200, 800);
    draw_mosaic_fixed(canvas_pedestal_median_pedestal_sigma, h2_pedestal_median_pedestal_sigma, topo_ped_median);
    output_root->cd();
    canvas_pedestal_median_pedestal_sigma.Modified();
    canvas_pedestal_median_pedestal_sigma.Update();
    canvas_pedestal_median_pedestal_sigma.Write();
    // save as png file too
    //canvas_pedestal_median_pedestal_sigma.SaveAs((script_output_folder + "/" + script_input_run_number + "_pedestal_median_pedestal_sigma.png").c_str());

    TCanvas canvas_adcmax_toa_correlation("canvas_adcmax_toa_correlation", "ADC Max vs ToA", 1200, 800);
    draw_mosaic_fixed(canvas_adcmax_toa_correlation, h2_adcmax_toa_correlation, topo_ped_median);
    output_root->cd();
    canvas_adcmax_toa_correlation.Modified();
    canvas_adcmax_toa_correlation.Update();
    canvas_adcmax_toa_correlation.Write();
    // save as png file too
    //canvas_adcmax_toa_correlation.SaveAs((script_output_folder + "/" + script_input_run_number + "_adcmax_toa_correlation.png").c_str());

    TCanvas canvas_hamming_code_errors("canvas_hamming_code_errors", "Hamming Code Errors per Event", 800, 600);
    h1_hamming_code_errors->Draw();
    output_root->cd();
    canvas_hamming_code_errors.Modified();
    canvas_hamming_code_errors.Update();
    canvas_hamming_code_errors.Write();

    TCanvas canvas_toa_distributions("canvas_toa_distributions", "ToA Distributions", 1200, 800);
    draw_mosaic_fixed(canvas_toa_distributions, h1_toa_distributions, topo_ped_median);
    output_root->cd();
    canvas_toa_distributions.Modified();
    canvas_toa_distributions.Update();
    canvas_toa_distributions.Write();

    TCanvas canvas_toa_correlation("canvas_toa_correlation", "ToA Code vs ToA ns Correlation", 1200, 800);
    draw_mosaic_fixed(canvas_toa_correlation, h2_toa_correlation, topo_ped_median);
    output_root->cd();
    canvas_toa_correlation.Modified();
    canvas_toa_correlation.Update();
    canvas_toa_correlation.Write();

    TCanvas canvas_waveform_toa_offset("canvas_waveform_toa_offset", "Waveform ToA Offset", 1200, 800);
    draw_mosaic_fixed(canvas_waveform_toa_offset, h2_waveform_toa_offset, topo_ped_median);
    output_root->cd();
    canvas_waveform_toa_offset.Modified();
    canvas_waveform_toa_offset.Update();
    canvas_waveform_toa_offset.Write();

    // write QA values to the root file as a TParameter
    TParameter<double> param_hamming_code_error_rate("hamming_code_error_rate", hamming_code_error_rate);
    param_hamming_code_error_rate.Write();
    TParameter<double> param_daqh_header_footer_error_rate("daqh_header_footer_error_rate", daqh_header_footer_error_rate);
    param_daqh_header_footer_error_rate.Write();
    TParameter<double> param_toa_valid_fraction("toa_valid_fraction", toa_valid_fraction);
    param_toa_valid_fraction.Write();

    output_root->Close();

    return 0;
}