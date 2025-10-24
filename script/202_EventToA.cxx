#include <fstream>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <TStyle.h>
#include <TParameter.h>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"

static double CrystalBallFcn(double *xx, double *pp) {
    const double x     = xx[0];
    const double N     = pp[0];
    const double alpha = pp[1];        // alpha>0 左尾；alpha<0 右尾
    const double n     = pp[2];
    const double mean  = pp[3];
    const double sigma = pp[4];

    if (sigma <= 0.0 || n <= 1.0) return 0.0;

    const double as   = std::fabs(alpha);
    const double t0   = (x - mean)/sigma * (alpha >= 0 ? 1.0 : -1.0); // 统一到“左尾”情况
    const double a2   = as*as;
    const double B    = n/as - as;
    const double A    = std::pow(n/as, n) * std::exp(-0.5*a2);

    double core, tail;
    if (t0 > -as) {
        core = std::exp(-0.5 * t0*t0);
        return N * core;
    } else {
        tail = A * std::pow(B - t0, -n);
        return N * tail;
    }
}

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

    const double toa_ns_min = -0.25 * (double(machine_gun_samples) * 25.0);
    const double toa_ns_max = 0.25 * (double(machine_gun_samples) * 25.0);

    std::vector<TH1D*> h_toa_distributions;
    // reserve for each event
    h_toa_distributions.reserve(entry_max);

    TH1I* h_max_sample_index = new TH1I("h_max_sample_index", "Max Sample Index Distribution;Sample Index;Counts", machine_gun_samples, -0.5, double(machine_gun_samples) - 0.5);

    // data buffers
    std::vector<double> toa_mean_values;
    std::vector<double> adc_sum_values;
    toa_mean_values.reserve(entry_max);
    adc_sum_values.reserve(entry_max);

    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        int hamming_code_errors_in_event = 0;
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            for (int asic_id = 0; asic_id < 2; asic_id++) {
                for (int half_id = 0; half_id < 2; half_id++) {
                    auto daqh_value = daqh_list_pools[vldb_id][0][asic_id * 2 + half_id];
                    auto hamming_code = (daqh_value & 0x00000070) >> 4;
                    if (hamming_code != 0) {
                        hamming_code_errors_in_event++;
                    }
                } // half_id
            } // asic_id
        } // vldb_id
        if (hamming_code_errors_in_event > 0) {
            spdlog::warn("Event {} has {} Hamming code errors", entry, hamming_code_errors_in_event);
            continue;
        }

        

        std::vector<double> toa_values_ns_list;
        std::vector<int> toa_value_channel_list;

        double adc_sum_event = 0.0;
        double adc_sum_valid_toa = 0.0;
        double event_toa_mean_weighted = 0.0;

        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int asic_id = channel / 76;
                int half_id = (channel % 76) / 38;

                UInt_t adc_max = 0;
                int adc_max_sample_index = -1;
                UInt_t toa_first = 0;
                int toa_first_index = -1;
                int non_zero_toa_count = 0;
                std::vector<UInt_t> adc_pedestal_samples;

                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    UInt_t adc_value = val0_list_pools[vldb_id][0][channel + sample * FPGA_CHANNEL_NUMBER];
                    UInt_t toa_value = val2_list_pools[vldb_id][0][channel + sample * FPGA_CHANNEL_NUMBER];
                    if (adc_value > 0 && (sample >= 5 && sample <= 8)) {
                        adc_max += adc_value;
                        adc_max_sample_index = sample;
                    }
                    if (toa_value != 0 && toa_first_index == -1) {
                        toa_first = toa_value;
                        toa_first_index = sample;
                    }
                    if (toa_value != 0) {
                        non_zero_toa_count++;
                    }
                    if (sample < 3) {
                        adc_pedestal_samples.push_back(adc_value);
                    }
                } // end of sample loop
                auto pedestal_median = pedestal_median_of_first3(adc_pedestal_samples);
                auto adc_max_minus_pedestal = double(adc_max) - pedestal_median * 4.0;
                if (adc_max_minus_pedestal < 0.0)
                    adc_max_minus_pedestal = 0.0;
                h_max_sample_index->Fill(adc_max_sample_index);
                adc_sum_event += adc_max_minus_pedestal;
                if (toa_first_index >= 0 && non_zero_toa_count == 1) {
                    auto toa_value_ns = decode_toa_value_ns(toa_first);
                    toa_value_ns += 25.0 * static_cast<double>(toa_first_index);
                    toa_manual_correction(vldb_id, asic_id, half_id, toa_first, toa_value_ns);
                    toa_value_ns -= find_closest_toa(timewalk_json, adc_max_minus_pedestal);
                    toa_values_ns_list.push_back(toa_value_ns);
                    toa_value_channel_list.push_back(vldb_id * FPGA_CHANNEL_NUMBER + channel);
                    adc_sum_valid_toa += adc_max_minus_pedestal;
                    event_toa_mean_weighted += toa_value_ns * adc_max_minus_pedestal;
                } // end of handling ToA
            } // end of channel loop
        } // end of vldb_id loop
        // print out all the valid ToA values in this event
        // spdlog::info("Event {}: Found {} valid ToA values", entry, toa_values_ns_list.size());
        // for (size_t i = 0; i < toa_values_ns_list.size(); i++) {
        //     auto channel = toa_value_channel_list[i] % FPGA_CHANNEL_NUMBER;
        //     auto vldb_id = toa_value_channel_list[i] / FPGA_CHANNEL_NUMBER;
        //     spdlog::info("  VLDB {} Channel {}: ToA = {:.2f} ns", vldb_id, channel, toa_values_ns_list[i]);
        // }
        double event_toa_mean = 0.0;
        if (adc_sum_valid_toa > 0.0)
            event_toa_mean = event_toa_mean_weighted / adc_sum_valid_toa;
        // h_toa_mean->Fill(event_toa_mean);
        // h_adc_sum_toa_mean->Fill(adc_sum_event, event_toa_mean);
        toa_mean_values.push_back(event_toa_mean);
        adc_sum_values.push_back(adc_sum_event);
    } // end of entry loop

    input_root->Close();

    output_root->cd();
    h_max_sample_index->Write();

    TH1D* h_toa_mean = new TH1D("h_toa_mean", "Mean ToA per Event", 100, toa_ns_min, toa_ns_max);
    for (auto toa_mean : toa_mean_values) {
        h_toa_mean->Fill(toa_mean);
    }
    // do a gaussian fit
    auto toa_mean_mean = h_toa_mean->GetMean();
    auto toa_mean_stddev = h_toa_mean->GetStdDev();
    TF1* fit_func = new TF1("fit_func", "gaus", toa_mean_mean - 2.0 * toa_mean_stddev, toa_mean_mean + 2.0 * toa_mean_stddev);
    h_toa_mean->Fit(fit_func, "RQN");
    fit_func->SetLineColor(kRed);

    TCanvas *c1 = new TCanvas("c1", "Mean ToA per Event", 800, 600);
    h_toa_mean->Draw();
    fit_func->Draw("same");
    c1->Write();

    // find the 90% largest adc sum
    std::vector<double> adc_sum_values_sorted = adc_sum_values;
    std::sort(adc_sum_values_sorted.begin(), adc_sum_values_sorted.end());
    double adc_sum_90pct = adc_sum_values_sorted[static_cast<size_t>(0.9 * adc_sum_values_sorted.size())];
    TH2D* h_adc_sum_toa_mean = new TH2D("h_adc_sum_toa_mean", "Mean ToA vs ADC Sum;ADC Sum;Mean ToA (ns)", 100, 0, adc_sum_90pct*1.3, 256, toa_ns_min, toa_ns_max);

    for (size_t i = 0; i < toa_mean_values.size(); i++) {
        if (adc_sum_values[i] > 0 && toa_mean_values[i] != 0.0) {
            h_adc_sum_toa_mean->Fill(adc_sum_values[i], toa_mean_values[i]);
        }
    }

    TCanvas *c2 = new TCanvas("c2", "Mean ToA vs ADC Sum", 800, 600);
    h_adc_sum_toa_mean->Draw("COLZ");
    c2->Write();

    const int k_adc_sum_bins = 200;
    TCanvas *canvas_adc_sums = new TCanvas("canvas_adc_sums", "ADC Sum Distribution", 800, 600);
    TH1D* h_all_adc_sums = new TH1D("h_all_adc_sums", "ADC Sum Distribution;ADC Sum;Counts", k_adc_sum_bins, 0, adc_sum_90pct*1.3);
    TH1D* h_toa_center_25ns_adc_sums = new TH1D("h_toa_center_25ns_adc_sums", "ADC Sum Distribution for ToA in [-25ns, 25ns];ADC Sum;Counts", k_adc_sum_bins, 0, adc_sum_90pct*1.3);
    TH1D* h_toa_center_10ns_adc_sums = new TH1D("h_toa_center_10ns_adc_sums", "ADC Sum Distribution for ToA in [-10ns, 10ns];ADC Sum;Counts", k_adc_sum_bins, 0, adc_sum_90pct*1.3);
    // remove statistics box
    h_all_adc_sums->SetStats(kFALSE);
    h_toa_center_25ns_adc_sums->SetStats(kFALSE);
    h_toa_center_10ns_adc_sums->SetStats(kFALSE);
    // remove title
    h_all_adc_sums->SetTitle("");
    h_toa_center_25ns_adc_sums->SetTitle("");
    h_toa_center_10ns_adc_sums->SetTitle("");
    size_t counter_all_events = 0;
    size_t counter_toa_25ns = 0;
    size_t counter_toa_10ns = 0;
    for (size_t i = 0; i < toa_mean_values.size(); i++) {
        h_all_adc_sums->Fill(adc_sum_values[i]);
        counter_all_events++;
        if (std::abs(toa_mean_values[i] - fit_func->GetParameter(1)) <= 12.5) {
            h_toa_center_25ns_adc_sums->Fill(adc_sum_values[i]);
            counter_toa_25ns++;
        }
        if (std::abs(toa_mean_values[i] - fit_func->GetParameter(1)) <= 2.5) {
            h_toa_center_10ns_adc_sums->Fill(adc_sum_values[i]);
            counter_toa_10ns++;
        }
    }
    h_all_adc_sums->SetLineColor(kBlack);
    h_toa_center_25ns_adc_sums->SetLineColor(kBlue);
    h_toa_center_10ns_adc_sums->SetLineColor(kRed);
    // normalize all the histograms
    h_all_adc_sums->Scale(1.0 / h_all_adc_sums->Integral());
    h_toa_center_25ns_adc_sums->Scale(1.0 / h_toa_center_25ns_adc_sums->Integral());
    h_toa_center_10ns_adc_sums->Scale(1.0 / h_toa_center_10ns_adc_sums->Integral());
    double max_y = std::max({h_all_adc_sums->GetMaximum(), h_toa_center_25ns_adc_sums->GetMaximum(), h_toa_center_10ns_adc_sums->GetMaximum()});
    h_all_adc_sums->SetMaximum(max_y * 1.3);

    // fit a gaussian to the ADC sum distributions
    double adc_sum_all_mean = h_all_adc_sums->GetMean();
    double adc_sum_all_stddev = h_all_adc_sums->GetStdDev();
    TF1* pre_pre_fit_adc_sum_all = new TF1("pre_pre_fit_adc_sum_all", "gaus", adc_sum_all_mean - 2.0 * adc_sum_all_stddev, adc_sum_all_mean + 2.0 * adc_sum_all_stddev);
    h_all_adc_sums->Fit(pre_pre_fit_adc_sum_all, "RQN");
    TF1* pre_fit_adc_sum_all = new TF1("pre_fit_adc_sum_all", "gaus", pre_pre_fit_adc_sum_all->GetParameter(1) - 1.0 * pre_pre_fit_adc_sum_all->GetParameter(2), pre_pre_fit_adc_sum_all->GetParameter(1) + 1.0 * pre_pre_fit_adc_sum_all->GetParameter(2));
    h_all_adc_sums->Fit(pre_fit_adc_sum_all, "RQN");

    double adc_sum_toa25_mean = h_toa_center_25ns_adc_sums->GetMean();
    double adc_sum_toa25_stddev = h_toa_center_25ns_adc_sums->GetStdDev();
    TF1* pre_pre_fit_adc_sum_toa25 = new TF1("pre_pre_fit_adc_sum_toa25", "gaus", adc_sum_toa25_mean - 2.0 * adc_sum_toa25_stddev, adc_sum_toa25_mean + 2.0 * adc_sum_toa25_stddev);
    h_toa_center_25ns_adc_sums->Fit(pre_pre_fit_adc_sum_toa25, "RQN");
    TF1* pre_fit_adc_sum_toa25 = new TF1("pre_fit_adc_sum_toa25", "gaus", pre_pre_fit_adc_sum_toa25->GetParameter(1) - 1.0 * pre_pre_fit_adc_sum_toa25->GetParameter(2), pre_pre_fit_adc_sum_toa25->GetParameter(1) + 1.0 * pre_pre_fit_adc_sum_toa25->GetParameter(2));
    h_toa_center_25ns_adc_sums->Fit(pre_fit_adc_sum_toa25, "RQN");

    double adc_sum_toa10_mean = h_toa_center_10ns_adc_sums->GetMean();
    double adc_sum_toa10_stddev = h_toa_center_10ns_adc_sums->GetStdDev();
    TF1* pre_pre_fit_adc_sum_toa10 = new TF1("pre_pre_fit_adc_sum_toa10", "gaus", adc_sum_toa10_mean - 2.0 * adc_sum_toa10_stddev, adc_sum_toa10_mean + 2.0 * adc_sum_toa10_stddev);
    h_toa_center_10ns_adc_sums->Fit(pre_pre_fit_adc_sum_toa10, "RQN");
    TF1* pre_fit_adc_sum_toa10 = new TF1("pre_fit_adc_sum_toa10", "gaus", pre_pre_fit_adc_sum_toa10->GetParameter(1) - 1.0 * pre_pre_fit_adc_sum_toa10->GetParameter(2), pre_pre_fit_adc_sum_toa10->GetParameter(1) + 1.0 * pre_pre_fit_adc_sum_toa10->GetParameter(2));
    h_toa_center_10ns_adc_sums->Fit(pre_fit_adc_sum_toa10, "RQN");

    h_all_adc_sums->Draw("hist");
    pre_fit_adc_sum_all->SetLineColor(kBlack);
    // pre_fit_adc_sum_all->Draw("same");
    h_toa_center_25ns_adc_sums->Draw("hist same");
    pre_fit_adc_sum_toa25->SetLineColor(kBlue);
    // pre_fit_adc_sum_toa25->Draw("same");
    h_toa_center_10ns_adc_sums->Draw("hist same");
    pre_fit_adc_sum_toa10->SetLineColor(kRed);
    // pre_fit_adc_sum_toa10->Draw("same");

    // do real fits
    // std::vector<double> fit_sigma_range_left = {0.7, 0.6, 0.5};
    std::vector<double> fit_sigma_range_left  = {4.0, 4.5, 3.5};
    std::vector<double> fit_sigma_range_right = {2.0, 2.25, 2.5};
    std::vector<double> adc_sum_all_fit_means;
    std::vector<double> adc_sum_toa25_fit_means;
    std::vector<double> adc_sum_toa10_fit_means;
    std::vector<double> adc_sum_all_fit_sigmas;
    std::vector<double> adc_sum_toa25_fit_sigmas;
    std::vector<double> adc_sum_toa10_fit_sigmas;
    for (auto left_sigma : fit_sigma_range_left) {
        for (auto right_sigma : fit_sigma_range_right) {
            const double xmin = pre_fit_adc_sum_all->GetParameter(1) - left_sigma  * pre_fit_adc_sum_all->GetParameter(2);
            const double xmax = pre_fit_adc_sum_all->GetParameter(1) + right_sigma * pre_fit_adc_sum_all->GetParameter(2);

            // 1) Tell TF1 there are 5 parameters
            // 2) Give each TF1 a unique name, or don't add to global list
            auto fname = Form("fit_adc_sum_all_%zu_%zu", static_cast<size_t>(left_sigma * 10), static_cast<size_t>(right_sigma * 10));
            TF1* fit_adc_sum_all = new TF1(fname, CrystalBallFcn, xmin, xmax, /*npar=*/5);

            fit_adc_sum_all->SetParNames("N","alpha","n","mean","sigma");
            fit_adc_sum_all->SetParameters(
                pre_fit_adc_sum_all->GetParameter(0), // N
                1.5,                                  // alpha
                3.0,                                  // n
                pre_fit_adc_sum_all->GetParameter(1), // mean
                pre_fit_adc_sum_all->GetParameter(2)  // sigma
            );
            // fit_adc_sum_all->SetParLimits(0, 0.0, 1e15);
            // fit_adc_sum_all->SetParLimits(2, 1.0, 80.0);
            // fit_adc_sum_all->SetParLimits(4, 1e-6, 1e6);
            fit_adc_sum_all->SetParLimits(0, pre_fit_adc_sum_all->GetParameter(0) * 0.9, pre_fit_adc_sum_all->GetParameter(0) * 1.2);
            fit_adc_sum_all->SetParLimits(3, pre_fit_adc_sum_all->GetParameter(1) - pre_fit_adc_sum_all->GetParameter(2) * 0.1,
                                          pre_fit_adc_sum_all->GetParameter(1) + pre_fit_adc_sum_all->GetParameter(2) * 0.1);

            // R: respect range, Q: quiet, N: do NOT store fit function with the histogram, +: add result to list
            h_all_adc_sums->Fit(fit_adc_sum_all, "RQN+");

            adc_sum_all_fit_means .push_back(fit_adc_sum_all->GetParameter(3));
            adc_sum_all_fit_sigmas.push_back(fit_adc_sum_all->GetParameter(4));

            fit_adc_sum_all->SetLineColorAlpha(kGray, 0.2);
            fit_adc_sum_all->Draw("same");

            const double xmin25 = pre_fit_adc_sum_toa25->GetParameter(1) - left_sigma  * pre_fit_adc_sum_toa25->GetParameter(2);
            const double xmax25 = pre_fit_adc_sum_toa25->GetParameter(1) + right_sigma * pre_fit_adc_sum_toa25->GetParameter(2);
            auto fname25 = Form("fit_adc_sum_toa25_%zu_%zu", static_cast<size_t>(left_sigma * 10), static_cast<size_t>(right_sigma * 10));
            TF1* fit_adc_sum_toa25 = new TF1(fname25, CrystalBallFcn, xmin25, xmax25, /*npar=*/5);
            fit_adc_sum_toa25->SetParNames("N","alpha","n","mean","sigma");
            fit_adc_sum_toa25->SetParameters(
                pre_fit_adc_sum_toa25->GetParameter(0), // N
                1.5,                                   // alpha
                3.0,                                   // n
                pre_fit_adc_sum_toa25->GetParameter(1), // mean
                pre_fit_adc_sum_toa25->GetParameter(2)  // sigma
            );
            fit_adc_sum_toa25->SetParLimits(0, pre_fit_adc_sum_toa25->GetParameter(0) * 0.9, pre_fit_adc_sum_toa25->GetParameter(0) * 1.2);
            fit_adc_sum_toa25->SetParLimits(3, pre_fit_adc_sum_toa25->GetParameter(1) - pre_fit_adc_sum_toa25->GetParameter(2) * 0.1,
                                            pre_fit_adc_sum_toa25->GetParameter(1) + pre_fit_adc_sum_toa25->GetParameter(2) * 0.1);
            h_toa_center_25ns_adc_sums->Fit(fit_adc_sum_toa25, "RQN+");
            adc_sum_toa25_fit_means .push_back(fit_adc_sum_toa25->GetParameter(3));
            adc_sum_toa25_fit_sigmas.push_back(fit_adc_sum_toa25->GetParameter(4));
            fit_adc_sum_toa25->SetLineColorAlpha(kBlue-10, 0.2);
            fit_adc_sum_toa25->Draw("same");
            
            const double xmin10 = pre_fit_adc_sum_toa10->GetParameter(1) - left_sigma  * pre_fit_adc_sum_toa10->GetParameter(2);;
            const double xmax10 = pre_fit_adc_sum_toa10->GetParameter(1) + right_sigma * pre_fit_adc_sum_toa10->GetParameter(2);;
            auto fname10 = Form("fit_adc_sum_toa10_%zu_%zu", static_cast<size_t>(left_sigma * 10), static_cast<size_t>(right_sigma * 10));
            TF1* fit_adc_sum_toa10 = new TF1(fname10, CrystalBallFcn, xmin10, xmax10, /*npar=*/5);
            fit_adc_sum_toa10->SetParNames("N","alpha","n","mean","sigma");
            fit_adc_sum_toa10->SetParameters(
                pre_fit_adc_sum_toa10->GetParameter(0), // N
                1.5,                                   // alpha
                3.0,                                   // n
                pre_fit_adc_sum_toa10->GetParameter(1), // mean
                pre_fit_adc_sum_toa10->GetParameter(2)  // sigma
            );
            fit_adc_sum_toa10->SetParLimits(0, pre_fit_adc_sum_toa10->GetParameter(0) * 0.9, pre_fit_adc_sum_toa10->GetParameter(0) * 1.2);
            fit_adc_sum_toa10->SetParLimits(3, pre_fit_adc_sum_toa10->GetParameter(1) - pre_fit_adc_sum_toa10->GetParameter(2) * 0.1,
                                            pre_fit_adc_sum_toa10->GetParameter(1) + pre_fit_adc_sum_toa10->GetParameter(2) * 0.1);
            h_toa_center_10ns_adc_sums->Fit(fit_adc_sum_toa10, "RQN+");
            adc_sum_toa10_fit_means .push_back(fit_adc_sum_toa10->GetParameter(3));
            adc_sum_toa10_fit_sigmas.push_back(fit_adc_sum_toa10->GetParameter(4));
            fit_adc_sum_toa10->SetLineColorAlpha(kRed-10, 0.2);
            fit_adc_sum_toa10->Draw("same");
        }
    }

    double adc_sum_all_fit_mean_avg = std::accumulate(adc_sum_all_fit_means.begin(), adc_sum_all_fit_means.end(), 0.0) / adc_sum_all_fit_means.size();
    double adc_sum_toa25_fit_mean_avg = std::accumulate(adc_sum_toa25_fit_means.begin(), adc_sum_toa25_fit_means.end(), 0.0) / adc_sum_toa25_fit_means.size();
    double adc_sum_toa10_fit_mean_avg = std::accumulate(adc_sum_toa10_fit_means.begin(), adc_sum_toa10_fit_means.end(), 0.0) / adc_sum_toa10_fit_means.size();
    double adc_sum_all_fit_sigma_avg = std::accumulate(adc_sum_all_fit_sigmas.begin(), adc_sum_all_fit_sigmas.end(), 0.0) / adc_sum_all_fit_sigmas.size();
    double adc_sum_toa25_fit_sigma_avg = std::accumulate(adc_sum_toa25_fit_sigmas.begin(), adc_sum_toa25_fit_sigmas.end(), 0.0) / adc_sum_toa25_fit_sigmas.size();
    double adc_sum_toa10_fit_sigma_avg = std::accumulate(adc_sum_toa10_fit_sigmas.begin(), adc_sum_toa10_fit_sigmas.end(), 0.0) / adc_sum_toa10_fit_sigmas.size();
    double adc_sum_all_resolution = adc_sum_all_fit_sigma_avg / adc_sum_all_fit_mean_avg * 100.0;
    double adc_sum_toa25_resolution = adc_sum_toa25_fit_sigma_avg / adc_sum_toa25_fit_mean_avg * 100.0;
    double adc_sum_toa10_resolution = adc_sum_toa10_fit_sigma_avg / adc_sum_toa10_fit_mean_avg * 100.0;
    
    // top left legend
    TLegend *legend = new TLegend(0.15, 0.65, 0.45, 0.90);
    legend->AddEntry(h_all_adc_sums, Form("All Events (%zu)", counter_all_events), "l");
    legend->AddEntry((TObject*)0, Form("All Events Fit Mean = %.2f", adc_sum_all_fit_mean_avg), "");
    legend->AddEntry((TObject*)0, Form("All Events Resolution = %.2f %%", adc_sum_all_resolution), "");

    legend->AddEntry(h_toa_center_25ns_adc_sums, Form("ToA in [-12.5ns, 12.5ns] (%zu)", counter_toa_25ns), "l");
    legend->AddEntry((TObject*)0, Form("ToA in [-12.5ns, 12.5ns] Fit Mean = %.2f", adc_sum_toa25_fit_mean_avg), "");
    legend->AddEntry((TObject*)0, Form("ToA in [-12.5ns, 12.5ns] Resolution = %.2f %%", adc_sum_toa25_resolution), "");

    legend->AddEntry(h_toa_center_10ns_adc_sums, Form("ToA in [-2.5ns, 2.5ns] (%zu)", counter_toa_10ns), "l");
    legend->AddEntry((TObject*)0, Form("ToA in [-2.5ns, 2.5ns] Fit Mean = %.2f", adc_sum_toa10_fit_mean_avg), "");
    legend->AddEntry((TObject*)0, Form("ToA in [-2.5ns, 2.5ns] Resolution = %.2f %%", adc_sum_toa10_resolution), "");

    legend->Draw();
    canvas_adc_sums->Write();
    // save as pdf
    canvas_adc_sums->SaveAs((script_output_folder + "/" + script_name + "Run" + script_input_run_number + "_adc_sum_distribution.pdf").c_str());


    output_root->Close();
    return 0;
}
