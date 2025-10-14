#include <fstream>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"

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

    TH2D* pedestal_distribution_th2d = new TH2D("pedestal_distribution", "Pedestal Distribution;Channel;ADC Value", FPGA_CHANNEL_NUMBER_VALID*vldb_number, 0, FPGA_CHANNEL_NUMBER_VALID*vldb_number, 512, 0, 512);
    pedestal_distribution_th2d->SetStats(0);
    pedestal_distribution_th2d->SetTitle("");

    TH2D* peak_distribution_th2d = new TH2D("peak_distribution", "Peak Distribution;Channel;ADC Value", FPGA_CHANNEL_NUMBER_VALID*vldb_number, 0, FPGA_CHANNEL_NUMBER_VALID*vldb_number, 512, 0, 1024);
    peak_distribution_th2d->SetStats(0);
    peak_distribution_th2d->SetTitle("");

    // if events are more than 10k, make 10 histograms each for 1/10 of the events
    // if events are more than 1k, make 5 histograms each for 1/5 of the events
    // other wise make 2 histograms
    int n_partial_hist = 2;
    if (entry_max > 100000) {
        n_partial_hist = 10;
    } else if (entry_max > 10000) {
        n_partial_hist = 5;
    } else {
        n_partial_hist = 2;
    }

    std::vector<double> event_adc_sum_list;
    event_adc_sum_list.reserve(entry_max);
    std::vector<double> event_adc_sum_list_copy;
    event_adc_sum_list.reserve(entry_max);
    // start event loop
    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        double adc_sum_all = 0;
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            // channel loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                int channel_index_in_h2g = channel % 76;
                int channel_index_in_h2g_half = channel_index_in_h2g % 38;
                int asic_id = channel / 76;
                int half_id = channel_index_in_h2g / 38;

                if (vldb_id == 0 && asic_id == 1)
                    continue;
                
                if (channel_index_in_h2g_half == 19)
                    continue;
                if (channel_index_in_h2g_half == 0)
                    continue;

                int channel_index_no_CM_Calib = channel_index_in_h2g_half;
                
                if (channel_index_in_h2g_half > 19)
                    channel_index_no_CM_Calib -= 2;
                else if (channel_index_in_h2g_half > 0)
                    channel_index_no_CM_Calib -= 1;

                int channel_number_valid = asic_id * 72 + half_id * 36 + (channel_index_no_CM_Calib);

                // pedestal to be the average of the smallest 2 of the first 4 samples
                double adc_pedestal = 0;
                std::vector<int> adc_samples;
                double adc_peak = 0;
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    auto &adc = val0_list_pools[vldb_id][0][channel + sample * FPGA_CHANNEL_NUMBER];
                    adc_samples.push_back(adc);
                    auto adc_double = static_cast<double>(adc);
                    if (adc_double > adc_peak) {
                        adc_peak = adc_double;
                    }
                }
                pedestal_subtraction_2minoffirst3(adc_samples, adc_pedestal);

                double adc_peak_subtracted = adc_peak - adc_pedestal;
                if (adc_peak_subtracted < 0) {
                    adc_peak_subtracted = 0;
                }

                pedestal_distribution_th2d->Fill(vldb_id * FPGA_CHANNEL_NUMBER_VALID + channel_number_valid, adc_pedestal);
                peak_distribution_th2d->Fill(vldb_id * FPGA_CHANNEL_NUMBER_VALID + channel_number_valid, adc_peak_subtracted);

                adc_sum_all += adc_peak_subtracted;
            } // end of channel loop
        } // end of vldb loop
        event_adc_sum_list.push_back(adc_sum_all);
        event_adc_sum_list_copy.push_back(adc_sum_all);
    } // end of event loop

    // calculate the 90% max of all adc sums
    // sort the adc sums
    if (event_adc_sum_list.empty()) {
        spdlog::warn("No events processed, using default ADC sum range.");
        event_adc_sum_list.push_back(0.0);
    }
    std::sort(event_adc_sum_list.begin(), event_adc_sum_list.end());
    // get the 90% max

    double adc_sum_max = event_adc_sum_list[std::min(static_cast<size_t>(event_adc_sum_list.size() * 0.9), event_adc_sum_list.size() - 1)];
    adc_sum_max = adc_sum_max * 1.3; // add 30% margin
    // make 10 histograms each for 1/10 of the events
    std::vector<TH1D*> run_adc_sum_partial_hist1d_list;

    spdlog::info("ADC sum max: {}", adc_sum_max);

    TH1D* run_adc_sum_all_hist1d = new TH1D("run_adc_sum_all", "Run ADC Sum All;ADC Sum;Counts", 256, 0, adc_sum_max);
    run_adc_sum_all_hist1d->SetStats(0);
    run_adc_sum_all_hist1d->SetTitle("");

    for (int i = 0; i < n_partial_hist; i++) {
        auto hist = new TH1D(("run_adc_sum_partial_" + std::to_string(i)).c_str(), ("Run ADC Sum Partial " + std::to_string(i) + ";""ADC Sum;Counts").c_str(), 256, 0, adc_sum_max);
        hist->SetStats(0);
        hist->SetTitle("");
        run_adc_sum_partial_hist1d_list.push_back(hist);
    }

    for (size_t i = 0; i < event_adc_sum_list.size(); i++) {
        auto& adc_sum = event_adc_sum_list_copy[i];
        run_adc_sum_all_hist1d->Fill(adc_sum);
        int partial_hist_index = i * n_partial_hist / event_adc_sum_list.size();
        run_adc_sum_partial_hist1d_list[partial_hist_index]->Fill(adc_sum);
    }

    double pre_fit_mean = run_adc_sum_all_hist1d->GetMean();
    double pre_fit_sigma = run_adc_sum_all_hist1d->GetRMS();
    double fit_min = std::max(0.0, pre_fit_mean - 2 * pre_fit_sigma);
    double fit_max = std::min(adc_sum_max, pre_fit_mean + 2 * pre_fit_sigma);

    TF1* gaus_fit_pre = new TF1("gaus_fit_pre", "gaus", fit_min, fit_max);
    gaus_fit_pre->SetLineColorAlpha(kWhite, 0.0);
    run_adc_sum_all_hist1d->Fit(gaus_fit_pre, "RQ");

    double fit_mean = gaus_fit_pre->GetParameter(1);
    double fit_sigma = gaus_fit_pre->GetParameter(2);

    // std::vector<double> fit_sigmas = {2.0, 2.5, 3.0};
    // std::vector<double> fit_offsets = {0.2, 0.4, 0.6};
    std::vector<double> fit_sigma_range_left = {0.3, 0.4, 0.5};
    std::vector<double> fit_sigma_range_right = {2.0, 2.5, 3.0};
    std::vector<double> fit_results_means;
    std::vector<double> fit_results_sigmas;
    // for (auto fit_sigma_multiplier : fit_sigmas) {
    //     for (auto fit_offset_multiplier : fit_offsets) {
    //         double fit_min_final = std::max(0.0, fit_mean - fit_sigma_multiplier * fit_sigma + fit_offset_multiplier * fit_sigma);
    //         double fit_max_final = std::min(adc_sum_max, fit_mean + fit_sigma_multiplier * fit_sigma + fit_offset_multiplier * fit_sigma);
    for (auto fit_sigma_multiplier : fit_sigma_range_left) {
        for (auto fit_offset_multiplier : fit_sigma_range_right) {
            double fit_min_final = std::max(0.0, fit_mean - fit_sigma_multiplier * fit_sigma);
            double fit_max_final = std::min(adc_sum_max, fit_mean + fit_offset_multiplier * fit_sigma);
            TF1* gaus_fit_final = new TF1(("gaus_fit_final_" + std::to_string(static_cast<int>(fit_sigma_multiplier * 10)) + "_" + std::to_string(static_cast<int>(fit_offset_multiplier * 10))).c_str(), "gaus", fit_min_final, fit_max_final);
            // set transparency
            gaus_fit_final->SetLineColorAlpha(kPink, 0.2);
            run_adc_sum_all_hist1d->Fit(gaus_fit_final, "RQ+");
            double final_fit_mean = gaus_fit_final->GetParameter(1);
            double final_fit_sigma = gaus_fit_final->GetParameter(2);
            fit_results_means.push_back(final_fit_mean);
            fit_results_sigmas.push_back(final_fit_sigma);
            spdlog::info("Final Fit with sigma range {} and offset {}: mean = {}, sigma = {}", fit_sigma_multiplier, fit_offset_multiplier, final_fit_mean, final_fit_sigma);
        }
    }

    // calculate the average of the fit results
    double final_mean = 0.0;
    double final_sigma = 0.0;
    double final_mean_error = 0.0;
    double final_sigma_error = 0.0;
    if (!fit_results_means.empty()) {
        final_mean = std::accumulate(fit_results_means.begin(), fit_results_means.end(), 0.0) / fit_results_means.size();
        final_sigma = std::accumulate(fit_results_sigmas.begin(), fit_results_sigmas.end(), 0.0) / fit_results_sigmas.size();
    } else {
        final_mean = fit_mean;
        final_sigma = fit_sigma;
    }
    if (fit_results_means.size() > 1) {
        double mean_variance = 0.0;
        double sigma_variance = 0.0;
        for (size_t i = 0; i < fit_results_means.size(); i++) {
            mean_variance += (fit_results_means[i] - final_mean) * (fit_results_means[i] - final_mean);
            sigma_variance += (fit_results_sigmas[i] - final_sigma) * (fit_results_sigmas[i] - final_sigma);
        }
        mean_variance /= (fit_results_means.size() - 1);
        sigma_variance /= (fit_results_sigmas.size() - 1);
        final_mean_error = std::sqrt(mean_variance / fit_results_means.size());
        final_sigma_error = std::sqrt(sigma_variance / fit_results_sigmas.size());
    } else {
        final_mean_error = gaus_fit_pre->GetParError(1);
        final_sigma_error = gaus_fit_pre->GetParError(2);
    }

    auto output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    output_root->cd();

    TCanvas *pedestal_distribution_th2d_canvas = new TCanvas("pedestal_distribution_canvas", "pedestal_distribution_canvas", 1200, 600);
    pedestal_distribution_th2d->Draw("COLZ");
    pedestal_distribution_th2d_canvas->SetLogz();

    TLatex latex_pedestal;
    latex_pedestal.SetTextColor(kGray+2);
    latex_pedestal.SetNDC();
    latex_pedestal.SetTextSize(0.05);
    latex_pedestal.SetTextFont(62);
    latex_pedestal.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());
    latex_pedestal.SetTextSize(0.035);
    latex_pedestal.SetTextFont(42);
    latex_pedestal.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
    // write run number, date time
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream date_stream;
    date_stream << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    latex_pedestal.DrawLatex(0.12, 0.76, (std::string("Run ") + script_input_run_number).c_str());
    latex_pedestal.DrawLatex(0.12, 0.72, "Pedestal Distribution");
    latex_pedestal.DrawLatex(0.12, 0.68, date_stream.str().c_str());    
    pedestal_distribution_th2d_canvas->Write();
    pedestal_distribution_th2d_canvas->Close();

    TCanvas *peak_distribution_th2d_canvas = new TCanvas("peak_distribution_canvas", "peak_distribution_canvas", 1200, 600);
    peak_distribution_th2d->Draw("COLZ");
    peak_distribution_th2d_canvas->SetLogz();
    TLatex latex_peak;
    latex_peak.SetTextColor(kGray+2);
    latex_peak.SetNDC();
    latex_peak.SetTextSize(0.05);
    latex_peak.SetTextFont(62);
    latex_peak.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());
    latex_peak.SetTextSize(0.035);
    latex_peak.SetTextFont(42);
    latex_peak.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
    // write run number, date time
    latex_peak.DrawLatex(0.12, 0.76, (std::string("Run ") + script_input_run_number).c_str());
    latex_peak.DrawLatex(0.12, 0.72, "Peak Distribution");
    latex_peak.DrawLatex(0.12, 0.68, date_stream.str().c_str());    
    peak_distribution_th2d_canvas->Write();
    peak_distribution_th2d_canvas->Close();

    TCanvas *run_adc_sum_canvas = new TCanvas("run_adc_sum_all_canvas", "run_adc_sum_all_canvas", 800, 600);
    // draw all the partial histograms in one canvas
    run_adc_sum_all_hist1d->SetLineColor(kBlack);
    run_adc_sum_all_hist1d->SetLineWidth(2);
    run_adc_sum_all_hist1d->Draw("HIST");
    TLegend *run_adc_sum_legend = new TLegend(0.80, 0.65, 0.89, 0.89);
    run_adc_sum_legend->SetBorderSize(0);
    run_adc_sum_legend->SetFillStyle(0);
    run_adc_sum_legend->SetTextFont(42);
    run_adc_sum_legend->SetTextSize(0.02);
    run_adc_sum_legend->AddEntry(run_adc_sum_all_hist1d, "All Events", "l");

    for (int i = 0; i < n_partial_hist; i++) {
        auto hist = run_adc_sum_partial_hist1d_list[i];
        hist->SetLineColor(kBlue + i * 2);
        hist->SetLineWidth(1);
        hist->Draw("HIST SAME");
        run_adc_sum_legend->AddEntry(hist, ("Partial " + std::to_string(i)).c_str(), "l");
    }
    run_adc_sum_legend->Draw();
    TLatex latex_fit;
    latex_fit.SetTextColor(kRed);
    latex_fit.SetNDC();
    latex_fit.SetTextSize(0.03);
    latex_fit.SetTextFont(42);
    latex_fit.DrawLatex(0.60, 0.60, ("Fit Mean: " + std::to_string(final_mean)).c_str());
    latex_fit.DrawLatex(
        0.60, 0.55, 
        ("Fit Sigma: " + std::to_string(final_sigma) + 
        " (" + std::to_string(static_cast<int>(final_sigma / final_mean * 10000) / 100.0) + "%)").c_str()
    );
    // write run number
    TLatex latex_run_adc_sum;
    latex_run_adc_sum.SetTextColor(kGray+2);
    latex_run_adc_sum.SetNDC();
    latex_run_adc_sum.SetTextSize(0.05);
    latex_run_adc_sum.SetTextFont(62);
    latex_run_adc_sum.DrawLatex(0.12, 0.85, annotation_canvas_title.c_str());
    latex_run_adc_sum.SetTextSize(0.035);
    latex_run_adc_sum.SetTextFont(42);
    latex_run_adc_sum.DrawLatex(0.12, 0.80, annotation_testbeam_title.c_str());
    // write run number, date time
    latex_run_adc_sum.DrawLatex(0.12, 0.76, (std::string("Run ") + script_input_run_number).c_str());
    latex_run_adc_sum.DrawLatex(0.12, 0.72, "Run ADC Sum");
    latex_run_adc_sum.DrawLatex(0.12, 0.68, date_stream.str().c_str());    
    run_adc_sum_canvas->Write();
    run_adc_sum_canvas->Close();

    output_root->Close();

    input_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}