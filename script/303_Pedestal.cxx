#include <fstream>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "TCanvas.h"
#include "TParameter.h"

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
    pedestal_distribution_th2d->SetDirectory(nullptr);

    TH2D* pedestal_distribution_average_th2d = new TH2D("pedestal_distribution_average", "Pedestal Distribution Average;Channel;ADC Value", FPGA_CHANNEL_NUMBER_VALID*vldb_number, 0, FPGA_CHANNEL_NUMBER_VALID*vldb_number, 512, 0, 512);
    pedestal_distribution_average_th2d->SetStats(0);
    pedestal_distribution_average_th2d->SetTitle("");
    pedestal_distribution_average_th2d->SetDirectory(nullptr);

    std::vector<double> event_adc_pedestal_list;
    event_adc_pedestal_list.reserve(entry_max);
    std::vector<double> event_adc_pedestal_mean_list;
    event_adc_pedestal_mean_list.reserve(entry_max);
    // start event loop
    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        double adc_pedestal_sum = 0.0;
        double adc_pedestal_mean_sum = 0.0;
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            // channel loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                std::vector<UInt_t> adc_pedestal_samples; // only take the first 3 samples
                adc_pedestal_samples.reserve(3);
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    UInt_t adc_value = val0_list_pools[vldb_id][0][idx];
                    UInt_t tot_value = val1_list_pools[vldb_id][0][idx];
                    UInt_t toa_value = val2_list_pools[vldb_id][0][idx];
                    if (sample < 3) {
                        adc_pedestal_samples.push_back(adc_value);
                    }
                } // end of sample loop
                // calculate the pedestal
                double adc_pedestal = pedestal_median_of_first3(adc_pedestal_samples);
                double adc_pedestal_mean = pedestal_average_of_first3(adc_pedestal_samples);
                pedestal_distribution_th2d->Fill(double(vldb_id)*double(FPGA_CHANNEL_NUMBER_VALID) + double(channel), adc_pedestal);
                pedestal_distribution_average_th2d->Fill(double(vldb_id)*double(FPGA_CHANNEL_NUMBER_VALID) + double(channel), adc_pedestal_mean);
                adc_pedestal_sum += adc_pedestal;
                adc_pedestal_mean_sum += adc_pedestal_mean;
            } // end of channel loop
        } // end of vldb loop
        event_adc_pedestal_list.push_back(adc_pedestal_sum);
        event_adc_pedestal_mean_list.push_back(adc_pedestal_mean_sum);
    } // end of event loop
    input_root->Close();

    spdlog::info("Finished processing {} events", entry_max);


    // * === Plot to output file =============================================================
    // * =====================================================================================
    auto output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    output_root->cd();

    auto pedestal_distribution_canvas = new TCanvas("pedestal_distribution_canvas", "Pedestal Distribution Canvas", 800, 600);
    pedestal_distribution_canvas->cd();
    pedestal_distribution_canvas->SetLogz();
    pedestal_distribution_th2d->Draw("COLZ");
    // pedestal_distribution_canvas->SetLogz();
    pedestal_distribution_canvas->Write();
    pedestal_distribution_canvas->Close();

    auto pedestal_distribution_average_canvas = new TCanvas("pedestal_distribution_average_canvas", "Pedestal Distribution Average Canvas", 800, 600);
    pedestal_distribution_average_canvas->cd();
    pedestal_distribution_average_canvas->SetLogz();
    pedestal_distribution_average_th2d->Draw("COLZ");
    // pedestal_distribution_average_canvas->SetLogz();
    pedestal_distribution_average_canvas->Write();
    pedestal_distribution_average_canvas->Close();

    TCanvas *pedestal_distribution_1d_canvas = new TCanvas("pedestal_distribution_1d_canvas", "Pedestal Distribution 1D Canvas", 800, 600);
    pedestal_distribution_1d_canvas->cd();
    // calculate the 90% max and 10% min for y axis
    auto sorted_event_adc_pedestal_list = event_adc_pedestal_list;
    std::sort(sorted_event_adc_pedestal_list.begin(), sorted_event_adc_pedestal_list.end());
    double y_min = sorted_event_adc_pedestal_list[static_cast<size_t>(0.1 * sorted_event_adc_pedestal_list.size())];
    double y_max = sorted_event_adc_pedestal_list[static_cast<size_t>(0.9 * sorted_event_adc_pedestal_list.size())];
    spdlog::info("Pedestal 10% min: {}, 90% max: {}", y_min, y_max);
    TH1D *pedestal_distribution_th1d = new TH1D("pedestal_distribution_1d", "Pedestal Distribution 1D", 200, y_min - 2.0*(y_max - y_min), y_max + 2.0*(y_max - y_min));
    pedestal_distribution_th1d->SetStats(0);
    pedestal_distribution_th1d->SetTitle("");
    pedestal_distribution_th1d->GetXaxis()->SetTitle("Sum of ADC Pedestal over all channels");
    pedestal_distribution_th1d->GetYaxis()->SetTitle("Entries");

    for (const auto& adc_pedestal : event_adc_pedestal_list) {
        pedestal_distribution_th1d->Fill(adc_pedestal);
    }
    pedestal_distribution_th1d->Draw();
    // do gaussian fit
    double fit_mean = pedestal_distribution_th1d->GetMean();
    double fit_sigma = pedestal_distribution_th1d->GetRMS();
    TF1 *gaus_fit = new TF1("gaus_fit", "gaus", fit_mean - 4.0*fit_sigma, fit_mean + 4.0*fit_sigma);
    pedestal_distribution_th1d->Fit(gaus_fit, "RQ");
    gaus_fit->SetLineColor(kRed);
    gaus_fit->Draw("same");
    pedestal_distribution_1d_canvas->Update();
    // write fit results
    double pedestal_fit_mean = gaus_fit->GetParameter(1);
    double pedestal_fit_sigma = gaus_fit->GetParameter(2);
    double pedestal_fit_mean_error = gaus_fit->GetParError(1);
    double pedestal_fit_sigma_error = gaus_fit->GetParError(2);
    TLatex fit_results_text;
    fit_results_text.SetNDC();
    fit_results_text.SetTextSize(0.03);
    fit_results_text.SetTextColor(kBlack);
    fit_results_text.DrawLatex(0.15, 0.85, Form("#mu = %.2f #pm %.2f", pedestal_fit_mean, pedestal_fit_mean_error));
    fit_results_text.DrawLatex(0.15, 0.80, Form("#sigma = %.2f #pm %.2f", pedestal_fit_sigma, pedestal_fit_sigma_error));

    // pedestal_distribution_1d_canvas->SetLogy();
    pedestal_distribution_1d_canvas->Write();
    pedestal_distribution_1d_canvas->Close();

    // write the fitting results to TParameter
    TParameter<double> pedestal_fit_mean_param("pedestal_fit_mean", pedestal_fit_mean);
    TParameter<double> pedestal_fit_sigma_param("pedestal_fit_sigma", pedestal_fit_sigma);
    pedestal_fit_mean_param.Write();
    pedestal_fit_sigma_param.Write();

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}