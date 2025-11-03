#include <fstream>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "TRandom3.h"

double quantile_from_tf1(TF1* f, double q, double xmin, double xmax)
{
    const double I_tot = f->Integral(xmin, xmax);
    double lo = xmin, hi = xmax;
    for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (lo + hi);
        double I   = f->Integral(xmin, mid);
        if (I / I_tot < q)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

struct MeanRMS90 {
    // central values
    double mean90;      
    double rms90;       
    double sigma_gauss; // = rms90 / k90
    double res_sigma_over_mean90; // (sigma_gauss / mean90)

    // Gaussian-approx stat errors
    double err_mean90;  
    double err_rms90;   
    double err_sigma_gauss;
    double err_res_sigma_over_mean90;

    // window info
    int    binLo, binHi;
    int    Neff;        // effective events inside the 90% window
};

// 如果 TH1 是“计数直方图”：useBinWidth=false（默认）；
// 若 TH1 是“密度直方图”（bin 已除以宽度），可把 useBinWidth=true 用于取均值/RMS，
// 但此时请用 nEventsOverride 指定总事件数，否则统计误差近似将不可靠。
MeanRMS90 computeMeanRMS90(const TH1* h, 
                           double frac = 0.90, 
                           bool useBinWidth = false,
                           long long nEventsOverride = -1)
{
    const int nb = h->GetNbinsX();
    MeanRMS90 out{};
    if (nb <= 0) { out.mean90 = out.rms90 = out.sigma_gauss = out.res_sigma_over_mean90 = NAN; return out; }

    // 预取边界、中心、权重（用于均值/RMS 计算）
    std::vector<double> edge(nb+1), x(nb+1), w(nb+1);
    for (int i=1; i<=nb; ++i) {
        edge[i-1] = h->GetXaxis()->GetBinLowEdge(i);
        x[i]      = h->GetXaxis()->GetBinCenter(i);
        double bw = h->GetXaxis()->GetBinWidth(i);
        double cnt= h->GetBinContent(i);
        w[i]      = useBinWidth ? (cnt * bw) : cnt; // 仅用于均值/RMS 计算
    }
    edge[nb] = h->GetXaxis()->GetBinUpEdge(nb);

    // 总计权重（若是计数直方图，就是总事件数）
    double totalW = 0.0;
    for (int i=1; i<=nb; ++i) totalW += w[i];
    if (totalW <= 0) { out.mean90 = out.rms90 = out.sigma_gauss = out.res_sigma_over_mean90 = NAN; return out; }

    const double targetW = frac * totalW;

    // 前缀和
    std::vector<double> csum(nb+1, 0.0);
    for (int i=1; i<=nb; ++i) csum[i] = csum[i-1] + w[i];

    // 双指针找“包含 frac 的最窄区间”
    int bestL = 1, bestR = nb;
    double bestWidth = edge[nb] - edge[0];

    int R = 1;
    for (int L = 1; L <= nb; ++L) {
        if (R < L) R = L;
        while (R <= nb && (csum[R] - csum[L-1]) < targetW) ++R;
        if (R > nb) break;
        double width = edge[R] - edge[L-1];
        if (width < bestWidth) { bestWidth = width; bestL = L; bestR = R; }
    }

    // 计算 Mean90 / RMS90
    double W = 0.0, Sx = 0.0, Sx2 = 0.0;
    for (int i=bestL; i<=bestR; ++i) {
        W   += w[i];
        Sx  += w[i] * x[i];
        Sx2 += w[i] * x[i] * x[i];
    }
    double mean90 = (W > 0) ? (Sx / W) : NAN;
    double var90  = (W > 0) ? (Sx2 / W - mean90*mean90) : NAN;
    double rms90  = (var90 > 0) ? std::sqrt(var90) : 0.0;

    // 高斯等效
    constexpr double k90 = 0.7893132990990863; // precise coefficient for central 90%
    double sigma_gauss = rms90 / k90;
    double res = (!std::isnan(mean90) && mean90 != 0.0) ? (sigma_gauss / mean90) : NAN;

    // -------- 高斯近似统计误差 ----------
    // 需要 N_eff：90% 区间内的“事件数”。对于计数直方图，W≈N_eff。
    long long Neff = -1;
    if (!useBinWidth && nEventsOverride < 0) {
        // 计数直方图：直接用该窗口内的计数总和
        double sumCounts = 0.0;
        for (int i=bestL; i<=bestR; ++i) sumCounts += h->GetBinContent(i);
        Neff = static_cast<long long>(std::llround(sumCounts));
    } else if (nEventsOverride >= 0) {
        // 提供了总事件数：用 0.9*N 作为 Neff 的近似
        Neff = static_cast<long long>(std::llround(frac * nEventsOverride));
    } else {
        // 无法确定 Neff：给出 NAN 提醒
        Neff = -1;
    }

    double err_rms90 = NAN, err_mean90 = NAN, err_sigma = NAN, err_res = NAN;
    if (Neff > 1 && std::isfinite(rms90)) {
        // ΔRMS90/RMS90 ≈ 1/sqrt(2*(Neff-1))
        double rel = 1.0 / std::sqrt(2.0 * (double)(Neff - 1));
        err_rms90  = rms90 * rel;

        // ΔMean90 ≈ RMS90 / sqrt(Neff)
        err_mean90 = rms90 / std::sqrt((double)Neff);

        // Δσ_Gauss = ΔRMS90 / k90
        err_sigma  = err_rms90 / k90;

        // Δ(R=σ/Mean) 误差传播（忽略两者相关性的近似）
        if (std::isfinite(res) && mean90 != 0.0 && sigma_gauss != 0.0) {
            double rel_sigma = err_sigma / sigma_gauss;
            double rel_mean  = err_mean90 / mean90;
            err_res = res * std::sqrt(rel_sigma*rel_sigma + rel_mean*rel_mean);
        }
    }

    out.mean90  = mean90;
    out.rms90   = rms90;
    out.sigma_gauss = sigma_gauss;
    out.res_sigma_over_mean90 = res;

    out.err_mean90 = err_mean90;
    out.err_rms90  = err_rms90;
    out.err_sigma_gauss = err_sigma;
    out.err_res_sigma_over_mean90 = err_res;

    out.binLo = bestL;
    out.binHi = bestR;
    out.Neff  = (int)Neff;
    return out;
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
    const int machine_gun_samples = 16;
    const int vldb_number = 2;

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

    // * === Define output data structures ===================================================
    // * =====================================================================================
    // * --- Accumulate waveform for each channel ---
    std::vector<TH2D*> h2_adc_waveforms(vldb_number * FPGA_CHANNEL_NUMBER, nullptr);
    for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
        for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
            std::string h2_name = "h2_adc_waveform_vldb" + std::to_string(vldb_id) + "_ch" + std::to_string(channel);
            std::string h2_title = "VLDB " + std::to_string(vldb_id) + " Channel " + std::to_string(channel) + " ADC Waveform;Sample Index;ADC Value";
            h2_adc_waveforms[vldb_id * FPGA_CHANNEL_NUMBER + channel] = new TH2D(
                h2_name.c_str(),
                h2_title.c_str(),
                machine_gun_samples, -0.5, machine_gun_samples - 0.5,
                512, 0, 1024
            );
            h2_adc_waveforms[vldb_id * FPGA_CHANNEL_NUMBER + channel]->SetDirectory(nullptr);
        }
    }

    // * --- ADC Peak Position Histogram ---
    TH1D* h1_adc_peak_position = new TH1D("h1_adc_peak_position", "ADC Peak Position;Sample Index;Counts", machine_gun_samples, -0.5, machine_gun_samples - 0.5);
    h1_adc_peak_position->SetDirectory(nullptr);

    // * --- ADC Peak Value - Channel Correlation Histogram ---
    TH2D *h2_adc_peak_channel_correlation = new TH2D(
        "h2_adc_peak_channel_correlation",
        "ADC Peak Value vs Channel Correlation;Channel;ADC Peak Value",
        FPGA_CHANNEL_NUMBER, -0.5, FPGA_CHANNEL_NUMBER - 0.5,
        512, 0, 1024
    );
    h2_adc_peak_channel_correlation->SetDirectory(nullptr);

    // * --- ADC Sum distribution Histogram ---
    std::vector<double> adc_sum_list;
    adc_sum_list.reserve(entry_max);

    const int adc_peak_min_index = 4;
    const int adc_peak_max_index = 8;

    // start event loop
    for (int entry = 0; entry < entry_max; entry++) {
        input_tree->GetEntry(entry);
        double adc_sum = 0.0;
        for (int vldb_id = 0; vldb_id < vldb_number; vldb_id++) {
            // channel loop
            for (int channel = 0; channel < FPGA_CHANNEL_NUMBER; channel++) {
                std::vector<UInt_t> adc_pedestal_samples; // only take the first 3 samples
                int adc_peak_index = -1;
                UInt_t adc_peak_value = 0;
                int adc_peak_ranged_index = -1;
                UInt_t adc_peak_ranged_value = 0;
                adc_pedestal_samples.reserve(3);
                for (int sample = 0; sample < machine_gun_samples; sample++) {
                    int idx = sample*FPGA_CHANNEL_NUMBER + channel;
                    UInt_t adc_value = val0_list_pools[vldb_id][0][idx];
                    UInt_t tot_value = val1_list_pools[vldb_id][0][idx];
                    UInt_t toa_value = val2_list_pools[vldb_id][0][idx];
                    // spdlog::info("Event {} VLDB {} Channel {} Sample {}: ADC {}, ToT {}, ToA {}", entry, vldb_id, channel, sample, adc_value, tot_value, toa_value);
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
                    // fill ADC waveform histogram
                    h2_adc_waveforms[vldb_id * FPGA_CHANNEL_NUMBER + channel]->Fill(sample, adc_value);
                    h1_adc_peak_position->Fill(adc_peak_index);
                } // end of sample loop
                // calculate the pedestal
                double adc_pedestal = pedestal_median_of_first3(adc_pedestal_samples);
                double adc_peak_value_pede_sub = static_cast<double>(adc_peak_ranged_value) - adc_pedestal;
                if (adc_peak_value_pede_sub > 0) {
                    h2_adc_peak_channel_correlation->Fill(channel, adc_peak_value_pede_sub);
                    adc_sum += adc_peak_value_pede_sub;
                }
            } // end of channel loop
        } // end of vldb loop
        adc_sum_list.push_back(adc_sum); // because we don't know the max adc sum value beforehand
    } // end of event loop
    input_root->Close();

    spdlog::info("Finished processing {} events", entry_max);

    // * === Plot to output file =============================================================
    // * =====================================================================================
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

    std::string annotation_canvas_title = CANVAS_TITLE;
    std::string annotation_testbeam_title = TESTBEAM_TITLE;
    output_root->cd();
    const std::string out_pdf = script_output_file.substr(0, script_output_file.find_last_of(".")) + ".pdf";

    // --- Write to output file ------------------------------------------------------------
    TCanvas adc_waveform_canvas("adc_waveform_canvas", "ADC Waveform Canvas", 1600, 1200);
    draw_mosaic_fixed(adc_waveform_canvas, h2_adc_waveforms, topo_wave);
    output_root->cd();
    adc_waveform_canvas.Modified();
    adc_waveform_canvas.Update();
    adc_waveform_canvas.Print((out_pdf + "[").c_str()); // begin of pdf
    adc_waveform_canvas.Write();
    adc_waveform_canvas.Close();

    TCanvas adc_peak_position_canvas("adc_peak_position_canvas", "ADC Peak Position Canvas", 800, 600);
    std::string canvas_info = "ADC Peak Position Distribution";
    format_1d_hist_canvas(&adc_peak_position_canvas, h1_adc_peak_position, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
    adc_peak_position_canvas.Print(out_pdf.c_str());
    adc_peak_position_canvas.Write();
    adc_peak_position_canvas.Close();

    TCanvas adc_peak_channel_correlation_canvas("adc_peak_channel_correlation_canvas", "ADC Peak Channel Correlation Canvas", 800, 600);
    canvas_info = "ADC Peak Value vs Channel Correlation";
    format_2d_hist_canvas(&adc_peak_channel_correlation_canvas, h2_adc_peak_channel_correlation, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
    adc_peak_channel_correlation_canvas.Print(out_pdf.c_str());
    adc_peak_channel_correlation_canvas.Write();
    adc_peak_channel_correlation_canvas.Close();

    // ! -- Raw ADC Sum distribution Histogram ---
    // ! =====================================================================================
    TCanvas adc_sum_distribution_canvas("adc_sum_distribution_canvas", "ADC Sum Distribution Canvas", 800, 600);
    canvas_info = "ADC Sum Distribution, Run " + script_input_run_number;
    // determine 90% percentile max value for better visualization
    auto sorted_adc_sum_list = adc_sum_list;
    std::sort(sorted_adc_sum_list.begin(), sorted_adc_sum_list.end());
    double adc_sum_90pct_max = sorted_adc_sum_list[static_cast<size_t>(0.9 * sorted_adc_sum_list.size()) - 1];
    // create histogram and fill
    TH1D* h1_adc_sum_distribution = new TH1D("h1_adc_sum_distribution", "ADC Sum Distribution;ADC Sum;Counts", 256, 0, adc_sum_90pct_max*1.8);
    for (const auto& adc_sum_value : adc_sum_list) {
        h1_adc_sum_distribution->Fill(adc_sum_value);
    }
    format_1d_hist_canvas(&adc_sum_distribution_canvas, h1_adc_sum_distribution, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
    // add additional statistics
    double adc_sum_mean = h1_adc_sum_distribution->GetMean();
    double adc_sum_entry = h1_adc_sum_distribution->GetEntries();
    double adc_sum_rms = h1_adc_sum_distribution->GetRMS();
    TPaveText *pave_stats = new TPaveText(0.6,0.7,0.88,0.88,"NDC");
    pave_stats->SetFillColorAlpha(kWhite,0.0);
    pave_stats->SetBorderSize(0);
    // right align text
    pave_stats->SetTextAlign(31);
    pave_stats->SetTextFont(42);
    pave_stats->SetTextSize(0.03);
    pave_stats->SetTextColor(kBlue+2);
    pave_stats->AddText(("Entries: " + std::to_string(static_cast<int>(adc_sum_entry))).c_str());
    // only keep two decimal points for mean
    pave_stats->AddText(("Mean: " + std::to_string(static_cast<int>(adc_sum_mean * 100)) .insert(std::to_string(static_cast<int>(adc_sum_mean * 100)).length() - 2, ".")).c_str());
    // only keep two decimal points for rms
    pave_stats->AddText(("RMS: " + std::to_string(static_cast<int>(adc_sum_rms * 100)) .insert(std::to_string(static_cast<int>(adc_sum_rms * 100)).length() - 2, ".")).c_str());
    pave_stats->Draw();

    adc_sum_distribution_canvas.Print(out_pdf.c_str());
    adc_sum_distribution_canvas.Write();
    adc_sum_distribution_canvas.Close();

    // ! -- Fitted ADC Sum distribution Histogram ---
    TCanvas adc_sum_distribution_fitted_canvas("adc_sum_distribution_fitted_canvas", "ADC Sum Distribution Fitted Canvas", 800, 600);
    canvas_info = "Fitted ADC Sum Distribution, Run " + script_input_run_number;
    TH1D* h1_adc_sum_distribution_fitted = (TH1D*) h1_adc_sum_distribution->Clone("h1_adc_sum_distribution_fitted");
    format_1d_hist_canvas(&adc_sum_distribution_fitted_canvas, h1_adc_sum_distribution_fitted, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);
    // pre-fit for initial values
    double fit_min = adc_sum_mean - 2 * adc_sum_rms;
    double fit_max = adc_sum_mean + 3 * adc_sum_rms;
    if (fit_min < 0) fit_min = 0;
    double fit_amp_init = h1_adc_sum_distribution_fitted->GetMaximum();
    TF1 *gaus_fit_pre = new TF1("gaus_fit_pre", "gaus", fit_min, fit_max);
    // set initial values
    gaus_fit_pre->SetParameter(0, fit_amp_init);
    gaus_fit_pre->SetParameter(1, adc_sum_mean);
    gaus_fit_pre->SetParameter(2, adc_sum_rms);
    h1_adc_sum_distribution_fitted->Fit(gaus_fit_pre, "RQ");
    double pre_fit_amp = gaus_fit_pre->GetParameter(0);
    double pre_fit_mean = gaus_fit_pre->GetParameter(1);
    double pre_fit_sigma = gaus_fit_pre->GetParameter(2);
    double pre_fit_chi2 = gaus_fit_pre->GetChisquare();
    double pre_fit_ndf = gaus_fit_pre->GetNDF();
    // draw fitted function
    gaus_fit_pre->SetLineColor(kRed);
    gaus_fit_pre->SetLineWidth(2);
    gaus_fit_pre->Draw("same");

    // do the second round of fit with initial values from pre-fit
    fit_min = pre_fit_mean - 1.5 * pre_fit_sigma;
    fit_max = pre_fit_mean + 2.5 * pre_fit_sigma;
    if (fit_min < 0) fit_min = 0;
    TF1 *gaus_fit_final = new TF1("gaus_fit_final", "gaus", fit_min, fit_max);
    // set initial values
    gaus_fit_final->SetParameter(0, pre_fit_amp);
    gaus_fit_final->SetParameter(1, pre_fit_mean);
    gaus_fit_final->SetParameter(2, pre_fit_sigma);
    h1_adc_sum_distribution_fitted->Fit(gaus_fit_final, "RQ+");

    double pre_fit_2_amp = gaus_fit_final->GetParameter(0);
    double pre_fit_2_mean = gaus_fit_final->GetParameter(1);
    double pre_fit_2_sigma = gaus_fit_final->GetParameter(2);
    double pre_fit_2_chi2 = gaus_fit_final->GetChisquare();
    double pre_fit_2_ndf = gaus_fit_final->GetNDF();
    // draw final fitted function
    gaus_fit_final->SetLineColor(kGreen+2);
    gaus_fit_final->SetLineWidth(2);
    gaus_fit_final->Draw("same");

    // add additional statistics
    TPaveText *pave_stats_pre = new TPaveText(0.55,0.73,0.90,0.88,"NDC");
    pave_stats_pre->SetFillColorAlpha(kWhite,0.0);
    pave_stats_pre->SetBorderSize(0);
    pave_stats_pre->SetTextAlign(31);
    pave_stats_pre->SetTextFont(42);
    pave_stats_pre->SetTextSize(0.03);
    pave_stats_pre->SetTextColor(kRed);
    pave_stats_pre->AddText("Pre-Fit:");
    pave_stats_pre->AddText(("  Mean: " + std::to_string(static_cast<int>(pre_fit_mean * 100)) .insert(std::to_string(static_cast<int>(pre_fit_mean * 100)).length() - 2, ".")).c_str());
    pave_stats_pre->AddText(("  Sigma: " + std::to_string(static_cast<int>(pre_fit_sigma * 100)) .insert(std::to_string(static_cast<int>(pre_fit_sigma * 100)).length() - 2, ".")).c_str());
    pave_stats_pre->AddText(("  Chi2/NDF: " + std::to_string(static_cast<int>(pre_fit_chi2 * 100)) .insert(std::to_string(static_cast<int>(pre_fit_chi2 * 100)).length() - 2, ".") + "/" + std::to_string(static_cast<int>(pre_fit_ndf))).c_str());
    pave_stats_pre->Draw();

    TPaveText *pave_stats_pre_2 = new TPaveText(0.55,0.58,0.90,0.73,"NDC");
    pave_stats_pre_2->SetFillColorAlpha(kWhite,0.0);
    pave_stats_pre_2->SetBorderSize(0);
    pave_stats_pre_2->SetTextAlign(31);
    pave_stats_pre_2->SetTextFont(42);
    pave_stats_pre_2->SetTextSize(0.03);
    pave_stats_pre_2->SetTextColor(kGreen+2);
    pave_stats_pre_2->AddText("Pre-Fit 2:");
    pave_stats_pre_2->AddText(("  Mean: " + std::to_string(static_cast<int>(pre_fit_2_mean * 100)) .insert(std::to_string(static_cast<int>(pre_fit_2_mean * 100)).length() - 2, ".")).c_str());
    pave_stats_pre_2->AddText(("  Sigma: " + std::to_string(static_cast<int>(pre_fit_2_sigma * 100)) .insert(std::to_string(static_cast<int>(pre_fit_2_sigma * 100)).length() - 2, ".")).c_str());
    pave_stats_pre_2->AddText(("  Chi2/NDF: " + std::to_string(static_cast<int>(pre_fit_2_chi2 * 100 + 0.5) / 100.0).substr(0, std::to_string(static_cast<int>(pre_fit_2_chi2 * 100 + 0.5) / 100.0).find('.') + 3) + "/" + std::to_string(static_cast<int>(pre_fit_2_ndf))).c_str());
    pave_stats_pre_2->Draw();

    adc_sum_distribution_fitted_canvas.Modified();
    adc_sum_distribution_fitted_canvas.Update();
    adc_sum_distribution_fitted_canvas.Print(out_pdf.c_str());
    adc_sum_distribution_fitted_canvas.Write();
    adc_sum_distribution_fitted_canvas.Close();

    // ! --- Do the real gaussian fit ---
    TCanvas adc_sum_distribution_gausfit_canvas("adc_sum_distribution_gausfit_canvas", "ADC Sum Distribution Gaus Fit Canvas", 800, 600);
    canvas_info = "Gaussian Fitted ADC, Run " + script_input_run_number;
    TH1D* h1_adc_sum_distribution_gausfit = (TH1D*) h1_adc_sum_distribution->Clone("h1_adc_sum_distribution_gausfit");
    format_1d_hist_canvas(&adc_sum_distribution_gausfit_canvas, h1_adc_sum_distribution_gausfit, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);

    // define the fitting range
    const std::vector<double> fit_range_sigmas = {2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5}; // in sigma
    const std::vector<double> fit_range_offsets = {-0.4, -0.2, 0.0, 0.2, 0.4}; // in sigma
    std::vector<double> fit_results_means;
    std::vector<double> fit_results_mean_errs;
    std::vector<double> fit_results_sigmas;
    std::vector<double> fit_results_sigma_errs;
    std::vector<double> fit_results_chi2s;
    std::vector<double> fit_results_ndfs;

    for (const auto& range_sigma : fit_range_sigmas) {
        for (const auto& range_offset : fit_range_offsets) {
            double fit_min = pre_fit_2_mean - range_sigma * pre_fit_2_sigma + range_offset * pre_fit_2_sigma;
            double fit_max = pre_fit_2_mean + range_sigma * pre_fit_2_sigma + range_offset * pre_fit_2_sigma;
            if (fit_min < 0) fit_min = 0;
            TF1 *gaus_fit_real = new TF1("gaus_fit_real", "gaus", fit_min, fit_max);
            // set initial values
            gaus_fit_real->SetParameter(0, pre_fit_2_amp);
            gaus_fit_real->SetParameter(1, pre_fit_2_mean);
            gaus_fit_real->SetParameter(2, pre_fit_2_sigma);
            h1_adc_sum_distribution_gausfit->Fit(gaus_fit_real, "RNQ+");

            double fit_mean = gaus_fit_real->GetParameter(1);
            double fit_mean_err = gaus_fit_real->GetParError(1);
            double fit_sigma = gaus_fit_real->GetParameter(2);
            double fit_sigma_err = gaus_fit_real->GetParError(2);

            double fit_chi2 = gaus_fit_real->GetChisquare();
            double fit_ndf = gaus_fit_real->GetNDF();

            fit_results_means.push_back(fit_mean);
            fit_results_mean_errs.push_back(fit_mean_err);
            fit_results_sigmas.push_back(fit_sigma);
            fit_results_sigma_errs.push_back(fit_sigma_err);
            fit_results_chi2s.push_back(fit_chi2);
            fit_results_ndfs.push_back(fit_ndf);

            gaus_fit_real->SetLineColorAlpha(kCyan+2, 0.2);
            gaus_fit_real->SetLineWidth(2);
            gaus_fit_real->Draw("same");

            spdlog::info("Gaus Fit Range Sigma: {}, Offset: {} => Mean: {}, Sigma: {}, Chi2/NDF: {}/{}", 
                range_sigma, range_offset, fit_mean, fit_sigma, fit_chi2, fit_ndf
            );
        }
    }

    double mean_avg=0.0;
    double mean_err_sys=0.0;
    double mean_err_stat=0.0;
    double sigma_avg=0.0;
    double sigma_err_sys=0.0;
    double sigma_err_stat=0.0;

    spdlog::info("Calculating mean and sigma from multiple fit results...");
    mean_sigma_list_calculator(
        fit_results_means,
        fit_results_mean_errs,
        fit_results_sigmas,
        fit_results_sigma_errs,
        mean_avg,
        mean_err_sys,
        mean_err_stat,
        sigma_avg,
        sigma_err_sys,
        sigma_err_stat
    );
    spdlog::info("Mean: {} +/- {} (stat) +/- {} (sys)", mean_avg, mean_err_stat, mean_err_sys);
    spdlog::info("Sigma: {} +/- {} (stat) +/- {} (sys)", sigma_avg, sigma_err_stat, sigma_err_sys);

    double resolution = (sigma_avg / mean_avg) * 100.0;
    double mean_err = std::sqrt(mean_err_stat * mean_err_stat + mean_err_sys * mean_err_sys);
    double sigma_err = std::sqrt(sigma_err_stat * sigma_err_stat + sigma_err_sys * sigma_err_sys);
    double resolution_err = resolution * std::sqrt( (sigma_err / sigma_avg) * (sigma_err / sigma_avg) + (mean_err / mean_avg) * (mean_err / mean_avg) );
    spdlog::info("Calculated Resolution: {} % +/- {} %", resolution, resolution_err);

    // add additional statistics
    TPaveText *pave_stats_gausfit = new TPaveText(0.55,0.6,0.90,0.88,"NDC");
    pave_stats_gausfit->SetFillColorAlpha(kWhite,0.0);
    pave_stats_gausfit->SetBorderSize(0);
    pave_stats_gausfit->SetTextAlign(31);
    pave_stats_gausfit->SetTextFont(42);
    pave_stats_gausfit->SetTextSize(0.03);
    pave_stats_gausfit->SetTextColor(kCyan+2);
    pave_stats_gausfit->AddText("Gaussian Fit Results:");
    
    auto format_to_2_decimals = [](double value) -> std::string {
        if (value == 0.0) {
            return "0.00";
        }
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << value;
        return oss.str();
    };
    
    std::string mean_str = format_to_2_decimals(mean_avg);
    std::string mean_err_stat_str = format_to_2_decimals(mean_err_stat);
    std::string mean_err_sys_str = format_to_2_decimals(mean_err_sys);
    std::string sigma_str = format_to_2_decimals(sigma_avg);
    std::string sigma_err_stat_str = format_to_2_decimals(sigma_err_stat);
    std::string sigma_err_sys_str = format_to_2_decimals(sigma_err_sys);
    
    std::string mean_text = "  Mean: " + mean_str + " #pm " + mean_err_stat_str + " (stat)";
    pave_stats_gausfit->AddText(mean_text.c_str());
    std::string mean_sys_text = "        #pm " + mean_err_sys_str + " (sys)";
    pave_stats_gausfit->AddText(mean_sys_text.c_str());

    std::string sigma_text = "  Sigma: " + sigma_str + " #pm " + sigma_err_stat_str + " (stat)";
    pave_stats_gausfit->AddText(sigma_text.c_str());
    std::string sigma_sys_text = "         #pm " + sigma_err_sys_str + " (sys)";
    pave_stats_gausfit->AddText(sigma_sys_text.c_str());

    std::string resolution_str = format_to_2_decimals(resolution);
    std::string resolution_err_str = format_to_2_decimals(resolution_err);
    std::string resolution_text = "  Resolution: " + resolution_str + " % #pm " + resolution_err_str + " %";
    pave_stats_gausfit->AddText(resolution_text.c_str());
    pave_stats_gausfit->Draw();

    adc_sum_distribution_gausfit_canvas.Modified();
    adc_sum_distribution_gausfit_canvas.Update();
    adc_sum_distribution_gausfit_canvas.Print(out_pdf.c_str());
    adc_sum_distribution_gausfit_canvas.Write();
    adc_sum_distribution_gausfit_canvas.Close();
    
     // ! --- Do the real crystalball fit ---
    TCanvas adc_sum_distribution_cbfit_canvas("adc_sum_distribution_cbfit_canvas", "ADC Sum Distribution CB Fit Canvas", 800, 600);
    canvas_info = "Crystal Ball Fitted ADC, Run " + script_input_run_number;
    TH1D* h1_adc_sum_distribution_cbfit = (TH1D*) h1_adc_sum_distribution->Clone("h1_adc_sum_distribution_cbfit");
    format_1d_hist_canvas(&adc_sum_distribution_cbfit_canvas, h1_adc_sum_distribution_cbfit, kBlue+2, annotation_canvas_title, annotation_testbeam_title, canvas_info);

    // use the same fitting range as gaussian fit
    std::vector<double> cb_fit_results_means;
    std::vector<double> cb_fit_results_mean_errs;
    std::vector<double> cb_fit_results_sigmas;
    std::vector<double> cb_fit_results_sigma_errs;
    std::vector<double> cb_fit_results_eff_sigmas;
    std::vector<double> cb_fit_results_eff_sigma_errs;
    std::vector<double> cb_fit_results_alphas;
    std::vector<double> cb_fit_results_nbs;
    std::vector<double> cb_fit_results_chi2s;
    std::vector<double> cb_fit_results_ndfs;

    for (const auto& range_sigma : fit_range_sigmas) {
        for (const auto& range_offset : fit_range_offsets) {
            double fit_min = pre_fit_2_mean - range_sigma * pre_fit_2_sigma + range_offset * pre_fit_2_sigma;
            double fit_max = pre_fit_2_mean + range_sigma * pre_fit_2_sigma + range_offset * pre_fit_2_sigma;
            if (fit_min < 0) fit_min = 0;
            TF1 *cb_fit_real = new TF1("cb_fit_real", crystallball_left, fit_min, fit_max, 5);
            // set initial values
            cb_fit_real->SetParameter(0, pre_fit_2_amp);      // N
            cb_fit_real->SetParameter(1, 1.5);                // alpha
            cb_fit_real->SetParameter(2, 5.0);                // n
            cb_fit_real->SetParameter(3, pre_fit_2_mean);     // mean
            cb_fit_real->SetParameter(4, pre_fit_2_sigma);    // sigma

            // cb_fit_real->SetParLimits(1, 1.5, 10.0); // alpha
            // cb_fit_real->SetParLimits(2, 5.0, 50.0); // n
            h1_adc_sum_distribution_cbfit->Fit(cb_fit_real, "RNQ+");

            double fit_mean = cb_fit_real->GetParameter(3);
            double fit_mean_err = cb_fit_real->GetParError(3);
            double fit_sigma = cb_fit_real->GetParameter(4);
            double fit_sigma_err = cb_fit_real->GetParError(4);
            double fit_alpha = cb_fit_real->GetParameter(1);
            double fit_n = cb_fit_real->GetParameter(2);

            double fit_chi2 = cb_fit_real->GetChisquare();
            double fit_ndf = cb_fit_real->GetNDF();

            cb_fit_results_means.push_back(fit_mean);
            cb_fit_results_mean_errs.push_back(fit_mean_err);
            cb_fit_results_sigmas.push_back(fit_sigma);
            cb_fit_results_sigma_errs.push_back(fit_sigma_err);
            cb_fit_results_alphas.push_back(fit_alpha);
            cb_fit_results_nbs.push_back(fit_n);
            cb_fit_results_chi2s.push_back(fit_chi2);
            cb_fit_results_ndfs.push_back(fit_ndf);

            // calculate effective sigma (68.3% quantile)
            // double q16 = quantile_from_tf1(cb_fit_real, 0.158655253931, h1_adc_sum_distribution_cbfit->GetXaxis()->GetXmin(), h1_adc_sum_distribution_cbfit->GetXaxis()->GetXmax());
            // double q84 = quantile_from_tf1(cb_fit_real, 0.841344746069, h1_adc_sum_distribution_cbfit->GetXaxis()->GetXmin(), h1_adc_sum_distribution_cbfit->GetXaxis()->GetXmax());
            // double effective_sigma = (q84 - q16) / 2.0;
            // double effective_sigma_err = fit_sigma_err;

            // define effective sigma as 67% area for better stability
            double prob[2] = {0.158655253931, 0.841344746069};
            double quantiles[2];
            h1_adc_sum_distribution_cbfit->GetQuantiles(2, quantiles, prob);
            double effective_sigma = (quantiles[1] - quantiles[0]) / 2.0;
            // give gaussian estimate for effective sigma error
            double effective_sigma_err = fit_sigma_err;

            // define effective sigma from RMS 90
            // auto r = computeMeanRMS90(h1_adc_sum_distribution_cbfit);
            // double effective_sigma = r.sigma_gauss;
            // double effective_sigma_err = r.err_sigma_gauss;

            cb_fit_results_eff_sigmas.push_back(effective_sigma);
            cb_fit_results_eff_sigma_errs.push_back(effective_sigma_err);

            cb_fit_real->SetLineColorAlpha(kMagenta+2, 0.2);
            cb_fit_real->SetLineWidth(2);
            cb_fit_real->Draw("same");

            spdlog::info("CB Fit Range Sigma: {}, Offset: {} => Mean: {}, Sigma: {}, Alpha: {}, N: {}, effective Sigma: {}, Chi2/NDF: {}/{}", 
                range_sigma, range_offset, fit_mean, fit_sigma, fit_alpha, fit_n, effective_sigma, fit_chi2, fit_ndf
            );
        }
    }

    double cb_mean_avg=0.0;
    double cb_mean_err_sys=0.0;
    double cb_mean_err_stat=0.0;
    double cb_sigma_avg=0.0;
    double cb_sigma_err_sys=0.0;
    double cb_sigma_err_stat=0.0;
    double cb_effective_sigma_avg=0.0;
    double cb_effective_sigma_err_sys=0.0;
    double cb_effective_sigma_err_stat=0.0;

    spdlog::info("Calculating mean and sigma from multiple CB fit results...");
    mean_sigma_list_calculator(
        cb_fit_results_means,
        cb_fit_results_mean_errs,
        cb_fit_results_sigmas,
        cb_fit_results_sigma_errs,
        cb_mean_avg,
        cb_mean_err_sys,
        cb_mean_err_stat,
        cb_sigma_avg,
        cb_sigma_err_sys,
        cb_sigma_err_stat
    );

    // do another round for effective sigma
    mean_sigma_list_calculator(
        cb_fit_results_means,
        cb_fit_results_mean_errs,
        cb_fit_results_eff_sigmas,
        cb_fit_results_eff_sigma_errs,
        cb_effective_sigma_avg,
        cb_effective_sigma_err_sys,
        cb_effective_sigma_err_stat,
        cb_effective_sigma_avg,
        cb_effective_sigma_err_sys,
        cb_effective_sigma_err_stat
    );
    spdlog::info("CB Mean: {} +/- {} (stat) +/- {} (sys)", cb_mean_avg, cb_mean_err_stat, cb_mean_err_sys);
    spdlog::info("CB Sigma: {} +/- {} (stat) +/- {} (sys)", cb_sigma_avg, cb_sigma_err_stat, cb_sigma_err_sys);
    spdlog::info("CB Effective Sigma: {} +/- {} (stat) +/- {} (sys)", cb_effective_sigma_avg, cb_effective_sigma_err_stat, cb_effective_sigma_err_sys);
    double cb_resolution = (cb_sigma_avg / cb_mean_avg) * 100.0;
    double cb_mean_err = std::sqrt(cb_mean_err_stat * cb_mean_err_stat + cb_mean_err_sys * cb_mean_err_sys);
    double cb_sigma_err = std::sqrt(cb_sigma_err_stat * cb_sigma_err_stat + cb_sigma_err_sys * cb_sigma_err_sys);
    double cb_resolution_err = cb_resolution * std::sqrt( (cb_sigma_err / cb_sigma_avg) * (cb_sigma_err / cb_sigma_avg) + (cb_mean_err / cb_mean_avg) * (cb_mean_err / cb_mean_avg) );
    spdlog::info("Calculated CB Resolution: {} % +/- {} %", cb_resolution, cb_resolution_err);
    double cb_effective_resolution = (cb_effective_sigma_avg / cb_mean_avg) * 100.0;
    double cb_effective_sigma_err = std::sqrt(cb_effective_sigma_err_stat * cb_effective_sigma_err_stat + cb_effective_sigma_err_sys * cb_effective_sigma_err_sys);
    double cb_effective_resolution_err = cb_effective_resolution * std::sqrt( (cb_effective_sigma_err / cb_effective_sigma_avg) * (cb_effective_sigma_err / cb_effective_sigma_avg) + (cb_mean_err / cb_mean_avg) * (cb_mean_err / cb_mean_avg) );
    spdlog::info("Calculated CB Effective Resolution: {} % +/- {} %", cb_effective_resolution, cb_effective_resolution_err);

    // add additional statistics
    TPaveText *pave_stats_cbfit = new TPaveText(0.55,0.4,0.90,0.88,"NDC");
    pave_stats_cbfit->SetFillColorAlpha(kWhite,0.0);
    pave_stats_cbfit->SetBorderSize(0);
    pave_stats_cbfit->SetTextAlign(31);
    pave_stats_cbfit->SetTextFont(42);
    pave_stats_cbfit->SetTextSize(0.03);
    pave_stats_cbfit->SetTextColor(kMagenta+2);
    pave_stats_cbfit->AddText("Crystal Ball Fit Results:");
    std::string cb_mean_str = format_to_2_decimals(cb_mean_avg);
    std::string cb_mean_err_stat_str = format_to_2_decimals(cb_mean_err_stat);
    std::string cb_mean_err_sys_str = format_to_2_decimals(cb_mean_err_sys);
    std::string cb_sigma_str = format_to_2_decimals(cb_sigma_avg);
    std::string cb_sigma_err_stat_str = format_to_2_decimals(cb_sigma_err_stat);
    std::string cb_sigma_err_sys_str = format_to_2_decimals(cb_sigma_err_sys);
    std::string cb_mean_text = "  Mean: " + cb_mean_str + " #pm " + cb_mean_err_stat_str + " (stat)";
    pave_stats_cbfit->AddText(cb_mean_text.c_str());
    std::string cb_mean_sys_text = "        #pm " + cb_mean_err_sys_str + " (sys)";
    pave_stats_cbfit->AddText(cb_mean_sys_text.c_str());
    std::string cb_sigma_text = "  Sigma: " + cb_sigma_str + " #pm " + cb_sigma_err_stat_str + " (stat)";
    pave_stats_cbfit->AddText(cb_sigma_text.c_str());
    std::string cb_sigma_sys_text = "         #pm " + cb_sigma_err_sys_str + " (sys)";
    pave_stats_cbfit->AddText(cb_sigma_sys_text.c_str());
    std::string cb_resolution_str = format_to_2_decimals(cb_resolution);
    std::string cb_resolution_err_str = format_to_2_decimals(cb_resolution_err);
    std::string cb_resolution_text = "  Core Resolution: " + cb_resolution_str + " % #pm " + cb_resolution_err_str + " %";
    pave_stats_cbfit->AddText(cb_resolution_text.c_str());
    std::string cb_effective_sigma_str = format_to_2_decimals(cb_effective_sigma_avg);
    std::string cb_effective_sigma_err_stat_str = format_to_2_decimals(cb_effective_sigma_err_stat);
    std::string cb_effective_sigma_err_sys_str = format_to_2_decimals(cb_effective_sigma_err_sys);
    std::string cb_effective_sigma_text = "  Eff Sigma: " + cb_effective_sigma_str;
    pave_stats_cbfit->AddText(cb_effective_sigma_text.c_str());
    std::string cb_effective_resolution_str = format_to_2_decimals(cb_effective_resolution);
    std::string cb_effective_resolution_text = "  Eff Resolution: " + cb_effective_resolution_str + " %";
    pave_stats_cbfit->AddText(cb_effective_resolution_text.c_str());
    pave_stats_cbfit->Draw();

    adc_sum_distribution_cbfit_canvas.Modified();
    adc_sum_distribution_cbfit_canvas.Update();
    adc_sum_distribution_cbfit_canvas.Print(out_pdf.c_str());
    adc_sum_distribution_cbfit_canvas.Write();
    adc_sum_distribution_cbfit_canvas.Close();


    // ! =====================================================================================


    // end of pdf, create dummy canvas
    TCanvas end_canvas("end_canvas", "End Canvas", 800, 600);
    end_canvas.Print((out_pdf + "]").c_str());
    end_canvas.Close();

    // write important values for comparison
    TParameter<double> param_gaus_fit_mean("gaus_fit_mean", mean_avg); 
    TParameter<double> param_gaus_fit_sigma("gaus_fit_sigma", sigma_avg);
    TParameter<double> param_gaus_fit_resolution("gaus_fit_resolution", resolution);
    TParameter<double> param_gaus_fit_mean_err_stat("gaus_fit_mean_err_stat", mean_err_stat);
    TParameter<double> param_gaus_fit_mean_err_sys("gaus_fit_mean_err_sys", mean_err_sys);
    TParameter<double> param_gaus_fit_sigma_err_stat("gaus_fit_sigma_err_stat", sigma_err);
    TParameter<double> param_gaus_fit_sigma_err_sys("gaus_fit_sigma_err_sys", sigma_err);
    TParameter<double> param_gaus_fit_resolution_err("gaus_fit_resolution_err", resolution_err);

    TParameter<double> param_cb_fit_mean("cb_fit_mean", cb_mean_avg);
    TParameter<double> param_cb_fit_sigma("cb_fit_sigma", cb_sigma_avg);
    TParameter<double> param_cb_fit_effective_sigma("cb_fit_effective_sigma", cb_effective_sigma_avg);
    TParameter<double> param_cb_fit_resolution("cb_fit_resolution", cb_resolution);
    TParameter<double> param_cb_fit_effective_resolution("cb_fit_effective_resolution", cb_effective_resolution);
    TParameter<double> param_cb_fit_mean_err_stat("cb_fit_mean_err_stat", cb_mean_err_stat);
    TParameter<double> param_cb_fit_mean_err_sys("cb_fit_mean_err_sys", cb_mean_err_sys);
    TParameter<double> param_cb_fit_sigma_err_stat("cb_fit_sigma_err_stat", cb_sigma_err);
    TParameter<double> param_cb_fit_sigma_err_sys("cb_fit_sigma_err_sys", cb_sigma_err);
    TParameter<double> param_cb_fit_effective_sigma_err_stat("cb_fit_effective_sigma_err_stat", cb_effective_sigma_err_stat);
    TParameter<double> param_cb_fit_effective_sigma_err_sys("cb_fit_effective_sigma_err_sys", cb_effective_sigma_err_sys);
    TParameter<double> param_cb_fit_resolution_err("cb_fit_resolution_err", cb_resolution_err);
    TParameter<double> param_cb_fit_effective_resolution_err("cb_fit_effective_resolution_err", cb_effective_resolution_err);

    output_root->cd();
    param_gaus_fit_mean.Write();
    param_gaus_fit_sigma.Write();
    param_gaus_fit_resolution.Write();
    param_gaus_fit_mean_err_stat.Write();
    param_gaus_fit_mean_err_sys.Write();
    param_gaus_fit_sigma_err_stat.Write();
    param_gaus_fit_sigma_err_sys.Write();
    param_gaus_fit_resolution_err.Write();

    param_cb_fit_mean.Write();
    param_cb_fit_sigma.Write();
    param_cb_fit_effective_sigma.Write();
    param_cb_fit_resolution.Write();
    param_cb_fit_effective_resolution.Write();
    param_cb_fit_mean_err_stat.Write();
    param_cb_fit_mean_err_sys.Write();
    param_cb_fit_sigma_err_stat.Write();
    param_cb_fit_sigma_err_sys.Write();
    param_cb_fit_effective_sigma_err_stat.Write();
    param_cb_fit_effective_sigma_err_sys.Write();
    param_cb_fit_resolution_err.Write();
    param_cb_fit_effective_resolution_err.Write();
    

    output_root->Close();

    spdlog::info("Output file {} has been saved.", script_output_file);

    return 0;
}

// inline double crystallball_left(double *x, double *par) {
//     // par[0] = A
//     // par[1] = a
//     // par[2] = n
//     // par[3] = mean
//     // par[4] = sigma
//     const double N      = par[0];
//     const double alpha  = par[1];
//     const double n      = par[2];
//     const double mean   = par[3];
//     const double sigma  = par[4];

//     double z = (x[0] - mean) / sigma;
//     if (z > -alpha) {
//         return N * exp(-0.5 * z * z);
//     } else {
//         double A = pow(n / fabs(alpha), n) * exp(-0.5 * alpha * alpha);
//         double B = n / fabs(alpha) - fabs(alpha);
//         return N * A * pow(B - z, -n);
//     }
// }

