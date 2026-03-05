#include <cmath>
#include <fstream>
#include <string>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "TKey.h"

int main(int argc, char **argv) {
    gROOT->SetBatch(kTRUE);
    std::string script_input_file, script_output_file;
    int script_n_events = -1;
    bool script_verbose = false;
    std::string script_name = __FILE__;
    const std::string script_version = "0.1";
    std::string script_output_folder;
    std::vector<int> color_list = {kRed, kBlue, kGreen+2, kMagenta+2, kCyan+2, kOrange+2, kViolet+2};

    script_name = script_name.substr(script_name.find_last_of("/\\") + 1).substr(0, script_name.find_last_of("."));

    configure_logger(false);

    cxxopts::Options options(script_name, "Generate heatmaps from machine gun data");

    options.add_options()
        ("f,file", "Input .json config file", cxxopts::value<std::string>())
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
    json config_json;
    std::ifstream config_ifs(script_input_file);
    if (!config_ifs.is_open()) {
        spdlog::error("Failed to open config json file {}", script_input_file);
        return 1;
    }
    config_ifs >> config_json;
    config_ifs.close();
    script_output_file = parsed["output"].as<std::string>();
    script_n_events    = parsed["events"].as<int>();
    script_verbose     = parsed["verbose"].as<bool>();

    configure_logger(script_verbose);

    std::vector<int> run_numbers;
    std::vector<double> run_energies;
    std::vector<std::string> run_labels;
    std::string beam_type;
    std::string config_description;
    
    bool use_peak_integral = false;
    bool enable_gaussian_fit = false;
    int target_event_number = -1;
    bool found_run_energies = false;

    try {
        run_numbers = config_json.at("run_numbers").get<std::vector<int>>();
        if (config_json.contains("run_energies")) {
            run_energies = config_json.at("run_energies").get<std::vector<double>>();
            found_run_energies = true;
        } else if (config_json.contains("run_labels")) {
            run_labels = config_json.at("run_labels").get<std::vector<std::string>>();
        } else {
            spdlog::error("Config json file {} must contain either run_energies or run_labels!", script_input_file);
            return 1;
        }
        if (config_json.contains("enable_gaussian_fit")) {
            enable_gaussian_fit = config_json.at("enable_gaussian_fit").get<bool>();
        }
        if (config_json.contains("use_peak_integral")) {
            use_peak_integral = config_json.at("use_peak_integral").get<bool>();
        }
        beam_type = config_json.at("beam_type").get<std::string>();
        config_description = config_json.at("description").get<std::string>();
        if (config_json.contains("n_events")) {
            target_event_number = config_json.at("n_events").get<int>();
        }
    } catch (const json::exception& e) {
        spdlog::error("Failed to parse config json file {}: {}", script_input_file, e.what());
        return 1;
    }

    std::vector<std::string> run_files;

    if (run_numbers.empty()) {
        spdlog::error("No run numbers specified in config json file {}", script_input_file);
        return 1;
    }
    const std::string run_file_folder = "dump/314_ToA_Analysis_LUT/beamtests/";
    for (const auto& run_number : run_numbers) {
        std::string run_number_str = std::to_string(run_number);
        run_number_str.insert(0, 4 - run_number_str.length(), '0');
        std::string run_file = run_file_folder + "Run" + run_number_str + ".root";
        if (access(run_file.c_str(), F_OK) == -1) {
            spdlog::error("Run file {} does not exist!", run_file);
            return 1;
        }
        run_files.push_back(run_file);
    }

    if (found_run_energies){
        std::vector<size_t> indices(run_files.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&run_energies](size_t i1, size_t i2) {
            return run_energies[i1] > run_energies[i2];
        });
        std::vector<std::string> sorted_run_files;
        std::vector<double> sorted_run_energies;
        std::vector<int> sorted_run_numbers;
        for (const auto& index : indices) {
            sorted_run_files.push_back(run_files[index]);
            sorted_run_energies.push_back(run_energies[index]);
            sorted_run_numbers.push_back(run_numbers[index]);
        }
        run_files = sorted_run_files;
        run_energies = sorted_run_energies;
        run_numbers = sorted_run_numbers;
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

    std::vector<TH1D*> hist_adc_25ns_list;
    std::vector<TH1D*> hist_adc_10ns_list;
    TKey *key = nullptr;
    TObject *primitive = nullptr;
    for (int _file_index = 0; _file_index < run_files.size(); ++_file_index) {
        const auto& run_file = run_files[_file_index];
        TFile* file = TFile::Open(run_file.c_str(), "READ");
        if (!file || file->IsZombie()) {
            spdlog::error("Failed to open run file {}", run_file);
            return 1;
        }
        
        TCanvas* toa_filtered_adc_canvas = nullptr;
        TH1D* toa_filtered_25ns_adc_hist = nullptr;
        TH1D* toa_filtered_10ns_adc_hist = nullptr;

        TIter next_toa_canvas(file->GetListOfKeys());
        while ((key = dynamic_cast<TKey*>(next_toa_canvas()))) {
            const char* key_name = key->GetName();
            if (key_name && std::string(key_name) == "canvas_adc_sum_comparison") {
                auto toa_canvas_in_file = dynamic_cast<TCanvas*>(key->ReadObj());
                if (!toa_canvas_in_file) {
                    continue;
                }
                toa_filtered_adc_canvas = (TCanvas*)toa_canvas_in_file->Clone(Form("toa_filtered_adc_canvas_Run%d", run_numbers[_file_index]));
                break;
            }
        }

        TIter next_toa_primitive(toa_filtered_adc_canvas->GetListOfPrimitives());
        while ((primitive = next_toa_primitive())) {
            if (auto hist = dynamic_cast<TH1D*>(primitive)) {
                if (hist->GetName() == std::string("h1d_event_adc_25ns_central_window")) {
                    toa_filtered_25ns_adc_hist = (TH1D*)hist->Clone("toa_filtered_25ns_adc_hist");
                    toa_filtered_25ns_adc_hist->SetDirectory(nullptr);
                    spdlog::info("Cloned 25ns histogram");
                } else if (hist->GetName() == std::string("h1d_event_adc_10ns_central_window")) {
                    toa_filtered_10ns_adc_hist = (TH1D*)hist->Clone("toa_filtered_10ns_adc_hist");
                    toa_filtered_10ns_adc_hist->SetDirectory(nullptr);
                spdlog::info("Cloned 10ns histogram");
                }
            }
        }

        if (!toa_filtered_25ns_adc_hist || !toa_filtered_10ns_adc_hist) {
            spdlog::error("Failed to find required histograms in run file {}", run_file);
            return 1;
        }
        hist_adc_25ns_list.push_back(toa_filtered_25ns_adc_hist);
        hist_adc_10ns_list.push_back(toa_filtered_10ns_adc_hist);

        file->Close();
    }

    auto output_file = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_file->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }

    // create a canvas for all 25ns histograms
    std::vector<double> cb_fit_mean_values;
    std::vector<double> cb_fit_mean_err_sys_values;
    std::vector<double> cb_fit_mean_err_stat_values;
    std::vector<double> cb_fit_sigma_values;
    std::vector<double> cb_fit_sigma_err_sys_values;
    std::vector<double> cb_fit_sigma_err_stat_values;
    std::vector<double> cb_fit_resolution_values;
    std::vector<double> cb_fit_resolution_err_values;

    std::vector<double> fit_range_sigma_values = {
        2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5
    };
    std::vector<double> fit_range_offset_values = {
        0.0, 0.5, 1.0
    };

    TCanvas* canvas_25ns = new TCanvas("canvas_25ns", "25ns ADC Comparison", 800, 600);
    canvas_25ns->SetLeftMargin(0.15);
    canvas_25ns->SetBottomMargin(0.10);
    for (size_t i = 0; i < hist_adc_25ns_list.size(); ++i) {
        double energy = found_run_energies ? run_energies[i] : static_cast<double>(i);
        hist_adc_25ns_list[i]->SetLineColor(i + 1);
        hist_adc_25ns_list[i]->SetMarkerColor(i + 1);
        hist_adc_25ns_list[i]->SetMarkerStyle(20 + i);
        hist_adc_25ns_list[i]->SetMarkerSize(1);
        if (i == 0) {
            hist_adc_25ns_list[i]->SetTitle("25ns ADC Comparison;ADC Value;Entries");
            hist_adc_25ns_list[i]->Draw("E");
        } else {
            hist_adc_25ns_list[i]->Draw("E SAME");
        }
        double fit_result_mean, fit_result_mean_err_sys, fit_result_mean_err_stat;
        double fit_result_sigma, fit_result_sigma_err_sys, fit_result_sigma_err_stat;
        double fit_result_resolution, fit_result_resolution_err;
        TCanvas fit_canvas = TCanvas(Form("fit_canvas_%dGeV", static_cast<int>(energy)), Form("CB Fit for %.0f GeV", energy), 800, 600);
        fit_canvas.cd();
        TH1D *fit_hist = (TH1D*)hist_adc_25ns_list[i]->Clone(Form("fit_hist_%dGeV", static_cast<int>(energy)));
        fit_hist->Rebin(2);
        fit_hist->Draw("HIST");
        crystalball_fit_th1d(fit_canvas, *fit_hist, fit_range_sigma_values, fit_range_offset_values, color_list[i % color_list.size()], fit_result_mean, fit_result_mean_err_sys, fit_result_mean_err_stat, fit_result_sigma, fit_result_sigma_err_sys, fit_result_sigma_err_stat, fit_result_resolution, fit_result_resolution_err);

        spdlog::info("Fit results for 25ns ADC at {} GeV: mean = {} +/- {} (sys) +/- {} (stat), sigma = {} +/- {} (sys) +/- {} (stat), resolution = {} +/- {}", energy, fit_result_mean, fit_result_mean_err_sys, fit_result_mean_err_stat, fit_result_sigma, fit_result_sigma_err_sys, fit_result_sigma_err_stat, fit_result_resolution, fit_result_resolution_err);

        cb_fit_mean_values.push_back(fit_result_mean);
        cb_fit_mean_err_sys_values.push_back(fit_result_mean_err_sys);
        cb_fit_mean_err_stat_values.push_back(fit_result_mean_err_stat);
        cb_fit_sigma_values.push_back(fit_result_sigma);
        cb_fit_sigma_err_sys_values.push_back(fit_result_sigma_err_sys);
        cb_fit_sigma_err_stat_values.push_back(fit_result_sigma_err_stat);
        cb_fit_resolution_values.push_back(fit_result_resolution);
        cb_fit_resolution_err_values.push_back(fit_result_resolution_err);

        fit_canvas.Write();
        fit_canvas.Close();
    }
    TLegend* legend_25ns = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t i = 0; i < hist_adc_25ns_list.size(); ++i) {
        std::string legend_entry = found_run_energies ? std::to_string(run_energies[i]) : run_labels[i];
        legend_25ns->AddEntry(hist_adc_25ns_list[i], legend_entry.c_str(), "lep");
    }
    legend_25ns->Draw();
    canvas_25ns->Write();
    canvas_25ns->Close();


    std::vector<double> cb_10ns_fit_mean_values;
    std::vector<double> cb_10ns_fit_mean_err_sys_values;
    std::vector<double> cb_10ns_fit_mean_err_stat_values;
    std::vector<double> cb_10ns_fit_sigma_values;
    std::vector<double> cb_10ns_fit_sigma_err_sys_values;
    std::vector<double> cb_10ns_fit_sigma_err_stat_values;
    std::vector<double> cb_10ns_fit_resolution_values;
    std::vector<double> cb_10ns_fit_resolution_err_values;

    TCanvas* canvas_10ns = new TCanvas("canvas_10ns", "10ns ADC Comparison", 800, 600);
    canvas_10ns->SetLeftMargin(0.15);
    canvas_10ns->SetBottomMargin(0.10);
    for (size_t i = 0; i < hist_adc_10ns_list.size(); ++i) {
        double energy = found_run_energies ? run_energies[i] : static_cast<double>(i);
        hist_adc_10ns_list[i]->SetLineColor(i + 1);
        hist_adc_10ns_list[i]->SetMarkerColor(i + 1);
        hist_adc_10ns_list[i]->SetMarkerStyle(20 + i);
        hist_adc_10ns_list[i]->SetMarkerSize(1);
        if (i == 0) {
            hist_adc_10ns_list[i]->SetTitle("10ns ADC Comparison;ADC Value;Entries");
            hist_adc_10ns_list[i]->Draw("E");
        } else {
            hist_adc_10ns_list[i]->Draw("E SAME");
        }
        double fit_10ns_result_mean, fit_10ns_result_mean_err_sys, fit_10ns_result_mean_err_stat;
        double fit_10ns_result_sigma, fit_10ns_result_sigma_err_sys, fit_10ns_result_sigma_err_stat;
        double fit_10ns_result_resolution, fit_10ns_result_resolution_err;
        TCanvas fit_canvas = TCanvas(Form("fit_canvas_10ns_%dGeV", static_cast<int>(energy)), Form("CB Fit for 10ns ADC at %.0f GeV", energy), 800, 600);
        fit_canvas.cd();
        TH1D *fit_hist = (TH1D*)hist_adc_10ns_list[i]->Clone(Form("fit_hist_10ns_%dGeV", static_cast<int>(energy)));
        fit_hist->Rebin(2);
        fit_hist->Draw("HIST");
        crystalball_fit_th1d(fit_canvas, *fit_hist, fit_range_sigma_values, fit_range_offset_values, color_list[i % color_list.size()], fit_10ns_result_mean, fit_10ns_result_mean_err_sys, fit_10ns_result_mean_err_stat, fit_10ns_result_sigma, fit_10ns_result_sigma_err_sys, fit_10ns_result_sigma_err_stat, fit_10ns_result_resolution, fit_10ns_result_resolution_err);

        spdlog::info("Fit results for 10ns ADC at {} GeV: mean = {} +/- {} (sys) +/- {} (stat), sigma = {} +/- {} (sys) +/- {} (stat), resolution = {} +/- {}", energy, fit_10ns_result_mean, fit_10ns_result_mean_err_sys, fit_10ns_result_mean_err_stat, fit_10ns_result_sigma, fit_10ns_result_sigma_err_sys, fit_10ns_result_sigma_err_stat, fit_10ns_result_resolution, fit_10ns_result_resolution_err);

        cb_10ns_fit_mean_values.push_back(fit_10ns_result_mean);
        cb_10ns_fit_mean_err_sys_values.push_back(fit_10ns_result_mean_err_sys);
        cb_10ns_fit_mean_err_stat_values.push_back(fit_10ns_result_mean_err_stat);
        cb_10ns_fit_sigma_values.push_back(fit_10ns_result_sigma);
        cb_10ns_fit_sigma_err_sys_values.push_back(fit_10ns_result_sigma_err_sys);
        cb_10ns_fit_sigma_err_stat_values.push_back(fit_10ns_result_sigma_err_stat);
        cb_10ns_fit_resolution_values.push_back(fit_10ns_result_resolution);
        cb_10ns_fit_resolution_err_values.push_back(fit_10ns_result_resolution_err);

        fit_canvas.Write();
        fit_canvas.Close();
    }
    TLegend* legend_10ns = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t i = 0; i < hist_adc_10ns_list.size(); ++i) {
        std::string legend_entry = found_run_energies ? std::to_string(run_energies[i]) : run_labels[i];
        legend_10ns->AddEntry(hist_adc_10ns_list[i], legend_entry.c_str(), "lep");
    }
    legend_10ns->Draw();
    canvas_10ns->Write();
    canvas_10ns->Close();

    // draw the resolution as a function of energy
    TCanvas* canvas_resolution = new TCanvas("canvas_resolution", "Energy Resolution", 800, 600);
    TGraphErrors* graph_resolution = new TGraphErrors("graph_resolution");
    graph_resolution->SetName("graph_resolution");
    for (size_t i = 0; i < cb_fit_mean_values.size(); ++i) {
        double x = found_run_energies ? run_energies[i] : static_cast<double>(i);
        double y = cb_fit_resolution_values[i];
        double xerr = 0.03 * x;
        double yerr = cb_fit_resolution_err_values[i];
        graph_resolution->SetPoint(i, x, y);
        graph_resolution->SetPointError(i, xerr, yerr);
    }
    // set axis range
    graph_resolution->GetXaxis()->SetLimits(40, 400);
    graph_resolution->GetYaxis()->SetRangeUser(5, 25);
    graph_resolution->SetMarkerStyle(20);
    graph_resolution->SetMarkerColor(kRed);
    graph_resolution->SetLineColor(kRed);
    graph_resolution->SetTitle("Energy Resolution;Energy (GeV);Resolution (%)");
    graph_resolution->Draw("AP");

    TGraphErrors* graph_resolution_10ns = new TGraphErrors("graph_resolution_10ns");
    graph_resolution_10ns->SetName("graph_resolution_10ns");
    for (size_t i = 0; i < cb_10ns_fit_mean_values.size(); ++i) {
        double x = found_run_energies ? run_energies[i] : static_cast<double>(i);
        double y = cb_10ns_fit_resolution_values[i];
        double xerr = 0.03 * x;
        double yerr = cb_10ns_fit_resolution_err_values[i];
        graph_resolution_10ns->SetPoint(i, x, y);
        graph_resolution_10ns->SetPointError(i, xerr, yerr);
    }
    graph_resolution_10ns->SetMarkerStyle(21);
    graph_resolution_10ns->SetMarkerColor(kBlue);
    graph_resolution_10ns->SetLineColor(kBlue);
    graph_resolution_10ns->Draw("P SAME");
    canvas_resolution->Write();
    // save as a separate pdf
    std::string resolution_pdf_name = script_output_folder + "/energy_resolution.pdf";
    canvas_resolution->SaveAs(resolution_pdf_name.c_str());
    canvas_resolution->Close();

    output_file->Close();
    return 0;
}