#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "H2GCROC_Common.hxx"
#include "H2GCROC_ADC_Analysis.hxx"
#include "H2GCROC_Toolbox.hxx"
#include "TRandom3.h"
#include "TKey.h"
#include "CommonParams.hxx"

double resolution_func(double *x, double *par) {
    // par[0] = p0 (constant term)
    // par[1] = p1 (stochastic term)
    // par[2] = p2 (noise term)
    double E = x[0];
    if (E == 0) return 0;
    double res = std::sqrt( std::pow(par[0], 2) + std::pow(par[1]/std::sqrt(E), 2) + std::pow(par[2]/E, 2) );
    return res;
}

int main(int argc, char **argv) {
    gROOT->SetBatch(kTRUE);

    int  script_n_events = -1;
    bool script_verbose = false;
    std::string         script_input_file;
    std::string         script_output_file;
    std::string         script_name = __FILE__;
    const std::string   script_version = "0.1";
    std::string         script_output_folder;
    script_name = script_name.substr(script_name.find_last_of("/\\") + 1).substr(0, script_name.find_last_of("."));
    configure_logger(true);  // Enable verbose logging for debugging

    const int example_channel = CommonParams::example_channel;

    std::vector<int> color_list = {kBlue, kRed, kGreen+2, kMagenta+2, kCyan+2, kOrange+2, kViolet+2, kPink+2, kAzure+2, kSpring+2, kTeal+2};
    int color_MC_curve = kGreen+2;
    int color_TB_curve = kPink+2;
    int color_previous_curve = kBrown+2;

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

    const double energy_axis_min = 40;
    const double energy_axis_max = 400;

    const double resolution_axis_min = 0;
    const double resolution_axis_max = 40;

    std::vector<double> previous_energies = {60, 80, 100, 150, 200, 250, 300, 350};
    std::vector<double> previous_resolutions = {0.202903,0.193606,0.180118,0.157828,0.143362,0.135984,0.128683,0.124429};
    std::vector<double> previous_resolution_errors = {0.0223685347754385,
                0.0276764080400618,
                0.0256984028686609,
                0.0235984254347615,
                0.016989855708628,
                0.0151618517338747,
                0.0109004483394033,
                0.0102402795860269
    };

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
    const std::string run_file_folder = "dump/501_ADC_Analysis/beamtests/";
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

    std::vector<TH1D*> hist_adc_sum_list;
    hist_adc_sum_list.reserve(run_files.size());
    TKey *key = nullptr;
    TObject *primitive = nullptr;
    for (int _file_index = 0; _file_index < run_files.size(); ++_file_index) {
        const auto& run_file = run_files[_file_index];
        TFile* file = TFile::Open(run_file.c_str(), "READ");
        if (!file || file->IsZombie()) {
            spdlog::error("Failed to open run file {}", run_file);
            return 1;
        }

        TCanvas* temp_canvas = nullptr;
        TH1D* hist_adc_sum = nullptr;

        TIter next_toa_canvas(file->GetListOfKeys());
        while ((key = dynamic_cast<TKey*>(next_toa_canvas()))) {
            const char* key_name = key->GetName();
             if (key_name && std::string(key_name) == "canvas") {
                auto toa_canvas_in_file = dynamic_cast<TCanvas*>(key->ReadObj());
                if (!toa_canvas_in_file) {
                    continue;
                }
                temp_canvas = (TCanvas*)toa_canvas_in_file->Clone(Form("canvas_adc_sum_run%d", run_numbers[_file_index]));
                break;
            }
        }
        TIter next_toa_primitive(temp_canvas->GetListOfPrimitives());
        while ((primitive = next_toa_primitive())) {
            if (auto hist = dynamic_cast<TH1D*>(primitive)) {
                if (hist->GetName() == std::string("adc_peak_sum_hist")){
                    hist_adc_sum = (TH1D*)hist->Clone(Form("hist_adc_sum_run%d", run_numbers[_file_index]));
                    hist_adc_sum->SetDirectory(nullptr); // Detach histogram from file
                    hist_adc_sum_list.push_back(hist_adc_sum);
                    break;
                }
            }
        }

        if (!hist_adc_sum) {
            spdlog::error("Failed to find histogram adc_peak_sum_hist in run file {}", run_file);
            return 1;
        }

        file->Close();

    } // end of file loop
        
    auto output_file = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_file->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }
    output_file->cd();

    const std::vector<double> fit_range_sigmas = {2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5}; // in sigma
    const std::vector<double> fit_range_offsets = {-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5}; // in sigma

    const std::string canvas_title = "FoCal-H Prototype 3";
    const std::string canvas_subtitle = "SPS H2 Test Beam 2025";

    std::vector<double> mix_mean_values;
    std::vector<double> mix_mean_errors_stat;
    std::vector<double> mix_mean_errors_sys;
    std::vector<double> mix_sigma_values;
    std::vector<double> mix_sigma_errors_stat;
    std::vector<double> mix_sigma_errors_sys;

    int cb_fit_color    = kPink+2;
    int gaus_fit_color  = kCyan+2;

    TCanvas* canvas_adc_sum = new TCanvas("canvas_adc_sum", "ADC Sum", 800, 600);
    TLegend* legend_adc_sum = new TLegend(0.75, 0.6, 0.9, 0.9);
    legend_adc_sum->SetBorderSize(0);
    legend_adc_sum->SetFillStyle(0);
    legend_adc_sum->SetTextSize(0.025);
    // canvas_adc_sum->SetLeftMargin(0.15);
    // canvas_adc_sum->SetBottomMargin(0.10);
    double max_y = 0;
    double max_x = 0;
    // go through the histograms once to find the global max for setting the same axis range
    for (const auto& hist_adc_sum : hist_adc_sum_list) {
        // normalize
        hist_adc_sum->Scale(1.0 / hist_adc_sum->Integral());
        double current_max_y = hist_adc_sum->GetMaximum();
        double current_max_x = hist_adc_sum->GetXaxis()->GetXmax();
        if (current_max_y > max_y) {
            max_y = current_max_y;
        }        if (current_max_x > max_x) {
            max_x = current_max_x;
        }
    }
    for (size_t i = 0; i < hist_adc_sum_list.size(); ++i) {
        TCanvas* canvas_hist = new TCanvas(Form("canvas_adc_sum_run%d", run_numbers[i]), Form("ADC Sum Run %d", run_numbers[i]), 800, 600);
        canvas_hist->cd();
        TH1D* hist_adc_sum = (TH1D*)hist_adc_sum_list[i]->Clone(Form("hist_adc_sum_run%d", run_numbers[i]));
        hist_adc_sum->SetDirectory(nullptr); // Detach histogram from file
        hist_adc_sum->SetMaximum(max_y*1.3);
        hist_adc_sum->SetAxisRange(0, max_x*0.8, "X");
        hist_adc_sum->Draw();  // Draw the histogram on the canvas
        int color = color_list[i % color_list.size()];

        double gaus_fit_mean_nom, gaus_fit_mean_err_stat, gaus_fit_mean_err_syst;
        double gaus_fit_sigma_nom, gaus_fit_sigma_err_stat, gaus_fit_sigma_err_syst;
        TF1* nominal_fit_func = nullptr;
        gaussian_fit_th1d_nom(
            *canvas_hist, *hist_adc_sum, fit_range_sigmas, fit_range_offsets, gaus_fit_color,
            gaus_fit_mean_nom, gaus_fit_sigma_nom, gaus_fit_mean_err_stat, gaus_fit_sigma_err_stat,
            gaus_fit_mean_err_syst, gaus_fit_sigma_err_syst, nominal_fit_func
        );
        
        double cb_fit_mean_nom, cb_fit_sigma_nom;
        double cb_fit_mean_err_stat, cb_fit_sigma_err_stat;
        double cb_fit_mean_err_syst, cb_fit_sigma_err_syst;
        TF1* cb_nominal_fit_func = nullptr;
        crystalball_fit_th1d_nom(
            *canvas_hist, *hist_adc_sum, fit_range_sigmas, fit_range_offsets, cb_fit_color,
            cb_fit_mean_nom, cb_fit_sigma_nom, cb_fit_mean_err_stat, cb_fit_sigma_err_stat,
            cb_fit_mean_err_syst, cb_fit_sigma_err_syst, cb_nominal_fit_func
        );

        double cb_gaus_mix_mean_nom, cb_gaus_mix_mean_err_stat, cb_gaus_mix_mean_err_syst;
        double cb_gaus_mix_sigma_nom, cb_gaus_mix_sigma_err_stat, cb_gaus_mix_sigma_err_syst;
        cb_gaus_combine(cb_fit_mean_nom, cb_fit_mean_err_stat, cb_fit_mean_err_syst, gaus_fit_mean_nom, gaus_fit_mean_err_stat, gaus_fit_mean_err_syst, cb_gaus_mix_mean_nom, cb_gaus_mix_mean_err_stat, cb_gaus_mix_mean_err_syst);
        cb_gaus_combine(cb_fit_sigma_nom, cb_fit_sigma_err_stat, cb_fit_sigma_err_syst, gaus_fit_sigma_nom, gaus_fit_sigma_err_stat, gaus_fit_sigma_err_syst, cb_gaus_mix_sigma_nom, cb_gaus_mix_sigma_err_stat, cb_gaus_mix_sigma_err_syst);
        
        mix_mean_values.push_back(cb_gaus_mix_mean_nom);
        mix_mean_errors_stat.push_back(cb_gaus_mix_mean_err_stat);
        mix_mean_errors_sys.push_back(cb_gaus_mix_mean_err_syst);
        mix_sigma_values.push_back(cb_gaus_mix_sigma_nom);
        mix_sigma_errors_stat.push_back(cb_gaus_mix_sigma_err_stat);
        mix_sigma_errors_sys.push_back(cb_gaus_mix_sigma_err_syst);

        canvas_hist->Modified();
        canvas_hist->Update();
        canvas_hist->Write();
        // make a clone for canvas_adc_sum
        TH1D* hist_adc_sum_clone = (TH1D*)hist_adc_sum->Clone(Form("hist_adc_sum_run%d", run_numbers[i]));
        hist_adc_sum_clone->SetDirectory(nullptr);
        // canvas_hist->Close();
        hist_adc_sum_clone->SetLineColor(color);
        hist_adc_sum_clone->SetMarkerColor(color);
        canvas_adc_sum->cd();
        if (i == 0) {
            hist_adc_sum_clone->GetYaxis()->SetTitle("Normalized Count");
            hist_adc_sum_clone->GetXaxis()->SetTitle("ADC Peak Sum");
            format_1d_hist_canvas(canvas_adc_sum, hist_adc_sum_clone, color, canvas_title.c_str(), canvas_subtitle.c_str(), "Hadron Energy Scan", false);
        } else {
            hist_adc_sum_clone->Draw("hist same");
        }
        legend_adc_sum->AddEntry(hist_adc_sum_clone, Form("%d GeV", (int)run_energies[i]), "l");
    }
    legend_adc_sum->Draw();
    canvas_adc_sum->Write();
    // save to a separate pdf file
    canvas_adc_sum->SaveAs((script_output_folder + "/" + script_name + "_adc_sum_comparison.pdf").c_str());
    canvas_adc_sum->Close();

    // draw the linearity plot
    TCanvas* canvas_linearity = new TCanvas("canvas_linearity", "Linearity", 800, 600);
    TLegend* legend_linearity = new TLegend(0.5, 0.1, 0.9, 0.3);
    legend_linearity->SetBorderSize(0);
    legend_linearity->SetFillStyle(0);
    legend_linearity->SetTextSize(0.025);
    // canvas_linearity->SetLeftMargin(0.15);
    // canvas_linearity->SetBottomMargin(0.10);
    TGraphErrors* graph_linearity = new TGraphErrors(run_energies.size());
    double linearity_y_scale_factor = 0.001;
    for (size_t i = 0; i < run_energies.size(); ++i) {
        double energy = run_energies[i];
        double mean = mix_mean_values[i];
        double mean_err_stat = mix_mean_errors_stat[i];
        double mean_err_syst = mix_mean_errors_sys[i];
        double mean_scaled = mean * linearity_y_scale_factor;
        double mean_err_stat_scaled = mean_err_stat * linearity_y_scale_factor;
        double mean_err_syst_scaled = mean_err_syst * linearity_y_scale_factor;
        double energy_err = 0.03 * energy; // assume 3% energy uncertainty
        // graph_linearity->SetPoint(i, energy, mean);
        // graph_linearity->SetPointError(i, energy_err, std::sqrt(mean_err_stat * mean_err_stat + mean_err_syst * mean_err_syst));
        graph_linearity->SetPoint(i, energy, mean_scaled);
        graph_linearity->SetPointError(i, energy_err, std::sqrt(mean_err_stat_scaled * mean_err_stat_scaled + mean_err_syst_scaled * mean_err_syst_scaled));
    }
    graph_linearity->SetMarkerStyle(20);
    graph_linearity->SetMarkerSize(1.0);
    graph_linearity->SetMarkerColor(color_TB_curve);
    graph_linearity->SetLineColor(color_TB_curve);
    if (linearity_y_scale_factor == 0.001) {
        graph_linearity->GetYaxis()->SetTitle("Mean ADC Sum (1000 ADC)");
    } else {
        graph_linearity->GetYaxis()->SetTitle("Mean ADC Sum");
    }
    graph_linearity->GetXaxis()->SetTitle("Beam Energy (GeV)");
    graph_linearity->GetXaxis()->SetRangeUser(energy_axis_min, energy_axis_max);
    format_tgrapherr_canvas(canvas_linearity, graph_linearity, color_TB_curve, canvas_title.c_str(), canvas_subtitle.c_str(), "Hadron Energy Scan");
    legend_linearity->AddEntry(graph_linearity, "Prototype 3", "pe");

    // do a linear fit to the linearity graph
    TF1* linear_fit = new TF1("linear_fit", "[0] + [1]*x", energy_axis_min, energy_axis_max);
    graph_linearity->Fit(linear_fit, "RNQ");
    double linear_fit_a = linear_fit->GetParameter(0);
    double linear_fit_b = linear_fit->GetParameter(1);
    double linear_fit_a_err = linear_fit->GetParError(0);
    double linear_fit_b_err = linear_fit->GetParError(1);
    spdlog::info("Linearity fit results: a = {} +/- {}, b = {} +/- {}", linear_fit_a, linear_fit_a_err, linear_fit_b, linear_fit_b_err);
    linear_fit->SetLineColor(color_TB_curve);
    linear_fit->Draw("same");

    // legend_linearity->AddEntry((TObject*)linear_fit, Form("Linear Fit: y = (%.2f #pm %.2f) + (%.2f #pm %.2f)x", linear_fit_a, linear_fit_a_err, linear_fit_b, linear_fit_b_err), "l");
    std::string linear_fit_label = Form("ADC Sum = %.2f + %.2f * E_{beam}", linear_fit_a, linear_fit_b);
    legend_linearity->AddEntry((TObject*)linear_fit, linear_fit_label.c_str(), "l");
    legend_linearity->Draw();
    canvas_linearity->Write();
    // save to a separate pdf file
    canvas_linearity->SaveAs((script_output_folder + "/" + script_name + "_linearity.pdf").c_str());
    canvas_linearity->Close();

    // draw the resolution plot
    TCanvas* canvas_resolution = new TCanvas("canvas_resolution", "Resolution", 800, 600);
    TLegend* legend_resolution = new TLegend(0.5, 0.55, 0.9, 0.9);
    legend_resolution->SetBorderSize(0);
    legend_resolution->SetFillStyle(0);
    legend_resolution->SetTextSize(0.025);
    // canvas_resolution->SetLeftMargin(0.15);
    // canvas_resolution->SetBottomMargin(0.10);
    TGraphErrors* graph_resolution = new TGraphErrors(run_energies.size());
    for (size_t i = 0; i < run_energies.size(); ++i) {
        double energy = run_energies[i];
        double sigma = mix_sigma_values[i];
        double sigma_err_stat = mix_sigma_errors_stat[i];
        double sigma_err_syst = mix_sigma_errors_sys[i];
        double mean = mix_mean_values[i];
        double mean_err_stat = mix_mean_errors_stat[i];
        double mean_err_syst = mix_mean_errors_sys[i];
        double resolution, resolution_err_stat, resolution_err_syst;
        resolution_calculator(mean, mean_err_stat, mean_err_syst, sigma, sigma_err_stat, sigma_err_syst, resolution, resolution_err_stat, resolution_err_syst);
        double energy_err = 0.03 * energy; // assume 3% energy uncertainty
        graph_resolution->SetPoint(i, energy, resolution);
        graph_resolution->SetPointError(i, energy_err, std::sqrt(resolution_err_stat * resolution_err_stat + resolution_err_syst * resolution_err_syst));
    }
    graph_resolution->SetMarkerStyle(20);
    graph_resolution->SetMarkerSize(1.0);
    graph_resolution->SetMarkerColor(color_TB_curve);
    graph_resolution->SetLineColor(color_TB_curve);
    graph_resolution->GetYaxis()->SetTitle("#sigma_{E} / E_{beam} (%)");
    graph_resolution->GetYaxis()->SetRangeUser(resolution_axis_min, resolution_axis_max);
    graph_resolution->GetXaxis()->SetTitle("Beam Energy (GeV)");
    graph_resolution->GetXaxis()->SetRangeUser(energy_axis_min, energy_axis_max);
    format_tgrapherr_canvas(canvas_resolution, graph_resolution, color_TB_curve, canvas_title.c_str(), canvas_subtitle.c_str(), "Hadron Energy Scan");
    legend_resolution->AddEntry(graph_resolution, "Prototype 3", "pe");

    // fit the resolution graph with the resolution function
    TF1* resolution_fit = new TF1("resolution_fit", resolution_func, energy_axis_min, energy_axis_max, 3);
    resolution_fit->SetParNames("p0", "p1", "p2");
    resolution_fit->SetLineColor(color_TB_curve);
    resolution_fit->SetParameters(5.0, 100.0, 0.0);
    // resolution_fit->FixParameter(2, 0.0); // fix noise term to 0 for better fit stability, can be released if needed
    graph_resolution->Fit(resolution_fit, "RNQ");
    double p0 = std::abs(resolution_fit->GetParameter(0));
    double p1 = std::abs(resolution_fit->GetParameter(1));
    double p2 = std::abs(resolution_fit->GetParameter(2));
    double p0_err = resolution_fit->GetParError(0);
    double p1_err = resolution_fit->GetParError(1);
    double p2_err = resolution_fit->GetParError(2);

    spdlog::info("Resolution fit results: p0 = {} +/- {}, p1 = {} +/- {}, p2 = {} +/- {}", p0, p0_err, p1, p1_err, p2, p2_err);

    resolution_fit->Draw("same");
    std::string resolution_fit_label = Form("#sigma/E = %.2f/#sqrt{E} #oplus %.2f / E #oplus %.2f", p1, p2, p0);
    legend_resolution->AddEntry((TObject*)resolution_fit, resolution_fit_label.c_str(), "l");

    TGraphErrors* graph_previous = new TGraphErrors(previous_energies.size());
    for (size_t i = 0; i < previous_energies.size(); ++i) {
        graph_previous->SetPoint(i, previous_energies[i], previous_resolutions[i]*100); // convert to percentage
        graph_previous->SetPointError(i, 0.03*previous_energies[i], previous_resolution_errors[i]*100); // convert to percentage
    }
    graph_previous->SetMarkerStyle(21);
    graph_previous->SetMarkerSize(1.0);
    graph_previous->SetMarkerColor(color_previous_curve);
    graph_previous->SetLineColor(color_previous_curve);
    graph_previous->Draw("PE same");
    legend_resolution->AddEntry(graph_previous, "Prototype 2", "pe");

    TF1* previous_fit = new TF1("previous_fit", resolution_func, energy_axis_min, energy_axis_max, 3);
    previous_fit->SetParNames("p0", "p1", "p2");
    previous_fit->SetLineColor(color_previous_curve);
    previous_fit->SetParameters(5.0, 100.0, 0.0);
    // previous_fit->FixParameter(2, 0.0); // fix noise term to 0 for better fit stability, can be released if needed
    graph_previous->Fit(previous_fit, "RNQ");
    double p0_prev = std::abs(previous_fit->GetParameter(0));
    double p1_prev = std::abs(previous_fit->GetParameter(1));
    double p2_prev = std::abs(previous_fit->GetParameter(2));
    double p0_prev_err = previous_fit->GetParError(0);
    double p1_prev_err = previous_fit->GetParError(1);
    double p2_prev_err = previous_fit->GetParError(2);

    spdlog::info("Previous resolution fit results: p0 = {} +/- {}, p1 = {} +/- {}, p2 = {} +/- {}", p0_prev, p0_prev_err, p1_prev, p1_prev_err, p2_prev, p2_prev_err);

    previous_fit->Draw("same");
    std::string resolution_fit_previous_label = Form("#sigma/E = %.2f/#sqrt{E} #oplus %.2f / E #oplus %.2f", p1_prev, p2_prev, p0_prev);
    legend_resolution->AddEntry((TObject*)previous_fit, resolution_fit_previous_label.c_str(), "l");

    legend_resolution->Draw();

    canvas_resolution->Write();
    // save to a separate pdf file
    canvas_resolution->SaveAs((script_output_folder + "/" + script_name + "_resolution.pdf").c_str());
    canvas_resolution->Close();


    output_file->Close();

    return 0;
}