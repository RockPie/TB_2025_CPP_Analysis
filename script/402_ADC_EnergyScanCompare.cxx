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
    std::vector<int> color_list = {kRed, kBlue, kGreen+2, kMagenta+2, kCyan+2, kOrange+2, kViolet+2};

    std::string output_file_name = "dump/402_ADC_EnergyScanCompare/ADC_comparison.root";
    std::string simulation_adc_folder = "data/simulation/Config24/";
    std::string tb_adc_file = "dump/305_ADC_Fit_Compare/scan_0_analysis.root";
    std::string adc_toa_file = "dump/322_ToA_ADC_Compare/scan_0_analysis.root";

    TFile *tb_adc_root_file = new TFile(tb_adc_file.c_str(), "READ");
    if (tb_adc_root_file->IsZombie()) {
        spdlog::error("Failed to open beam test ADC file {}", tb_adc_file);
        return 1;
    }
    // search for TCanvas with the name canvas_resolution_to_previous
    TCanvas *tb_adc_canvas = nullptr;
    TIter next(tb_adc_root_file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        if (std::string(key->GetClassName()) == "TCanvas" && std::string(key->GetName()) == "canvas_resolution_to_previous") {
            tb_adc_canvas = (TCanvas*)key->ReadObj();
            break;
        }
    }
    if (!tb_adc_canvas) {
        spdlog::error("Failed to find canvas named canvas_resolution_to_previous in beam test ADC file {}", tb_adc_file);
        return 1;
    }
    // print all the primitives in the canvas
    spdlog::info("Primitives in beam test ADC canvas:");
    TIter next_primitive(tb_adc_canvas->GetListOfPrimitives());
    TObject *primitive;
    TGraphErrors *cb_fit_graph = nullptr;
    TGraphErrors *previous_prototype_graph = nullptr;
    while ((primitive = (TObject*)next_primitive())) {
        spdlog::info("Primitive: {}, Class: {}", primitive->GetName(), primitive->ClassName());
        // look for graph_resolution_cb and graph_resolution_previous
        if (std::string(primitive->GetName()) == "graph_resolution_cb") {
            auto cb_fit_graph_input = (TGraphErrors*)primitive;
            cb_fit_graph = (TGraphErrors*)cb_fit_graph_input->Clone("cb_fit_graph");
            spdlog::info("Found graph named graph_resolution_cb");
        } else if (std::string(primitive->GetName()) == "graph_resolution_previous") {
            auto previous_prototype_graph_input = (TGraphErrors*)primitive;
            previous_prototype_graph = (TGraphErrors*)previous_prototype_graph_input->Clone("previous_prototype_graph");
            spdlog::info("Found graph named graph_resolution_previous");
        }
    }

    TCanvas *tb_adc_linear_canvas = nullptr;
    TIter next_canvas(tb_adc_root_file->GetListOfKeys());
    while ((key = (TKey*)next_canvas())) {
        if (std::string(key->GetClassName()) == "TCanvas" && std::string(key->GetName()) == "canvas_linearity") {
            tb_adc_linear_canvas = (TCanvas*)key->ReadObj();
            break;
        }
    }
    if (!tb_adc_linear_canvas) {
        spdlog::error("Failed to find canvas named canvas_linearity in beam test ADC file {}", tb_adc_file);
        return 1;
    }

    TGraphErrors *linearity_graph_cb = nullptr;
    TIter next_linear_primitive(tb_adc_linear_canvas->GetListOfPrimitives());
    while ((primitive = (TObject*)next_linear_primitive())) {
        spdlog::info("Primitive: {}, Class: {}", primitive->GetName(), primitive->ClassName());
        if (std::string(primitive->GetName()) == "graph_linearity_cb") {
            auto linearity_graph_cb_input = (TGraphErrors*)primitive;
            linearity_graph_cb = (TGraphErrors*)linearity_graph_cb_input->Clone("linearity_graph_cb");
            spdlog::info("Found graph named graph_linearity_cb");
        }
    }

    tb_adc_root_file->Close();

    TFile *adc_toa_root_file = new TFile(adc_toa_file.c_str(), "READ");
    if (adc_toa_root_file->IsZombie()) {
        spdlog::error("Failed to open ADC ToA file {}", adc_toa_file);
        return 1;
    }
    TCanvas *adc_toa_canvas = nullptr;
    spdlog::info("Keys in ADC ToA file:");
    TIter next_adc_toa(adc_toa_root_file->GetListOfKeys());
    while ((key = (TKey*)next_adc_toa())) {
        spdlog::info("Key: {}, Class: {}", key->GetName(), key->GetClassName());
        if (std::string(key->GetClassName()) == "TCanvas" && std::string(key->GetName()) == "canvas_resolution") {
            auto adc_toa_canvas_input = (TCanvas*)key->ReadObj();
            adc_toa_canvas = (TCanvas*)adc_toa_canvas_input->Clone("adc_toa_canvas");
            break;
        }
    }
    if (!adc_toa_canvas) {
        spdlog::error("Failed to find canvas named canvas_resolution in ADC ToA file {}", adc_toa_file);
        return 1;
    }
    TGraphErrors *toa_cb_25ns_fit_graph = nullptr;
    TGraphErrors *toa_cb_10ns_fit_graph = nullptr;
    TIter next_toa_primitive(adc_toa_canvas->GetListOfPrimitives());
    while ((primitive = (TObject*)next_toa_primitive())) {
        spdlog::info("Primitive: {}, Class: {}", primitive->GetName(), primitive->ClassName());
        if (std::string(primitive->GetName()) == "graph_resolution") {
            auto toa_cb_25ns_fit_graph_input = (TGraphErrors*)primitive;
            toa_cb_25ns_fit_graph = (TGraphErrors*)toa_cb_25ns_fit_graph_input->Clone("toa_cb_25ns_fit_graph");
            spdlog::info("Found graph named graph_resolution_cb in ADC ToA file, cloned to toa_cb_25ns_fit_graph");
        } else if (std::string(primitive->GetName()) == "graph_resolution_10ns") {
            auto toa_cb_10ns_fit_graph_input = (TGraphErrors*)primitive;
            toa_cb_10ns_fit_graph = (TGraphErrors*)toa_cb_10ns_fit_graph_input->Clone("toa_cb_10ns_fit_graph");
            spdlog::info("Found graph named graph_resolution_cb_10ns in ADC ToA file, cloned to toa_cb_10ns_fit_graph");
        }
    }
    adc_toa_root_file->Close();

    // get all files ending with .root in the simulation_adc_folder
    std::vector<std::string> simulation_adc_files;
    std::vector<double> simulation_file_beam_energies;
    // file name is like data/simulation/Config40/10_1000_250-Config40.root
    for (const auto& entry : std::filesystem::directory_iterator(simulation_adc_folder)) {
        if (entry.path().extension() == ".root") {
            simulation_adc_files.push_back(entry.path().string());
            std::string filename = entry.path().filename().string();
            // extract the energy label from the file name, which is after the second _ and before the -
            std::string energy_label = "";
            size_t first_underscore = filename.find_first_of("_");
            size_t second_underscore = filename.find_first_of("_", first_underscore + 1);
            size_t dash = filename.find_first_of("-");
            if (first_underscore != std::string::npos && second_underscore != std::string::npos && dash != std::string::npos && second_underscore < dash) {
                energy_label = filename.substr(second_underscore + 1, dash - second_underscore - 1);
            } else {
                spdlog::warn("Failed to extract energy label from file name {}, expected format like 10_1000_250-Config40.root", filename);
            }
            try {
                double energy = std::stod(energy_label);
                simulation_file_beam_energies.push_back(energy);
                spdlog::info("Found simulation ADC file {} with extracted energy {}", entry.path().string(), energy);
            } catch (const std::exception& e) {
                spdlog::warn("Failed to extract energy from file name {}, error: {}", filename, e.what());
                simulation_file_beam_energies.push_back(-1); // Use -1 to indicate unknown energy
            }
        }
    }
    if (simulation_adc_files.empty()) {
        spdlog::error("No .root files found in simulation ADC folder {}", simulation_adc_folder);
        return 1;
    } else {
        spdlog::info("Found {} .root files in simulation ADC folder {}", simulation_adc_files.size(), simulation_adc_folder);
        // print all the energies
        spdlog::info("Extracted energies from simulation file names:");
        for (size_t i = 0; i < simulation_file_beam_energies.size(); ++i) {
            spdlog::info("File: {}, Energy: {}", simulation_adc_files[i], simulation_file_beam_energies[i]);
        }
    }

    std::vector<TH1D*> simulation_adc_hists;
    for (size_t i = 0; i < simulation_adc_files.size(); ++i) {
        const std::string& sim_file = simulation_adc_files[i];
        double energy = simulation_file_beam_energies[i];
        spdlog::info("Processing simulation ADC file {} with energy {}", sim_file, energy);
        // Here you can add code to read the ADC histogram from the simulation file and compare it with the beam test data
        TFile *sim_adc_file = new TFile(sim_file.c_str(), "READ");
        if (sim_adc_file->IsZombie()) {
            spdlog::error("Failed to open simulation ADC file {}", sim_file);
            continue; // Skip this file and move to the next one
        }
        TH1D *sim_adc_hist = nullptr;
        if (sim_adc_file->cd("Photons")){
            TObject *adc_obj = gDirectory->Get("ADC");
            if (adc_obj){
                TH1D *hist = dynamic_cast<TH1D*>(adc_obj);
                sim_adc_hist = (TH1D*)hist->Clone(Form("sim_adc_hist_%dGeV", static_cast<int>(energy)));
                sim_adc_hist->SetDirectory(nullptr); // Detach histogram from file
                simulation_adc_hists.push_back(sim_adc_hist);
            }
        }
        sim_adc_file->Close();
    }

    std::string folder_path = simulation_adc_folder;
    // Remove trailing slash if present
    if (!folder_path.empty() && (folder_path.back() == '/' || folder_path.back() == '\\')) {
        folder_path.pop_back();
    }

    // create the output file and save the histograms
    TFile *output_file = new TFile(output_file_name.c_str(), "RECREATE");
    if (output_file->IsZombie()) {
        spdlog::error("Failed to create output file {}", output_file_name);
        return 1;
    }
    TCanvas *comparison_canvas = new TCanvas("comparison_canvas", "ADC Comparison", 800, 600);
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

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

    for (size_t i = 0; i < simulation_adc_hists.size(); ++i) {
        TH1D *sim_adc_hist = simulation_adc_hists[i];
        double energy = simulation_file_beam_energies[i];
        // normalize
        sim_adc_hist->Scale(1.0 / sim_adc_hist->Integral());
        // set max y
        sim_adc_hist->SetMaximum(sim_adc_hist->GetMaximum() * 1.4);
        if (i == 0)
            format_1d_hist_canvas(comparison_canvas, sim_adc_hist, color_list[i % color_list.size()], "H2GCROC Beam Test", "MC Simulation Comparison", Form("Simulated ADC at %.0f GeV", energy), false);
        else {
            sim_adc_hist->SetLineColor(color_list[i % color_list.size()]); // Use different colors for different energies
            sim_adc_hist->SetLineWidth(2);
            sim_adc_hist->Draw("HIST SAME");
        }
        // ! === Do the fitting here ===
        double fit_result_mean, fit_result_mean_err_sys, fit_result_mean_err_stat;
        double fit_result_sigma, fit_result_sigma_err_sys, fit_result_sigma_err_stat;
        double fit_result_resolution, fit_result_resolution_err;
        TCanvas fit_canvas = TCanvas(Form("fit_canvas_%dGeV", static_cast<int>(energy)), Form("CB Fit for %.0f GeV", energy), 800, 600);
        fit_canvas.cd();
        // copy the histogram to the fit canvas
        TH1D *fit_hist = (TH1D*)sim_adc_hist->Clone(Form("fit_hist_%dGeV", static_cast<int>(energy)));
        fit_hist->Rebin(2);
        fit_hist->Draw("HIST");
        crystalball_fit_th1d(fit_canvas, *fit_hist, fit_range_sigma_values, fit_range_offset_values, color_list[i % color_list.size()], fit_result_mean, fit_result_mean_err_sys, fit_result_mean_err_stat, fit_result_sigma, fit_result_sigma_err_sys, fit_result_sigma_err_stat, fit_result_resolution, fit_result_resolution_err);

        spdlog::info("Fit results for simulation ADC at {} GeV: mean = {} +/- {} (sys) +/- {} (stat), sigma = {} +/- {} (sys) +/- {} (stat), resolution = {} +/- {}", energy, fit_result_mean, fit_result_mean_err_sys, fit_result_mean_err_stat, fit_result_sigma, fit_result_sigma_err_sys, fit_result_sigma_err_stat, fit_result_resolution, fit_result_resolution_err);

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

        legend->AddEntry(sim_adc_hist, Form("Simulation ADC (%.0f GeV)", energy), "l");
        sim_adc_hist->Write();
    }
    legend->Draw();
    comparison_canvas->Modified();
    comparison_canvas->Update();
    comparison_canvas->Write();
    comparison_canvas->SaveAs("dump/402_ADC_EnergyScanCompare/ADC_comparison.pdf");
    comparison_canvas->Close();

    // draw the mean as a function of energy
    TCanvas mean_canvas("mean_canvas", "ADC Mean", 800, 600);
    TLegend legend_mean(0.6, 0.7, 0.9, 0.9);
    legend_mean.SetBorderSize(0);
    legend_mean.SetFillStyle(0);
    TGraphErrors mean_graph(cb_fit_mean_values.size());
    for (size_t i = 0; i < cb_fit_mean_values.size(); ++i) {
        double beam_energy = simulation_file_beam_energies[i];
        mean_graph.SetPoint(i, beam_energy, cb_fit_mean_values[i]);
        mean_graph.SetPointError(i, 0.03 * beam_energy, cb_fit_mean_err_stat_values[i]); // Use statistical error for the error bars
    }
    mean_graph.SetMarkerStyle(20);
    mean_graph.SetMarkerSize(1);
    mean_graph.SetLineColor(kBlue);
    // set axis range    mean_graph.GetXaxis()->SetLimits(40, 400);
    mean_graph.GetYaxis()->SetRangeUser(0, *std::max_element(cb_fit_mean_values.begin(), cb_fit_mean_values.end()) * 1.4);
    mean_graph.SetTitle("ADC Mean;Energy (GeV);ADC Mean");
    // set y axis range
    mean_graph.GetYaxis()->SetRangeUser(0, 650000);
    mean_graph.Draw("AP");
    std::string sim_data_label2 = "Simulation " + folder_path.substr(folder_path.find_last_of("/\\") + 1);
    legend_mean.AddEntry(&mean_graph, sim_data_label2.c_str(), "p");
    if (linearity_graph_cb) {
        linearity_graph_cb->SetMarkerStyle(21);
        linearity_graph_cb->SetMarkerSize(1);
        linearity_graph_cb->SetLineColor(kMagenta+2);
        linearity_graph_cb->SetMarkerColor(kMagenta+2);
        linearity_graph_cb->Draw("P SAME");
        legend_mean.AddEntry(linearity_graph_cb, "Beam Test CB Fit Mean", "p");
    }
    legend_mean.Draw();
    mean_canvas.Modified();
    mean_canvas.Update();
    mean_canvas.Write();
    mean_canvas.SaveAs("dump/402_ADC_EnergyScanCompare/adc_mean_comparison.pdf");
    mean_canvas.Close();

    // draw the resolution as a function of energy
    TCanvas resolution_canvas("resolution_canvas", "Energy Resolution", 800, 600);
    TLegend legend_resolution(0.6, 0.7, 0.9, 0.9);
    legend_resolution.SetBorderSize(0);
    legend_resolution.SetFillStyle(0);
    TGraphErrors resolution_graph(cb_fit_mean_values.size());
    for (size_t i = 0; i < cb_fit_mean_values.size(); ++i) {
        double beam_energy = simulation_file_beam_energies[i];
        resolution_graph.SetPoint(i, beam_energy, cb_fit_resolution_values[i]);
        resolution_graph.SetPointError(i, 0.03 * beam_energy, cb_fit_resolution_err_values[i]);
    }
    resolution_graph.SetMarkerStyle(20);
    resolution_graph.SetMarkerSize(1);
    resolution_graph.SetLineColor(kBlue);
    // set axis range
    resolution_graph.GetXaxis()->SetLimits(40, 400);
    resolution_graph.GetYaxis()->SetRangeUser(5, 25);
    resolution_graph.SetTitle("Energy Resolution;Energy (GeV);Resolution (%)");
    resolution_graph.Draw("AP");
    // add the base folder name as the label for the simulation data in the legend

    std::string sim_data_label = "Simulation " + folder_path.substr(folder_path.find_last_of("/\\") + 1);
    legend_resolution.AddEntry(&resolution_graph, sim_data_label.c_str(), "p");
    if (cb_fit_graph) {
        cb_fit_graph->SetMarkerStyle(21);
        cb_fit_graph->SetMarkerSize(1);
        cb_fit_graph->SetLineColor(kMagenta+2);
        cb_fit_graph->SetMarkerColor(kMagenta+2);
        cb_fit_graph->Draw("P SAME");
        legend_resolution.AddEntry(cb_fit_graph, "Beam Test CB Fit", "p");
    }
    if (previous_prototype_graph) {
        previous_prototype_graph->SetMarkerStyle(22);
        previous_prototype_graph->SetMarkerSize(1);
        previous_prototype_graph->SetLineColor(kRed);
        previous_prototype_graph->SetMarkerColor(kRed);
        previous_prototype_graph->Draw("P SAME");
        legend_resolution.AddEntry(previous_prototype_graph, "Previous Prototype", "p");
    }
    if (toa_cb_25ns_fit_graph) {
        toa_cb_25ns_fit_graph->SetMarkerStyle(23);
        toa_cb_25ns_fit_graph->SetMarkerSize(1);
        toa_cb_25ns_fit_graph->SetLineColor(kGreen+2);
        toa_cb_25ns_fit_graph->SetMarkerColor(kGreen+2);
        toa_cb_25ns_fit_graph->Draw("P SAME");
        legend_resolution.AddEntry(toa_cb_25ns_fit_graph, "ToA 25ns CB Fit", "p");
    }
    if (toa_cb_10ns_fit_graph) {
        toa_cb_10ns_fit_graph->SetMarkerStyle(24);
        toa_cb_10ns_fit_graph->SetMarkerSize(1);
        toa_cb_10ns_fit_graph->SetLineColor(kCyan+2);
        toa_cb_10ns_fit_graph->SetMarkerColor(kCyan+2);
        toa_cb_10ns_fit_graph->Draw("P SAME");
        legend_resolution.AddEntry(toa_cb_10ns_fit_graph, "ToA 10ns CB Fit", "p");
    }
    legend_resolution.Draw();
    resolution_canvas.Modified();
    resolution_canvas.Update();
    resolution_canvas.Write();
    resolution_canvas.SaveAs("dump/402_ADC_EnergyScanCompare/energy_resolution.pdf");
    resolution_canvas.Write();

    output_file->Close();


    return 0;
}
