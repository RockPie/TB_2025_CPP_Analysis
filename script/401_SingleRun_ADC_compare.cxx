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

    std::string raw_adc_analysis_root_file = "dump/304_RawADC/beamtests/Run0452.root";
    // get the energy label from the run number in the file name
    std::string run_info_from_files = raw_adc_analysis_root_file.substr(raw_adc_analysis_root_file.find_last_of("/\\") + 1);
    std::string toa_filtered_adc_analysis_root_file_prefix = "dump/314_ToA_Analysis_LUT/beamtests/";
    std::string output_file_prefix = "dump/401_SingleRun_ADC_compare/";

    std::string toa_filtered_adc_analysis_root_file = toa_filtered_adc_analysis_root_file_prefix + run_info_from_files;
    std::string output_root_file = output_file_prefix + run_info_from_files.substr(0, run_info_from_files.find_last_of(".")) + "_ADC_comparison.root";


    // std::string toa_filtered_adc_analysis_root_file = "dump/314_ToA_Analysis_LUT/beamtests/Run0452.root";
    // std::string output_root_file = "dump/401_SingleRun_ADC_compare/Run0452_ADC_comparison.root";

    std::map<std::string, std::string> run_energies = {
        {"Run0452", "_60"},
        {"Run0453", "_350"},
        {"Run0447", "_250"}
    };


    std::string energy_label = "";
    for (const auto& [run_number, label] : run_energies) {
        if (run_info_from_files.find(run_number) != std::string::npos) {
            energy_label = label;
            break;
        }
    }
    if (energy_label.empty()) {
        spdlog::warn("Could not determine energy label from file name {}, defaulting to empty string", raw_adc_analysis_root_file);
        return 1;
    } else {
        spdlog::info("Determined energy label {} from file name {}", energy_label, raw_adc_analysis_root_file);
    }
    std::string simulation_adc_analysis_root_file = "data/simulation/Config86/10_200" + energy_label + "-Config86-M.root";
    // std::string simulation_adc_analysis_root_file = "data/simulation/Config63/10_1000_60-Proton.root";
    // std::string simulation_adc_analysis_root_file = "data/simulation/Config63/10_1000" + energy_label + "-Proton.root";

    // std::string simulation_adc_analysis_root_file = "data/simulation/Config40/10_1000_60-Config40.root";

    TFile *raw_adc_file = new TFile(raw_adc_analysis_root_file.c_str(), "READ");
    if (raw_adc_file->IsZombie()) {
        spdlog::error("Failed to open raw ADC analysis file {}", raw_adc_analysis_root_file);
        return 1;
    }
    // find canvas with name adc_sum_distribution_canvas
    TCanvas *raw_adc_canvas = nullptr;
    TIter next(raw_adc_file->GetListOfKeys());
    TKey *key = nullptr;
    while ((key = dynamic_cast<TKey*>(next()))) {
        const char *key_name = key->GetName();
        if (key_name && std::string(key_name) == "adc_sum_distribution_canvas") {
            auto raw_adc_canvas_in_file = dynamic_cast<TCanvas*>(key->ReadObj());
            if (!raw_adc_canvas_in_file) {
                continue;
            }
            raw_adc_canvas = (TCanvas*)raw_adc_canvas_in_file->Clone("raw_adc_canvas");
            break;
        }
    }
    if (!raw_adc_canvas) {
        spdlog::error("Failed to find canvas adc_sum_distribution_canvas in file {}", raw_adc_analysis_root_file);
        return 1;
    }

    TH1D *raw_adc_hist = nullptr;
    // find the first TH1D in the canvas (before closing the file!)
    TIter next_primitive(raw_adc_canvas->GetListOfPrimitives());
    TObject *primitive = nullptr;
    while ((primitive = next_primitive())) {
        auto hist = dynamic_cast<TH1D*>(primitive);
        if (hist) {
            raw_adc_hist = (TH1D*)hist->Clone("raw_adc_hist");
            raw_adc_hist->SetDirectory(nullptr); // Detach from file to prevent it from being deleted when file is closed
            break;
        }
        // print the class name of the primitive for debugging
        spdlog::debug("Found primitive of class {}", primitive->ClassName());
    }
    if (!raw_adc_hist) {
        spdlog::error("Failed to find TH1D in canvas adc_sum_distribution_canvas in file {}", raw_adc_analysis_root_file);
        return 1;
    }
    
    // Close the file after extracting and cloning the histogram
    raw_adc_file->Close();

    spdlog::info("Opening simulation file: {}", simulation_adc_analysis_root_file);
    
    // Suppress ROOT warnings about missing dictionaries for objects we don't need
    Int_t oldIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    
    TFile *simulation_adc_file = new TFile(simulation_adc_analysis_root_file.c_str(), "READ");
    spdlog::info("File opened, checking if zombie...");
    
    if (simulation_adc_file->IsZombie()) {
        gErrorIgnoreLevel = oldIgnoreLevel;
        spdlog::error("Failed to open simulation ADC analysis file {}", simulation_adc_analysis_root_file);
        return 1;
    }
    
    spdlog::info("File OK, navigating to Photons directory...");
    
    // find th1d under Photons folder with name ADC
    // Navigate directly using cd() to avoid triggering reads of other objects
    TH1D *simulation_adc_hist = nullptr;
    if (simulation_adc_file->cd("Photons")) {
        spdlog::info("In Photons directory, getting ADC histogram...");
        TObject *adc_obj = gDirectory->Get("ADC");
        if (adc_obj) {
            spdlog::info("Found ADC object, casting to TH1D...");
            TH1D *hist_from_file = dynamic_cast<TH1D*>(adc_obj);
            if (hist_from_file) {
                spdlog::info("Cloning histogram...");
                simulation_adc_hist = (TH1D*)hist_from_file->Clone("simulation_adc_hist");
                simulation_adc_hist->SetDirectory(nullptr); // Detach from file to prevent it from being deleted when file is closed
                spdlog::info("Clone complete");
            }
        }
    }
    if (!simulation_adc_hist) {
        gErrorIgnoreLevel = oldIgnoreLevel;
        spdlog::error("Failed to find histogram ADC in folder Photons in file {}", simulation_adc_analysis_root_file);
        return 1;
    }
    
    spdlog::info("Closing simulation file...");
    // Close the file after cloning the histogram
    simulation_adc_file->Close();
    spdlog::info("Simulation file closed");

    TFile *toa_filtered_adc_file = new TFile(toa_filtered_adc_analysis_root_file.c_str(), "READ");
    if (toa_filtered_adc_file->IsZombie()) {
        gErrorIgnoreLevel = oldIgnoreLevel;
        spdlog::error("Failed to open ToA filtered ADC analysis file {}", toa_filtered_adc_analysis_root_file);
        return 1;
    }
    // find canvas with name canvas_adc_sum_comparison
    TCanvas *toa_filtered_adc_canvas = nullptr;
    TIter next_toa_canvas(toa_filtered_adc_file->GetListOfKeys());
    while ((key = dynamic_cast<TKey*>(next_toa_canvas()))) {
        const char *key_name = key->GetName();
        if (key_name && std::string(key_name) == "canvas_adc_sum_comparison") {
            auto toa_canvas_in_file = dynamic_cast<TCanvas*>(key->ReadObj());
            if (!toa_canvas_in_file) {
                continue;
            }
            toa_filtered_adc_canvas = (TCanvas*)toa_canvas_in_file->Clone("toa_filtered_adc_canvas");
            break;
        }
    }
    if (!toa_filtered_adc_canvas) {
        gErrorIgnoreLevel = oldIgnoreLevel;
        spdlog::error("Failed to find canvas canvas_adc_sum_comparison in file {}", toa_filtered_adc_analysis_root_file);
        return 1;
    }

    // list the primitives in the canvas for debugging
    TH1D *toa_filtered_25ns_adc_hist = nullptr;
    TH1D *toa_filtered_10ns_adc_hist = nullptr;
    spdlog::info("Listing primitives in ToA filtered ADC canvas:");
    TIter next_toa_primitive(toa_filtered_adc_canvas->GetListOfPrimitives());
    while ((primitive = next_toa_primitive())) {
        spdlog::info("Found primitive of class {}", primitive->ClassName());
        // if it is TH1D, print its name
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
    
    // Restore error level
    gErrorIgnoreLevel = oldIgnoreLevel;

    TFile *output_root = new TFile(output_root_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", output_root_file);
        return 1;
    }
    // draw the simulation and raw ADC histograms on the same canvas for comparison
    TCanvas *comparison_canvas = new TCanvas("comparison_canvas", "ADC Comparison", 800, 600);
    TLegend *legend = new TLegend(0.7, 0.7, 0.89, 0.89);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    // normalize both histograms to unit area for better shape comparison
    raw_adc_hist->Scale(1.0 / raw_adc_hist->Integral());
    std::string energy_label_for_title = energy_label.empty() ? "" : (" at " + energy_label.substr(1) + " GeV");
    format_1d_hist_canvas(comparison_canvas, raw_adc_hist, kRed, "H2GCROC Beam Test", "MC Simulation Comparison", energy_label_for_title.c_str(), false);
    legend->AddEntry(raw_adc_hist, "Beam Test ADC", "l");

    simulation_adc_hist->Rebin(2);
    simulation_adc_hist->SetLineColor(kBlue);
    simulation_adc_hist->SetLineWidth(2);
    simulation_adc_hist->Scale(1.0 / simulation_adc_hist->Integral());
    simulation_adc_hist->Draw("HIST SAME");
    legend->AddEntry(simulation_adc_hist, "Simulation ADC", "l");

    // toa_filtered_10ns_adc_hist->Rebin(1);
    // toa_filtered_10ns_adc_hist->SetLineColor(kGreen+2);
    // toa_filtered_10ns_adc_hist->SetLineWidth(2);
    // toa_filtered_10ns_adc_hist->Scale(0.8 / toa_filtered_10ns_adc_hist->Integral());
    // toa_filtered_10ns_adc_hist->Draw("HIST SAME");
    // legend->AddEntry(toa_filtered_10ns_adc_hist, "ToA Filtered ADC (10ns)", "l");

    toa_filtered_25ns_adc_hist->Rebin(2);
    toa_filtered_25ns_adc_hist->SetLineColor(kMagenta+2);
    toa_filtered_25ns_adc_hist->SetLineWidth(2);
    toa_filtered_25ns_adc_hist->Scale(1.0 / toa_filtered_25ns_adc_hist->Integral());
    toa_filtered_25ns_adc_hist->Draw("HIST SAME");
    legend->AddEntry(toa_filtered_25ns_adc_hist, "ToA Filtered ADC (25ns)", "l");

    legend->Draw();
    comparison_canvas->Modified();
    comparison_canvas->Update();
    comparison_canvas->Write();
    comparison_canvas->SaveAs((output_root_file.substr(0, output_root_file.find_last_of(".")) + "_adc_comparison.pdf").c_str());
    comparison_canvas->Close();

    output_root->Close();

    return 0;
}