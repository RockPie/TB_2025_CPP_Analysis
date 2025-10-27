#include <iostream>
#include <exception>
#include <unistd.h>
#include <sys/stat.h>
#include "TCanvas.h" 
#include "TVectorD.h"
#include "TVector.h"
#include "TF1.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TParameter.h"
// #include "easylogging++.h"
// #include "argparse/argparse.hpp"

#include "cxxopts.hpp"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#ifndef FPGA_CHANNEL_NUMBER
#define FPGA_CHANNEL_NUMBER 152
#endif


int main(int argc, char **argv){

    // * --- Initial settings -----------------------------------------------------------
    // * --------------------------------------------------------------------------------
    std::string script_input_file, script_output_file;
    int script_n_events;
    std::string script_name = __FILE__;
    std::string script_version = "0.1";
    std::string script_output_folder;

    // START_EASYLOGGINGPP(argc, argv);
    // set_easylogger();

    script_name = script_name.substr(script_name.find_last_of("/\\") + 1).substr(0, script_name.find_last_of("."));

    cxxopts::Options options(script_name, "Event reconstruction utility");
    options.add_options()
        ("f,file", "Input .root file", cxxopts::value<std::string>())
        ("o,output", "Output .root file", cxxopts::value<std::string>())
        ("e,events", "Number of events to process", cxxopts::value<int>()->default_value("-1"))
        ("h,help", "Print help");

    cxxopts::ParseResult parsed_opts;
    try {
        parsed_opts = options.parse(argc, argv);
        if (parsed_opts.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }
        if (!parsed_opts.count("file") || !parsed_opts.count("output")) {
            spdlog::error("Both --file and --output options are required.");
            std::cout << options.help() << std::endl;
            return 1;
        }
        script_input_file  = parsed_opts["file"].as<std::string>();
        script_output_file = parsed_opts["output"].as<std::string>();
        script_n_events    = parsed_opts["events"].as<int>();
    } catch (const std::exception& err) {
        spdlog::error("{}", err.what());
        std::cout << options.help() << std::endl;
        return 1;
    }

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

    // * --- Read input file ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    int machine_gun_samples = 0;
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
    TNamed *legal_fpga_id_list_tnamed = (TNamed*) input_root->Get("Rootifier_legal_fpga_id_list");
    if (legal_fpga_id_list_tnamed == nullptr) {
        spdlog::error("Failed to get legal fpga id list from input file {}", script_input_file);
        return 1;
    }
    std::string legal_fpga_id_list_str = legal_fpga_id_list_tnamed->GetTitle();
    std::vector <UShort_t> legal_fpga_id_list;
    std::istringstream legal_fpga_id_list_stream(legal_fpga_id_list_str);
    UShort_t legal_fpga_id;
    while (legal_fpga_id_list_stream >> legal_fpga_id) {
        legal_fpga_id_list.push_back(legal_fpga_id);
    }
    if (!legal_fpga_id_list.empty() && script_n_events > 0) {
        script_n_events = legal_fpga_id_list.size() * script_n_events;
    }
    spdlog::info("Number of events: {}", script_n_events);
    spdlog::info("Legal FPGA ID list: {}", legal_fpga_id_list_str);
    int total_input_entries = input_tree->GetEntries();
    auto entries_to_process = total_input_entries;
    spdlog::info("Total entries in input file: {}", total_input_entries);
    if (script_n_events > 0 && script_n_events < total_input_entries) {
        entries_to_process = script_n_events;
    }
    spdlog::info("Number of events to process: {}", entries_to_process);

    UShort_t  input_fpga_id;
    ULong64_t input_timestamp;
    UInt_t    input_daqh_list[4];
    Bool_t    input_tc_list[FPGA_CHANNEL_NUMBER];
    Bool_t    input_tp_list[FPGA_CHANNEL_NUMBER];
    UInt_t    input_val0_list[FPGA_CHANNEL_NUMBER];
    UInt_t    input_val1_list[FPGA_CHANNEL_NUMBER];
    UInt_t    input_val2_list[FPGA_CHANNEL_NUMBER];
    UInt_t    input_crc32_list[4];
    UInt_t    input_last_heartbeat;

    input_tree->SetBranchAddress("fpga_id", &input_fpga_id);
    input_tree->SetBranchAddress("timestamp", &input_timestamp);
    input_tree->SetBranchAddress("daqh_list", input_daqh_list);
    input_tree->SetBranchAddress("tc_list", input_tc_list);
    input_tree->SetBranchAddress("tp_list", input_tp_list);
    input_tree->SetBranchAddress("val0_list", input_val0_list);
    input_tree->SetBranchAddress("val1_list", input_val1_list);
    input_tree->SetBranchAddress("val2_list", input_val2_list);
    input_tree->SetBranchAddress("crc32_list", input_crc32_list);
    input_tree->SetBranchAddress("last_heartbeat", &input_last_heartbeat);

    int SWMA_window_size = 200;
    int SWMA_core_size   = 100;
    int SWMA_threshold   = 1;

    std::vector <UShort_t> fpga_id_pool;
    std::vector <ULong64_t> timestamp_pool;
    std::vector <UInt_t*> daqh_list_pool;
    std::vector <Bool_t*> tc_list_pool;
    std::vector <Bool_t*> tp_list_pool;
    std::vector <UInt_t*> val0_list_pool;
    std::vector <UInt_t*> val1_list_pool;
    std::vector <UInt_t*> val2_list_pool;
    std::vector <UInt_t*> crc32_list_pool;
    std::vector <UInt_t> last_heartbeat_pool;
    std::vector <bool> entry_matched;
    fpga_id_pool.reserve(SWMA_window_size);
    timestamp_pool.reserve(SWMA_window_size);
    entry_matched.reserve(SWMA_window_size);

    std::vector <int> window_matched_count;

    for (int _entry = 0; _entry < entries_to_process; _entry++) {
        input_tree->GetEntry(_entry);
        // only accept the first fpga ids
        // if (!legal_fpga_id_list.empty()) {
        //     // debug output
        //     if (input_fpga_id != legal_fpga_id_list[1]) {
        //         continue;
        //     }
        // }
        if (fpga_id_pool.size() < SWMA_window_size) {
            fpga_id_pool.push_back(input_fpga_id);
            timestamp_pool.push_back(input_timestamp);
            entry_matched.push_back(false);
        } else {
            for (int _pool_index = 0; _pool_index < SWMA_core_size; _pool_index++) {
                // spdlog::info("SWMA Pool Index {}: FPGA ID {}, Timestamp {}", _pool_index, fpga_id_pool[_pool_index], timestamp_pool[_pool_index]);
                if (entry_matched[_pool_index]) {
                    continue;
                }
                auto _timestamp_seed = timestamp_pool[_pool_index];
                auto _fpga_id_seed = fpga_id_pool[_pool_index];
                int _matched_count = 1;
                for (int _match_search_index = _pool_index + 1; _match_search_index < timestamp_pool.size(); _match_search_index++) {
                    if (entry_matched[_match_search_index]) {
                        continue;
                    }
                    uint32_t diff = (timestamp_pool[_match_search_index] > _timestamp_seed) ? (timestamp_pool[_match_search_index] - _timestamp_seed) : (_timestamp_seed - timestamp_pool[_match_search_index]);
                    if (diff > (UINT32_MAX / 2)) { 
                        diff = UINT32_MAX - diff; // Correct for overflow
                    }
                    if (diff < SWMA_threshold && fpga_id_pool[_match_search_index] == _fpga_id_seed) {
                        _matched_count++;
                        entry_matched[_match_search_index] = true;
                        // spdlog::info("Matched entry {} with seed entry {}: FPGA ID {}, Timestamp {}, diff {}", _match_search_index, _pool_index, _fpga_id_seed, timestamp_pool[_match_search_index], diff);
                    }
                    
                }
                if (_matched_count > 1) {
                    // LOG(INFO) << "Matched " << _matched_count << " entries with seed " << _fpga_id_seed << " at " << _timestamp_seed;
                    window_matched_count.push_back(_matched_count);
                    // if (_matched_count > machine_gun_samples) {
                    //     machine_gun_samples = _matched_count;
                        
                    // }
                }
            }
            // delete the core size
            fpga_id_pool.erase(fpga_id_pool.begin(), fpga_id_pool.begin() + SWMA_core_size);
            timestamp_pool.erase(timestamp_pool.begin(), timestamp_pool.begin() + SWMA_core_size);
            entry_matched.erase(entry_matched.begin(), entry_matched.begin() + SWMA_core_size);
        }
    }
    spdlog::info("Machine gun samples: {}", machine_gun_samples);
    // print all the matched counts
    // for (auto count : window_matched_count) {
    //     spdlog::info("Window matched count: {}", count);
    // }
    // set machine_gun_samples to the most often matched count
    std::map<int, int> count_map;
    for (auto count : window_matched_count) {
        count_map[count]++;
    }
    int max_count = 0;
    int most_often_count = 0;;
    for (auto const& [count, freq] : count_map) {
        if (freq > max_count) {
            max_count = freq;
            most_often_count = count;
        }
    }
    machine_gun_samples = most_often_count;
    spdlog::info("Determined machine gun samples: {}", machine_gun_samples);
    // print the percentage of the most often count
    double percentage = static_cast<double>(max_count) / static_cast<double>(window_matched_count.size()) * 100.0;
    spdlog::info("Most often count percentage: {:.2f}%", percentage);
    if (machine_gun_samples == 0) {
        spdlog::info("No machine gun samples found.");
        return 0;
    }

    // * --- Create output file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", script_output_file);
        return 1;
    }
    TTree* output_tree = new TTree("data_tree", "Data tree");
    output_tree->SetDirectory(output_root);

    auto *branch_fpga_id    = new UShort_t();   // 0 - 8
    auto *branch_timestamps = new ULong64_t[machine_gun_samples];  // 64 bits
    auto *branch_daqh_list  = new UInt_t[4 * machine_gun_samples];    // 32 bits
    auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
    auto *branch_crc32_list = new UInt_t[4 * machine_gun_samples];
    auto *branch_last_heartbeat = new UInt_t[machine_gun_samples];

    output_tree->Branch("fpga_id", branch_fpga_id, "fpga_id/s");
    output_tree->Branch("timestamps", branch_timestamps, ("timestamps[" + std::to_string(machine_gun_samples) + "]/l").c_str());
    output_tree->Branch("daqh_list", branch_daqh_list, ("daqh_list[" + std::to_string(4 * machine_gun_samples) + "]/i").c_str());
    output_tree->Branch("tc_list", branch_tc_list, ("tc_list[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/O").c_str());
    output_tree->Branch("tp_list", branch_tp_list, ("tp_list[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/O").c_str());
    output_tree->Branch("val0_list", branch_val0_list, ("val0_list[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/i").c_str());
    output_tree->Branch("val1_list", branch_val1_list, ("val1_list[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/i").c_str());
    output_tree->Branch("val2_list", branch_val2_list, ("val2_list[" + std::to_string(FPGA_CHANNEL_NUMBER * machine_gun_samples) + "]/i").c_str());
    output_tree->Branch("crc32_list", branch_crc32_list, ("crc32_list[" + std::to_string(4 * machine_gun_samples) + "]/i").c_str());
    output_tree->Branch("last_heartbeat", branch_last_heartbeat, ("last_heartbeat[" + std::to_string(machine_gun_samples) + "]/i").c_str());

    // * --- Fill output file -----------------------------------------------------------
    // * --------------------------------------------------------------------------------
    // clear the pools
    fpga_id_pool.clear();
    timestamp_pool.clear();
    entry_matched.clear();
    Long64_t good_machine_gun_events = 0;
    Long64_t bad_machine_gun_events = 0;
    bool flag_last_entry = false;
    for (int _entry = 0; _entry < entries_to_process; _entry++) {
        input_tree->GetEntry(_entry);
        if (_entry == entries_to_process - 1) {
            flag_last_entry = true;
        }
        // print out the timestamp and fpga_id
        // spdlog::info("Entry {}: FPGA ID: {}, Timestamp: {}", _entry, input_fpga_id, input_timestamp);
        if (fpga_id_pool.size() < SWMA_window_size || flag_last_entry) {
            fpga_id_pool.push_back(input_fpga_id);
            timestamp_pool.push_back(input_timestamp);
            // print out the timestamp and fpga_id
            UInt_t *daqh_list = new UInt_t[4];
            Bool_t *tc_list = new Bool_t[FPGA_CHANNEL_NUMBER];
            Bool_t *tp_list = new Bool_t[FPGA_CHANNEL_NUMBER];
            UInt_t *val0_list = new UInt_t[FPGA_CHANNEL_NUMBER];
            UInt_t *val1_list = new UInt_t[FPGA_CHANNEL_NUMBER];
            UInt_t *val2_list = new UInt_t[FPGA_CHANNEL_NUMBER];
            UInt_t *crc32_list = new UInt_t[4];

            for (int _daqh_index = 0; _daqh_index < 4; _daqh_index++) {
                daqh_list[_daqh_index] = input_daqh_list[_daqh_index];
                crc32_list[_daqh_index] = input_crc32_list[_daqh_index];
            }
            for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                // if (input_val0_list[_channel_index] > 200) {
                //     LOG(WARNING) << "Channel " << _channel_index << " has value " << input_val0_list[_channel_index] << ", which is greater than 200";
                // }
                tc_list[_channel_index] = input_tc_list[_channel_index];
                tp_list[_channel_index] = input_tp_list[_channel_index];
                val0_list[_channel_index] = input_val0_list[_channel_index];
                // if (_channel_index > 76 && _channel_index < 76 + 76) {
                //     LOG(INFO) << "FPGA" << input_fpga_id << " Chn" << _channel_index << " has value: " << input_val0_list[_channel_index];
                // }
                val1_list[_channel_index] = input_val1_list[_channel_index];
                val2_list[_channel_index] = input_val2_list[_channel_index];
            }

            // LOG(DEBUG) << "Entry " << _entry << " FPGA ID: " << input_fpga_id << " Timestamp: " << input_timestamp;

            daqh_list_pool.push_back(daqh_list);
            tc_list_pool.push_back(tc_list);
            tp_list_pool.push_back(tp_list);
            val0_list_pool.push_back(val0_list);
            val1_list_pool.push_back(val1_list);
            val2_list_pool.push_back(val2_list);
            crc32_list_pool.push_back(crc32_list);
            last_heartbeat_pool.push_back(input_last_heartbeat);

            entry_matched.push_back(false);
        }
        if (fpga_id_pool.size() >= SWMA_window_size || flag_last_entry) {
            // LOG(DEBUG) << "Processing window " << _entry;
            int _pool_index_max = SWMA_core_size;
            if (flag_last_entry) {
                _pool_index_max = fpga_id_pool.size();
            }
            for (int _pool_index = 0; _pool_index < _pool_index_max; _pool_index++) {
                if (entry_matched[_pool_index]) {
                    continue;
                }
                auto _timestamp_seed = timestamp_pool[_pool_index];
                auto _fpga_id_seed = fpga_id_pool[_pool_index];
                int _matched_count = 1;
                std::vector <int> _matched_index;
                _matched_index.push_back(_pool_index);
                for (int _match_search_index = _pool_index + 1; _match_search_index < timestamp_pool.size(); _match_search_index++) {
                    if (entry_matched[_match_search_index]) {
                        continue;
                    }
                    uint32_t diff = (timestamp_pool[_match_search_index] > _timestamp_seed) ? (timestamp_pool[_match_search_index] - _timestamp_seed) : (_timestamp_seed - timestamp_pool[_match_search_index]);
                    if (diff > (UINT32_MAX / 2)) { 
                        diff = UINT32_MAX - diff; // Correct for overflow
                    }
                    if (diff < SWMA_threshold && fpga_id_pool[_match_search_index] == _fpga_id_seed) {
                        _matched_count++;
                        entry_matched[_match_search_index] = true;
                        _matched_index.push_back(_match_search_index);
                    }
                }
                if (_matched_count == machine_gun_samples) {
                    good_machine_gun_events++;
                    // -- Put to output tree --------------------------------------------
                    *branch_fpga_id = _fpga_id_seed;
                    for (int _matched_index_index = 0; _matched_index_index < _matched_index.size(); _matched_index_index++) {
                        branch_timestamps[_matched_index_index] = timestamp_pool[_matched_index[_matched_index_index]];
                        for (int _daqh_index = 0; _daqh_index < 4; _daqh_index++) {
                            branch_daqh_list[_matched_index_index * 4 + _daqh_index] = daqh_list_pool[_matched_index[_matched_index_index]][_daqh_index];
                            branch_crc32_list[_matched_index_index * 4 + _daqh_index] = crc32_list_pool[_matched_index[_matched_index_index]][_daqh_index];
                        }
                        for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                            // if (_channel_index == 10) {
                            //     LOG(DEBUG) << "Channel 10: " << input_val0_list[_channel_index];
                            // }
                            branch_tc_list[_matched_index_index * FPGA_CHANNEL_NUMBER + _channel_index] = tc_list_pool[_matched_index[_matched_index_index]][_channel_index];
                            branch_tp_list[_matched_index_index * FPGA_CHANNEL_NUMBER + _channel_index] = tp_list_pool[_matched_index[_matched_index_index]][_channel_index];
                            branch_val0_list[_matched_index_index * FPGA_CHANNEL_NUMBER + _channel_index] = val0_list_pool[_matched_index[_matched_index_index]][_channel_index];
                            branch_val1_list[_matched_index_index * FPGA_CHANNEL_NUMBER + _channel_index] = val1_list_pool[_matched_index[_matched_index_index]][_channel_index];
                            branch_val2_list[_matched_index_index * FPGA_CHANNEL_NUMBER + _channel_index] = val2_list_pool[_matched_index[_matched_index_index]][_channel_index];
                        }
                        branch_last_heartbeat[_matched_index_index] = last_heartbeat_pool[_matched_index[_matched_index_index]];
                    }
                    output_tree->Fill();
                } else {
                    bad_machine_gun_events++;
                }
            }
            // delete the core size
            fpga_id_pool.erase(fpga_id_pool.begin(), fpga_id_pool.begin() + SWMA_core_size);
            timestamp_pool.erase(timestamp_pool.begin(), timestamp_pool.begin() + SWMA_core_size);
            for (int _pool_index = 0; _pool_index < SWMA_core_size; _pool_index++) {
                delete daqh_list_pool[_pool_index];
                delete tc_list_pool[_pool_index];
                delete tp_list_pool[_pool_index];
                delete val0_list_pool[_pool_index];
                delete val1_list_pool[_pool_index];
                delete val2_list_pool[_pool_index];
                delete crc32_list_pool[_pool_index];
            }
            daqh_list_pool.erase(daqh_list_pool.begin(), daqh_list_pool.begin() + SWMA_core_size);
            tc_list_pool.erase(tc_list_pool.begin(), tc_list_pool.begin() + SWMA_core_size);
            tp_list_pool.erase(tp_list_pool.begin(), tp_list_pool.begin() + SWMA_core_size);
            val0_list_pool.erase(val0_list_pool.begin(), val0_list_pool.begin() + SWMA_core_size);
            val1_list_pool.erase(val1_list_pool.begin(), val1_list_pool.begin() + SWMA_core_size);
            val2_list_pool.erase(val2_list_pool.begin(), val2_list_pool.begin() + SWMA_core_size);
            crc32_list_pool.erase(crc32_list_pool.begin(), crc32_list_pool.begin() + SWMA_core_size);
            last_heartbeat_pool.erase(last_heartbeat_pool.begin(), last_heartbeat_pool.begin() + SWMA_core_size);
            entry_matched.erase(entry_matched.begin(), entry_matched.begin() + SWMA_core_size);
        }
    }

    if ((good_machine_gun_events + bad_machine_gun_events) == 0) {
        spdlog::error("No machine gun events found!");
    } else {
        double good_ratio = static_cast<double>(good_machine_gun_events) /
                            static_cast<double>(good_machine_gun_events + bad_machine_gun_events) * 100.0;
        spdlog::info("Good machine gun events: {} ({:.2f}%)", good_machine_gun_events, good_ratio);
    }

    // Read the meta data from the input file
    TNamed *input_script_name_tnamed = (TNamed*) input_root->Get("Rootifier_script_name");
    TNamed *input_script_version_tnamed = (TNamed*) input_root->Get("Rootifier_script_version");
    TNamed *input_script_input_file_tnamed = (TNamed*) input_root->Get("Rootifier_script_input_file");
    TNamed *input_script_output_file_tnamed = (TNamed*) input_root->Get("Rootifier_script_output_file");
    TNamed *input_script_n_events_tnamed = (TNamed*) input_root->Get("Rootifier_script_n_events");
    TNamed *input_script_output_folder_tnamed = (TNamed*) input_root->Get("Rootifier_script_output_folder");

    // Read the running info from the input file
    TNamed *input_counter_header_line_tnamed = (TNamed*) input_root->Get("Rootifier_counter_header_line");
    TNamed *input_counter_heartbeat_packet_tnamed = (TNamed*) input_root->Get("Rootifier_counter_heartbeat_packet");
    TNamed *input_counter_complete_20_lines_tnamed = (TNamed*) input_root->Get("Rootifier_counter_complete_20_lines");
    TNamed *input_counter_not_complete_20_lines_tnamed = (TNamed*) input_root->Get("Rootifier_counter_not_complete_20_lines");
    TNamed *input_counter_not_complete_20_lines_total_tnamed = (TNamed*) input_root->Get("Rootifier_counter_not_complete_20_lines_total");
    TNamed *input_counter_invalid_line_tnamed = (TNamed*) input_root->Get("Rootifier_counter_invalid_line");
    TNamed *input_counter_valid_line_tnamed = (TNamed*) input_root->Get("Rootifier_counter_valid_line");
    TNamed *input_legal_fpga_id_list_tnamed = (TNamed*) input_root->Get("Rootifier_legal_fpga_id_list");

    input_root->Close();

    output_root->cd();
    output_tree->Write();

    // -- Write the TNamed objects from the input file
    if (input_script_name_tnamed != nullptr) {
        input_script_name_tnamed->Write();
    }
    if (input_script_version_tnamed != nullptr) {
        input_script_version_tnamed->Write();
    }
    if (input_script_input_file_tnamed != nullptr) {
        input_script_input_file_tnamed->Write();
    }
    if (input_script_output_file_tnamed != nullptr) {
        input_script_output_file_tnamed->Write();
    }
    if (input_script_n_events_tnamed != nullptr) {
        input_script_n_events_tnamed->Write();
    }
    if (input_script_output_folder_tnamed != nullptr) {
        input_script_output_folder_tnamed->Write();
    }
    if (input_counter_header_line_tnamed != nullptr) {
        input_counter_header_line_tnamed->Write();
    }
    if (input_counter_heartbeat_packet_tnamed != nullptr) {
        input_counter_heartbeat_packet_tnamed->Write();
    }
    if (input_counter_complete_20_lines_tnamed != nullptr) {
        input_counter_complete_20_lines_tnamed->Write();
    }
    if (input_counter_not_complete_20_lines_tnamed != nullptr) {
        input_counter_not_complete_20_lines_tnamed->Write();
    }
    if (input_counter_not_complete_20_lines_total_tnamed != nullptr) {
        input_counter_not_complete_20_lines_total_tnamed->Write();
    }
    if (input_counter_invalid_line_tnamed != nullptr) {
        input_counter_invalid_line_tnamed->Write();
    }
    if (input_counter_valid_line_tnamed != nullptr) {
        input_counter_valid_line_tnamed->Write();
    }
    if (input_legal_fpga_id_list_tnamed != nullptr) {
        input_legal_fpga_id_list_tnamed->Write();
    }

    // -- generate machinegun 2D histogram for each channel
    // -- X axis: machinegun index
    // -- Y axis: ADC/ToT/ToA
    // int channel_mg_y_bins = 256 ;
    // for (int _fpga_index = 0; _fpga_index < legal_fpga_id_list.size(); _fpga_index++) {
    //     TCanvas *canvas = new TCanvas(
    //         ("fpga" + std::to_string(legal_fpga_id_list[_fpga_index]) + "_canvas").c_str(),
    //         ("FPGA" + std::to_string(legal_fpga_id_list[_fpga_index]) + " Canvas").c_str(),
    //         1200, 800);
    //     canvas->Divide(int(std::ceil(std::sqrt(FPGA_CHANNEL_NUMBER))), int(std::ceil(std::sqrt(FPGA_CHANNEL_NUMBER))));
    //     auto _fpga_id = legal_fpga_id_list[_fpga_index];
    //     for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
    //         TH2D *hist2d_val0 = new TH2D(
    //             ("fpga" + std::to_string(_fpga_id) + "_chn" + std::to_string(_channel_index) + "_val0").c_str(),
    //             ("FPGA" + std::to_string(_fpga_id) + " Chn" + std::to_string(_channel_index) + " Val0").c_str(),
    //             machine_gun_samples, 0, machine_gun_samples,
    //             channel_mg_y_bins, 0, 1024);
    //         TH2D *hist2d_val1 = new TH2D(
    //             ("fpga" + std::to_string(_fpga_id) + "_chn" + std::to_string(_channel_index) + "_val1").c_str(),
    //             ("FPGA" + std::to_string(_fpga_id) + " Chn" + std::to_string(_channel_index) + " Val1").c_str(),
    //             machine_gun_samples, 0, machine_gun_samples,
    //             channel_mg_y_bins, 0, 1024);
    //         TH2D *hist2d_val2 = new TH2D(
    //             ("fpga" + std::to_string(_fpga_id) + "_chn" + std::to_string(_channel_index) + "_val2").c_str(),
    //             ("FPGA" + std::to_string(_fpga_id) + " Chn" + std::to_string(_channel_index) + " Val2").c_str(),
    //             machine_gun_samples, 0, machine_gun_samples,
    //             channel_mg_y_bins, 0, 1024);
    //         output_root->cd();
    //         int mg_index = 0;
    //         for (int _entry = 0; _entry < output_tree->GetEntries(); _entry++) {
    //             output_tree->GetEntry(_entry);
    //             if (*branch_fpga_id != _fpga_id) {
    //                 continue;
    //             }
    //             for (int _mg_sample_index = 0; _mg_sample_index < machine_gun_samples; _mg_sample_index++) {
    //                 hist2d_val0->Fill(mg_index, branch_val0_list[_mg_sample_index * FPGA_CHANNEL_NUMBER + _channel_index]);
    //                 hist2d_val1->Fill(mg_index, branch_val1_list[_mg_sample_index * FPGA_CHANNEL_NUMBER + _channel_index]);
    //                 hist2d_val2->Fill(mg_index, branch_val2_list[_mg_sample_index * FPGA_CHANNEL_NUMBER + _channel_index]);
    //                 mg_index++;
    //             }
    //         }
    //         canvas->cd(_channel_index + 1);
    //         hist2d_val0->GetXaxis()->SetTitle("Machine Gun Index");
    //         hist2d_val0->GetYaxis()->SetTitle("ADC Value");
    //         hist2d_val0->Draw("COLZ");
    //     }
    //     output_root->cd();
    //     canvas->Write();
    // }

    

    // -- Write the meta data
    TNamed("EventRecon_script_name", script_name.c_str()).Write();
    TNamed("EventRecon_script_version", script_version.c_str()).Write();
    TNamed("EventRecon_script_input_file", script_input_file.c_str()).Write();
    TNamed("EventRecon_script_output_file", script_output_file.c_str()).Write();
    TNamed("EventRecon_script_n_events", std::to_string(script_n_events).c_str()).Write();
    TNamed("EventRecon_script_output_folder", script_output_folder.c_str()).Write();
    // -- Write the running info
    TNamed("EventRecon_good_machine_gun_events", std::to_string(good_machine_gun_events).c_str()).Write();
    TNamed("EventRecon_bad_machine_gun_events", std::to_string(bad_machine_gun_events).c_str()).Write();
    TNamed("EventRecon_machine_gun_samples", std::to_string(machine_gun_samples).c_str()).Write();
    // --- Write with TParameter
    TParameter<double>("EventRecon_good_machine_gun_events_ratio", 
        static_cast<double>(good_machine_gun_events) /
        static_cast<double>(good_machine_gun_events + bad_machine_gun_events) * 100.0).Write();
    TParameter<int>("EventRecon_machine_gun_samples", machine_gun_samples).Write();
    
    output_root->Close();
    return 0;
}

// void set_easylogger(){
//     el::Configurations defaultConf;
//     defaultConf.setToDefault();
//     defaultConf.setGlobally(el::ConfigurationType::Format, "%datetime{%H:%m:%s}[%levshort] (%fbase) %msg");
//     defaultConf.set(el::Level::Info,    el::ConfigurationType::Format, 
//         "%datetime{%H:%m:%s}[\033[1;34m%levshort\033[0m] (%fbase) %msg");
//     defaultConf.set(el::Level::Warning, el::ConfigurationType::Format, 
//         "%datetime{%H:%m:%s}[\033[1;33m%levshort\033[0m] (%fbase) %msg");
//     defaultConf.set(el::Level::Error,   el::ConfigurationType::Format, 
//         "%datetime{%H:%m:%s}[\033[1;31m%levshort\033[0m] (%fbase) %msg");
//     el::Loggers::reconfigureLogger("default", defaultConf);
// }