#include <binparse/parser.hpp>
#include <binparse/tail.hpp>
#include <binparse/bytecursor.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <array>
#include <unordered_set>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <THttpServer.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TVectorD.h>
#include <TLatex.h>

#include "cxxopts.hpp"

#define FPGA_CHANNEL_NUMBER 152
#define BX_PER_ORBIT 3564
#define FPGA_NUMBER 2

#define VERBOSE_LEVEL_SILENT    0
#define VERBOSE_LEVEL_ERROR     1
#define VERBOSE_LEVEL_WARN      2
#define VERBOSE_LEVEL_INFO      3
#define VERBOSE_LEVEL_DEBUG     4

static inline uint32_t crc32_msb_table(uint8_t idx) {
    // One-time table generator (lazy, thread-safe in C++11 with function-local statics)
    static std::array<uint32_t, 256> T = []{
        std::array<uint32_t, 256> t{};
        const uint32_t poly = 0x04C11DB7u;
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i << 24;
            for (int k = 0; k < 8; ++k)
                c = (c & 0x80000000u) ? ((c << 1) ^ poly) : (c << 1);
            t[i] = c;
        }
        return t;
    }();
    return T[idx];
}

inline uint32_t crc32_msb_update(uint32_t crc, const uint8_t* p, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        uint8_t idx = static_cast<uint8_t>((crc >> 24) ^ p[i]);
        crc = (crc << 8) ^ crc32_msb_table(idx);
    }
    return crc; // xorOut = 0
}

// Feed raw bytes
inline uint32_t crc32_msb(const uint8_t* data, size_t len,
                          uint32_t init = 0x00000000u) {
    return crc32_msb_update(init, data, len);
}

// Feed a vector<uint32_t> as BIG-ENDIAN bytes (matches your example)
inline uint32_t crc32_u32vec_be(const std::vector<uint32_t>& v,
                                uint32_t init = 0x00000000u) {
    uint32_t crc = init;
    for (uint32_t x : v) {
        uint8_t b[4] = {
            static_cast<uint8_t>((x >> 24) & 0xFF),
            static_cast<uint8_t>((x >> 16) & 0xFF),
            static_cast<uint8_t>((x >>  8) & 0xFF),
            static_cast<uint8_t>( x        & 0xFF)
        };
        crc = crc32_msb_update(crc, b, 4);
    }
    return crc; // xorOut = 0
}

void parse_word(uint32_t& _word, Bool_t& _tc, Bool_t& _tp, UInt_t& _val0, UInt_t& _val1, UInt_t& _val2) {
    _tc   = (_word >> 31) & 0x1;
    _tp   = (_word >> 30) & 0x1;
    _val0 = (_word >> 20) & 0x3FF;
    _val1 = (_word >> 10) & 0x3FF;
    _val2 = (_word      ) & 0x3FF;
    return;
}

bool parse_160byte_frame(std::vector<std::vector<uint32_t>>& _frame_words, std::vector<UInt_t>& _val0_array, std::vector<UInt_t>& _val1_array, std::vector<UInt_t>& _val2_array, std::vector<Bool_t>& _tc_array, std::vector<Bool_t>& _tp_array, std::vector<UInt_t>& _daqh_list, std::vector<UInt_t>& _crc_list, bool _check_crc, int _verbose=0) {
    if (_frame_words.empty() || _frame_words[0].size() != 40) {
        spdlog::error("Frame size is not 160 bytes (40 words). Actual size: {}", _frame_words.size() * 4);
        return false;
    }
    _val0_array.clear();    // FPGA_CHANNEL_NUMBER
    _val1_array.clear();    // FPGA_CHANNEL_NUMBER
    _val2_array.clear();    // FPGA_CHANNEL_NUMBER

    _tc_array.clear();  // FPGA_CHANNEL_NUMBER
    _tp_array.clear();  // FPGA_CHANNEL_NUMBER

    _daqh_list.clear(); // 4 words
    _crc_list.clear();  // 4 words

    _crc_list.push_back(_frame_words[0].back());
    _crc_list.push_back(_frame_words[1].back());
    _crc_list.push_back(_frame_words[2].back());
    _crc_list.push_back(_frame_words[3].back());

    if (_check_crc) {
        for (int i = 0; i < 4; i++) {
            auto crc_39words = crc32_u32vec_be(std::vector<uint32_t>(_frame_words[i].begin(), _frame_words[i].end() - 1));
            if (crc_39words != _crc_list[i]) {
                if (_verbose >= VERBOSE_LEVEL_ERROR) {
                    spdlog::error("CRC32 mismatch in part {}! C: 0x{:08X}, E: 0x{:08X}", i, crc_39words, _crc_list[i]);
                }
                return false;
            }
        }
    }

    _daqh_list.push_back(_frame_words[0].front());
    _daqh_list.push_back(_frame_words[1].front());
    _daqh_list.push_back(_frame_words[2].front());
    _daqh_list.push_back(_frame_words[3].front());

    for (int _asic_half_index = 0; _asic_half_index < 4; _asic_half_index++) {
        for (int _word_index = 1; _word_index < 39; _word_index++) {
            auto &_word = _frame_words[_asic_half_index][_word_index];
            Bool_t tc, tp;
            UInt_t val0, val1, val2;
            // print the word in hex
            // std::cout << "Word " << (_asic_half_index * 39 + _word_index - 1) << ": 0x" << std::hex << _word << std::dec << std::endl;
            parse_word(_word, tc, tp, val0, val1, val2);
            // print the parsed values
            // if (true) {
            //     spdlog::info("  Word {}: TC={}, TP={}, Val0={}, Val1={}, Val2={}", (_asic_half_index * 39 + _word_index - 1), tc, tp, val0, val1, val2);
            // }
            _tc_array.push_back(tc);
            _tp_array.push_back(tp);
            _val0_array.push_back(val0);
            _val1_array.push_back(val1);
            _val2_array.push_back(val2);
        }
    }

    return true;

}

void trigger_data_processing(bp::TrgLine& _line, std::vector<std::vector<bp::DataLine>>& _vldb_data_lines, std::size_t _trg_index, UShort_t* _branch_fpga_id, ULong64_t* _branch_timestamp, UInt_t* _branch_daqh_list, Bool_t* _branch_tc_list, Bool_t* _branch_tp_list, UInt_t* _branch_val0_list, UInt_t* _branch_val1_list, UInt_t* _branch_val2_list, UInt_t* _branch_crc32_list, UInt_t* _branch_last_heartbeat, TTree* _output_tree, std::size_t &_n_trigger_40line_error, std::size_t &_n_trigger_160byte_error, std::size_t &_n_trigger_header_ts_error, std::size_t &_n_legal_samples, int _verbose=0) {
    // if (true) {
    //     spdlog::info("Trg {}: bx=0x{:03X}, ob=0x{:03X}", (int)_trg_index, _line.bx_cnt, _line.ob_cnt);
    // }
    auto& trg_bx = _line.bx_cnt;
    auto& trg_ob = _line.ob_cnt;
    const int timestamp_diff_mg               = 41;
    const int timestamp_diff_mg_tolerance     = 2;
    const int timestamp_diff_sample           = 39;
    const int timestamp_diff_sample_tolerance = 0;
    uint64_t trg_timestamp = (static_cast<uint64_t>(trg_ob) * BX_PER_ORBIT) + trg_bx;

    std::vector<UInt_t> val0_array_min;
    for (int i = 0; i < FPGA_CHANNEL_NUMBER; i++) {
        val0_array_min.push_back(1023);
    }

    

    std::vector<std::vector<uint32_t>> frame_words(6);

    
    for (auto &_data_lines : _vldb_data_lines) {
        uint64_t idle_line_cnt = 0;
        uint64_t line_skipped_cnt = 0;
        uint64_t first_valid_data_timestamp = 0;
        uint64_t current_lead_data_timestamp = 0;
        int max_machinegun_index = 0;
        int L1A_cmd_counter = 0;
        bool first_data_found   = false;
        bool first_mg_flag       = false;
        bool looking_for_header = true;
        bool skip_line_mask = false;
        int line_candidate_counter = 0;
        int machinegun_index = -1;
        int candidate_header_line_index = -1;
        for (auto data_line_index = 0; data_line_index < _data_lines.size(); data_line_index++) {
            auto& data_line = _data_lines[data_line_index];
            auto& data_bx = data_line.bx_cnt;
            auto& data_ob = data_line.ob_cnt;
            uint64_t data_timestamp = (static_cast<uint64_t>(data_ob) * BX_PER_ORBIT) + data_bx;
            // print the bx and ob in hex
            // if (data_bx > 0xf00 || data_bx < 0x100) {
            //    spdlog::info("  Data line: vldb_id={}, bx=0x{:03X}, ob=0x{:03X}, timestamp=0x{:08X}", (int)data_line.header_vldb_id, data_bx, data_ob, data_timestamp);
            // }
            bool is_idle_line = false;
            if ((data_line.data_word0 == 0xACCCCCCC) || (data_line.data_word1 == 0xACCCCCCC) || (data_line.data_word2 == 0xACCCCCCC) || (data_line.data_word3 == 0xACCCCCCC)) {
                if (skip_line_mask == false){
                    idle_line_cnt++;
                }
                is_idle_line = true;
            }
            if (data_line.data_word4 == 0x4B00004B){
                if (skip_line_mask == false){
                    L1A_cmd_counter++;
                    // print line info and data
                    if (_verbose >= VERBOSE_LEVEL_INFO) {
                        spdlog::info("  L1A cmd: vldb_id={}, bx=0x{:03X}, ob=0x{:03X}, line_index={}",
                                     (int)data_line.header_vldb_id, data_bx, data_ob, data_line_index);
                        spdlog::info("    data_word0=0x{:08X}, data_word1=0x{:08X}, data_word2=0x{:08X}, data_word3=0x{:08X}, data_word4=0x{:08X}, data_word5=0x{:08X}",
                                     data_line.data_word0, data_line.data_word1, data_line.data_word2, data_line.data_word3, data_line.data_word4, data_line.data_word5);
                    }
                }
            }
            
            if (looking_for_header) {
                // check if the word start with 0x5 and ends with 0x5
                // if ((data_line.data_word0 >> 28 == 0x5 && (data_line.data_word0 & 0xF) == 0x5) || (data_line.data_word1 >> 28 == 0x5 && (data_line.data_word1 & 0xF) == 0x5) || (data_line.data_word2 >> 28 == 0x5 && (data_line.data_word2 & 0xF) == 0x5) || (data_line.data_word3 >> 28 == 0x5 && (data_line.data_word3 & 0xF) == 0x5)) {
                // check if the word start with 0xF and ends with 0x5 or 0x2
                if (
                    (data_line.data_word0 >> 28 == 0xF && ((data_line.data_word0 & 0xF) == 0x5 || (data_line.data_word0 & 0xF) == 0x2)) && 
                    (data_line.data_word1 >> 28 == 0xF && ((data_line.data_word1 & 0xF) == 0x5 || (data_line.data_word1 & 0xF) == 0x2)) &&
                    (data_line.data_word2 >> 28 == 0xF && ((data_line.data_word2 & 0xF) == 0x5 || (data_line.data_word2 & 0xF) == 0x2)) &&
                    (data_line.data_word3 >> 28 == 0xF && ((data_line.data_word3 & 0xF) == 0x5 || (data_line.data_word3 & 0xF) == 0x2))
                ) {
                    if (!first_data_found) {    // if this is the first data frame for the trigger, assume it is mg index 0
                        first_data_found = true;
                        first_valid_data_timestamp = data_timestamp;
                        current_lead_data_timestamp = data_timestamp;
                        auto first_data_shift = data_timestamp - trg_timestamp;
                        first_mg_flag = true;
                        // if (true) {
                        //     spdlog::info("  First: vldb_id={}, bx=0x{:03X}, ob=0x{:03X}, first_data_shift={}",
                        //                  (int)data_line.header_vldb_id, data_bx, data_ob, first_data_shift);
                        // }
                    } else {
                        first_mg_flag = false;
                        skip_line_mask = false;
                        // to find out which machinegun index this data belongs to
                        auto header_shift = data_timestamp - first_valid_data_timestamp;
                        current_lead_data_timestamp = data_timestamp;
                        int diff_mg = header_shift % timestamp_diff_mg;
                        if (diff_mg > timestamp_diff_mg_tolerance || diff_mg < -timestamp_diff_mg_tolerance) {
                            if (_verbose >= VERBOSE_LEVEL_WARN) {
                                spdlog::warn("    Warning: header_shift={}, current_lead_data_timestamp=0x{:08X}, first_valid_data_timestamp=0x{:08X} differs from expected mg multiple (mg {})",
                                             header_shift, current_lead_data_timestamp, first_valid_data_timestamp, machinegun_index);
                                spdlog::warn("      vldb_id={}, bx=0x{:03X}, ob=0x{:03X}, line_index={}",
                                             (int)data_line.header_vldb_id, data_bx, data_ob, data_line_index);
                            }
                            _n_trigger_header_ts_error++;
                            // print the bx and ob in hex
                            // if (true) {
                            //     spdlog::info("  (b)Data: vldb_id={}, bx=0x{:03X}, ob=0x{:03X}, header_shift={}",
                            //                  (int)data_line.header_vldb_id, data_bx, data_ob, header_shift);
                            // }
                            continue;   
                        } else {
                            machinegun_index = std::round(static_cast<double>(header_shift) / timestamp_diff_mg);
                            candidate_header_line_index = data_line_index;
                            if (machinegun_index > max_machinegun_index) {
                                max_machinegun_index = 16;
                            }
                            // if (true) {
                            //     spdlog::info("  (g)Data: vldb_id={}, bx=0x{:03X}, ob=0x{:03X}, header_shift={}",
                            //                  (int)data_line.header_vldb_id, data_bx, data_ob, header_shift);
                            // }
                            if (_verbose >= VERBOSE_LEVEL_INFO) {
                                spdlog::info("    Info: header_shift={} corresponds to machinegun index {}",
                                             header_shift, machinegun_index);
                            }
                        }
                        if (_verbose >= VERBOSE_LEVEL_INFO) {
                            spdlog::info("  Data: vldb_id={}, bx=0x{:03X}, ob=0x{:03X}, header_shift={}",
                                         (int)data_line.header_vldb_id, data_bx, data_ob, header_shift);
                        }
                    } // end of mg index assignment
                    looking_for_header = false;
                    line_candidate_counter = 1;

                    frame_words[0].clear();
                    frame_words[1].clear();
                    frame_words[2].clear();
                    frame_words[3].clear();
                    frame_words[4].clear();
                    frame_words[5].clear();

                    frame_words[0].push_back(data_line.data_word0);
                    frame_words[1].push_back(data_line.data_word1);
                    frame_words[2].push_back(data_line.data_word2);
                    frame_words[3].push_back(data_line.data_word3);
                    frame_words[4].push_back(data_line.data_word4);
                    frame_words[5].push_back(data_line.data_word5);
                }
                else {
                    line_skipped_cnt++;
                    if (is_idle_line==false && skip_line_mask==false) {
                        // if ((int)data_line.header_vldb_id == 1) {
                        //     spdlog::info("    Debug: skipping line without header, vldb_id={}, bx=0x{:03X}, ob=0x{:03X}, line_index={}",
                        //                  (int)data_line.header_vldb_id, data_bx, data_ob, data_line_index);
                        //     spdlog::info("      data_word0=0x{:08X}, data_word1=0x{:08X}, data_word2=0x{:08X}, data_word3=0x{:08X}",
                        //                  data_line.data_word0, data_line.data_word1, data_line.data_word2, data_line.data_word3);
                        // }
                    }
                }
            } // end of looking for header
            else {
                line_candidate_counter++;
                frame_words[0].push_back(data_line.data_word0);
                frame_words[1].push_back(data_line.data_word1);
                frame_words[2].push_back(data_line.data_word2);
                frame_words[3].push_back(data_line.data_word3);
                frame_words[4].push_back(data_line.data_word4);
                frame_words[5].push_back(data_line.data_word5);
                    if (line_candidate_counter >= 40) {
                        // give up looking for header
                        auto last_data_shift = data_timestamp - current_lead_data_timestamp;
                        auto diff_sample = std::abs(static_cast<int>(last_data_shift) - timestamp_diff_sample);
                        if (diff_sample > timestamp_diff_sample_tolerance) {
                            if (first_mg_flag){
                                // if this is the first mg, do not count as error
                                if (true) {
                                    spdlog::warn("    Warning: first mg, last_data_shift={} differs from expected {} (mg {})",
                                                 last_data_shift, timestamp_diff_sample, machinegun_index);
                                }
                            }
                            if (true) {
                                spdlog::warn("    Warning: last_data_shift={} differs from expected {} (mg {})",
                                             last_data_shift, timestamp_diff_sample, machinegun_index);
                            }
                            _n_trigger_40line_error++;
                            // roll back the data lines
                            // spdlog::info(".   lineindex back from {} to {}", data_line_index, candidate_header_line_index);
                            data_line_index = candidate_header_line_index; // -1 because of the for loop increment
                            skip_line_mask = true;
                        } else {
                            // * handle the 160-byte frame
                            UShort_t fpga_id = static_cast<UShort_t>(data_line.header_vldb_id);
                            std::vector<UInt_t> val0_array;
                            std::vector<UInt_t> val1_array;
                            std::vector<UInt_t> val2_array;
                            std::vector<Bool_t> tc_array;
                            std::vector<Bool_t> tp_array;
                            std::vector<UInt_t> daqh_list;
                            std::vector<UInt_t> crc_list;
                            auto parse_success = parse_160byte_frame(frame_words, val0_array, val1_array, val2_array, tc_array, tp_array, daqh_list, crc_list, true, _verbose);
                            if (parse_success) {
                                _n_legal_samples++;
                                // if (true) {
                                //     spdlog::info("    Frame {}: fpga_id={}", _trg_index, (int)fpga_id);
                                //     spdlog::info("    Frame {}: Val0: {}", _trg_index, fmt::join(val0_array, ", "));
                                //     spdlog::info("    Frame {}: Val1: {}", _trg_index, fmt::join(val1_array, ", "));
                                //     spdlog::info("    Frame {}: Val2: {}", _trg_index, fmt::join(val2_array, ", "));
                                //     spdlog::info("    Frame {}: TC: {}", _trg_index, fmt::join(tc_array, ", "));
                                //     spdlog::info("    Frame {}: TP: {}", _trg_index, fmt::join(tp_array, ", "));
                                //     spdlog::info("    Frame {}: DAQH: 0x{:08X}, 0x{:08X}, 0x{:08X}, 0x{:08X}", _trg_index, daqh_list[0], daqh_list[1], daqh_list[2], daqh_list[3]);
                                //     spdlog::info("    Frame {}: CRC: 0x{:08X}, 0x{:08X}, 0x{:08X}, 0x{:08X}", _trg_index, crc_list[0], crc_list[1], crc_list[2], crc_list[3]);
                                // }
                                // fill branches
                                *(_branch_fpga_id) = fpga_id;
                                *(_branch_timestamp) = data_timestamp;
                                for (int i = 0; i < 4; i++) {
                                    *(_branch_daqh_list + i) = daqh_list[i];
                                    *(_branch_crc32_list + i) = crc_list[i];
                                }
                                const int n_chan = std::min<int>(FPGA_CHANNEL_NUMBER, (int)val0_array.size());
                                for (int i = 0; i < n_chan; ++i) {
                                    _branch_tc_list[i]   = tc_array[i];
                                    _branch_tp_list[i]   = tp_array[i];
                                    // if (val0_array[i] < val0_array_min[i]) {
                                    //     val0_array_min[i] = val0_array[i];
                                    // }
                                    _branch_val0_list[i] = val0_array[i];
                                    _branch_val1_list[i] = val1_array[i];
                                    _branch_val2_list[i] = val2_array[i];
                                }
                                // for (int i = 0; i < FPGA_CHANNEL_NUMBER; i++) {
                                //     *(_branch_tc_list + i) = tc_array[i];
                                //     *(_branch_tp_list + i) = tp_array[i];
                                //     *(_branch_val0_list + i) = val0_array[i];
                                //     *(_branch_val1_list + i) = val1_array[i];
                                //     *(_branch_val2_list + i) = val2_array[i];
                                // }
                                *(_branch_last_heartbeat) = static_cast<UInt_t>(last_data_shift);
                                _output_tree->Fill();
                            } else {
                                _n_trigger_160byte_error++;
                            } // end of parse fail
                            if (_verbose >= VERBOSE_LEVEL_INFO) {
                                spdlog::info("  Data: vldb_id={}, bx=0x{:03X}, ob=0x{:03X}, last_data_shift={}",
                                            (int)data_line.header_vldb_id, data_bx, data_ob, last_data_shift);
                            }
                        } // end of legal frame
                        looking_for_header = true;
                    } // end of if line_candidate_counter >= 40
            } // end of processing one data line
        } // end of for each data line loop
        
        // spdlog::info("  Finished processing trigger, total data lines: {}, legal samples: {}, 40-line errors: {}, 160-byte errors: {}, header timestamp errors: {}, max mg index: {}, idle lines: {}, skipped lines: {}, L1A sent: {}",
        //              _data_lines.size(), _n_legal_samples, _n_trigger_40line_error, _n_trigger_160byte_error, _n_trigger_header_ts_error, max_machinegun_index, idle_line_cnt, line_skipped_cnt, L1A_cmd_counter);
        _data_lines.clear();
    } // for each vldb+
}

int main(int argc, char** argv) {
    // * Parse command line arguments
    // * ----------------------------------------------------------------------
    std::string script_full_name = argv[0];
    auto script_name = script_full_name.substr(
        script_full_name.find_last_of("/\\") + 1,
        script_full_name.find_last_of('.') - script_full_name.find_last_of("/\\") - 1
    );
    spdlog::info("Running script: {}", script_name);
    spdlog::info("----------------------------------------");
    cxxopts::Options options(script_name, "Convert raw binary data to ROOT format");
    options.add_options()
        ("i,input", "Input file", cxxopts::value<std::string>())
        ("o,output", "Output file", cxxopts::value<std::string>())
        ("n,number", "Number of vldb", cxxopts::value<int>()->default_value("2"))
        ("v,verbose", "Verbose mode", cxxopts::value<int>()->default_value("2"))
        ("h,help", "Print help");
    auto parsed_opts = options.parse(argc, argv);
    if (parsed_opts.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }
    // check if the input file is provided
    std::string input_file;
    if (!parsed_opts.count("input")) {
        // give default input file name
        spdlog::warn("Using the example input file.");
        input_file = "data/internal_examples/2025_09_25_cru_dma_log_example";
    } else {
        input_file = parsed_opts["input"].as<std::string>();
    }
    int n_vldb = parsed_opts["number"].as<int>();
    spdlog::info("Number of VLDB: {}", n_vldb);
    std::string input_file_name = input_file.substr(
        input_file.find_last_of("/\\") + 1,
        input_file.find_last_of('.') - input_file.find_last_of("/\\") - 1
    );
    // check if the output file is provided
    std::string output_file;
    if (!parsed_opts.count("output")) {
        // give default output file name
        spdlog::warn("Using the default output file.");
        output_file = "dump/" + script_name + "/" + input_file_name + ".root";
    } else {
        output_file = parsed_opts["output"].as<std::string>();
    }
    auto global_verbose = parsed_opts["verbose"].as<int>();
    // Create output file
    TFile *output_root = TFile::Open(output_file.c_str(), "RECREATE");
    if (output_root == nullptr || output_root->IsZombie()) {
        spdlog::error("Failed to create output file {}", output_file);
        return 1;
    }
    TTree* output_tree = new TTree("data_tree", "Data tree");
    output_tree->SetDirectory(output_root);
    auto *branch_fpga_id    = new UShort_t();   // 0 - 8
    auto *branch_timestamp  = new ULong64_t();  // 64 bits
    auto *branch_daqh_list  = new UInt_t[4];    // 32 bits
    auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER];
    auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER];
    auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER];
    auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER];
    auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER];
    auto *branch_crc32_list = new UInt_t[4];
    auto *branch_last_heartbeat = new UInt_t();

    // output_tree->Branch("fpga_id", branch_fpga_id, "fpga_id/s");
    // output_tree->Branch("timestamp", branch_timestamp, "timestamp/l");
    // output_tree->Branch("daqh_list", branch_daqh_list, "daqh_list[4]/i");
    // output_tree->Branch("tc_list", branch_tc_list, ("tc_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/O").c_str());
    // output_tree->Branch("tp_list", branch_tp_list, ("tp_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/O").c_str());
    // output_tree->Branch("val0_list", branch_val0_list, ("val0_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str());
    // output_tree->Branch("val1_list", branch_val1_list, ("val1_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str());
    // output_tree->Branch("val2_list", branch_val2_list, ("val2_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str());
    // output_tree->Branch("crc32_list", branch_crc32_list, "crc32_list[4]/i");
    // output_tree->Branch("last_heartbeat", branch_last_heartbeat, "last_heartbeat/i");

    output_tree->Branch("fpga_id",          branch_fpga_id,          "fpga_id/s");               // was /s
    output_tree->Branch("timestamp",        branch_timestamp,        "timestamp/l");             // was /l
    output_tree->Branch("daqh_list",        branch_daqh_list,        "daqh_list[4]/i");          // was /i
    output_tree->Branch("tc_list",          branch_tc_list,          ("tc_list["  + std::to_string(FPGA_CHANNEL_NUMBER) + "]/O").c_str());
    output_tree->Branch("tp_list",          branch_tp_list,          ("tp_list["  + std::to_string(FPGA_CHANNEL_NUMBER) + "]/O").c_str());
    output_tree->Branch("val0_list",        branch_val0_list,        ("val0_list["+ std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str()); // was /i
    output_tree->Branch("val1_list",        branch_val1_list,        ("val1_list["+ std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str()); // was /i
    output_tree->Branch("val2_list",        branch_val2_list,        ("val2_list["+ std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str()); // was /i
    output_tree->Branch("crc32_list",       branch_crc32_list,       "crc32_list[4]/i");         // was /i
    output_tree->Branch("last_heartbeat",   branch_last_heartbeat,   "last_heartbeat/i");        // was /i (this one is fine either way)

    // * ----------------------------------------------------------------------

    bp::TailOptions opts;
    opts.poll_ms = 1;                 // check for new data every 1 ms
    opts.read_chunk = 1u << 20;        // 1 MB read chunk
    opts.inactivity_timeout_ms = 10000; // exit if no new data for 1 second

    // --- data qc variables ---
    bool qc_rdh0_found = false;

    std::vector<std::byte> stash;
    std::size_t total_bytes       = 0;
    std::size_t total_lines       = 0;
    std::size_t n_packets         = 0;
    std::size_t n_heartbeats      = 0;
    std::size_t n_syncs           = 0;

    std::size_t n_rdh_l0           = 0;
    std::size_t n_rdh_l1           = 0;
    std::size_t n_data_lines       = 0;
    std::size_t n_trg_lines        = 0;

    std::size_t n_trigger_processed = 0;
    std::size_t n_trigger_illegal_40_timestamps = 0;
    std::size_t n_trigger_160byte_parse_failures = 0;
    std::size_t n_trigger_header_ts_errors = 0;
    std::size_t n_legal_samples = 0;

    const uint32_t idle_pattern = 0xACCCCCCC;
    const uint32_t L1A_pattern  = 0x4B00004B;
    auto t_start = std::chrono::steady_clock::now();

    // std::vector<std::vector<uint32_t>> vldb_timestamps(n_vldb);
    std::vector<std::vector<bp::DataLine>> vldb_data_lines_buffer(n_vldb);
    // set of vldb id
    std::unordered_set<uint8_t> vldb_ids;
    bp::TrgLine last_trg_line;

    bp::StreamParser parser(
        /* on_packet */
        [&](const bp::Packet& pkt) {
            (void)pkt;
            n_packets++;
        },
        /* on_heartbeat */
        [&](const bp::Heartbeat&) {
            n_heartbeats++;
        },
        /* on_sync */
        [&](std::span<const std::byte>) {
        },
        /* on_rdh_l0 */
        [&](const bp::RDH_L0& rdh, std::span<const std::byte> raw){
            // use rdh / raw if needed
            (void)rdh; (void)raw;
            n_rdh_l0++;
            if (qc_rdh0_found == false) {
                qc_rdh0_found = true;
            } else {
                spdlog::warn("[L{}]: Consecutive RDH L0 headers", total_lines);
            }
        },
        /* on_rdh_l1 */
        [&](const bp::RDH_L1& rdh, std::span<const std::byte> raw){
            (void)rdh; (void)raw;
            n_rdh_l1++;
            if (qc_rdh0_found) {
                qc_rdh0_found = false;
            } else {
                spdlog::warn("[L{}]: RDH L1 header found without RDH L0.", total_lines);
            }
        },
        /* on_data_line */
        [&](const bp::DataLine& line, std::span<const std::byte> raw){
            // (void)line; (void)raw;
            n_data_lines++;
            auto &orbit_id  = line.ob_cnt;
            auto &bc_id     = line.bx_cnt;
            auto &vldb_id   = line.header_vldb_id;
            auto unique_timestamp = (static_cast<uint64_t>(orbit_id) * BX_PER_ORBIT) + bc_id;
            // check if vldb_id is in the set
            int vldb_index = -1;
            if (vldb_ids.find(vldb_id) == vldb_ids.end()) {
                vldb_ids.insert(vldb_id);
                spdlog::info("Found new VLDB ID: {} (total {})", (int)vldb_id, vldb_ids.size());
                if (vldb_ids.size() > static_cast<size_t>(n_vldb)) {
                    spdlog::warn("Number of unique VLDB IDs ({}) exceeds the specified number of VLDBs ({}).",
                                 vldb_ids.size(), n_vldb);
                }
                vldb_index = vldb_ids.size() - 1;
            } else {
                vldb_index = std::distance(vldb_ids.begin(), vldb_ids.find(vldb_id));
            }
            // append to timestamp_array
            // vldb_timestamps[vldb_index].push_back(static_cast<uint32_t>(unique_timestamp & 0xFFFFFFFF));
            vldb_data_lines_buffer[vldb_index].push_back(line);
        },
        /* on_trg_line */
        [&](const bp::TrgLine& line, std::span<const std::byte> raw){
            (void)line; (void)raw;
            n_trg_lines++;
            
            // auto &trg_bx_id = line.bx_cnt;
            // auto &trg_orbit_id = line.ob_cnt;
            // auto trg_unique_timestamp = (static_cast<uint64_t>(trg_orbit_id) << 12) | trg_bx_id;
            
            // for (size_t i = 0; i < vldb_timestamps.size(); i++) {
            //     if (!vldb_timestamps[i].empty()) {
            //         auto min_timestamp = *std::min_element(vldb_timestamps[i].begin(), vldb_timestamps[i].end());
            //         auto adjusted_timestamp = vldb_timestamps[i].back() - min_timestamp;
            //         spdlog::info("  VLDB {}: trg timestamp: {}, adjusted: {}", i, vldb_timestamps[i].back(), adjusted_timestamp);
            //         vldb_timestamps[i].clear();
            //     } else {
            //         spdlog::info("  VLDB {}: no timestamps recorded", i);
            //     }
            // }

            // process trigger data
            if (n_trg_lines > 1){
                n_trigger_processed++;
                trigger_data_processing(const_cast<bp::TrgLine&>(last_trg_line), vldb_data_lines_buffer, n_trg_lines, 
                                        branch_fpga_id, branch_timestamp, branch_daqh_list, branch_tc_list, branch_tp_list,
                                        branch_val0_list, branch_val1_list, branch_val2_list, branch_crc32_list, branch_last_heartbeat,
                                        output_tree, n_trigger_illegal_40_timestamps, n_trigger_160byte_parse_failures, n_trigger_header_ts_errors, n_legal_samples);
            }
            last_trg_line = line;

            // spdlog::info("  trg bx_id: {}, orbit_id: {}, unique_timestamp: {}",
            //              trg_bx_id, trg_orbit_id, trg_unique_timestamp);
        }
    );

    bp::tail_growing_file(input_file, opts, [&](std::span<const std::byte> chunk) {
        total_bytes += chunk.size();

        // merge with stash to ensure full 32-byte lines
        if (!stash.empty()) {
            const std::size_t need = bp::ByteCursor::kLineSize - stash.size();
            if (chunk.size() >= need) {
                std::vector<std::byte> one(stash.begin(), stash.end());
                one.insert(one.end(), chunk.begin(), chunk.begin() + need);
                parser.feed(one);
                chunk = chunk.subspan(need);
                stash.clear();
            } else {
                stash.insert(stash.end(), chunk.begin(), chunk.end());
                return;
            }
        }

        const std::size_t remainder = chunk.size() % bp::ByteCursor::kLineSize;
        const auto main_part = chunk.first(chunk.size() - remainder);
        if (!main_part.empty()) parser.feed(main_part);
        if (remainder) {
            stash.assign(chunk.end() - remainder, chunk.end());
        }

        // show simple progress every ~1 MB
        if (total_bytes % (1 << 20) < bp::ByteCursor::kLineSize) {
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - t_start);
            std::cout << "[Progress] "
                      << total_bytes / 1e6 << " MB read, "
                      << total_lines << " lines parsed, "
                      << "time elapsed: " << elapsed.count() << " ms\r"
                      << std::flush;
        }
    });

    auto t_end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);

    spdlog::info("=== Parsing completed =============================");
    spdlog::info("Total bytes read: {} ({:.2f} MB)", total_bytes, total_bytes / 1e6);
    spdlog::info("Total lines parsed: {}", total_lines);
    spdlog::info("Total packets: {}", n_packets);
    spdlog::info("Total heartbeats: {}", n_heartbeats);
    spdlog::info("Total RDH L0: {}", n_rdh_l0);
    spdlog::info("Total RDH L1: {}", n_rdh_l1);
    spdlog::info("Total data lines: {}", n_data_lines);
    spdlog::info("Total trg lines: {}", n_trg_lines);
    spdlog::info("Total syncs: {}", n_syncs);
    spdlog::info("Time elapsed: {} ms", elapsed.count());
    spdlog::info("Average speed: {:.2f} MB/s", total_bytes / 1e6 / (elapsed.count() / 1000.0));
    spdlog::info("Total triggers processed: {}", n_trigger_processed);
    spdlog::info("Illegal 40 timestamps: {}", n_trigger_illegal_40_timestamps);
    spdlog::info("160-byte parse failures: {}", n_trigger_160byte_parse_failures);
    spdlog::info("Header timestamp errors: {}", n_trigger_header_ts_errors);
    const double avg_per_trigger = n_trigger_processed > 0
        ? static_cast<double>(n_legal_samples) / static_cast<double>(n_trigger_processed)
        : 0.0;
    spdlog::info("Legal samples: {} ({:.4g} per trigger)", n_legal_samples, avg_per_trigger);
    spdlog::info("===================================================");

    spdlog::info("Writing output to {}", output_file);

    int chn_map_y_bin = 512;
    // x: channel 0-152
    // y: adc value 0-1024
    TCanvas* c_chn_map = new TCanvas("c_chn_map", "Channel Map", 1200, 2400);
    TH2D* adc_chn_map_0 = new TH2D("adc_chn_map", "ADC Channel Map;Channel;ADC Value", FPGA_CHANNEL_NUMBER, 0, FPGA_CHANNEL_NUMBER, chn_map_y_bin, 0, 1024);
    TH2D* tot_chn_map_0 = new TH2D("tot_chn_map", "ToT Channel Map;Channel;ToT Value", FPGA_CHANNEL_NUMBER, 0, FPGA_CHANNEL_NUMBER, chn_map_y_bin, 0, 1024);
    TH2D* toa_chn_map_0 = new TH2D("toa_chn_map", "ToA Channel Map;Channel;ToA Value", FPGA_CHANNEL_NUMBER, 0, FPGA_CHANNEL_NUMBER, chn_map_y_bin, 0, 1024);
    TH2D* adc_chn_map_1 = new TH2D("adc_chn_map", "ADC Channel Map;Channel;ADC Value", FPGA_CHANNEL_NUMBER, 0, FPGA_CHANNEL_NUMBER, chn_map_y_bin, 0, 1024);
    TH2D* tot_chn_map_1 = new TH2D("tot_chn_map", "ToT Channel Map;Channel;ToT Value", FPGA_CHANNEL_NUMBER, 0, FPGA_CHANNEL_NUMBER, chn_map_y_bin, 0, 1024);
    TH2D* toa_chn_map_1 = new TH2D("toa_chn_map", "ToA Channel Map;Channel;ToA Value", FPGA_CHANNEL_NUMBER, 0, FPGA_CHANNEL_NUMBER, chn_map_y_bin, 0, 1024);


    // Fill the histograms
    Long64_t total_entries = output_tree->GetEntries();
    for (Long64_t entry = 0; entry < total_entries; ++entry)
    {
        output_tree->GetEntry(entry);
        if (branch_fpga_id[0] == 0)
        {
            for (int ch = 0; ch < FPGA_CHANNEL_NUMBER; ++ch)
            {
                adc_chn_map_0->Fill(ch, branch_val0_list[ch]);
                tot_chn_map_0->Fill(ch, branch_val1_list[ch]);
                toa_chn_map_0->Fill(ch, branch_val2_list[ch]);
            }
            continue;
        }
        if (branch_fpga_id[0] == 1)
        {
            for (int ch = 0; ch < FPGA_CHANNEL_NUMBER; ++ch)
            {
                adc_chn_map_1->Fill(ch, branch_val0_list[ch]);
                tot_chn_map_1->Fill(ch, branch_val1_list[ch]);
                toa_chn_map_1->Fill(ch, branch_val2_list[ch]);
            }
            continue;
        }
    }

    adc_chn_map_0->SetStats(0);
    tot_chn_map_0->SetStats(0);
    toa_chn_map_0->SetStats(0);
    adc_chn_map_0->SetTitle("");
    tot_chn_map_0->SetTitle("");
    toa_chn_map_0->SetTitle("");

    adc_chn_map_1->SetStats(0);
    tot_chn_map_1->SetStats(0);
    toa_chn_map_1->SetStats(0);
    adc_chn_map_1->SetTitle("");
    tot_chn_map_1->SetTitle("");
    toa_chn_map_1->SetTitle("");

    c_chn_map->Divide(2, 3);
    c_chn_map->cd(1);
    adc_chn_map_0->Draw("COLZ");
    c_chn_map->cd(3);
    tot_chn_map_0->Draw("COLZ");
    c_chn_map->cd(5);
    toa_chn_map_0->Draw("COLZ");
    c_chn_map->cd(2);
    adc_chn_map_1->Draw("COLZ");
    c_chn_map->cd(4);
    tot_chn_map_1->Draw("COLZ");
    c_chn_map->cd(6);
    toa_chn_map_1->Draw("COLZ");

    // add annotation to the plots
    c_chn_map->cd(1);
    TLatex latex_adc;
    latex_adc.SetNDC();
    latex_adc.SetTextSize(0.05);
    latex_adc.DrawLatex(0.10, 0.92, ("ADC (VLDB 0), Input: " + input_file_name).c_str());
    c_chn_map->cd(3);
    TLatex latex_tot;
    latex_tot.SetNDC();
    latex_tot.SetTextSize(0.05);
    latex_tot.DrawLatex(0.10, 0.92, ("ToT (VLDB 0), Input: " + input_file_name).c_str());
    c_chn_map->cd(5);
    TLatex latex_toa;
    latex_toa.SetNDC();
    latex_toa.SetTextSize(0.05);
    latex_toa.DrawLatex(0.10, 0.92, ("ToA (VLDB 0), Input: " + input_file_name).c_str());
    c_chn_map->cd(2);
    TLatex latex_adc1;
    latex_adc1.SetNDC();
    latex_adc1.SetTextSize(0.05);
    latex_adc1.DrawLatex(0.10, 0.92, ("ADC (VLDB 1), Input: " + input_file_name).c_str());
    c_chn_map->cd(4);
    TLatex latex_tot1;
    latex_tot1.SetNDC();
    latex_tot1.SetTextSize(0.05);
    latex_tot1.DrawLatex(0.10, 0.92, ("ToT (VLDB 1), Input: " + input_file_name).c_str());
    c_chn_map->cd(6);
    TLatex latex_toa1;
    latex_toa1.SetNDC();
    latex_toa1.SetTextSize(0.05);
    latex_toa1.DrawLatex(0.10, 0.92, ("ToA (VLDB 1), Input: " + input_file_name).c_str());

    output_root->cd();
    c_chn_map->Write();

    std::string _legal_fpga_id_list_str = "";
    for (const auto& id : vldb_ids) {
        _legal_fpga_id_list_str += std::to_string(static_cast<int>(id)) + " ";
    }
    TNamed("Rootifier_legal_fpga_id_list", _legal_fpga_id_list_str.c_str()).Write();
    output_root->Write();
    output_root->Close();

    c_chn_map->Close();

    return 0;
}