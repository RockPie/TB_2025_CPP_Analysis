#include <binparse/parser.hpp>
#include <binparse/tail.hpp>
#include <binparse/bytecursor.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include <iostream>
#include <vector>
#include <chrono>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include "cxxopts.hpp"

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
        ("o,output", "Output file", cxxopts::value<std::string>()->default_value("out.root"))
        ("n,number", "Number of vldb", cxxopts::value<int>()->default_value("2"))
        ("v,verbose", "Verbose mode", cxxopts::value<bool>()->default_value("false"))
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
    // * ----------------------------------------------------------------------

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

    std::size_t n_idle_words_0     = 0;
    std::size_t n_idle_words_1     = 0;
    std::size_t n_idle_words_2     = 0;
    std::size_t n_idle_words_3     = 0;
    std::size_t n_idle_words_4     = 0;
    std::size_t n_idle_words_5     = 0;

    std::size_t n_l1a_commands     = 0;

    const uint32_t idle_pattern = 0xACCCCCCC;
    const uint32_t L1A_pattern  = 0x4B00004B;
    std::size_t n_data_lines_between_triggers = 0;
    auto t_start = std::chrono::steady_clock::now();

    std::vector<std::vector<uint32_t>> vldb_timestamps(n_vldb);
    // set of vldb id
    std::unordered_set<uint8_t> vldb_ids;

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
        },
        /* on_rdh_l1 */
        [&](const bp::RDH_L1& rdh, std::span<const std::byte> raw){
            (void)rdh; (void)raw;
            n_rdh_l1++;
        },
        /* on_data_line */
        [&](const bp::DataLine& line, std::span<const std::byte> raw){
            // (void)line; (void)raw;
            n_data_lines++;
            n_data_lines_between_triggers++;
            auto &orbit_id  = line.ob_cnt;
            auto &bc_id     = line.bx_cnt;
            auto &vldb_id   = line.header_vldb_id;
            auto unique_timestamp = (static_cast<uint64_t>(orbit_id) << 12) | bc_id;
            if (line.data_word0 == idle_pattern)
               n_idle_words_0++;
            if (line.data_word1 == idle_pattern)
                n_idle_words_1++;
            if (line.data_word2 == idle_pattern)
                n_idle_words_2++;
            if (line.data_word3 == idle_pattern)
                n_idle_words_3++;
            if (line.data_word4 == idle_pattern)
                n_idle_words_4++;
            if (line.data_word5 == idle_pattern)
                n_idle_words_5++;

            if (line.data_word4 == L1A_pattern)
                n_l1a_commands++;
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
            vldb_timestamps[vldb_index].push_back(static_cast<uint32_t>(unique_timestamp & 0xFFFFFFFF));
        },
        /* on_trg_line */
        [&](const bp::TrgLine& line, std::span<const std::byte> raw){
            (void)line; (void)raw;
            n_trg_lines++;
            spdlog::info("Trigger received after {} data lines", n_data_lines_between_triggers);
            n_data_lines_between_triggers = 0;
            auto &trg_bx_id = line.bx_cnt;
            auto &trg_orbit_id = line.ob_cnt;
            auto trg_unique_timestamp = (static_cast<uint64_t>(trg_orbit_id) << 12) | trg_bx_id;
            
            for (size_t i = 0; i < vldb_timestamps.size(); i++) {
                if (!vldb_timestamps[i].empty()) {
                    auto min_timestamp = *std::min_element(vldb_timestamps[i].begin(), vldb_timestamps[i].end());
                    auto adjusted_timestamp = vldb_timestamps[i].back() - min_timestamp;
                    spdlog::info("  VLDB {}: trg timestamp: {}, adjusted: {}",
                                 i, vldb_timestamps[i].back(), adjusted_timestamp);
                    vldb_timestamps[i].clear();
                } else {
                    spdlog::info("  VLDB {}: no timestamps recorded", i);
                }
            }

            spdlog::info("  trg bx_id: {}, orbit_id: {}, unique_timestamp: {}",
                         trg_bx_id, trg_orbit_id, trg_unique_timestamp);
            // print idle words count
            spdlog::info("  Idle words between triggers: 0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}",
                         n_idle_words_0, n_idle_words_1, n_idle_words_2,
                         n_idle_words_3, n_idle_words_4, n_idle_words_5);
            n_idle_words_0 = 0;
            n_idle_words_1 = 0;
            n_idle_words_2 = 0;
            n_idle_words_3 = 0;
            n_idle_words_4 = 0;
            n_idle_words_5 = 0;

            // print L1A commands count
            spdlog::info("  L1A commands between triggers: {}", n_l1a_commands);
            n_l1a_commands = 0;
            // for each vldb, find the closest timestamp
        }
    );

    bp::TailOptions opts;
    opts.poll_ms = 50;                 // check for new data every 50 ms
    opts.read_chunk = 1u << 20;        // 1 MB read chunk
    opts.inactivity_timeout_ms = 1000; // exit if no new data for 1 second
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

    spdlog::info("=== Parsing completed ===");
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

    return 0;
}