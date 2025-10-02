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
        spdlog::warn("No input file provided. Using the example input file.");
        input_file = "data/internal_examples/2025_09_25_cru_dma_log_example";
    } else {
        input_file = parsed_opts["input"].as<std::string>();
    }

    // * ----------------------------------------------------------------------

    std::vector<std::byte> stash;
    std::size_t total_bytes       = 0;
    std::size_t total_lines       = 0;
    std::size_t n_packets         = 0;
    std::size_t n_heartbeats      = 0;
    std::size_t n_syncs           = 0;
    auto t_start = std::chrono::steady_clock::now();

    bp::StreamParser parser(
        /* on_packet */
        [&](const bp::Packet& pkt) {
            ++n_packets;
            total_lines += pkt.block.size() / bp::ByteCursor::kLineSize;
        },
        /* on_heartbeat */
        [&](const bp::Heartbeat&) {
            ++n_heartbeats;
            total_lines += 2;
        },
        /* on_sync */
        [&](std::span<const std::byte>) {
            ++n_syncs;
            total_lines += 1;
        }
    );

    bp::tail_growing_file(input_file, {}, [&](std::span<const std::byte> chunk) {
        total_bytes += chunk.size();

        // merge with stash to ensure full 40-byte lines
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

    spdlog::info("\nParsing completed.");
    spdlog::info("Total bytes read: {} ({:.2f} MB)", total_bytes, total_bytes / 1e6);
    spdlog::info("Total lines parsed: {}", total_lines);
    spdlog::info("Total packets: {}", n_packets);
    spdlog::info("Total heartbeats: {}", n_heartbeats);
    spdlog::info("Total syncs: {}", n_syncs);
    spdlog::info("Time elapsed: {} ms", elapsed.count());
    spdlog::info("Average speed: {:.2f} MB/s", total_bytes / 1e6 / (elapsed.count() / 1000.0));

    return 0;
}