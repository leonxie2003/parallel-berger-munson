/**
 * Sequential Berger-Munson program
 *
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include "bm_seq.h"
#include "parse_fasta.h"
#include "align.h"
#include "bm_utils.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <string>
#include <vector>

#include <unistd.h>

int main(int argc, char *argv[]) {
    const auto start_time = std::chrono::steady_clock::now();

    // Parse cmd line args
    std::string input_filename;
    std::string output_filename;
    int random_mode = DEVICERANDOM;

    int opt;
    while((opt = getopt(argc, argv, "i:o:r:")) != -1) {
        switch (opt) {
            case 'i':
                input_filename = optarg;
                break;
            case 'o':
                output_filename = optarg;
                break;
            case 'r':
                if (optarg[0] == 'R')
                    random_mode = DEVICERANDOM;
                else if (optarg[0] == 'P')
                    random_mode = PSEUDORANDOM;
                break;
        default:
            std::cerr << "Usage: " << argv[0] << " -i input_filename -o output_filename\n";
            exit(EXIT_FAILURE);
        }
    }

    if (empty(input_filename) || empty(output_filename)) {
        std::cerr << "Usage: " << argv[0] << " -i input_filename -o output_filename\n";
        exit(EXIT_FAILURE);
    }

    // Parse FASTA file
    std::cout << "Input file: " << input_filename << "\n";
    std::vector<fasta_seq_t> fasta_seqs = parse_fasta(input_filename);

    // Initialize program state
    align_params_t params{};

    int glbl_idx = 0; // Berger-Munson iteration number

    seq_group_t cur_alnmt = naiive_alnmt(fasta_seqs);
    int best_score = INT_MIN;
    int best_glbl_idx = -1;

    int num_seqs = fasta_seqs.size();
    int num_partns = num_seqs + (num_seqs * (num_seqs - 1)) / 2;

    std::string accept_reject_chain = "";

    while (glbl_idx - (best_glbl_idx + 1) < num_partns) {
        // Partition into two groups
        seq_group_t group1{};
        seq_group_t group2{};
        select_partn(cur_alnmt, glbl_idx, random_mode, group1, group2);

        remove_glbl_gaps(group1);
        remove_glbl_gaps(group2);

        // Compute alignment between two groups
        gap_pos_t gap_pos{};
        int cur_score = align_groups(group1, group2, params, gap_pos);

        if (cur_score > best_score) {
            // Update program state
            best_score = cur_score;
            best_glbl_idx = glbl_idx;
            cur_alnmt = update_alnmt(group1, group2, gap_pos);
            accept_reject_chain += 'A';
        } else {
            accept_reject_chain += 'R';
        }

        glbl_idx++;
    }

    const auto end_time = std::chrono::steady_clock::now();
    const double runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    std::cout << "Ran for " << glbl_idx << " iterations.\n";
    std::cout << "Runtime (sec): " << runtime << "\n";
    std::cout << "Alignment score: " << best_score << "\n";
    std::cout << "Accepts and rejects: " << accept_reject_chain << "\n";

    std::ofstream fout(output_filename);

    fout << "Final alignment (score = " << best_score << "):\n";
    fout << "Accepts and rejects: " << accept_reject_chain << "\n";
    fout << "\n\n";
    for (seq_t seq : cur_alnmt) {
        fout << "seq " << std::setw(3) << seq.id << ": ";
        fout << seq.data << "\n";
    }
}
