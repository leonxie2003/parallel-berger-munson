/**
 * Sequential Berger-Munson program
 *
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include "bm_seq.h"
#include "parse_fasta.h"
#include "align.h"
#include "bm_utils.h"

#include <iomanip>
#include <iostream>
#include <limits.h>
#include <string>
#include <vector>

#include <unistd.h>

#define ACCEPT true
#define REJECT false

int main(int argc, char *argv[]) {
    /* --- parse cmd line args --- */
    std::string input_filename;
    std::string output_filename;

    int opt;
    while((opt = getopt(argc, argv, "i:o:")) != -1) {
        switch (opt) {
            case 'i':
                input_filename = optarg;
                break;
            case 'o':
                output_filename = optarg;
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

    /* --- parse FASTA file --- */
    std::cout << "Input file: " << input_filename << "\n";
    std::vector<fasta_seq_t> fasta_seqs = parse_fasta(input_filename); // TODO use a struct for seqs to maintain id & desc

    /* TODO remove debug prints
    for (fasta_seq_t fasta_seq : fasta_seqs) {
        print_fasta_seq(fasta_seq);
    }
    */

    /* --- Berger-Munson algorithm --- */

    // construct naiive alignment (add gaps until all same length)
    seq_group_t cur_alnmt = naiive_alnmt(fasta_seqs);
    /* TODO remove debug prints
    std::cout << "Naiive alignment:\n";
    for (seq_t seq : cur_alnmt) {
        std::cout << seq.data << "\n";
    }
    std::cout << "\n";
    */

    // iteratively improve alignment

    // TODO Replace with "q reject in a row", right now is just constant num of
    // iterations
    int glbl_idx = 0;
    int best_score = INT_MIN;
    int best_glbl_idx = -1;
    int num_seqs = fasta_seqs.size();
    int num_partns = num_seqs + (num_seqs * (num_seqs - 1)) / 2;

    std::string accept_reject_chain = "";

    while (glbl_idx - (best_glbl_idx + 1) < num_partns) {
        // partition into two groups
        seq_group_t group1{};
        seq_group_t group2{};
        select_partn(cur_alnmt, glbl_idx, group1, group2);

        remove_glbl_gaps(group1);
        remove_glbl_gaps(group2);

        /* TODO remove debug prints
        std::cout << "===== glbl_indx: " << glbl_idx << " =====\n";
        std::cout << "Group 1: \n";
        for (seq_t group1_seq: group1) {
            std::cout << group1_seq.data << "\n";
        }
        std::cout << "\n";

        std::cout << "Group 2: \n";
        for (seq_t group2_seq: group2) {
            std::cout << group2_seq.data << "\n";
        }
        std::cout << "\n";
        */

        gap_pos_t gap_pos{};
        align_params_t params{};
        params.match_reward = 1;
        params.gap_penalty = -1;
        params.sub_penalty = 0;

        int cur_score = align_groups(group1, group2,  params, gap_pos);

        if (cur_score > best_score) {
            best_score = cur_score;
            best_glbl_idx = glbl_idx;
            cur_alnmt = update_alnmt(group1, group2, gap_pos);
            accept_reject_chain += 'A';
        } else {
            accept_reject_chain += 'R';
        }

        glbl_idx++;
    }

    std::cout << "Ran for " << glbl_idx << " iterations.\n";
    std::cout << "Accepts and rejects: " << accept_reject_chain << "\n";
    std::cout << "Final alignment (score = " << best_score << "):\n";
    for (seq_t seq : cur_alnmt) {
        std::cout << "seq " << std::setw(3) << seq.id << " :";
        std::cout << seq.data << "\n";
    }
}
