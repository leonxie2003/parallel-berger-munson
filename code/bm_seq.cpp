/**
 * Sequential Berger-Munson program
 *
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include "bm_seq.h"
#include "parse_fasta.h"
#include "align.h"
#include "bm_utils.h"

#include <iostream>
#include <vector>
#include <string>

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

    for (fasta_seq_t seq : fasta_seqs) {
        print_fasta_seq(seq);
    }

    /* --- Berger-Munson algorithm --- */

    // construct naiive alignment (add gaps until all same length)
    seq_group_t curr_alnmt = naiive_alnmt(fasta_seqs);
    std::cout << "Naiive alignment:\n";
    for (seq_t seq : curr_alnmt) {
        std::cout << seq << "\n\n";
    }

    // iteratively improve alignment

    // TODO Replace with "q reject in a row", right now is just constant num of
    // iterations
    int glbl_idx = 0;
    while (glbl_idx < 1) {
        // partition into two groups
        seq_group_t group1{};
        seq_group_t group2{};
        select_partn(curr_alnmt, glbl_idx, group1, group2);

        std::cout << "Group 1:\n";
        for (seq_t seq : group1) {
            std::cout << seq << "\n\n";
        }

        std::cout << "Group 2:\n";
        for (seq_t seq : group2) {
            std::cout << seq << "\n\n";
        }

        remove_glbl_gaps(group1);
        remove_glbl_gaps(group2);

        glbl_idx++;
    }
}
