/**
 * Sequential Berger-Munson program
 *
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include "bm_seq.h"
#include "parse_fasta.h"

#include <iostream>
#include <string>

#include <unistd.h>

int main(int argc, char *argv[]) {
    /* --- parse cmd line args --- */
    std::string input_filename;
    std::string output_filename;
    int num_procs;

    int opt;
    while((opt = getopt(argc, argv, "i:o:n:")) != -1) {
        switch (opt) {
            case 'i':
                input_filename = optarg;
                break;
            case 'o':
                output_filename = optarg;
                break;
            case 'n':
                num_procs = atoi(optarg);
                break;
        default:
            std::cerr << "Usage: " << argv[0] << " -i input_filename -o output_filename -n num_processors\n";
            exit(EXIT_FAILURE);
        }
    }

    if (empty(input_filename) || empty(output_filename) || num_procs < 1) {
        std::cerr << "Usage: " << argv[0] << " -i input_filename -o output_filename -n num_processors\n";
        exit(EXIT_FAILURE);
    }

    std::cout << "Input file: " << input_filename << "\n";
    std::vector<fasta_seq_t> seqs = parse_fasta(input_filename); // TODO use a struct for seqs to maintain id & desc

    // parse FASTA file
    std::cout << "Input file: " << input_filename << "\n";

    for (fasta_seq_t seq : seqs) {
        print_fasta_seq(seq);
    }

    // begin alignment
}
