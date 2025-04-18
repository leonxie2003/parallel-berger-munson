/**
 * Sequential Berger-Munson program
 *
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include "bm_seq.h"
#include "parse_fasta.h"

#include <iostream>
#include <vector>
#include <string>

#include <unistd.h>

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

    // parse FASTA file
    std::cout << "Input file: " << input_filename << "\n";
    std::vector<std::string> seqs = parse_fasta(input_filename);

    for (std::string seq : seqs) {
        std::cout << seq << "\n";
    }

    // begin alignment
}
