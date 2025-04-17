/**
 * Sequential Berger-Munson program
 *
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include <iostream>
#include <string>

#include <unistd.h>

int main(int argc, char *argv[]) {
    /* --- parse cmd line args --- */
    std::string in_filename;
    std::string out_filename;

    int opt;
    while((opt = getopt(argc, argv, "i:o:")) != -1) {
        switch (opt) {
            case 'i':
                in_filename = optarg;
                break;
            case 'o':
                out_filename = optarg;
                break;
        default:
            std::cerr << "Usage: " << argv[0] << " -i input_filename -o output_filename\n";
            exit(EXIT_FAILURE);
        }
    }

    if (empty(in_filename) || empty(out_filename)) {
        std::cerr << "Usage: " << argv[0] << " -i input_filename -o output_filename\n";
        exit(EXIT_FAILURE);
    }

    // parse FASTA file
    std::cout << "Input file: " << in_filename << "\n";

    // begin alignment
}
