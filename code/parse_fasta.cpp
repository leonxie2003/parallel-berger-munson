/**
 * Parse FASTA files.
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>

std::vector<std::string> parse_fasta(std::string filename) {
    std::vector<std::string> seqs{};

    std::ifstream fin(filename);

    if (!fin) {
        std::cerr << "Unable to open file: " << filename << ".\n";
        exit(EXIT_FAILURE);
    }

    std::string input_line{};
    std::string cur_seq;
    bool seen_seq = false;

    while (std::getline(fin, input_line)) {
        if (input_line[0] == '>') {
            if (seen_seq) {
                seqs.push_back(cur_seq);
            } else {
                seen_seq = true;
            }

            cur_seq = "";
        } else {
            cur_seq.append(input_line);
        }
    }
    seqs.push_back(cur_seq); // flush last seq

    return seqs;
}
