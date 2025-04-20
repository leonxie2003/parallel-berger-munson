/**
 * Parse FASTA files.
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include <parse_fasta.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>

std::vector<fasta_seq_t> parse_fasta(std::string filename) {
    std::vector<fasta_seq_t> seqs{};

    std::ifstream fin(filename);

    if (!fin) {
        std::cerr << "Unable to open file: " << filename << ".\n";
        exit(EXIT_FAILURE);
    }

    std::string input_line{};
    fasta_seq_t cur_seq;
    bool seen_first_seq = false;

    while (std::getline(fin, input_line)) {
        if (input_line[0] == '>') {
            if (seen_first_seq) {
                seqs.push_back(cur_seq);
            } else {
                seen_first_seq = true;
            }

            size_t space_pos = input_line.find(" ");
            if (space_pos != std::string::npos) {
                cur_seq.id = input_line.substr(1, space_pos - 1);
                cur_seq.desc = input_line.substr(space_pos);
            } else {
                cur_seq.id = input_line.substr(1);
                cur_seq.desc = "";
            }
            cur_seq.seq = "";
        } else {
            cur_seq.seq.append(input_line);
        }
    }
    seqs.push_back(cur_seq); // flush last seq

    return seqs;
}

void print_fasta_seq(fasta_seq_t seq) {
    std::cout << "FASTA sequence ID: " << seq.id << "\n";
    std::cout << "Description: " << seq.desc << "\n";
    std::cout << "Sequence: " << seq.seq << "\n";
    std::cout << "\n";
}
