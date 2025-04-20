/**
 * Parse FASTA files.
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include <string>
#include <vector>

typedef struct fasta_seq {
    std::string id;
    std::string desc;
    std::string seq;
} fasta_seq_t;

std::vector<fasta_seq_t> parse_fasta(std::string filename);

void print_fasta_seq(fasta_seq_t seq);
