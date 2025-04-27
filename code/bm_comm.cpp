#include "bm_comm.h"

#include "align.h"
#include "parse_fasta.h"

#include <cstring>
#include <cassert>
#include <iostream>

#define FASTA_IDENT 0
#define FASTA_DESC 1
#define FASTA_SEQ 2

char *serialize_fasta_seqs(std::vector<fasta_seq_t>& fasta_seqs, size_t& num_bytes) {
    num_bytes = 0;
    for (fasta_seq_t fasta_seq : fasta_seqs) {
        num_bytes += fasta_seq.ident.size() + 1;
        num_bytes += fasta_seq.desc.size() + 1;
        num_bytes += fasta_seq.seq.size() + 1;
    }

    char *bytes = (char *) malloc(num_bytes);
    size_t pos = 0;

    for (fasta_seq_t fasta_seq : fasta_seqs) {
        memcpy(&bytes[pos], fasta_seq.ident.c_str(), fasta_seq.ident.size() + 1);
        pos += fasta_seq.ident.size() + 1;
        memcpy(&bytes[pos], fasta_seq.desc.c_str(), fasta_seq.desc.size() + 1);
        pos += fasta_seq.desc.size() + 1;
        memcpy(&bytes[pos], fasta_seq.seq.c_str(), fasta_seq.seq.size() + 1);
        pos += fasta_seq.seq.size() + 1;
    }
    assert(pos == num_bytes);

    return bytes;
}

std::vector<fasta_seq_t> deserialize_fasta_seqs(char *bytes, size_t num_bytes) {
    std::vector<fasta_seq_t> fasta_seqs{};

    int fasta_parse_state = FASTA_IDENT;
    static const fasta_seq_t empty_seq;
    fasta_seq_t curr_seq{};
    char *curr_str_start = bytes;
    for (size_t pos = 0; pos < num_bytes; pos++) {
        if (bytes[pos] == '\0') {
            if (fasta_parse_state == FASTA_IDENT) {
                curr_seq.ident = std::string(curr_str_start);
                fasta_parse_state = FASTA_DESC;
            } else if (fasta_parse_state == FASTA_DESC) {
                curr_seq.desc = std::string(curr_str_start);
                fasta_parse_state = FASTA_SEQ;
            } else if (fasta_parse_state == FASTA_SEQ) {
                curr_seq.seq = std::string(curr_str_start);
                fasta_seqs.push_back(curr_seq);
                curr_seq = empty_seq;
                fasta_parse_state = FASTA_IDENT;
            }

            curr_str_start = &bytes[pos+1];
        }
    }

    return fasta_seqs;
}
