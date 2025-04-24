#include "bm_utils.h"
#include "align.h"
#include "parse_fasta.h"

#include <string>
#include <vector>

seq_group_t naiive_alnmt(std::vector<fasta_seq_t> fasta_seqs) {
    seq_group_t naiive_alnmt{};
    size_t longest = 0;
    for (fasta_seq_t fasta_seq : fasta_seqs) {
        naiive_alnmt.push_back(fasta_seq.seq);

        if (fasta_seq.seq.size() > longest) {
            longest = fasta_seq.seq.size();
        }

    }

    for (seq_t& seq : naiive_alnmt) {
        size_t seq_size = seq.size();
        for (size_t i = 0; i < longest - seq_size; i++)
            seq.append("-");
    }

    return naiive_alnmt;
}

void select_partn(seq_group_t seqs, int glbl_idx, seq_group_t& group1, seq_group_t& group2) {
    // TODO actually implement
    group1.push_back(seqs[0]);
    group1.push_back(seqs[1]);
    for (size_t i = 2; i < seqs.size(); i++) {
        group2.push_back(seqs[i]);
    }
}

void remove_glbl_gaps(seq_group_t& group) {
    // TODO implement
    return;
}
