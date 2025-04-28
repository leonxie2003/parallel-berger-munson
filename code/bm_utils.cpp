#include "bm_utils.h"
#include "align.h"
#include "parse_fasta.h"

#include <string>
#include <vector>
#include <random>

seq_group_t naiive_alnmt(std::vector<fasta_seq_t> fasta_seqs) {
    seq_group_t naiive_alnmt{};
    size_t longest = 0;
    for (size_t i = 0; i < fasta_seqs.size(); i++) {
        seq_t seq;
        seq.id = i;
        seq.data = fasta_seqs[i].seq;
        naiive_alnmt.push_back(seq);

        if (fasta_seqs[i].seq.size() > longest) {
            longest = fasta_seqs[i].seq.size();
        }

    }

    for (seq_t& seq : naiive_alnmt) {
        size_t seq_len = seq.data.size();
        for (size_t i = 0; i < longest - seq_len; i++)
            seq.data.append("-");
    }

    return naiive_alnmt;
}

void select_partn(seq_group_t seqs, int glbl_idx, int random_mode, seq_group_t& group1, seq_group_t& group2) {
    int num_seqs = seqs.size();

    std::mt19937 gen{};
    if (random_mode == DEVICERANDOM) {
        std::random_device rd;
        gen.seed(rd());
    } else if (random_mode == PSEUDORANDOM) {
        gen.seed(glbl_idx);
    }

    int num_partns = num_seqs + (num_seqs * (num_seqs - 1)) / 2;
    std::uniform_int_distribution<> distr(0, num_partns - 1);

    int partn_num = distr(gen);

    if (partn_num < num_seqs) {
        // group1 has 1 sequence
        group1.push_back(seqs[partn_num]);
        for (int i = 0; i < num_seqs; i++) {
            if (i != partn_num)
                group2.push_back(seqs[i]);
        }
        return;
    } else {
        // group1 has 2 sequences
        partn_num -= num_seqs;
        int group1_first = 0;
        while (partn_num >= num_seqs - (group1_first + 1)) {
            partn_num -= num_seqs - (group1_first + 1);
            group1_first++;
        }
        int group1_second = group1_first + 1 + partn_num;

        group1.push_back(seqs[group1_first]);
        group1.push_back(seqs[group1_second]);
        for (int i = 0; i < num_seqs; i++)
            if (i != group1_first && i != group1_second)
                group2.push_back(seqs[i]);
    }
}

void remove_glbl_gaps(seq_group_t& group) {
    size_t seq_len = group[0].data.size();

    for (int i = seq_len - 1; i >= 0; i--) {

        bool all_gap = true;
        for (seq_t& seq : group) {
            if (seq.data[i] != '-') {
                all_gap = false;
                break;
            }
        }

        if (all_gap) {
            for (seq_t& seq : group)
                seq.data.erase(i, 1);
        }
    }

    return;
}
