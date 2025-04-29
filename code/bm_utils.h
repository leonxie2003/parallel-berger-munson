#ifndef __BM_UTILS_H__
#define __BM_UTILS_H__

#include "align.h"
#include "parse_fasta.h"

#include <vector>

#define ACCEPT 1
#define REJECT 2

#define DEVICERANDOM 1
#define PSEUDORANDOM 2

#define CLOCK_NOW (std::chrono::steady_clock::now())
#define TIME_SEC(START, END) (std::chrono::duration_cast<std::chrono::duration<double>>((END) - (START)).count())

seq_group_t naiive_alnmt(std::vector<fasta_seq_t> fasta_seqs);

void select_partn(seq_group_t seqs, int glbl_idx, int random_mode, seq_group_t& group1, seq_group_t& group2);

void remove_glbl_gaps(seq_group_t& group);

#endif
