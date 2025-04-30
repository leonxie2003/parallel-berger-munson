/** @file bm_utils.h
 *  Provides utility functions for Berger-Munson.
 *  @author Leon Xie (leonx)
 *  @author Taekseung Kim (taekseuk)
 */

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

/**
 * Constructs a naiive alignment by adding gaps to the end.
 */
seq_group_t naiive_alnmt(std::vector<fasta_seq_t> fasta_seqs);

/**
 * Randomly (or pseudorandomly) selects a partition of the sequence.
 *
 * Constructs the partition into group1 and group2.
 */
void select_partn(seq_group_t seqs, int glbl_idx, int random_mode, seq_group_t& group1, seq_group_t& group2);

/**
 * Removes global gaps (gaps that exist in every sequence of a group).
 */
void remove_glbl_gaps(seq_group_t& group);

#endif
