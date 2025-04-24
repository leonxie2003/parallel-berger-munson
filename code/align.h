/** @file align.h
 *  @brief Alignment code for Berger-Munson.
 *  @author Leon Xie (leonx)
 *  @author Taekseung Kim (taekseuk)
 */

#ifndef __ALIGN_H__
#define __ALIGN_H__

#include <string>
#include <vector>

typedef std::string seq_t;
typedef std::vector<seq_t> seq_group_t;



/**
 * Paramters that are used for alignment process
 * gapPenalty: Penalty for using a gap. 
 * diff: 
 */
int gapPenalty;
int diff;
std::vector<std::vector<int>> scoreMatrix;
std::vector<std::vector<int>> backTrackMatrix;


/**
 * Data structure to represent gap positions. After alignment, each position
 * in the new alignment between groups has a gap in at most one of the two
 * groups.
 */
typedef struct gap_option {
    bool group1_gap;
    bool group2_gap;
} gap_option_t;
typedef std::vector<gap_option_t> gap_pos_t;

/**
 * Aligns two sequence groups, saving the new gaps.
 *
 *
 * @param group1
 * @param group2
 * @param gap_pos
 * @return Score of the resulting alignment.
 */
int align_groups(seq_group_t group1, seq_group_t group2, gap_pos_t& gap_pos);

void update_alnmt(seq_group_t& alnmt, gap_pos_t& gap_pos);

#endif
