/** @file align.h
 *  @brief Alignment code for Berger-Munson.
 *  @author Leon Xie (leonx)
 *  @author Taekseung Kim (taekseuk)
 */

#ifndef __ALIGN_H__
#define __ALIGN_H__

#include <string>
#include <vector>

/**
 * Data structure to represent sequences.
 */
typedef struct seq {
    int id;
    std::string data;
} seq_t;
typedef std::vector<seq_t> seq_group_t;

/**
 * Represents 2D matrices of integers.
 */
typedef std::vector<std::vector<int>> matrix_t;

/**
 * Represents aligment parameters
 */
typedef struct align_params {
    int match_reward = 1;
    int gap_penalty = -1;
    int sub_penalty = 0; // Replace with subst matrix for better aligments.
} align_params_t;

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
 * Aligns two sequence groups, saving the new gaps in gap_pos.
 *
 * @param group1
 * @param group2
 * @param gap_pos
 * @return Score of the resulting alignment.
 */
int align_groups(seq_group_t& group1, seq_group_t& group2, align_params_t& params, gap_pos_t& gap_pos);

/**
 * Updates an alignment with new gap positons. group1 and group2 represent a partition of the alignment.
 * @param group1
 * @param group2
 * @param gap_pos
 */
seq_group_t update_alnmt(seq_group_t& group1, seq_group_t& group2, gap_pos_t& gap_pos);

#endif
