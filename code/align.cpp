#include "align.h"

#include <algorithm>
#include <string>
#include <vector>

#define HORIZONTAL 0
#define VERTICAL 1
#define DIAGONAL 2

// Substituion score for one residue against another.
int sub_residue(char res1, char res2, align_params_t& params) {
    if (res1 == '-' && res2 == '-') {
        return 0;
    } else if (res1 == res2) {
        return params.match_reward;
    } else if (res1 == '-' || res2 == '-') {
        return params.gap_penalty;
    }  else {
        return params.sub_penalty;
    }
}

// Score of inserting a gap to the *other* group, against index i of group
int gap_score(int num_gaps, seq_group_t& group, int i, align_params_t& params) {
    int score = 0;

    // Calculate score within group
    for (size_t k = 0; k < group.size() - 1; k++) {
        for (size_t m = k + 1; m < group.size(); m++) {
            score += sub_residue(group[k].data[i], group[m].data[i], params);
        }
    }

    // Calculate score within gaps (should be all 0, gap-against-gap is 0)
    // for (size_t l = 0; l < num_gaps - 1; l++) {
    //    for (size_t m = l + 1; m < num_gaps, m++) {
    //        score += sub_residue('-', '-', params);
    //    }
    // }

    // Calcualte score between gaps and group
    for (size_t k = 0; k < group.size(); k++) {
        // for (size_t l = 0; l < num_gaps; l++) {
        //    score += sub_residue(group[k].data[i], '-');
        // }
        score += num_gaps * sub_residue(group[k].data[i], '-', params);
    }

    return score;
}

// Score calculation between two sequence groups.
int sub_score(seq_group_t& group1, seq_group_t& group2, int i, int j, align_params_t& params){
    int score = 0;

    // Calculate score within group 1
    for (size_t k = 0; k < group1.size() - 1; k++) {
        for (size_t m = k + 1; m < group1.size(); m++) {
            score += sub_residue(group1[k].data[i], group1[m].data[i], params);
        }
    }

    // Calculate score within group 2
    for (size_t l = 0; l < group2.size() - 1; l++) {
        for (size_t m = l + 1; m < group2.size(); m++) {
            score += sub_residue(group2[l].data[j], group2[m].data[j], params);
        }
    }

    // Calculate score between groups
    for(size_t k = 0; k < group1.size(); k++){
        for(size_t l = 0; l < group2.size(); l++){
            score += sub_residue(group1[k].data[i], group2[l].data[j], params);
        }
    }

    return score;
}

/*
 * group1 is represented along the vertical axis, and group2 is on the
 * horizontal axis.
 *
 * @pre group1 and group2 each have uniform length sequences
 * @return score of resulting alignment
 */
int forward_pass(seq_group_t& group1, seq_group_t& group2, align_params_t& params, matrix_t& score, matrix_t& backtrack){
    /**
     * S[i,j] = max {
     *   S[i,j-1] + group1.size() * gap,
     *   S[i-1,j] + group2.size() * gap,
     *   S[i-1,j-1] + sub'(group1, group2, i, j)
     * }
     */

    int num_rows = score.size();
    int num_cols = score[0].size();

    score[0][0] = 0;

    // Initialize first row and column with gap penalties
    for (int i = 1; i < num_rows; i++) {
        score[i][0] = score[i-1][0] + gap_score(group2.size(), group1, i-1, params);
        backtrack[i][0] = 1; // Vertical movement
    }

    for (int j = 1; j < num_cols; j++) {
        score[0][j] = score[0][j-1] + gap_score(group1.size(), group2, j-1, params);
        backtrack[0][j] = 0; // Horizontal movement
    }

    // Fill the matrices
    for (int i = 1; i < num_rows; i++) {
        for (int j = 1; j < num_cols; j++) {
            // Calculate scores for three possible moves
            int horizontal = score[i][j-1] + gap_score(group1.size(), group2, j-1, params);  // Gap in group1
            int vertical = score[i-1][j] + gap_score(group2.size(), group1, i-1, params);    // Gap in group2
            int diagonal = score[i-1][j-1] + sub_score(group1, group2, i-1, j-1, params);

            // Find the maximum score
            int maxScore = horizontal;
            int direction = HORIZONTAL; // 0 for horizontal, 1 for vertical, 2 for diagonal

            if (vertical > maxScore) {
                maxScore = vertical;
                direction = VERTICAL;
            }

            if (diagonal > maxScore) {
                maxScore = diagonal;
                direction = DIAGONAL;
            }

            // Update matrices
            score[i][j] = maxScore;
            backtrack[i][j] = direction;
        }
    }

    // Store best score
    int alnmt_score = score[num_rows-1][num_cols-1];
    return alnmt_score;
}

// Backtracking pass. Builds gap positions based on result of forward pass
void backward_pass(matrix_t& backtrack, gap_pos_t& gap_pos){
    int i = backtrack.size() - 1;
    int j = backtrack[0].size() - 1;

    // Clear previous gap positions
    gap_pos.clear();

    // Backtrack from bottom-right to top-left
    while (i > 0 || j > 0) {
        gap_option_t gap;
        gap.group1_gap = false;
        gap.group2_gap = false;

        if (i > 0 && j > 0 && backtrack[i][j] == 2) {
            // Diagonal move - no gaps
            i--;
            j--;
        } else if (j > 0 && (i == 0 || backtrack[i][j] == 0)) {
            // Horizontal move - gap in group1
            gap.group1_gap = true;
            j--;
        } else if (i > 0 && (j == 0 || backtrack[i][j] == 1)) {
            // Vertical move - gap in group2
            gap.group2_gap = true;
            i--;
        }

        gap_pos.push_back(gap);
    }

    // Reverse the gap positions to get them in correct order (from left to right)
    std::reverse(gap_pos.begin(), gap_pos.end());
}

// Implements align_groups, described in align.h
int align_groups(seq_group_t& group1, seq_group_t& group2, align_params_t& params, gap_pos_t& gap_pos) {
    // Initialize matrices
    int num_rows = group1[0].data.length() + 1;
    int num_cols = group2[0].data.length() + 1;

    matrix_t score{};
    matrix_t backtrack{};
    score.resize(num_rows, std::vector<int>(num_cols, 0));
    backtrack.resize(num_rows, std::vector<int>(num_cols, 0));

    int alnmt_score = forward_pass(group1, group2, params, score, backtrack);
    backward_pass(backtrack, gap_pos);

    return alnmt_score;
}

// Implements update_alnmt, described in align.h
seq_group_t update_alnmt(seq_group_t& group1, seq_group_t& group2, gap_pos_t& gap_pos) {
    int num_seqs = group1.size() + group2.size();

    // Initialize new alignment
    seq_group_t new_alnmt{};
    for (int i = 0; i < num_seqs; i++) {
        seq_t seq_i;
        seq_i.id = i;
        seq_i.data = "";
        new_alnmt.push_back(seq_i);
    }

    // Build new alignment
    int group1_pos = 0;
    int group2_pos = 0;
    for (size_t i = 0; i < gap_pos.size(); i++) {
        bool group1_gap = gap_pos[i].group1_gap;
        bool group2_gap = gap_pos[i].group2_gap;

        if (group1_gap) {
            for (seq_t& group1_seq : group1)
                new_alnmt[group1_seq.id].data += '-';
        } else {
            for (seq_t& group1_seq : group1) {
                char residue = group1_seq.data[group1_pos];
                new_alnmt[group1_seq.id].data += residue;
            }
            group1_pos++;
        }

        if (group2_gap) {
            for (seq_t& group2_seq : group2)
                new_alnmt[group2_seq.id].data += '-';
        } else {
            for (seq_t& group2_seq : group2) {
                char residue = group2_seq.data[group2_pos];
                new_alnmt[group2_seq.id].data += residue;
            }
            group2_pos++;
        }
    }

    return new_alnmt;
}
