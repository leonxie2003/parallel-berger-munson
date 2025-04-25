#include "align.h"

#include <algorithm>
#include <string>
#include <vector>

#define HORIZONTAL 0
#define VERTICAL 1
#define DIAGONAL 2


// Score calculation between two sequence groups.
int score_group_subst(const seq_group_t& group1, const seq_group_t& group2, int i, int j, align_params_t params){
    int result = 0;

    for(int k = 0; k < group1.size(); k++){
        for(int l = 0; l < group2.size(); l++){
            if(group1[k][i] != group2[l][j]){
                result += params.sub_penalty;
            }
        }
    }

    return result;
}

/**
 * group1 is represented along the vertical axis, and group2 is on the
 * horizontal axis.
 *
 * @pre group1 and group2 each have uniform length sequences
 * @return score of resulting alignment
 */
int forward_pass(seq_group_t& group1, seq_group_t& group2, align_params_t params, matrix_t& score, matrix_t& backtrack){
    /**
     * S[i,j] = max {
     *   S[i,j-1] + group1.size() * gap,
     *   S[i-1,j] + group2.size() * gap,
     *   S[i-1,j-1] + sub'(group1, group2, i, j)
     * }
     */

    int num_rows = score.size();
    int num_cols = score[0].size();

    // Initialize first row and column with gap penalties
    for (int i = 1; i < num_rows; i++) {
        score[i][0] = i * params.gap_penalty * group2.size();
        backtrack[i][0] = 1; // Vertical movement
    }

    for (int j = 1; j < num_cols; j++) {
        score[0][j] = j * params.gap_penalty * group1.size();
        backtrack[0][j] = 0; // Horizontal movement
    }

    // Fill the matrices
    for (int i = 1; i < num_rows; i++) {
        for (int j = 1; j < num_cols; j++) {
            // Calculate scores for three possible moves
            int horizontal = score[i][j-1] + params.gap_penalty * group1.size();  // Gap in group1
            int vertical = score[i-1][j] + params.gap_penalty * group2.size();    // Gap in group2
            int diagonal = score[i-1][j-1] + score_group_subst(group1, group2, i-1, j-1, params);

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

// Back tracking phase
// Takes in the backTracking matrix, based on the backtracking Matrix update the gapPositions
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



int align_groups(seq_group_t& group1, seq_group_t& group2, align_params_t params, gap_pos_t& gap_pos) {
    // Initialize matrices
    int num_rows = group1[0].length() + 1;
    int num_cols = group2[0].length() + 1;

    matrix_t score{};
    matrix_t backtrack{};
    score.resize(num_rows, std::vector<int>(num_cols, 0));
    backtrack.resize(num_rows, std::vector<int>(num_cols, 0));

    forward_pass(group1, group2, params, score, backtrack);
    backward_pass(backtrack, gap_pos);

    return 0;
}

void update_alnmt(seq_group_t& alnmt, gap_pos_t& gap_pos) {
    // TODO implement
    return;
}
