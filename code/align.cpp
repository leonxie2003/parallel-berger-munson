#include "align.h"



void forwardPass(seq_group_t group1, seq_group_t group2, gap_pos_t& gap_pos){
    // Initialize matrix 

    // Assume that group 1 and group 2 has same length sequences only in the group

    //S[i,j] = max(S[i,j-1]) + sizeof(g1(vertical)) * gap, S[i-1,j]) + sizeof(g2(horizontal)) * gap,
    //S[i-1,j-1] + scorebtwseqgroup(group1, group2, i, j))

    // We should update the backTrack matrix also, from the max value chosen. If we choose something 
    // from horizontal, store 0 in the backtrack matrix. Vertical store 1, diagonal store 2
    
    // We would probably want to fill DP matrix by filling horizontally first, and then adding up the
    // other parts.


     // Initialize matrices
    int rows = group1[0].length() + 1;
    int cols = group2[0].length() + 1;
    
    // Resize matrices
    scoreMatrix.resize(rows, std::vector<int>(cols, 0));
    backTrackMatrix.resize(rows, std::vector<int>(cols, 0));
    
    // Initialize first row and column with gap penalties
    for (int i = 1; i < rows; i++) {
        scoreMatrix[i][0] = i * gapPenalty * group2.size();
        backTrackMatrix[i][0] = 1; // Vertical movement
    }
    
    for (int j = 1; j < cols; j++) {
        scoreMatrix[0][j] = j * gapPenalty * group1.size();
        backTrackMatrix[0][j] = 0; // Horizontal movement
    }
    
    // Fill the matrices
    for (int i = 1; i < rows; i++) {
        for (int j = 1; j < cols; j++) {
            // Calculate scores for three possible moves
            int horizontal = scoreMatrix[i][j-1] + gapPenalty * group1.size();  // Gap in group1
            int vertical = scoreMatrix[i-1][j] + gapPenalty * group2.size();    // Gap in group2
            int diagonal = scoreMatrix[i-1][j-1] + scoreBtwGroupSeq(group1, group2, i-1, j-1);
            
            // Find the maximum score
            int maxScore = horizontal;
            int direction = 0; // 0 for horizontal, 1 for vertical, 2 for diagonal
            
            if (vertical > maxScore) {
                maxScore = vertical;
                direction = 1;
            }
            
            if (diagonal > maxScore) {
                maxScore = diagonal;
                direction = 2;
            }
            
            // Update matrices
            scoreMatrix[i][j] = maxScore;
            backTrackMatrix[i][j] = direction;
        }
    }
    
    // Store best score
    int bestScore = scoreMatrix[rows-1][cols-1];

}

// Back tracking phase
// Takes in the backTracking matrix, based on the backtracking Matrix update the gapPositions
void reversePass(seq_group_t group1, seq_group_t group2, gap_pos_t& gap_pos){
    int i = group1[0].length();
    int j = group2[0].length();

    // Clear previous gap positions
    gap_pos.clear();

    // Backtrack from bottom-right to top-left
    while (i > 0 || j > 0) {
        gap_option_t gap;
        gap.group1_gap = false;
        gap.group2_gap = false;
        
        if (i > 0 && j > 0 && backTrackMatrix[i][j] == 2) {
            // Diagonal move - no gaps
            i--;
            j--;
        } else if (j > 0 && (i == 0 || backTrackMatrix[i][j] == 0)) {
            // Horizontal move - gap in group1
            gap.group1_gap = true;
            j--;
        } else if (i > 0 && (j == 0 || backTrackMatrix[i][j] == 1)) {
            // Vertical move - gap in group2
            gap.group2_gap = true;
            i--;
        }
    
    gap_pos.push_back(gap);
}

// Reverse the gap positions to get them in correct order (from left to right)
std::reverse(gap_pos.begin(), gap_pos.end());
}


// Score calculation between two sequence groups. 
int scoreBtwGroupSeq(const seq_group_t& group1, const seq_group_t& group2, int i, int j){
    int result = 0;

    
    //We have the group1 first element, as standard? 

    for(int k = 0; k < group1.size(); k++){
        for(int l = 0; l < group2.size(); l++){
            if(group1[k][i] != group2[l][j]){
                result += diff;
            }
        }
    }

    return result;
}


int align_groups(seq_group_t group1, seq_group_t group2, gap_pos_t& gap_pos) {
    // TODO implement
    return 0;
}

void update_alnmt(seq_group_t& alnmt, gap_pos_t& gap_pos) {
    // TODO implement
    return;
}
