#include <vector>
#include <string>
#include <stdio.h>


class multipleAlignmentObject{

    public:
        typedef struct{
            std::string seq;
        }seq;

        typedef struct{
            bool group1_gap;
            bool group2_gap;
        }gapPosition;

        int gapPenalty;
        // Value for difference when we calculate score
        int diff;
        std::vector<std::vector<int>> scoreMatrix;
        std::vector<std::vector<int>> backTrackMatrix;
        
        std::vector<seq> group1;
        std::vector<seq> group2;

        // From the backtracking matrix, when we do backtracking and we go horizontal we put gap on group1,
        // And when we go vertical put gap on group2, and then when we go diagonal we do not put gap.
        std::vector<gapPosition> gapPositions;

        int bestScore = 0;

        // Takes in group of sequences, divide it into two
        void partitionGroup(const std::vector<seq>& group){
            //Randomize grouping
            //Add one or two sequences to group1
            
            //Everything else in group2
            return;
        }

        // Fill in the score matrix and backtracking matrix

        void forwardPass(){
            // Initialize matrix 

            // Assume that group 1 and group 2 has same length sequences only in the group

            //S[i,j] = max(S[i,j-1]) + sizeof(g1(vertical)) * gap, S[i-1,j]) + sizeof(g2(horizontal)) * gap,
            //S[i-1,j-1] + scorebtwseqgroup(group1, group2, i, j))

            // We should update the backTrack matrix also, from the max value chosen. If we choose something 
            // from horizontal, store 0 in the backtrack matrix. Vertical store 1, diagonal store 2
            
            // We would probably want to fill DP matrix by filling horizontally first, and then adding up the
            // other parts.


             // Initialize matrices
            int rows = group1[0].seq.length() + 1;
            int cols = group2[0].seq.length() + 1;
            
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
            bestScore = scoreMatrix[rows-1][cols-1];

        }

        // Back tracking phase
        // Takes in the backTracking matrix, based on the backtracking Matrix update the gapPositions
        void reversePass(){
            int i = group1[0].seq.length();
        int j = group2[0].seq.length();
        
        // Clear previous gap positions
        gapPositions.clear();
        
        // Backtrack from bottom-right to top-left
        while (i > 0 || j > 0) {
            gapPosition gap;
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
            
            gapPositions.push_back(gap);
        }
        
        // Reverse the gap positions to get them in correct order (from left to right)
        std::reverse(gapPositions.begin(), gapPositions.end());
        }


        // Score calculation between two sequence groups. 
        int scoreBtwGroupSeq(const std::vector<seq>& group1, const std::vector<seq>& group2, int i, int j){
            int result = 0;

            
            //We have the group1 first element, as standard? 

            seq firstSeq = group1[0];
            char standard = firstSeq.seq[i];

            for(int k = 0; k < group1.size(); k++){
                for(int l = 0; l < group2.size(); l++){
                    if((group1[k].seq)[i] != (group2[l].seq)[j]){
                        result += diff;
                    }
                }
            }

            return result;
        }


        
        
        
        // Takes in group of sequences, gapPositions, and then based on that gap position and group number
        // Modify the orignal group
        // TODO: Let's figure
        void modify(std::vector<seq>& group, std::vector<gapPosition>, int groupNumber){

        }


        
        
    
};