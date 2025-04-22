#include <vector>
#include <string>
#include <stdio.h>

#include 

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
        std::vector<std::vector<int>> scoreMatrix;
        std::vector<std::vector<int>> backTrackMatrix;
        
        std::vector<seq> group1;
        std::vector<seq> group2;

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

            // matrix size would be the maximum length of each group(Q. Should we start from sequences
            // that have same length?)

            //S[i,j] = max(S[i,j-1]) + sizeof(g1(vertical)) * gap, S[i,j-1]) + sizeof(g2(horizontal)) * gap,
            //S[i-1,j-1] + scorebtwseqgroup(group1, group2, i, j))
            throw("forwardPass NOT IMPLEMENTED");
        }

        // Back tracking phase
        // Takes in the backTracking matrix, based on the backtracking Matrix update the gapPositions
        void reversePass(){
            throw("reversePass NOT IMPLEMENTED");
        }

        // Score calculation between two sequences. Helper function for scorebtwgroupseq
        int scoreBtwSeq(const seq& seq1, const seq& seq2){
            throw("scoreBtwSeq NOT IMPLEMENTED");
        }

        // Score calculation between two sequence groups. 
        int scoreBtwGroupSeq(const std::vector<seq>& group1, const std::vector<seq>& group2){
            throw("scoreBtwGroupSeq NOT IMPLEMENTED");
        }



        
        
        // Takes in group of sequences, gapPositions, and then based on that gap position and group number
        // Modify the orignal group
        // TODO: Let's figure
        void modify(std::vector<seq>& group, std::vector<gapPosition>, int groupNumber){

        }


        
        
    
};