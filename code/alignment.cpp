#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <ctime>
#include <sstream>

using namespace std;

class MultiGroupAlignment {
private:
    int gapOpen;
    int gapExtend;
    vector<vector<int>> matrix; // Substitution matrix
    
    // Dynamic programming matrices
    vector<int> HH; // Main score matrix
    vector<int> DD; // Gap in sequence 2
    vector<int> displ; // For alignment path
    
    // For backtracking and alignment
    int maxScore;
    int sb1, sb2; // Start positions for backtracking
    int se1, se2; // End positions (highest score)
    int printPtr, lastPrint;
    
    // Group information
    vector<int> group1; // Will contain exactly 2 sequences
    vector<int> group2; // Will contain all other sequences
    
    // Score matrix between all sequence pairs
    vector<vector<int>> scoreMatrix;
    
    // Best pair indices for alignment
    int bestSeq1Index;
    int bestSeq2Index;
    
    // Adds a diagonal match (0) or horizontal gap (positive value) to alignment path
    void add(int v) {
        if (lastPrint < 0) {
            displ[printPtr - 1] = v;
            displ[printPtr++] = lastPrint;
        }
        else {
            lastPrint = displ[printPtr++] = v;
        }
    }
    
    // Adds a vertical gap (negative value) to alignment path
    void del(int k) {
        if (lastPrint < 0) {
            lastPrint = displ[printPtr - 1] -= k;
        }
        else {
            lastPrint = displ[printPtr++] = -(k);
        }
    }
    
public:
    MultiGroupAlignment(int gapOpenPenalty, int gapExtendPenalty) 
        : gapOpen(gapOpenPenalty), 
          gapExtend(gapExtendPenalty),
          maxScore(0),
          sb1(0), sb2(0),
          se1(0), se2(0),
          printPtr(0), lastPrint(0),
          bestSeq1Index(-1),
          bestSeq2Index(-1) {
        // Initialize substitution matrix (simplified)
        matrix.resize(256, vector<int>(256, -5)); // Default mismatch penalty
        
        // Matching scores (this is simplified - normally would vary by amino acid)
        for (int i = 0; i < 256; i++) {
            matrix[i][i] = 10; // Match score
        }
    }
    
    // Divide sequences into two groups: one with two sequences, one with the rest
    void divideIntoGroups(const vector<vector<int>>& sequences) {
        // Clear previous groups
        group1.clear();
        group2.clear();
        
        int numSeq = sequences.size();
        if (numSeq < 3) {
            cout << "Need at least 3 sequences for this grouping approach." << endl;
            return;
        }
        
        // Create vector of indices
        vector<int> indices(numSeq);
        for (size_t i = 0; i < numSeq; i++) {
            indices[i] = i;
        }
        
        // Randomly shuffle indices
        random_device rd;
        mt19937 g(rd());
        shuffle(indices.begin(), indices.end(), g);
        
        // Put first two in group1, rest in group2
        group1.push_back(indices[0]);
        group1.push_back(indices[1]);
        
        for (size_t i = 2; i < numSeq; i++) {
            group2.push_back(indices[i]);
        }
        
        cout << "Group 1 contains sequences: ";
        for (int idx : group1) cout << idx << " ";
        cout << endl;
        
        cout << "Group 2 contains sequences: ";
        for (int idx : group2) cout << idx << " ";
        cout << endl;
    }
    
    // Forward pass to fill the scoring matrix
    void forwardPass(const vector<int>& seq1, const vector<int>& seq2) {
        int n = seq1.size() - 1; // Assuming 1-indexed sequences as in original code
        int m = seq2.size() - 1;
        
        // Resize matrices
        HH.resize(m + 1, 0);
        DD.resize(m + 1, 0);
        
        // Initialize first row
        for (int i = 0; i <= m; i++) {
            HH[i] = 0;
            DD[i] = -gapOpen;
        }
        
        maxScore = 0;
        se1 = se2 = 0;
        
        // Fill the matrix
        for (int i = 1; i <= n; i++) {
            int hh = 0;
            int p = 0;
            int f = -gapOpen; // Gap in sequence 1
            
            for (int j = 1; j <= m; j++) {
                // Update F (horizontal gap) - gap in sequence 1
                f -= gapExtend;
                int t = hh - gapOpen - gapExtend;
                if (f < t) {
                    f = t;
                }
                
                // Update DD (vertical gap) - gap in sequence 2
                DD[j] -= gapExtend;
                t = HH[j] - gapOpen - gapExtend;
                if (DD[j] < t) {
                    DD[j] = t;
                }
                
                // Calculate match/mismatch score
                hh = p + matrix[seq1[i]][seq2[j]];
                
                // Take maximum of match, horizontal gap, vertical gap
                if (hh < f) {
                    hh = f;
                }
                if (hh < DD[j]) {
                    hh = DD[j];
                }
                if (hh < 0) {
                    hh = 0; // Local alignment: scores can't be negative
                }
                
                // Save current cell for next iteration
                p = HH[j];
                HH[j] = hh;
                
                // Keep track of maximum score for local alignment
                if (hh > maxScore) {
                    maxScore = hh;
                    se1 = i;
                    se2 = j;
                }
            }
        }
    }
    
    // Calculate all pairwise scores between groups
    void calculateScoreMatrix(const vector<vector<int>>& sequences) {
        int numSeq = sequences.size();
        scoreMatrix.resize(numSeq, vector<int>(numSeq, 0));
        
        cout << "Calculating pairwise scores between groups..." << endl;
        
        // Calculate scores between each sequence in group1 and each sequence in group2
        for (int i : group1) {
            for (int j : group2) {
                // Run forward pass to calculate alignment score
                forwardPass(sequences[i], sequences[j]);
                
                // Store the maximum score
                scoreMatrix[i][j] = maxScore;
                scoreMatrix[j][i] = maxScore; // Matrix is symmetric
                
                cout << "Score between seq " << i << " (Group 1) and seq " << j 
                     << " (Group 2): " << maxScore << endl;
            }
        }
    }
    
    // Find best pair of sequences for alignment (one from each group)
    pair<int, int> findBestPair() {
        int bestScore = -1;
        pair<int, int> bestPair(-1, -1);
        
        // Find the best scoring pair across the two groups
        for (int i : group1) {
            for (int j : group2) {
                if (scoreMatrix[i][j] > bestScore) {
                    bestScore = scoreMatrix[i][j];
                    bestPair = make_pair(i, j);
                }
            }
        }
        
        return bestPair;
    }
    
    // Reverse pass (backtracking) to find starting positions
    void reversePass(const vector<int>& seq1, const vector<int>& seq2) {
        int cost = 0;
        sb1 = sb2 = 1; // Default starting positions
        
        vector<int> HH_rev(se2 + 1, -1);
        vector<int> DD_rev(se2 + 1, -1);
        
        for (int i = se1; i > 0; i--) {
            int hh = -1;
            int f = -1;
            int p = (i == se1) ? 0 : -1;
            
            for (int j = se2; j > 0; j--) {
                // Update F (gap in sequence 1)
                f -= gapExtend;
                int t = hh - gapOpen - gapExtend;
                if (f < t) {
                    f = t;
                }
                
                // Update DD (gap in sequence 2)
                DD_rev[j] -= gapExtend;
                t = HH_rev[j] - gapOpen - gapExtend;
                if (DD_rev[j] < t) {
                    DD_rev[j] = t;
                }
                
                // Calculate match score
                hh = p + matrix[seq1[i]][seq2[j]];
                
                // Take maximum
                if (hh < f) {
                    hh = f;
                }
                if (hh < DD_rev[j]) {
                    hh = DD_rev[j];
                }
                
                p = HH_rev[j];
                HH_rev[j] = hh;
                
                // Track best starting point
                if (hh > cost) {
                    cost = hh;
                    sb1 = i;
                    sb2 = j;
                    if (cost >= maxScore) {
                        break;
                    }
                }
            }
            
            if (cost >= maxScore) {
                break;
            }
        }
        
        cout << "Starting positions for alignment: (" << sb1 << ", " << sb2 << ")" << endl;
    }
    
    // Backtracking to generate alignment path
    void backtrack(const vector<int>& seq1, const vector<int>& seq2) {
        int n = se1 - sb1 + 1;
        int m = se2 - sb2 + 1;
        
        // Prepare for alignment path
        displ.resize(2 * (n + m) + 1);
        printPtr = 1;
        lastPrint = 0;
        
        // For a simplified version, we'll use a basic backtracking approach
        vector<vector<int>> score(n + 1, vector<int>(m + 1, 0));
        vector<vector<char>> direction(n + 1, vector<char>(m + 1, 0));
        
        // Reconstruct score matrix for backtracking
        for (int ii = 1; ii <= n; ii++) {
            for (int jj = 1; jj <= m; jj++) {
                // Calculate diagonal score (match/mismatch)
                int diag = score[ii-1][jj-1] + matrix[seq1[sb1+ii-1]][seq2[sb2+jj-1]];
                
                // Calculate horizontal score (gap in seq2)
                int horiz = score[ii][jj-1] - gapOpen - gapExtend;
                
                // Calculate vertical score (gap in seq1)
                int vert = score[ii-1][jj] - gapOpen - gapExtend;
                
                // Select best score
                if (diag >= horiz && diag >= vert) {
                    score[ii][jj] = diag;
                    direction[ii][jj] = 'd'; // Diagonal
                } else if (horiz >= vert) {
                    score[ii][jj] = horiz;
                    direction[ii][jj] = 'h'; // Horizontal
                } else {
                    score[ii][jj] = vert;
                    direction[ii][jj] = 'v'; // Vertical
                }
            }
        }
        
        // Backtrack to generate alignment
        int i = n;
        int j = m;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && direction[i][j] == 'd') {
                // Diagonal move (match/mismatch)
                add(0);
                i--;
                j--;
            } else if (j > 0 && (i == 0 || direction[i][j] == 'h')) {
                // Horizontal move (gap in seq1)
                add(1);
                j--;
            } else {
                // Vertical move (gap in seq2)
                del(1);
                i--;
            }
        }
    }
    
    // Trace path to generate alignment strings
    string tracePath(const vector<int>& seq1, const vector<int>& seq2) {
        int toDo = printPtr - 1;
        int i1 = sb1;
        int i2 = sb2;
        
        string aligned1, aligned2, middle;
        
        for (int i = 1; i <= toDo; ++i) {
            if (displ[i] == 0) {
                // Match or mismatch
                aligned1 += char(seq1[i1]);
                aligned2 += char(seq2[i2]);
                
                if (seq1[i1] == seq2[i2]) {
                    middle += '|'; // Match
                } else {
                    middle += ' '; // Mismatch
                }
                
                ++i1;
                ++i2;
            } else if (displ[i] > 0) {
                // Gap in sequence 1
                for (int k = 0; k < displ[i]; k++) {
                    aligned1 += '-';
                    aligned2 += char(seq2[i2]);
                    middle += ' ';
                    i2++;
                }
            } else {
                // Gap in sequence 2
                for (int k = 0; k < -displ[i]; k++) {
                    aligned1 += char(seq1[i1]);
                    aligned2 += '-';
                    middle += ' ';
                    i1++;
                }
            }
        }
        
        string result = "Alignment Score: " + to_string(maxScore) + "\n";
        result += "Sequence " + to_string(bestSeq1Index) + ": " + aligned1 + "\n";
        result += "                   " + middle + "\n";
        result += "Sequence " + to_string(bestSeq2Index) + ": " + aligned2 + "\n";
        
        return result;
    }
    
    // Main alignment function
    string alignSequences(const vector<vector<int>>& sequences) {
        if (sequences.size() < 3) {
            return "At least 3 sequences are required for this grouping approach.";
        }
        
        // Divide into two groups (2 sequences in group1, rest in group2)
        divideIntoGroups(sequences);
        
        // Calculate pairwise scores between groups
        calculateScoreMatrix(sequences);
        
        // Find best pair across groups
        pair<int, int> bestPair = findBestPair();
        
        if (bestPair.first == -1 || bestPair.second == -1) {
            return "Failed to find a good pair for alignment.";
        }
        
        bestSeq1Index = bestPair.first;
        bestSeq2Index = bestPair.second;
        
        cout << "Best pair for alignment: Sequence " << bestSeq1Index 
             << " and Sequence " << bestSeq2Index 
             << " (Score: " << scoreMatrix[bestSeq1Index][bestSeq2Index] << ")" << endl;
        
        // Forward pass on best pair to fill score matrix
        forwardPass(sequences[bestSeq1Index], sequences[bestSeq2Index]);
        
        // Reverse pass to find optimal starting positions
        reversePass(sequences[bestSeq1Index], sequences[bestSeq2Index]);
        
        // Backtrack to generate alignment path
        backtrack(sequences[bestSeq1Index], sequences[bestSeq2Index]);
        
        // Trace path to generate aligned sequences
        return tracePath(sequences[bestSeq1Index], sequences[bestSeq2Index]);
    }
    
    // Utility function to convert string to int vector (1-indexed like original code)
    static vector<int> stringToIntVector(const string& seq) {
        vector<int> result(seq.size() + 1);
        // Start from index 1 as in the original code
        for (size_t i = 0; i < seq.size(); i++) {
            result[i + 1] = seq[i];
        }
        return result;
    }
};

// Example usage
int main() {
    // Example parameters
    int gapOpen = 10;
    int gapExtend = 1;
    
    // Create aligner
    MultiGroupAlignment aligner(gapOpen, gapExtend);
    
    // Example sequences (4 or more sequences)
    vector<string> seq_strings = {
        "ACTGACGTAA",
        "ACGTACGTAC",
        "ACGTGACGTG",
        "ACTGACTGAC"
    };
    
    // Convert to int vectors (1-indexed)
    vector<vector<int>> sequences;
    for (const string& seq_str : seq_strings) {
        sequences.push_back(MultiGroupAlignment::stringToIntVector(seq_str));
    }
    
    // Perform alignment between groups
    string alignment_result = aligner.alignSequences(sequences);
    
    // Display result
    cout << "\nAlignment Result:" << endl;
    cout << alignment_result << endl;
    
    return 0;
}