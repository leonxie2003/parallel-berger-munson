/**
 * Sequential Berger-Munson program
 *
 * Leon Xie (leonx), Taekseung Kim (taekseuk)
 */

#include "bm_seq.h"
#include "parse_fasta.h"
#include "align.h"
#include "bm_utils.h"
#include "bm_comm.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <string>
#include <vector>

#include <unistd.h>
#include <mpi.h>

#define ACCEPT true
#define REJECT false

int main(int argc, char *argv[]) {
    const auto start_time = std::chrono::steady_clock::now();
    int pid;
    int nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    /* --- parse cmd line args --- */

    std::string input_filename;
    std::string output_filename;

    int opt;
    while((opt = getopt(argc, argv, "i:o:")) != -1) {
        switch (opt) {
            case 'i':
                input_filename = optarg;
                break;
            case 'o':
                output_filename = optarg;
                break;
        default:
            std::cerr << "Usage: " << argv[0] << " -i input_filename -o output_filename\n";
            exit(EXIT_FAILURE);
        }
    }

    if (empty(input_filename) || empty(output_filename)) {
        std::cerr << "Usage: " << argv[0] << " -i input_filename -o output_filename\n";
        exit(EXIT_FAILURE);
    }

    /* --- parse FASTA file --- */
    std::vector<fasta_seq_t> fasta_seqs{};
    size_t num_bytes;
    char *fasta_seqs_buf;
    if (pid == 0) {
        std::cout << "Input file: " << input_filename << "\n";
        fasta_seqs = parse_fasta(input_filename); // TODO use a struct for seqs to maintain id & desc
        fasta_seqs_buf = serialize_fasta_seqs(fasta_seqs, num_bytes);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&num_bytes, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (pid > 0) {
        fasta_seqs_buf = (char *) malloc(num_bytes);
    }

    MPI_Bcast(fasta_seqs_buf, num_bytes, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (pid > 0) {
        fasta_seqs = deserialize_fasta_seqs(fasta_seqs_buf, num_bytes);
    }

    seq_group_t cur_alnmt = naiive_alnmt(fasta_seqs);
    /* TODO remove debug prints
    for (fasta_seq_t seq : seqs) {
        print_fasta_seq(seq);
    }
    */


    // iterate sequentially until first reject
    align_params_t params{};
    params.match_reward = 1;
    params.gap_penalty = -1;
    params.sub_penalty = 0;

    int glbl_idx = 0;
    int best_score = INT_MIN;
    int best_glbl_idx = -1;
    int num_seqs = fasta_seqs.size();
    int num_partns = num_seqs + (num_seqs * (num_seqs - 1)) / 2;

    std::string accept_reject_chain = "";
    bool flag = ACCEPT;
    do {
        // partition into two group
        seq_group_t group1{};
        seq_group_t group2{};
        select_partn(cur_alnmt, glbl_idx, group1, group2);

        remove_glbl_gaps(group1);
        remove_glbl_gaps(group2);

        gap_pos_t gap_pos{};
        int cur_score = align_groups(group1, group2, params, gap_pos);

        if (cur_score > best_score) {
            best_score = cur_score;
            best_glbl_idx = glbl_idx;
            cur_alnmt = update_alnmt(group1, group2, gap_pos);
            accept_reject_chain += 'A';
            flag = ACCEPT;
        } else {
            accept_reject_chain = 'R';
            flag = REJECT;
        }
    } while (flag == ACCEPT);

    // begin speculation

    const auto end_time = std::chrono::steady_clock::now();
    const double runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    if (pid == 0) {
        std::cout << "Ran for " << glbl_idx << " iterations.\n";
        std::cout << "Runtime (sec): " << runtime << "\n";
        std::cout << "Accepts and rejects: " << accept_reject_chain << "\n";

        std::ofstream fout(output_filename);

        fout << "Final alignment (score = " << best_score << "):\n";
        for (seq_t seq : cur_alnmt) {
            fout << "seq " << std::setw(3) << seq.id << ": ";
            fout << seq.data << "\n";
        }
    }

    MPI_Finalize();
}
