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

int main(int argc, char *argv[]) {
    const auto start_time = std::chrono::steady_clock::now();

    // Initialize MPI
    int pid;
    int nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Parse cmd line args
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

    // P0 parses and serializes FASTA file
    std::vector<fasta_seq_t> fasta_seqs{};
    size_t num_bytes = -1;
    char *fasta_seqs_buf = NULL;
    if (pid == 0) {
        std::cout << "Input file: " << input_filename << "\n";
        fasta_seqs = parse_fasta(input_filename);
        fasta_seqs_buf = serialize_fasta_seqs(fasta_seqs, num_bytes);
    }

    // P0 broadcasts FASTA input to other procs
    MPI_Bcast(&num_bytes, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (pid > 0)
        fasta_seqs_buf = (char *) malloc(num_bytes);
    MPI_Bcast(fasta_seqs_buf, num_bytes, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (pid > 0)
        fasta_seqs = deserialize_fasta_seqs(fasta_seqs_buf, num_bytes);

    // Initialize program state
    align_params_t params{};

    int glbl_idx = 0; // Berger-Munson iteration number
    int seq_step = 0; // Sequential step count

    seq_group_t cur_alnmt = naiive_alnmt(fasta_seqs);
    int best_score = INT_MIN;
    int best_glbl_idx = -1;

    int num_seqs = fasta_seqs.size();
    int num_partns = num_seqs + (num_seqs * (num_seqs - 1)) / 2;

    int flag;
    std::string accept_reject_chain = "";

    // Register custom reduction op with MPI
    MPI_Op MPI_accept_op;
    MPI_Datatype MPI_pid_flag_t;
    MPI_Type_contiguous(2, MPI_INT, &MPI_pid_flag_t);
    MPI_Type_commit(&MPI_pid_flag_t);
    MPI_Op_create(accept_op, true, &MPI_accept_op);

    // Begin speculative computation
    while (glbl_idx - (best_glbl_idx + 1) < num_partns) {
        // Partition into two groups
        seq_group_t group1{};
        seq_group_t group2{};
        select_partn(cur_alnmt, glbl_idx + pid, group1, group2);

        remove_glbl_gaps(group1);
        remove_glbl_gaps(group2);

        // Compute alignment between two groups
        gap_pos_t gap_pos{};
        int cur_score = align_groups(group1, group2, params, gap_pos);
        if (cur_score > best_score)
            flag = ACCEPT;
        else
            flag = REJECT;

        // Check if any processors accepted, and take the one with lowest pid
        pid_flag_t send_pid_flag{};
        send_pid_flag.pid = pid;
        send_pid_flag.flag = flag;
        pid_flag_t recv_pid_flag{};
        MPI_Allreduce(&send_pid_flag, &recv_pid_flag, 1, MPI_pid_flag_t, MPI_accept_op, MPI_COMM_WORLD);

        if (recv_pid_flag.flag == ACCEPT) {
            int accepted_pid = recv_pid_flag.pid;

            // Broadcast data from accepted processor to others
            // index 0 --> size of group1 of partition (1 or 2)
            // index 1 --> first seq id of group1
            // index 2 --> second seq id of group1 (if there is one)
            // index 3 --> score of resulting alignment
            // index 4 --> length of resulting alignment
            int accepted_data[5];
            if (pid == accepted_pid) {
                accepted_data[0] = static_cast<int>(group1.size());
                accepted_data[1] = group1[0].id;
                if (group1.size() == 2)
                    accepted_data[2] = group1[1].id;
                else
                    accepted_data[2] = -1;
                accepted_data[3] = cur_score;
                accepted_data[4] = static_cast<int>(gap_pos.size());
            }
            MPI_Bcast(accepted_data, 5, MPI_INT, accepted_pid, MPI_COMM_WORLD);

            if (pid != accepted_pid) {
                // Reconstruct partition of accepted processor
                group1.clear();
                group2.clear();

                int group1_size = accepted_data[0];
                if (group1_size == 1) {
                    int group1_first = accepted_data[1];
                    group1.push_back(cur_alnmt[group1_first]);
                    for (int i = 0; i < num_seqs; i++) {
                        if (i != group1_first)
                            group2.push_back(cur_alnmt[i]);
                    }
                } else if (group1_size == 2) {
                    int group1_first = accepted_data[1];
                    int group1_second = accepted_data[2];
                    group1.push_back(cur_alnmt[group1_first]);
                    group1.push_back(cur_alnmt[group1_second]);
                    for (int i = 0; i < num_seqs; i++) {
                        if (i != group1_first && i != group1_second)
                            group2.push_back(cur_alnmt[i]);
                    }
                }

                remove_glbl_gaps(group1);
                remove_glbl_gaps(group2);
            }

            // Broadcast gap positions from accepted processor
            int gap_pos_len = accepted_data[4];
            char *gap_pos_bytes = (char *) malloc(gap_pos_len * 2);

            // Accepted processor serializes and broadcasts gap positions
            if (pid == accepted_pid) {
                for (int i = 0; i < gap_pos_len; i++) {
                    gap_pos_bytes[2*i] = 1 ? gap_pos[i].group1_gap : 0;
                    gap_pos_bytes[2*i+1] = 1 ? gap_pos[i].group2_gap : 0;
                }
            }
            MPI_Bcast(gap_pos_bytes, gap_pos_len * 2, MPI_CHAR, accepted_pid, MPI_COMM_WORLD);

            // Other processors deserialize gap positions
            if (pid != accepted_pid) {
                gap_pos.clear();
                for (int i = 0; i < gap_pos_len; i++) {
                    gap_option_t gap_opt;
                    gap_opt.group1_gap = true ? gap_pos_bytes[2*i] == 1 : false;
                    gap_opt.group2_gap = true ? gap_pos_bytes[2*i+1] == 1 : false;
                    gap_pos.push_back(gap_opt);
                }
            }

            // Update program state for next iteration
            int accepted_score = accepted_data[3];
            best_score = accepted_score;
            best_glbl_idx = glbl_idx + accepted_pid;
            cur_alnmt = update_alnmt(group1, group2, gap_pos);

            // Extend the accept-reject chain
            for (int i = 0; i < accepted_pid; i++)
                accept_reject_chain += 'R';
            accept_reject_chain += 'A';

            glbl_idx += accepted_pid + 1;
        } else if (recv_pid_flag.flag == REJECT) {
            // All processors have rejected
            for (int i = 0; i < nproc; i++)
                accept_reject_chain += 'R';
            glbl_idx += nproc;
        }

        seq_step++;
    }

    const auto end_time = std::chrono::steady_clock::now();
    const double runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    if (pid == 0) {
        std::cout << "Ran for " << glbl_idx << " iterations.\n";
        std::cout << "Took " << seq_step << " sequential steps.\n";
        std::cout << "Runtime (sec): " << runtime << "\n";
        std::cout << "Alignment score: " << best_score << "\n";
        std::cout << "Accepts and rejects: " << accept_reject_chain << "\n";

        std::ofstream fout(output_filename);

        fout << "Final alignment (score = " << best_score << "):\n";
        fout << "Accepts and rejects: " << accept_reject_chain << "\n";
        fout << "\n\n";
        for (seq_t seq : cur_alnmt) {
            fout << "seq " << std::setw(3) << seq.id << ": ";
            fout << seq.data << "\n";
        }
    }

    MPI_Op_free(&MPI_accept_op);
    MPI_Finalize();
}
