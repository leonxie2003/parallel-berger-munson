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
    char *fasta_seqs_buf = NULL;
    if (pid == 0) {
        std::cout << "Input file: " << input_filename << "\n";
        fasta_seqs = parse_fasta(input_filename); // TODO use a struct for seqs to maintain id & desc
        fasta_seqs_buf = serialize_fasta_seqs(fasta_seqs, num_bytes);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&num_bytes, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (pid > 0)
        fasta_seqs_buf = (char *) malloc(num_bytes);

    MPI_Bcast(fasta_seqs_buf, num_bytes, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (pid > 0)
        fasta_seqs = deserialize_fasta_seqs(fasta_seqs_buf, num_bytes);

    seq_group_t cur_alnmt = naiive_alnmt(fasta_seqs);

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
    int flag = ACCEPT;
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

        glbl_idx++;
    } while (flag == ACCEPT);

    // begin speculation
    MPI_Op MPI_accept_op;
    MPI_Datatype MPI_pid_flag_t;
    MPI_Type_contiguous(2, MPI_INT, &MPI_pid_flag_t);
    MPI_Type_commit(&MPI_pid_flag_t);

    MPI_Op_create(accept_op, true, &MPI_accept_op);

    while (glbl_idx - (best_glbl_idx + 1) < num_partns) {
        seq_group_t group1{};
        seq_group_t group2{};
        select_partn(cur_alnmt, glbl_idx + pid, group1, group2);

        remove_glbl_gaps(group1);
        remove_glbl_gaps(group2);

        gap_pos_t gap_pos{};
        int cur_score = align_groups(group1, group2, params, gap_pos);
        if (cur_score > best_score) {
            best_score = cur_score;
            flag = ACCEPT;
        }
        else {
            flag = REJECT;
        }

        // Reduce with other procs
        pid_flag_t send_pid_flag{};
        send_pid_flag.pid = pid;
        send_pid_flag.flag = flag;
        pid_flag_t recv_pid_flag{};

        // Implicit barrier here?
        MPI_Allreduce(&send_pid_flag, &recv_pid_flag, 1, MPI_pid_flag_t, MPI_accept_op, MPI_COMM_WORLD);

        if (recv_pid_flag.flag == ACCEPT) {
            int accepted_pid = recv_pid_flag.pid;
            // copy from accepted pid
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
                    int group1_first = accepted_data[2];
                    int group1_second = accepted_data[3];
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

            // Broadcast gap_pos
            if (pid == accepted_pid) {

            }

            //Bcast

            if (pid != accepted_pid) {

            }

            // TODO update best_score, update best_glbl_idx, update cur_alnmt
            // TODO update accept-reject chain...

            glbl_idx += accepted_pid + 1;
        } else if (recv_pid_flag.flag == REJECT) {
            glbl_idx += nproc;
        }
    }

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

    MPI_Op_free(&MPI_accept_op);
    MPI_Finalize();
}
