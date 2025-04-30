/** @file bm_comm.h
 *  Communication subroutines for parallel Berger-Munson.
 *  @author Leon Xie (leonx)
 *  @author Taekseung Kim (taekseuk)
 */

#ifndef __BM_COMM_H__
#define __BM_COMM_H__

#include "align.h"
#include "parse_fasta.h"

#include <vector>
#include <mpi.h>

/**
 * Represents a process ID, and a accept-reject flag.
 */
typedef struct pid_flag {
    int pid;
    int flag;
} pid_flag_t;

/**
 * Serializes FASTA sequences into bytes, for communication.
 */
char *serialize_fasta_seqs(std::vector<fasta_seq_t>& fasta_seqs, size_t& num_bytes);

/**
 * Deserializes bytes into FASTA sequences.
 */
std::vector<fasta_seq_t> deserialize_fasta_seqs(char *bytes, size_t num_bytes);

/**
 * Custom operation for reducing accepts and flags. Has type MPI_User_function.
 */
void accept_op(void *in, void *inout, int *len, MPI_Datatype *dptr);

#endif
