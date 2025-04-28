#ifndef __BM_COMM_H__
#define __BM_COMM_H__

#include "align.h"
#include "parse_fasta.h"

#include <vector>
#include <mpi.h>

typedef struct pid_flag {
    int pid;
    int flag;
} pid_flag_t;

char *serialize_fasta_seqs(std::vector<fasta_seq_t>& fasta_seqs, size_t& num_bytes);

std::vector<fasta_seq_t> deserialize_fasta_seqs(char *bytes, size_t num_bytes);

void accept_op(void *in, void *inout, int *len, MPI_Datatype *dptr);

#endif
