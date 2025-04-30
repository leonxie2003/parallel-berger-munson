/** @file parse_fasta.h
 *  Parse FASTA files.
 *  @author Leon Xie (leonx)
 *  @author Taekseung Kim (taekseuk)
 */

#ifndef __PARSE_FASTA_H__
#define __PARSE_FASTA_H__

#include <string>
#include <vector>

/**
 * Represents a sequence parsed from a FASTA file. Includes an identifier,
 * description, and the actual sequence.
 */
typedef struct fasta_seq {
    std::string ident;
    std::string desc;
    std::string seq;
} fasta_seq_t;


/**
 * Parses a FASTA file into memory.
 */
std::vector<fasta_seq_t> parse_fasta(std::string filename);

/**
 * Prints a FASTA seq's data to std::cout.
 */
void print_fasta_seq(fasta_seq_t seq);

#endif
