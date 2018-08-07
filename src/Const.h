/*
 * Const.h
 *
 *  Created on: Oct 21, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#ifndef CONST_H_
#define CONST_H_
#include <string>
#include <boost/unordered_map.hpp>

// Codes for MPI and speed advantages.
static const int ACTION_EXIT=1;
static const int ACTION_SPLIT=2;
static const int ACTION_PRE_PROCESSING=3;
static const int ACTION_ALIGNMENT=4;
static const int ACTION_ASSEMBLY=5;
static const int ACTION_TOTAL_PARTS=6;
static const int ACTION_RETURN=7;
static const int ACTION_LOAD_PREVIOUS=8;
static const int ACTION_MEMDIR=9;
static const int ACTION_SAVE=10;
static const int ACTION_CLEAN=11;
static const int FASTQ_JOINED=1;
static const int FASTQ_INTERLEAVED=2;
static const int FORMAT_FASTQ=0;
static const int FORMAT_FASTA=1;
static const int DIRECTION_FR=0;
static const int DIRECTION_RF=1;
static const int LEFT_READ=0;
static const int RIGHT_READ=1;
static const std::string TYPE_PROTEIN = "protein";
static const std::string TYPE_DNA = "dna";

// Default parameters.
static const int CLEAN_ROUND=4;
static const int CONTIG_LIMIT=500;
static const int SPLICED_ALIGNMENT_PROGRAM = 1;
static const int ASSEMBLER_PROGRAM = 0;
static const int GENE_FINDING_PROGRAM = 0;
static const int ASSEMBLY_ROUND=1;
static const int READS_PER_FILE=500000;
static const int START_K=15;
static const int END_K=45;
static const int STEP_K=10;
static const int END_SEARCH_LENGTH=0;
static const int END_SEARCH_ROUND=2;
static const int MASKING_ROUND=1;
static const int NUM_ROUNDS=10;
static const int INSERT_SIZE=300;
static const int INIT_MATCH_LENGTH_PROTEIN=10;
static const int INIT_MATCH_LENGTH_DNA=30;
static const int MISMATCH_ALLOWED_PROTEIN=1;
static const int MISMATCH_ALLOWED_DNA=2;
static const int RECUR_MATCH_LENGTH=30;
static const double MIN_SCORE=0.5;
static const double MIN_COVERAGE=0.8;
static const unsigned int QUERY_CONTIG_MIN=200;
static const unsigned int MIN_CONTIG_LGTH=200;
static const unsigned int MAX_CONTIG_LGTH=10000;
static const bool VERBOSE=false;
static const bool PREPROCESSING_ONLY=false;
static const bool OVER_WRITE = false;
static const bool CHECK_GENE_ASSEMBLED = true;
static const std::string READS_DATA = "processed_reads";
static const std::string QUERY_TYPE = TYPE_PROTEIN;
static const std::string DEFAULT_SPECIES = "arabidopsis";
static const std::string OUT_DIR = "./SRAssembler_output";
static const std::string TMP_LOC = "/tmp";



typedef boost::unordered_map<std::string,std::string> Params;

#endif /* CONST_H_ */
