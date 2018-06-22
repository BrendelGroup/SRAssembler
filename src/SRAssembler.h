/*
 * SRAssembler.h
 *
 *  Created on: Oct 12, 2011
 *      Author: hchou
 */

#ifndef SRASSEMBLER_H_
#define SRASSEMBLER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include "Utility.h"
#include "Const.h"
#include "Library.h"
#include "Aligner.h"
#include "Assembler.h"
#include "SplicedAligner.h"
#include "GeneFinder.h"
#include "MPIWrapper.h"
#include "Logger.h"
#include <boost/unordered_set.hpp>

using namespace std;

typedef struct {
	int total_contig;
	int n50;
	int n90;
	unsigned int longest_contig;
} Assembly_stats;

class SRAssembler {
public:
	SRAssembler();
	virtual ~SRAssembler();
	virtual int init(int argc, char * argv[], int rank, int mpiSize);
	virtual void do_preprocessing()=0;
	virtual void do_walking()=0;
	virtual void show_usage()=0;
	virtual void print_message(const std::string&)=0;
	static SRAssembler* getInstance(int pid);
	Logger* get_logger();
	//void signalHandler(int signum);
protected:
	int do_alignment(int round, int lib_idx, int idx);
	void do_spliced_alignment();
	string_map do_spliced_alignment(int);
	void do_gene_finding();
	void do_assembly(int, int, int threads);
	void create_index(int round);
	std::string get_contigs_index_name(int round);
	std::string get_query_fasta_file_name(int round);
	void mask_contigs(int round);
	std::string get_query_fasta_file_name_masked(int round);
	std::string get_contig_file_name(int round);
	std::string get_matched_reads_file_name(int round);
	std::string get_vmatch_output_filename(int round, int lib_idx, int idx);
	std::string get_type(int round);
	int get_match_length(int round);
	int get_mismatch_allowed(int round);
	std::string get_assembly_file_name(int round, int k);
	std::string get_assembled_scaf_file_name(int round, int k);
	void do_split_files(string read_file);
	void preprocess_read_part(int lib_idx, int file_part);
	int get_file_count(std::string);
	int count_preprocessed_reads(int lib_idx);
	void merge_mapped_files(int round);
	Assembly_stats get_assembly_stats(int round, int k);
	void save_found_reads(int round);
	void load_found_reads(int round);
	void remove_unmapped_reads(unsigned int lib_idx, int round);
	//void prepare_contig_file(int round, int k);
	void keep_long_contigs(string in_file, string out_file, unsigned int min_length);
	long get_total_read_count(int round);
	void send_code(const int& to, const int& action, const int& value1, const int& value2, const int& value3);
	void broadcast_code(const int& action, const int& value1, const int& value2, const int& value3);
	Aligner* get_aligner(int round);
	Assembler* get_assembler();
	SplicedAligner* get_spliced_aligner();
	GeneFinder* get_gene_finder();
	boost::unordered_map<std::string,Params> read_param_file();
	Params get_parameters(string program_name);
	std::string query_file, species, type, out_dir;
	int init_match_length;
	int recur_match_length;
	int mismatch_allowed;
	int num_rounds;
	int verbose;
	int preprocessing_only;
	int assembly_round;
	int clean_round;
	int contig_limit;
	int over_write;
	int check_gene_assembled;
	int reads_per_file;
	int start_k;
	int end_k;
	int step_k;
	int mpiSize;
	int rank;
	int fastq_format;
	int start_round;
	int spliced_alignment_program;
	int gene_finding_program;
	int assembler_program;
	int merge_factor;
	int edge_cov_cutoff;
	int masking;
	double min_score;
	double min_coverage;
	// A dictionary for tracking the best contigs found between rounds.
	tuple_map best_hits;
	unsigned int ini_contig_size;
	unsigned int min_contig_lgth;
	unsigned int max_contig_lgth;
	bool preprocessed_exist;
	bool ignore_contig_explosion;
	std::string library_file;
	std::string tmp_dir;
	std::string mem_dir;
	std::string intermediate_dir;
	std::string data_dir;
	std::string results_dir;
	std::string log_file;
	std::string param_file;
	std::string spliced_alignment_output_file;
	std::string gene_finding_output_file;
	std::string gene_finding_output_protein_file;
	std::string usage;
	std::string final_scaf_file;
	std::string final_contigs_file;
	std::string hit_contigs_file;
	std::string final_long_contig_file;
	std::string summary_file;
	std::string mapped_readnumbers_file;
	boost::unordered_set<std::string> found_reads;
	std::vector<Library> libraries;
	Logger* logger;
	boost::unordered_map<std::string,Params> parameters_dict;
private:
	static SRAssembler* _srassembler;
	bool read_library_file();
};


#endif /* SRASSEMBLER_H_ */
