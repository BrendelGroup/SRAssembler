/*
 * VmatchAligner.h
 *
 *  Created on: Oct 15, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#ifndef VMATCHALIGNER_H_
#define VMATCHALIGNER_H_

#include "Aligner.h"

using namespace boost;

class VmatchAligner: public Aligner {
public:
	VmatchAligner(int,string);
	void create_index(const string& index_name, const string& db_type, const string& fasta_file);
	void do_alignment(const string& index_name, const string& alignment_type, int match_length, int mismatch_allowed, const string& reads_file, const Params& params, const string& output_file);
	int parse_output(const string& output_file, unordered_set<string>& found_reads, const int lib_idx, const int read_chunk, const string& left_read_index, const string& right_read_index, const string& out_left_read, const string& out_right_read);
	int ignore_output(const string& output_file, unordered_set<string>& found_reads, const int lib_idx, const int read_chunk);
	void align_long_contigs(const string& long_contig_candidate_file, const string& aux_dir, const string& contig_file, const int max_contig_size, boost::unordered_set<string>& candidate_ids, boost::unordered_set<string>& long_contig_ids);
	bool is_available();
	string get_program_name();
	int get_format();
	virtual ~VmatchAligner();
};

#endif /* VMATCHALIGNER_H_ */
