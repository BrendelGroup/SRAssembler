/*
 * VmatchBWAfastmapAligner.h
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#ifndef VMATCHBWAFASTMAPALIGNER_H_
#define VMATCHBWAFASTMAPALIGNER_H_

#include "Aligner.h"

using namespace boost;

class VmatchBWAfastmapAligner: public Aligner {
public:
	VmatchBWAfastmapAligner(int,string);
	void create_index(const string& index_name, const string& type, const string& fasta_file);
	void do_alignment(const string& fasta_file, const string& index_name, const string& type, int match_length, int mismatch_allowed, const string& reads_file, const Params& params, const string& output_file);
	int parse_output(const string& output_file, unordered_set<string>& mapped_reads, const string& source_read, const string& out_left_read, const string& out_right_read, int fastq_format, int format);
	unordered_set<string> get_hit_list(const string& output_file);
	void align_long_contigs(const string& long_contig_candidate_file, const string& tmp_dir, const string& contig_file, const int max_contig_size, boost::unordered_set<string>& candicate_ids, boost::unordered_set<string>& long_contig_ids);
	bool is_available();
	string get_program_name();
	int get_format();
	virtual ~VmatchBWAfastmapAligner();
};

#endif /* VMATCHBWAFASTMAPALIGNER_H_ */
