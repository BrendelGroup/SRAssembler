/*
 * GTHAligner.h
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#ifndef GTHALIGNER_H_
#define GTHALIGNER_H_

#include "SplicedAligner.h"

class GTHAligner: public SplicedAligner {
public:
	GTHAligner(int, string);
	virtual ~GTHAligner();
	bool is_available();
	void do_spliced_alignment(const string& genomic_file, const string& type, const string& query_file, const string& species, const Params& params, const string& output_file);
	string_map get_aligned_contigs(const double& min_score, const double& min_coverage, const unsigned int& min_contig_lgth, const string& all_contig_file, const string& hit_contig_file, const string& alignment_file);
	string get_program_name();
	void get_hit_contigs(const double& min_score, const double& min_coverage, const unsigned int& min_contig_lgth, const string& final_contigs_file, const string& hit_contig_file, const string& alignment_file);
	void clean_files(const string&);
};

#endif /* GTHALIGNER_H_ */
