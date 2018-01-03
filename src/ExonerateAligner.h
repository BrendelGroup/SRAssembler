/*
 * ExonerateAligner
 *
 *  Created on: Nov 14, 2012
 *      Author: hchou
 */

#ifndef EXONERATEALIGNER_H_
#define EXONERATEALIGNER_H_

#include "SplicedAligner.h"

class ExonerateAligner: public SplicedAligner {
public:
	ExonerateAligner(int, string);
	virtual ~ExonerateAligner();
	bool is_available();
	void do_spliced_alignment(const string& genomic_file, const string& type, const string& query_file, const string& species, const Params& params, const string& output_file);
	string_map get_aligned_contigs(const double& min_score, const double& min_coverage, const unsigned int& min_contig_lgth, const string& all_contig_file, const string& hit_contig_file, const string& alignment_file);
	string get_program_name();
	void clean_files(const string&);
};

#endif /* EXONERATEALIGNER_H_ */
