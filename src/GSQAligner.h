/*
 * GSQAligner.h
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#ifndef GSQALIGNER_H_
#define GSQALIGNER_H_

#include "SplicedAligner.h"
#include <boost/unordered_map.hpp>

class GSQAligner: public SplicedAligner {
public:
	GSQAligner(int, string);
	virtual ~GSQAligner();
	bool is_available();
	string_map get_aligned_contigs(const double& min_score, const double& min_coverage, const unsigned int& min_contig_lgth, const string& all_contig_file, const string& hit_contig_file, const string& alignment_file);
	void do_spliced_alignment(const string& genomic_file, const string& type, const string& query_file, const string& species, const Params& params, const string& output_file);
	string get_program_name();
	void clean_files(const string&);
private:
	boost::unordered_map<string, string> species_names;
};

#endif /* GSQALIGNER_H_ */
