/*
 * Snap.h
 *
 *  Created on: 2012-11-14
 *      Author: hchou
 */

#ifndef SNAPGENEFINDER_H_
#define SNAPGENEFINDER_H_

#include "GeneFinder.h"

class SnapGeneFinder: public GeneFinder {
public:
	SnapGeneFinder(int, string);
	virtual ~SnapGeneFinder();
	void do_gene_finding(const string& genomic_file, const string& species, const Params& params, const string& output_file, const string& protein_output_file);
	string get_output_summary();
	string get_program_name();
	bool is_available();
};

#endif /* SNAPGENEFINDER_H_ */
