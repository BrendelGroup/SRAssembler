/*
 * AbyssAssembler.h
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#ifndef ABYSSASSEMBLER_H_
#define ABYSSASSEMBLER_H_

#include "Assembler.h"

class AbyssAssembler: public Assembler {
public:
	AbyssAssembler(int, string);
	void do_assembly(int kmer, const vector<Library>& libraries, const string& output_file);
	bool is_available();
	void clean_files(const string& dir);
	string get_output_contig_file_name(string prefix);
	string get_output_scaffold_file_name(string prefix);
	virtual ~AbyssAssembler();
};

#endif /* ABYSSASSEMBLER_H_ */
