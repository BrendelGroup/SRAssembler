/*
 * SOAPDenovoAssembler.h
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#ifndef SOAPDENOVOASSEMBLER_H_
#define SOAPDENOVOASSEMBLER_H_

#include "Assembler.h"

class SOAPDenovoAssembler: public Assembler {
public:
	SOAPDenovoAssembler(int, string);
	void do_assembly(int kmer, const vector<Library>& libraries, const string& output_file);
	bool is_available();
	void clean_files(const string& dir);
	string get_output_contig_file_name(string prefix);
	string get_output_scaffold_file_name(string prefix);
	virtual ~SOAPDenovoAssembler();
};

#endif /* SOAPDENOVOASSEMBLER_H_ */
