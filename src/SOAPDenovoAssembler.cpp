/*
 * SOAPDenovoAssembler.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 *  Updated for use with SOAPdenovo2: Aug 21, 2013 (VB)
 */

#include "SOAPDenovoAssembler.h"
#include <iostream>
#include <fstream>

SOAPDenovoAssembler::SOAPDenovoAssembler(int log_level, string log_file):Assembler(log_level, log_file){
}

SOAPDenovoAssembler::~SOAPDenovoAssembler() {
	// TODO Auto-generated destructor stub
}
bool SOAPDenovoAssembler::is_available(){
	int ret_127 = system("SOAPdenovo-127mer > /dev/null 2>&1");
	if (WEXITSTATUS(ret_127) != 0 && WEXITSTATUS(ret_127) != 1) {
		cout << "Cannot find SOAPdenovo-127mer, check your PATH variable!" << endl;
		return false;
	}
	int ret_63 = system("SOAPdenovo-63mer > /dev/null 2>&1");
	if (WEXITSTATUS(ret_63) != 0 && WEXITSTATUS(ret_63) != 1) {
		cout << "Cannot find SOAPdenovo-63mer, check your PATH variable!" << endl;
		return false;
	}
	return true;
}
void SOAPDenovoAssembler::do_assembly(int kmer, const vector<Library>& libraries, const string& output_file)
{
	string config_file = output_file + ".conf";
	ofstream outFile(config_file.c_str());

	for (unsigned int i=0; i<libraries.size();i++){
		Library lib = libraries[i];
		if (get_file_size(lib.get_matched_left_read_name()) == 0) continue;
		outFile << "[LIB]" << endl;
		outFile << "asm_flags=3" << endl;
		outFile << "rank=" << (i+1) << endl;
		if (lib.get_paired_end()) {
			outFile << "avg_ins=" << lib.get_insert_size() << endl;
			if (lib.get_format() == FORMAT_FASTQ) {
				outFile << "q1=" << lib.get_matched_left_read_name() << endl;
				outFile << "q2=" << lib.get_matched_right_read_name() << endl;
			} else {
				outFile << "f1=" << lib.get_matched_left_read_name() << endl;
				outFile << "f2=" << lib.get_matched_right_read_name() << endl;
			}
			if (lib.get_reversed()) {
				outFile << "reverse_seq=1" << endl;
			}
		} else {
			if (lib.get_format() == FORMAT_FASTQ)
				outFile << "q=" << lib.get_matched_left_read_name() << endl;
			else
				outFile << "f=" << lib.get_matched_left_read_name() << endl;
		}
	}
	outFile.close();
	string program = "SOAPdenovo-127mer";
	if (kmer <= 63)
		program = "SOAPdenovo-63mer";
	string cmd = program + " pregraph -s " + config_file + " -o " + output_file + " -p 1 -K " + int2str(kmer) + " >> " + logger->get_log_file() + " 2>&1";
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = program + " contig -g " + output_file + " >> " + logger->get_log_file() + " 2>&1";
	logger->debug(cmd);
	run_shell_command(cmd);
}

void SOAPDenovoAssembler::clean_files(const string& dir){
	string cmd = "rm -rf " + dir + "/assembly_*";
	logger->debug(cmd);
	run_shell_command(cmd);
}

string SOAPDenovoAssembler::get_output_contig_file_name(string prefix){
	return prefix + ".contig";
}

string SOAPDenovoAssembler::get_output_scaffold_file_name(string prefix){
	return prefix + ".scafSeq";
}
