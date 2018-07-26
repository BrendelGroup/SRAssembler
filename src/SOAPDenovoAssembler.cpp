/*
 * SOAPDenovoAssembler.cpp
 *
 *  Created on: Oct 15, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 *  Updated for use with SOAPdenovo2: Aug 21, 2013 (VB)
 */

#include "SOAPDenovoAssembler.h"
#include <iostream>
#include <fstream>

SOAPDenovoAssembler::SOAPDenovoAssembler(int log_level, string log_file):Assembler(log_level, log_file){
}

SOAPDenovoAssembler::~SOAPDenovoAssembler() {
	// Auto-generated destructor stub
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
void SOAPDenovoAssembler::do_assembly(int kmer, const vector<Library>& libraries, const string& output_file, int threads, boost::unordered_map<std::string,Params> parameters_dict){
	Params pregraph_params;
	Params contig_params;
	string pregraph_param_list = "";
	string contig_param_list = "";
	string config_file = output_file + ".conf";
	ofstream outFile(config_file.c_str());
	// Prepare the SOAPdenovo configuration file specific for this run.
	outFile << "max_rd_len=5000" << '\n';
	for (unsigned int i=0; i<libraries.size();i++){
		Library lib = libraries[i];
		if (get_file_size(lib.get_matched_left_reads_filename()) == 0) continue;
		outFile << "[LIB]" << '\n';
		outFile << "asm_flags=3" << '\n';
		outFile << "rank=" << (i+1) << '\n';
		if (lib.get_paired_end()) {
			outFile << "avg_ins=" << lib.get_insert_size() << '\n';
			outFile << "f1=" << lib.get_matched_left_reads_filename() << '\n';
			outFile << "f2=" << lib.get_matched_right_reads_filename() << '\n';
			if (lib.get_reversed()) {
				outFile << "reverse_seq=1" << '\n';
			}
		} else {
				outFile << "f=" << lib.get_matched_left_reads_filename() << '\n';
		}
	}
	outFile.close();
	// Import the parameters for the SOAPdenovo pregraph command.
	pregraph_params = parameters_dict["SOAPdenovo_pregraph"];
	for ( Params::const_iterator it = pregraph_params.begin(); it != pregraph_params.end(); ++it ){
		pregraph_param_list += " -" + it->first + " " + it->second;
	}
	// Import the parameters for the SOAPdenovo contig command.
	contig_params = parameters_dict["SOAPdenovo_contig"];
	for ( Params::const_iterator it = contig_params.begin(); it != contig_params.end(); ++it ){
		contig_param_list += " -" + it->first + " " + it->second;
	}
	string program = "SOAPdenovo-127mer";
	if (kmer <= 63)
		program = "SOAPdenovo-63mer";
	string cmd = program + " pregraph -s " + config_file + " -o " + output_file + " -p " + int2str(threads) + " -K " + int2str(kmer) + pregraph_param_list + " >> " + logger->get_log_file() + " 2>&1";
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = program + " contig -g " + output_file + contig_param_list + " >> " + logger->get_log_file() + " 2>&1";
	logger->debug(cmd);
	run_shell_command(cmd);
	standardize_contigs(get_output_contig_file_name(output_file));
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
