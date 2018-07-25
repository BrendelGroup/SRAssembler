/*
 * AbyssAssembler.cpp
 *
 *  Created on: Oct 15, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#include "AbyssAssembler.h"


AbyssAssembler::AbyssAssembler(int log_level, string log_string):Assembler(log_level, log_string) {
}

AbyssAssembler::~AbyssAssembler() {
}

bool AbyssAssembler::is_available(){
	int ret = system("abyss-pe > /dev/null 2>&1");
	if (WEXITSTATUS(ret) != 0 && WEXITSTATUS(ret) != 2) {
		cout << "Cannot find ABySS, check your PATH variable!" << WEXITSTATUS(ret) << endl;
		return false;
	}
	return true;
}

void AbyssAssembler::do_assembly(int kmer, const vector<Library>& libraries, const string& output_file, int threads, boost::unordered_map<std::string,Params> parameters_dict){
	Params contig_params;
	string contig_param_list = "";
	string lib_list = "";
	string paired_files = "";
	string single_files = "";
	for (unsigned int i=0; i<libraries.size();i++){
		Library lib = libraries[i];
		if (get_file_size(lib.get_matched_left_reads_filename()) == 0) continue;
		if (lib.get_paired_end()) {
			string lib_name = "pe" + int2str(lib.get_insert_size());
			lib_list += lib_name + " ";
			paired_files += lib_name + "='" + lib.get_matched_left_reads_filename() + " " + lib.get_matched_right_reads_filename() + "' ";
		} else
			single_files += lib.get_matched_left_reads_filename() + " ";
	}
	if (lib_list != ""){
		lib_list = "lib='" + lib_list + "'";
	}
	// Import the parameters for the Abyss "abyss-pe contig" command.
	contig_params = parameters_dict["Abyss_contig"];
	for ( Params::const_iterator it = contig_params.begin(); it != contig_params.end(); ++it ){
		contig_param_list += " -" + it->first + " " + it->second;
	}
	single_files = "se='" + single_files + "'";
	string cmd = "abyss-pe contigs k=" + int2str(kmer) + " name=" + output_file + " " + lib_list + " " + paired_files + " " + single_files + contig_param_list + ">> " + logger->get_log_file() + " 2>&1";
	logger->debug(cmd);
	run_shell_command(cmd);
	standardize_contigs(get_output_contig_file_name(output_file));
}

void AbyssAssembler::clean_files(const string& dir){
	string cmd = "rm -rf " + dir + "/assembly_*;" + " rm -rf coverage.hist pe*hist pe*dist pe*sam.gz";
	logger->debug(cmd);
	run_shell_command(cmd);
}

string AbyssAssembler::get_output_contig_file_name(string prefix){
	int i=1;
	while (file_exists(prefix + "-" + int2str(i) + ".fa")){
		i++;
	}
	return prefix + "-" + int2str(i-1) + ".fa";
}

string AbyssAssembler::get_output_scaffold_file_name(string prefix){
	return get_output_contig_file_name(prefix);
}
