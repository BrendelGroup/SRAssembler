/*
 * SnapGeneFinder.cpp
 *
 *  Created on: 2012-11-14
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#include "SnapGeneFinder.h"

SnapGeneFinder::SnapGeneFinder(int log_level, string log_file):GeneFinder(log_level, log_file){
	// Auto-generated constructor stub
}

SnapGeneFinder::~SnapGeneFinder() {
	// Auto-generated destructor stub
}

void SnapGeneFinder::do_gene_finding(const string& genomic_file, const string& species, const Params& params, const string& output_file, const string& protein_output_file){
	string param_list = "";
	string hmm = "";
	for ( Params::const_iterator it = params.begin(); it != params.end(); ++it ){
		if (it->first == "snaphmm")
			hmm = it->second;
		param_list += " -" + it->first + " " + it->second;
	}
	string cmd = "snap " + hmm + " " + genomic_file + " -aa " + protein_output_file + " > " + output_file + " 2>&1";
	logger->debug(cmd);
	run_shell_command(cmd);

}

string SnapGeneFinder::get_output_summary(){
	return "";
}
string SnapGeneFinder::get_program_name(){
	return "Snap";
}
bool SnapGeneFinder::is_available(){
	int ret = system("snap > /dev/null 2>&1");
	if (WEXITSTATUS(ret) != 1) {
		logger->error("Cannot find snap, check your PATH variable!");
		return false;
	}
	return true;
}
