/*
 * ExonerateAligner.cpp
 *
 *  Created on: Nov 14, 2012
 *      Author: hchou
 */

#include "ExonerateAligner.h"
#include <iostream>

ExonerateAligner::ExonerateAligner(int log_level, string log_file):SplicedAligner(log_level, log_file) {
}
ExonerateAligner::~ExonerateAligner() {
}

bool ExonerateAligner::is_available(){
	int ret = system("exonerate > /dev/null 2>&1");
	if (WEXITSTATUS(ret) != 1) {
		cout << "Cannot find exonerate, check your PATH variable!" << endl;
		return false;
	}
	return true;
}

void ExonerateAligner::do_spliced_alignment(const string& genomic_file, const string& type, const string& query_file, const string& species, const Params& params, const string& output_file, const string& hit_contig_file){
	string param_list = "";
	for ( Params::const_iterator it = params.begin(); it != params.end(); ++it ){
		string param = it->first;
		string value = it->second;
		string dash = (param.length() == 1)? "-" : "--";
		param_list += " "  + dash + param + " " + value;
	}
	string type_str = "est2genome";
	if (type == "protein")
		type_str = "protein2genome";
	string cmd = "exonerate --model " + type_str + param_list + " " + query_file + " " + genomic_file + " > " + output_file;
	logger->debug(cmd);
	run_shell_command(cmd);
	//get_aligned_contigs(genomic_file, hit_contig_file, output_file);
}

string_map ExonerateAligner::get_aligned_contigs(const double& min_score, const double& min_coverage, const int& min_contig_lgth, const string& all_contig_file, const string& hit_contig_file, const string& alignment_file){
	ifstream old_contig_fs(all_contig_file.c_str());
	ifstream alignment_fs(alignment_file.c_str());
	ofstream new_contig_fs(hit_contig_file.c_str());
	string line;
	vector<string> contig_list;
	logger->debug("Finding the aligned contigs");
	num_matches = 0;
	output_string = "<B>" + string_format("%-15s %-8s %-30s %-10s %-15s %-15s","Contig","Strand","Query","Score","Query_Range","Target_range") + "</B>\n";
	output_string += "-----------------------------------------------------------------------------------------------\n";
	while (getline(alignment_fs, line)) {
		if (line.substr(0,7) == "vulgar:"){
			vector<string> tokens;
			tokenize(line, tokens, " ");
			string contig_id = tokens[5];
			contig_list.push_back(contig_id);
			string contig = tokens[5];
			string strand = tokens[8];
			string query = tokens[1];
			string score = tokens[9];
			string query_range = tokens[2] + " - " + tokens[3];
			string contig_range = tokens[6] + " - " + tokens[7];
			int range = str2int(tokens[3]) -str2int(tokens[2]) + 1;
			output_string += string_format("%-15s %-8s %-30s %-10s %-15s %-15s",contig.c_str(),strand.c_str(),query.c_str(),score.c_str(),query_range.c_str(),contig_range.c_str()) + "\n";
			num_matches++;
			if (range > min_coverage && str2int(score) > min_score){
				aligned_query_list[query] = contig_id;
			}
		}
	}
	alignment_fs.close();

	while (getline(old_contig_fs, line)) {
		if (line.substr(0,1) == ">"){
			vector<string> tokens;
			tokenize(line, tokens, " ");
			string contig_id = tokens[0].substr(1, tokens[0].length()-1);
			if (std::find(contig_list.begin(), contig_list.end(), contig_id)!=contig_list.end()){
				new_contig_fs << line << endl;
				getline(old_contig_fs, line);
				new_contig_fs << line << endl;
			}
		}
	}
	old_contig_fs.close();
	new_contig_fs.close();
	return aligned_query_list;
}

string ExonerateAligner::get_program_name(){
	return "Exonerate";
}
void ExonerateAligner::clean_files(const string& file){

}
