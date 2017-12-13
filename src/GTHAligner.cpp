/*
 * GTHAligner.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#include "GTHAligner.h"
#include <iostream>
#include <fstream>

GTHAligner::GTHAligner(int log_level, string log_file):SplicedAligner(log_level, log_file) {
}
GTHAligner::~GTHAligner() {
}

bool GTHAligner::is_available(){
	int ret = system("gth > /dev/null 2>&1");
	if (WEXITSTATUS(ret) != 1) {
		cout << "Cannot find gth, check your PATH variable!" << endl;
		return false;
	}
	char* env = getenv ("BSSMDIR");
	if (env==NULL) {
		cout << "Cannot find $BSSMDIR variable!" << endl;
		return false;
	}
	env = getenv ("GTHDATADIR");
	if (env==NULL) {
		cout << "Cannot find $GTHDATADIR variable!" << endl;
		return false;
	}
	return true;
}

void GTHAligner::do_spliced_alignment(const string& genomic_file, const string& type, const string& query_file, const string& species, const Params& params, const string& output_file, const string& hit_contig_file){
	string param_list = "";
	for ( Params::const_iterator it = params.begin(); it != params.end(); ++it ){
		param_list += " -" + it->first + " " + it->second;
	}
	string cmd = "gth -genomic " + genomic_file + " -" + type + " " + query_file + " -species " + species + param_list + " > " + output_file;
	logger->debug(cmd);
	run_shell_command(cmd);
	//get_aligned_contigs(genomic_file, hit_contig_file, output_file);
}


string_map GTHAligner::get_aligned_contigs(const double& min_score, const double& min_coverage, const unsigned int& min_contig_lgth, const string& all_contig_file, const string& hit_contig_file, const string& alignment_file){
	ifstream alignment_fs(alignment_file.c_str());
	ofstream new_contig_fs(hit_contig_file.c_str());
	string line;
	vector<string> contig_list;
	char currentid[20];
	unsigned int contig_length;
	aligned_query_list.clear();
	logger->debug("Finding the aligned contigs");
	num_matches = 0;
	output_string = "<B>Contig\tStrand\tQuery\tStrand\tScore\tLength\tCov\tG/P/C</B>\n";
	output_string = "<B>" + string_format("%-15s %-8s %-30s %-8s %-10s %-15s %-15s %-15s","Contig","Strand","Query","Strand","Score","Length","Coverage","G/P/C") + "</B>\n";
	output_string += "--------------------------------------------------------------------------------------------------------------------\n";
	while (getline(alignment_fs, line)) {
		if (line.substr(0,16) == "Genomic Template"){
			sscanf(line.c_str(),"Genomic Template: %*s %*s %*s %*s description=%s length %d %*s",currentid,&contig_length);
			logger->info("... checking contig:\t" + std::string(currentid) + "\tof length:\t" + int2str(contig_length));
		}
		if (line.substr(0,5) == "MATCH"){
			vector<string> tokens;
			tokenize(line, tokens, "\t");
			string contig_id = tokens[1].substr(0, tokens[1].length()-1);
			contig_list.push_back(contig_id);
			string strand = tokens[1].substr(tokens[1].length()-1,1);
			string query = tokens[2].substr(0,tokens[2].length()-1);;
			string qstrand = tokens[2].substr(tokens[2].length()-1,1);
			string score = tokens[3];
			string length = tokens[4];
			string cov = tokens[5];
			string type = tokens[6];
			output_string += string_format("%-15s %-8s %-30s %-8s %-10s %-15s %-15s %-15s",contig_id.c_str(),strand.c_str(),query.c_str(),qstrand.c_str(),score.c_str(),length.c_str(),cov.c_str(),type.c_str()) + "\n";
			num_matches++;
			output_string += "\n";
			logger->info("   ... MATCH found with coverage:\t" + cov + "\tscore:\t" + score + "\tlength:\t" + int2str(contig_length));
			if (type == "P" || type == "C"){
				if (str2double(cov) > min_coverage && str2double(score) > min_score && contig_length >= min_contig_lgth)
					aligned_query_list[query] = contig_id;
			}
		}
	}
	output_string += "\nLength: cumulative length of scored exons\nCov G/P/C: coverage of contig (G) or cDNA (C) or protein (P), whichever is highest";
	alignment_fs.close();

	ifstream old_contig_fs(all_contig_file.c_str());
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

string GTHAligner::get_program_name(){
	return "GenomeThreader";
}
void GTHAligner::clean_files(const string& file){
	string cmd = "rm " + file + ".dna.*; rm " + file + ".md5";
	logger->debug(cmd);
	run_shell_command(cmd);
}
