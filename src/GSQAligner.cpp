/*
 * GSQAligner.cpp
 *
 *  Created on: Oct 15, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#include "GSQAligner.h"
#include <iostream>
#include <tuple>

GSQAligner::GSQAligner(int log_level, string log_file):SplicedAligner(log_level, log_file) {
	species_names["arabidopsis"] = "Arabidopsis";
	species_names["drosophila"] = "Drosophila";
	species_names["human"] = "human";
	species_names["mouse"] = "mouse";
	species_names["rat"] = "rat";
	species_names["chicken"] = "chicken";
	species_names["nematode"] = "nematode";
	species_names["fission_yeast"] = "yeast";
	species_names["aspergillus"] = "Aspergillus";
	species_names["maize"] = "maize";
	species_names["rice"] = "rice";
	species_names["medicago"] = "Medicago";
}
GSQAligner::~GSQAligner() {
}

bool GSQAligner::is_available(){
	int ret = system("GeneSeqer > /dev/null 2>&1");
	if (WEXITSTATUS(ret) != 0) {
		cout << "Cannot find GeneSeqer, check your PATH variable!" << endl;
		return false;
	}
	return true;
}

void GSQAligner::do_spliced_alignment(const string& genomic_file, const string& query_type, const string& query_file, const string& species, const Params& params, const string& output_file){
	string param_list = "";
	for ( Params::const_iterator it = params.begin(); it != params.end(); ++it ){
		param_list += " -" + it->first + " " + it->second;
	}
	string species_str = species_names[species];
	if (species_str == "")
		species_str = "generic";
	string type_str = "E";
	if (query_type == "protein")
		type_str = "Q";
	string cmd = "GeneSeqer -L " + genomic_file + " -" + type_str + " " + query_file + " -species " + species_str + " " + param_list + " -o " + output_file + " >> /dev/null 2>> " + logger->get_log_file();
	logger->debug(cmd);
	run_shell_command(cmd);
}

// This is for each round, to see if the ending criteria have been met.
string_map GSQAligner::get_aligned_contigs(const double& min_score, const double& min_coverage, const unsigned int& min_contig_lgth, const string& all_contig_file, const string& hit_contig_file, const string& alignment_file, const int round, tuple_map& best_hits){
	double best_coverage = std::get<1>(best_hits["coverage"]);
	ifstream old_contig_fs(all_contig_file.c_str());
	ifstream alignment_fs(alignment_file.c_str());
	ofstream new_contig_fs(hit_contig_file.c_str());
	string line;
	char currentid[20];
	unsigned int contig_length;
	vector<string> contig_list;
	logger->running("Finding the aligned contigs");
	num_matches = 0;
	output_string = "<B>" + string_format("%-15s %-8s %-30s %-10s %-15s %-15s %-15s","Contig","Strand","Query","Score","Length","Coverage","G/P/C") + "</B>\n";
	output_string += "----------------------------------------------------------------------------------------------------------\n";
	while (getline(alignment_fs, line)) {
		if (line.substr(0,8) == "Sequence"){
			sscanf(line.c_str(),"Sequence %*d: %80[^,], from %*d to %d,%*s",currentid,&contig_length);
			logger->debug("... checking contig:\t" + std::string(currentid) + "\tof length:\t" + int2str(contig_length));
		}
		if (line.substr(0,5) == "MATCH"){
			vector<string> tokens;
			tokenize(line, tokens, "\t");
			string contig_id = tokens[1].substr(0, tokens[1].length()-1);
			contig_list.push_back(contig_id);
			string contig = tokens[1].substr(0,tokens[1].length()-1);
			string strand = tokens[1].substr(tokens[1].length()-1,1);
			string query = tokens[2];
			string score = tokens[3];
			string length = tokens[4];
			string cov = tokens[5];
			string type = tokens[6];
			output_string += string_format("%-15s %-8s %-30s %-10s %-15s %-15s %-15s",contig.c_str(),strand.c_str(),query.c_str(),score.c_str(),length.c_str(),cov.c_str(),type.c_str()) + "\n";
			num_matches++;
			output_string += "\n";
			logger->info("   ... MATCH found with coverage:\t" + cov + " " + type + "\tscore:\t" + score + "\tcontig length:\t" + int2str(contig_length));
			if (type == "P" || type == "C"){
				if (str2double(score) > min_score && str2double(cov) > best_coverage) {
					get<0>(best_hits["coverage"]) = round;
					get<1>(best_hits["coverage"]) = str2double(cov);
					best_coverage = str2double(cov);
				}
				if (str2double(score) > min_score && str2double(cov) > min_coverage && contig_length >= min_contig_lgth)
					aligned_query_list[query] = contig_id;
			}
		}
	}
	output_string += "\nLength: cumulative length of scored exons\nCov G/P/C: coverage of contig (G) or cDNA (C) or protein (P), whichever is highest";
	alignment_fs.close();

	while (getline(old_contig_fs, line)) {
		if (line.substr(0,1) == ">"){
			vector<string> tokens;
			tokenize(line, tokens, " ");
			string contig_id = tokens[0].substr(1, tokens[0].length()-1);
			if (std::find(contig_list.begin(), contig_list.end(), contig_id)!=contig_list.end()){
				new_contig_fs << line << '\n';
				getline(old_contig_fs, line);
				new_contig_fs << line << '\n';
			}
		}
	}
	old_contig_fs.close();
	new_contig_fs.close();
	return aligned_query_list;
}

// This is for the final round. The hit contigs are identified and put in the final hit_contigs.fasta file.
void GSQAligner::get_hit_contigs(const double& min_score, const double& min_coverage, const unsigned int& min_contig_lgth, const string& final_contigs_file, const string& hit_contig_file, const string& alignment_file, tuple_map& best_hits){
	double best_coverage = std::get<1>(best_hits["coverage"]);
	double final_high_coverage = 0.0;
	ifstream old_contig_fs(final_contigs_file.c_str());
	ifstream alignment_fs(alignment_file.c_str());
	ofstream new_contig_fs(hit_contig_file.c_str());
	string line;
	string query_sequence_line;
	char currentid[20];
	unsigned int contig_length;
	unsigned int query_length;
	double min_match_length;
	char query_type[10];
	vector<string> contig_list;
	logger->running("Finding the hit contigs");
	num_matches = 0;
	output_string = "<B>" + string_format("%-15s %-8s %-30s %-10s %-15s %-15s %-15s","Contig","Strand","Query","Score","Length","Coverage","G/P/C") + "</B>\n";
	output_string += "----------------------------------------------------------------------------------------------------------\n";
	while (getline(alignment_fs, line)) {
		// Identify the contig being checked
		if (line.substr(0,8) == "Sequence"){
			sscanf(line.c_str(),"Sequence %*d: %80[^,], from %*d to %d,%*s",currentid,&contig_length);
			logger->info("... checking contig:\t" + std::string(currentid) + "\tof length:\t" + int2str(contig_length));
			continue;
		}
		// Measure the length of the query
		if (line.substr(0,5) == "Query") {
			query_length = 0;
			sscanf(line.c_str(), "Query %s sequence %*s", query_type);
			while (getline(alignment_fs, query_sequence_line)) {
				if (query_sequence_line.substr(0,9) == "Predicted") {
					if (std::string(query_type) == "protein") {
						min_match_length = query_length * 3 * min_coverage;
					} else {
						min_match_length = query_length * min_coverage;
					}
					break;
				}
				query_length += count_letters(query_sequence_line);
			}
			continue;
		}
		// Parse the MATCH line
		if (line.substr(0,5) == "MATCH"){
			vector<string> tokens;
			tokenize(line, tokens, "\t");
			string contig_id = tokens[1].substr(0, tokens[1].length()-1);
			string strand = tokens[1].substr(tokens[1].length()-1,1);
			string query = tokens[2];
			string score = tokens[3];
			string length = tokens[4];
			string cov = tokens[5];
			string type = tokens[6];

			if (str2double(cov) > final_high_coverage) {
				final_high_coverage = str2double(cov);
			}

			// It annoys me that this can actually decide that a match with high coverage of the CONTIG meets requirements, but better than a false negative.
			if (str2double(score) > min_score && str2double(cov) > min_coverage && contig_length >= min_contig_lgth && str2double(length) >= min_match_length) {
				output_string += string_format("%-15s %-8s %-30s %-10s %-15s %-15s %-15s",contig_id.c_str(),strand.c_str(),query.c_str(),score.c_str(),length.c_str(),cov.c_str(),type.c_str()) + "\n";
				num_matches++;
				output_string += "\n";
				logger->info("   ... HIT MATCH found with coverage:\t" + cov + " " + type + "\tscore:\t" + score + "\tcontig length:\t" + int2str(contig_length));
				contig_list.push_back(contig_id);
			}
		}
	}

	if (best_coverage > final_high_coverage) {
		logger->warn("Contig with better coverage found in round " + int2str(std::get<0>(best_hits["coverage"])));
	}
	output_string += "\nLength: cumulative length of scored exons\nCov G/P/C: coverage of contig (G) or cDNA (C) or protein (P), whichever is highest";
	alignment_fs.close();

	while (getline(old_contig_fs, line)) {
		if (line.substr(0,1) == ">"){
			vector<string> tokens;
			tokenize(line, tokens, " ");
			string contig_id = tokens[0].substr(1, tokens[0].length()-1);
			if (std::find(contig_list.begin(), contig_list.end(), contig_id)!=contig_list.end()){
				new_contig_fs << line << '\n';
				getline(old_contig_fs, line);
				new_contig_fs << line << '\n';
			}
		}
	}
	old_contig_fs.close();
	new_contig_fs.close();
}

int GSQAligner::get_longest_match(int round, int k, const double& min_score, const unsigned int& query_contig_min, const string& alignment_file, tuple_map& best_hits){
	double best_coverage = std::get<1>(best_hits["coverage"]);
	int best_length = 0;
	ifstream alignment_fs(alignment_file.c_str());
	string line;
	unsigned int contig_length;
	vector<string> contig_list;
	while (getline(alignment_fs, line)) {
		// Get the length of the contig being checked.
		if (line.substr(0,8) == "Sequence"){
			sscanf(line.c_str(),"Sequence %*d: %*80[^,], from %*d to %d,%*s",&contig_length);
			continue;
		}
		// Parse the MATCH line
		if (line.substr(0,5) == "MATCH"){
			vector<string> tokens;
			tokenize(line, tokens, "\t");
			string contig_id = tokens[1].substr(0, tokens[1].length()-1);
			string strand = tokens[1].substr(tokens[1].length()-1,1);
			string query = tokens[2];
			string score = tokens[3];
			string length = tokens[4];
			string cov = tokens[5];
			string type = tokens[6];
			logger->debug("   ... " + int2str(k) + "-mer check MATCH found with coverage:\t" + cov + " " + type + "\tscore:\t" + score + "\tlength:\t" + length);
			// Keep track of better coverage if it is seen. This should correspond to the aligned length.
			if (type == "P" || type == "C"){
				if (str2double(score) > min_score && str2double(cov) > best_coverage) {
					get<0>(best_hits["coverage"]) = round;
					get<1>(best_hits["coverage"]) = str2double(cov);
					best_coverage = str2double(cov);
				}
			}
			// If this has a longer alignment length than seen before, update best_length.
			if (str2double(score) > min_score && contig_length >= query_contig_min && str2int(length) > best_length) {
				best_length = str2int(length);
			}
		}
	}
	alignment_fs.close();
	return best_length;
}

string GSQAligner::get_program_name(){
	return "GeneSeqer";
}
void GSQAligner::clean_files(const string& file){

}
