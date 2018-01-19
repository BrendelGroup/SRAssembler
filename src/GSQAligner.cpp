/*
 * GSQAligner.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#include "GSQAligner.h"
#include <iostream>

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

void GSQAligner::do_spliced_alignment(const string& genomic_file, const string& type, const string& query_file, const string& species, const Params& params, const string& output_file){
	string param_list = "";
	for ( Params::const_iterator it = params.begin(); it != params.end(); ++it ){
		param_list += " -" + it->first + " " + it->second;
	}
	string type_str = "E";
	string species_str = species_names[species];
	if (species_str == "")
		species_str = "generic";
	if (type == "protein")
		type_str = "Q";
	string cmd = "GeneSeqer -L " + genomic_file + " -" + type_str + " " + query_file + " -species " + species_str + " " + param_list + " -o " + output_file + " >> " + logger->get_log_file() + " 2>&1";
	logger->debug(cmd);
	run_shell_command(cmd);
}

string_map GSQAligner::get_aligned_contigs(const double& min_score, const double& min_coverage, const unsigned int& min_contig_lgth, const string& all_contig_file, const string& hit_contig_file, const string& alignment_file){
	ifstream old_contig_fs(all_contig_file.c_str());
	ifstream alignment_fs(alignment_file.c_str());
	ofstream new_contig_fs(hit_contig_file.c_str());
	string line;
	char currentid[20];
	unsigned int contig_length;
	vector<string> contig_list;
	logger->debug("Finding the aligned contigs");
	num_matches = 0;
	//output_string = "<B>Contig\tStrand\tQuery\tScore\tLength\tCov\tG/P/C</B>\n";
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
			logger->debug("   ... MATCH found with coverage:\t" + cov + " " + type + "\tscore:\t" + score + "\tlength:\t" + int2str(contig_length));
			if (type == "P" || type == "C"){
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

void GSQAligner::get_hit_contigs(const double& min_score, const double& min_coverage, const unsigned int& min_contig_lgth, const string& final_contigs_file, const string& hit_contig_file, const string& alignment_file){
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
	logger->debug("Finding the aligned contigs");
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
			//cerr << "Query is " + line << endl;
			query_length = 0;
			sscanf(line.c_str(), "Query %s sequence %*s", query_type);
			while (getline(alignment_fs, query_sequence_line)) {
				if (query_sequence_line.substr(0,9) == "Predicted") {
					if (std::string(query_type) == "protein") {
						min_match_length = query_length * 3 * min_coverage;
					} else {
						min_match_length = query_length * min_coverage;
					}
					//cerr << "Query length is " + int2str(query_length) << endl;
					//cerr << "min_coverage is " << min_coverage << endl;
					//cerr << "Min match length is " + int2str(min_match_length) << endl;
					break;
				}
				//cerr << query_sequence_line << endl;
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

			if (str2double(score) > min_score && str2double(cov) > min_coverage && contig_length >= min_contig_lgth && str2double(length) >= min_match_length) {
				output_string += string_format("%-15s %-8s %-30s %-10s %-15s %-15s %-15s",contig_id.c_str(),strand.c_str(),query.c_str(),score.c_str(),length.c_str(),cov.c_str(),type.c_str()) + "\n";
				num_matches++;
				output_string += "\n";
				logger->info("   ... HIT found with coverage:\t" + cov + " " + type + "\tscore:\t" + score + "\tlength:\t" + int2str(contig_length));
				contig_list.push_back(contig_id);
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
}

string GSQAligner::get_program_name(){
	return "GeneSeqer";
}
void GSQAligner::clean_files(const string& file){

}
