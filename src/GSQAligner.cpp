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

void GSQAligner::do_spliced_alignment(const string& genomic_file, const string& type, const string& query_file, const string& species, const Params& params, const string& output_file, const string& hit_contig_file){
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
	string cmd = "GeneSeqer -L " + genomic_file + " -" + type_str + " " + query_file + " -species " + species_str + " " + param_list + " -O " + output_file + " >> " + logger->get_log_file() + " 2>&1";
	logger->debug(cmd);
	run_shell_command(cmd);
	//get_aligned_contigs(genomic_file, hit_contig_file, output_file);
}

string_map GSQAligner::get_aligned_contigs(const double& min_score, const double& min_coverage, const int& min_contig_lgth, const string& all_contig_file, const string& hit_contig_file, const string& alignment_file){
	ifstream old_contig_fs(all_contig_file.c_str());
	ifstream alignment_fs(alignment_file.c_str());
	ofstream new_contig_fs(hit_contig_file.c_str());
	string line;
	vector<string> contig_list;
	logger->debug("Finding the aligned contigs");
	num_matches = 0;
	output_string = "<B>Contig\tStrand\tQuery\tScore\tLength\tCov\tG/P/C</B>\n";
	output_string = "<B>" + string_format("%-15s %-8s %-30s %-10s %-15s %-15s %-15s","Contig","Strand","Query","Score","Length","Coverage","G/P/C") + "</B>\n";
	output_string += "----------------------------------------------------------------------------------------------------------\n";
	while (getline(alignment_fs, line)) {
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
			if (type == "P" || type == "C"){
				if (str2double(cov) > min_coverage && str2double(score) > min_score)
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

string GSQAligner::get_program_name(){
	return "GeneSeqer";
}
void GSQAligner::clean_files(const string& file){

}
