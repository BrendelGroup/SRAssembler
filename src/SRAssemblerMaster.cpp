/*
 * SRAssemblerMaster.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: hchou
 */
#include "SRAssemblerMaster.h"


SRAssemblerMaster::SRAssemblerMaster() {
	// TODO Auto-generated constructor stub
}

int SRAssemblerMaster::init(int argc, char * argv[], int rank, int mpiSize) {
	string cmd;
	int ret = SRAssembler::init(argc, argv, rank, mpiSize);
	if (ret == -1)
		return ret;
	if (!get_aligner(Aligner::PROTEIN_ALIGNER)->is_available()) return -1;
	if (!get_aligner(Aligner::DNA_ALIGNER)->is_available()) return -1;
	if (!get_assembler()->is_available()) return -1;
	if (!get_spliced_aligner()->is_available()) return -1;
	this->start_round = (over_write)?1:get_start_round();
	if (this->start_round == 1)
		create_folders();
	cmd = "Command: ";
	for (int i=0; i<argc;i++){
		cmd.append(argv[i]).append(" ");
	}
	logger->info(cmd);
	output_header();
	output_libraries();
	get_query_list();
	load_saved_contigs();
	return 0;
}

void SRAssemblerMaster::get_query_list(){
	if (start_round > 1)
		load_query_list();
	else {
		ifstream query_file(this->query_file.c_str());
		string line;
		while (getline(query_file, line)){
			if (line.substr(0,1) == ">"){
				vector<string> tokens;
				tokenize(line.substr(1), tokens, " ");
				query_list.push_back(tokens[0]);
			}
		}
		query_file.close();
	}
}

void SRAssemblerMaster::output_header(){
	output_content += "<HTML>\n<HEAD>\n<TITLE>SRAssembler Output</TITLE>\n</HEAD>\n<BODY>\n<A NAME=top><A NAME=file1>\n<H3>Libraries summary</H3><HR>";
}

void SRAssemblerMaster::output_libraries(){
	if (libraries.size() == 1)
		output_content += "<PRE>\n\n" + int2str(libraries.size()) + " library is tested\n\n";
	else
		output_content += "<PRE>\n\n" + int2str(libraries.size()) + " libraries are tested\n\n";
	output_content += "<B>Library\t\tInsert size\tPaired\treads</B>\n";
	output_content += "------------------------------------------------------------------------------\n";
	output_content += "<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"22\">";
	logger->info("We got " + int2str(libraries.size()) + " libraries");
	for (unsigned int i=0;i<this->libraries.size();i++){
		logger->info("library " + int2str(i+1) + " :");
		logger->info("insert size: " + int2str(libraries[i].get_insert_size()));
		logger->info("left read: " + libraries[i].get_left_read());
		logger->info("right read: " + libraries[i].get_right_read());
		logger->info("reversed: " + int2str(libraries[i].get_reversed()));
		logger->info("Paired-end: " + int2str(libraries[i].get_paired_end()));
		//output_content += "<tr><td bgcolor=\"#000000\">" + int2str(i+1) + "</td><td bgcolor=\"#000000\">" + int2str(libraries[i].get_insert_size()) + "</td><td bgcolor=\"#000000\">" + libraries[i].get_left_read() + "</td><td bgcolor=\"#000000\">" + libraries[i].get_right_read() + "</td></tr>";
		output_content += "<tr><td bgcolor=\"#ffffff\"><pre>\n";
		output_content += int2str(i+1) + "\t\t" + int2str(libraries[i].get_insert_size()) + "\t\t" + int2str(libraries[i].get_paired_end()) + "\t" + libraries[i].get_left_read() + "," + libraries[i].get_right_read() + "\n";
		output_content += "</td></tr>";
	}
	output_content += "</table>";
	output_content += "Paired: 1 : paired-end, 0: single";
}
void SRAssemblerMaster::output_summary_header(){
	output_content += "</PRE><H3>Assembly summary</H3><HR>";
	summary_header = "<B>Round\t";
	for (int k=start_k;k<=end_k;k+=step_k) {
		summary_header += "k=" + int2str(k) + "\t";
	}
	summary_header += "</B>\n";
}
void SRAssemblerMaster::output_summary(int round){
	output_content += "<PRE>\n\nThe final contigs file(<a href=\"" + get_file_name(this->final_contig_file) + "\">" + get_file_name(this->final_contig_file) + "</a>) are generated from round " + int2str(round) + " with k=" + int2str(best_k) + "\n";
	output_content += "\n</PRE><H5>The best k and the number of matched reads:</H5></PRE>\n<PRE>";
	output_content += "<B>Round\tBest_k\tMatched_reads</B>\n";
	output_content += "----------------------------------------------------\n";
	output_content += "<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"22\">";
	output_content += "<tr><td bgcolor=\"#ffffff\"><pre>\n";
	output_content += summary_best;
	output_content += "</td></tr></table></PRE>";
	output_content += "\n\n<H5>N50:</H5>\n<PRE>";
	output_content += summary_header;
	output_content += "----------------------------------------------------\n";
	output_content += "<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"22\">";
	output_content += "<tr><td bgcolor=\"#ffffff\"><pre>\n";
	output_content += summary_n50;
	output_content += "</td></tr></table></PRE>";
	output_content += "\n\n<H5>N90:</H5>\n<PRE>";
	output_content += summary_header;
	output_content += "----------------------------------------------------\n";
	output_content += "<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"22\">";
	output_content += "<tr><td bgcolor=\"#ffffff\"><pre>\n";
	output_content += summary_n90;
	output_content += "</td></tr></table></PRE>";

	output_content += "\n\n<H5>Longest contig length:</H5>\n<PRE>";
	output_content += summary_header;
	output_content += "----------------------------------------------------\n";
	output_content += "<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"22\">";
	output_content += "<tr><td bgcolor=\"#ffffff\"><pre>\n";
	output_content += summary_max;
	output_content += "</td></tr></table></PRE>";

	output_content += "\n\n<H5>Total contigs:</H5>\n<PRE>";
	output_content += summary_header;
	output_content += "----------------------------------------------------\n";
	output_content += "<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"22\">";
	output_content += "<tr><td bgcolor=\"#ffffff\"><pre>\n";
	output_content += summary_total;
	output_content += "</td></tr></table>";
	output_content += "\nThe red number refers to \"Best k\"</PRE>";
}

void SRAssemblerMaster::output_spliced_alignment(){
	//output_content += summary_best + "\nN50:\n\n" + title + summary_n50 + "\n\nN90:\n\n" + title + summary_n90 + "\n\nLongest contig length:\n\n" + title + summary_max + "\n\nTotal contigs:\n\n" + title + summary_total + "\n\n";
	//output_content += "\n<Spliced alignment summary>\n";
	SplicedAligner* spliced_aligner = get_spliced_aligner();
	output_content += "</PRE><H3>Spliced alignment - " + spliced_aligner->get_program_name() + "</H3><HR>";
	output_content += "<PRE>\n";
	string spliced_filename = get_file_name(this->spliced_alignment_output_file);
	string hit_contigs_name = get_file_name(this->aligned_contig_file);
	int num_matches = spliced_aligner->get_match_num();
	if (num_matches == 0)
		output_content += "\nNo matches found with " + spliced_aligner->get_program_name() + "\n";
	else {
		output_content += "\nFound " + int2str(num_matches) + " matches\n";
		output_content += "The detailed " + spliced_aligner->get_program_name() + " alignment report is <a href=\"" + spliced_filename + "\">" + spliced_filename + "</a>\n";
		output_content += "The contigs aligned to query sequences only can be found in <a href=\"" + hit_contigs_name + "\">" + hit_contigs_name + "</a>\n";
		output_content += "\n\n" + spliced_aligner->get_output_summary() + "\n";
	}

}

void SRAssemblerMaster::output_gene_finding(){
	GeneFinder* gene_finder = get_gene_finder();
	output_content += "</PRE><H3>Gene Finding - " + gene_finder->get_program_name() + "</H3><HR>";
	output_content += "<PRE>\n";
	string gene_finding_filename = get_file_name(this->gene_finding_output_file);
	output_content += "The detailed " + gene_finder->get_program_name() + " report is <a href=\"" + gene_finding_filename + "\">" + gene_finding_filename + "</a>\n";

}

void SRAssemblerMaster::print_message(const std::string& msg){
	cout << msg << endl;

}

void SRAssemblerMaster::show_usage(){
	cout << usage <<endl;
}

void SRAssemblerMaster::do_preprocessing(){
	logger->info("Now pre-processing the reads files ...");
	string cmd;
	logger->debug("Checking the existence of split reads data ...");
	for (unsigned lib_index=0;lib_index<this->libraries.size();lib_index++) {
		Library* lib = &this->libraries[lib_index];
		lib->set_num_parts(1);
		// test if files are generated before
		if (file_exists(lib->get_split_file_name(1, LEFT_READ))){
			//long library_read_count = get_read_count(lib->get_left_read(), lib->get_format()) + get_read_count(lib->get_right_read(), lib->get_format());
			long split_read_count = count_preprocessed_reads(lib_index);
			split_read_count /= 2; // Don't count header lines
			logger->debug("split_read_count: " + int2str(split_read_count));
			lib->set_num_parts(get_file_count(data_dir + "/lib" + int2str(lib_index+1) + "/" + get_file_base_name(lib->get_left_read()) + "_*." + "fasta"));
			// Test if split reads have been indexed
			if (file_exists(lib->get_read_part_index_name(lib->get_num_parts(), LEFT_READ) + ".skp")){
				logger->info("Using previously split files for read library " + int2str(lib_index+1));
				broadcast_code(ACTION_TOTAL_PARTS, lib_index, lib->get_num_parts(), 0);
				//broadcast_code(ACTION_EXIT, 0, 0);
				continue;
			}
		}
		logger->info("Splitting read library " + int2str(lib_index+1) + " ...");
		cmd = "rm -f " + data_dir + "/lib" + int2str(lib_index+1) + "/" + get_file_base_name(lib->get_left_read()) + "*";      //delete old files
		logger->debug(cmd);
		run_shell_command(cmd);
		cmd = "rm -f " + data_dir + "/lib" + int2str(lib_index+1) + "/" + get_file_base_name(lib->get_right_read()) + "*";      //delete old files
		logger->debug(cmd);
		run_shell_command(cmd);
		// Why would you name a variable 'from'?
		int from;
		int part = 0;
		int code_value;
		mpi_code code;
		// File splitting is handled by the actual library. It does not seem to take file type into acount.
		if (lib->get_paired_end() && mpiSize > 2){
			send_code(1, ACTION_SPLIT, lib_index, 1, 0);
			system("sleep 0.5");
			send_code(2, ACTION_SPLIT, lib_index, 2, 0);
			mpi_receive(code_value, from);
			mpi_receive(code_value, from);
		}
		else {
			//'LEFT_READ' and 'RIGHT_READ' are constants
			lib->do_split_files(LEFT_READ, this->reads_per_file);
			if (lib->get_paired_end())
				lib->do_split_files(RIGHT_READ, this->reads_per_file);
		}
		// We need to count all files because we are not interleaving.
		//if (lib->get_paired_end()) {
			//lib->set_num_parts(get_file_count(lib->get_split_read_prefix(lib->get_left_read()) + "*") + get_file_count(lib->get_split_read_prefix(lib->get_right_read()) + "*"));
		//} else {
			//lib->set_num_parts(get_file_count(lib->get_split_read_prefix(lib->get_left_read()) + "*"));
		//}
		lib->set_num_parts(get_file_count(lib->get_split_read_prefix(lib->get_left_read()) + "*.fasta"));
		logger->info("We have " + int2str(lib->get_num_parts()) +" split files");

		broadcast_code(ACTION_TOTAL_PARTS, lib_index, lib->get_num_parts(), 0);
		int completed = 0;

		if (mpiSize == 1){
			for (part=1; part<=lib->get_num_parts(); part++)
				SRAssembler::preprocess_read_part(lib_index, part);
		} else {
			if (lib->get_num_parts() < mpiSize){
				for (part=1; part<=lib->get_num_parts(); part++){
					//system("sleep 0.5");
					send_code(part, ACTION_PRE_PROCESSING, lib_index, part, 0);
				}
				while(completed < lib->get_num_parts()){
					mpi_receive(code_value, from);
					completed++;
				}
			}
			else {
				for (part=1;part<mpiSize;part++){
					//system("sleep 0.5");
					send_code(part, ACTION_PRE_PROCESSING, lib_index, part, 0);
				}
				while (completed < lib->get_num_parts()){
					mpi_receive(code_value, from);
					code = get_mpi_code(code_value);
					completed++;
					if (part <= lib->get_num_parts()){
						//cout << "seneding to " << from << ". lib: " << lib_index << ". part:" << part << endl;
						send_code(from, ACTION_PRE_PROCESSING, lib_index, part, 0);
						part++;
					}
				}
			}
		}
	}
	logger->info("Preprocessing done.");
}

int SRAssemblerMaster::get_start_round(){
	int start_round = 1;
	if (file_exists(tmp_dir)){
		for (int i=this->num_rounds;i>1;i--){
			bool found_previous = true;
			int mpiSize = (this->mpiSize == 0)? 1: this->mpiSize;
			// Why is this a for loop?
			//for (int j=0;j<mpiSize;j++) {
				// if reads file and next round index file exist (means assembled), we can continue from here.
				for (unsigned int lib_idx=0;lib_idx<this->libraries.size();lib_idx++){
					Library lib = this->libraries[lib_idx];
					// Why do we keep re-assigning found_previous?
					//found_previous = file_exists(lib.get_matched_left_read_filename(i));
					//if (lib.get_paired_end())
						//found_previous = file_exists(lib.get_matched_right_read_filename(i));
					found_previous = file_exists(get_query_fasta_file_name(i+1));
				}
			//}
			if (found_previous) {
				start_round = i+1;
				if (start_round > this->num_rounds)
					return start_round;
				logger->info("Previous results found. SRAssembler starts from round " + int2str(start_round));
				//clean the temp results if it is not complete.
				run_shell_command("rm -f " + tmp_dir + "/matched_reads_left_" + "r" + int2str(start_round) + "*");
				run_shell_command("rm -f " + tmp_dir + "/matched_reads_right_" + "r" + int2str(start_round) + "*");
				for (unsigned int lib_idx=0;lib_idx<this->libraries.size();lib_idx++){
					Library lib = this->libraries[lib_idx];
					run_shell_command("cp " + lib.get_matched_left_read_filename(i) + " " + lib.get_matched_left_read_filename());
					if(lib.get_paired_end()){
						run_shell_command("cp " + lib.get_matched_right_read_filename(i) + " " + lib.get_matched_right_read_filename());
					}
				}
				if (mpiSize > 1){
					for (int i=1;i<mpiSize;i++)
						 send_code(i, ACTION_LOAD_PREVIOUS, start_round - 1, 0, 0);
					int completed = 0;
					int code_value = 0;
					int from = 0;
					while(completed < mpiSize - 1){
						  mpi_receive(code_value, from);
						  completed++;
					}
				}
				else
					load_mapped_reads(start_round - 1);
				load_long_contigs();
				break;
			}
		}
	}
	return start_round;
}

void SRAssemblerMaster::do_walking(){
	if (preprocessing_only) {
		logger->info("Do pre-processing of reads only. The chromosome walking is skipped.");
		broadcast_code(ACTION_EXIT, 0, 0, 0);
		return;
	}
	logger->info("Start chromosome walking ...");
	logger->info("Total processors: " + int2str(mpiSize));
	int from;
	int read_part = 0;
	int code_value;
	mpi_code code;
	int round = this->start_round;
	if (round > this->num_rounds){
		logger->info("The previous assembly has been completed");
		broadcast_code(ACTION_EXIT, 0, 0, 0);
		return;
	}
	output_summary_header();
	while(true){
		logger->info("Starting round " + int2str(round) + " ...");
		int new_reads_count = 0;

		// Index the contigs in round 1 for the dnavsprot vmatch step
		if (round == 1)
			create_index(round);

		// For each library
		for (unsigned lib_idx=0;lib_idx<this->libraries.size();lib_idx++){
			int completed = 0;
			Library lib = this->libraries[lib_idx];
			// If not parallelized, start a new alignment every 1/2 second?
			if (mpiSize == 1){
				for (read_part=1; read_part<=lib.get_num_parts(); read_part++){
					sleep(0.5);
					new_reads_count += do_alignment(round, lib_idx, read_part);
				}
			} else {
				// If there are more split read files than processors
				if (lib.get_num_parts() < mpiSize){
					for (read_part=1; read_part<=lib.get_num_parts(); read_part++){
						sleep(0.5);
						send_code(read_part, ACTION_ALIGNMENT, round, read_part, lib_idx);
					}
					while(completed < lib.get_num_parts()){
						mpi_receive(code_value, from);
						code = get_mpi_code(code_value);
						int found_new_reads = code.value1;
						new_reads_count += found_new_reads;
						completed++;
					}
				// If there are fewer split read files than processors
				} else {
					for (read_part=1;read_part<mpiSize;read_part++){
						sleep(0.5);
						send_code(read_part, ACTION_ALIGNMENT, round, read_part, lib_idx);
					}
					while (completed < lib.get_num_parts()){
						mpi_receive(code_value, from);
						code = get_mpi_code(code_value);
						int found_new_reads = code.value1;
						int file_idx = code.value2;
						new_reads_count += found_new_reads;
						completed++;
						int next_file_idx = file_idx + mpiSize - 1;
						if (next_file_idx <= lib.get_num_parts())
							send_code(from, ACTION_ALIGNMENT, round, next_file_idx, lib_idx);
					}
				}
			}
		}

		if (new_reads_count == 0) {
			logger->info("The walking is terminated: No new reads found.");
			break;
		}
		merge_mapped_files(round);
		int read_count = get_total_read_count(round);
		logger->info("Found matched reads: " + int2str(read_count));
		if (assembly_round <= round){
			int longest_contig = do_assembly(round);
			summary_best += int2str(read_count) + "\n";
			bool no_reads = true;
			for (unsigned int lib_idx=0; lib_idx < this->libraries.size(); lib_idx++){
				if (get_file_size(libraries[lib_idx].get_matched_left_read_filename()) > 0) {
					no_reads = false;
					break;
				}
			}
			// if no reads found, stop
			if (no_reads){
				logger->info("The walking is terminated: No new reads found after removing reads associated with the assembled contigs.");
				break;
			}
			// if no contigs found, stop
			if (get_file_size(get_contig_file_name(round)) == 0) {
				logger->info("The walking is terminated: No contigs found.");
				break;
			}
			logger->info("Round " + int2str(round) + " is done!");
			if (round == num_rounds) {
				logger->info("The walking is terminated: The maximum round (" + int2str(num_rounds) + ") has been reached.");
				break;
			}
			if (longest_contig == 0) {
				assembly_round = round + 1;
				//RM HERE
				clean_tmp_files(round-1);
				round++;
				continue;
			}
			//if reach the max round, stop

			bool cleaned = false;
			//do spliced alignment and remove the query sequences already assembled
			if (round > 1 && check_gene_assembled){
				string_map query_map = do_spliced_alignment(round);
				//string_map query_map = this->get_spliced_aligner()->get_aligned_query_list();
				vector<string> contig_list;
				BOOST_FOREACH(string_map::value_type item, query_map) {
					contig_list.push_back(item.second);
					for (int i=this->query_list.size()-1;i>=0;i--){
						  if (query_list[i] == item.first){
							  query_list.erase(query_list.begin()+i);
						}
					}
				}
				save_query_list();
				if (query_list.size() == 0){
					logger->info("The walking is terminated: All homologous sequences have been assembled.");
					break;
				}
				if (contig_list.size() > 0) {
					remove_hit_contigs(contig_list, round);
					if (get_file_size(get_contig_file_name(round)) == 0){
						logger->info("The walking is terminated: All contigs have enough coverage and score.");
						break;
					}
					clean_unmapped_reads(round);
					cleaned = true;
				}
			}
			if (clean_round > 2 && !cleaned)
				if (round % clean_round == 0)
					clean_unmapped_reads(round);
		} else {
			logger->info("Round " + int2str(round) + " is done!");
			if (round == num_rounds) {
				logger->info("The walking is terminated: The maximum round (" + int2str(num_rounds) + ") has been reached.");
				break;
			}
		}
		if (round > 1) {
			//RM HERE
			clean_tmp_files(round-1);
		}
		round++;
	}
	if (round > 1) {
		//RM HERE
		clean_tmp_files(round-1);
	}
	// notify all slaves to stop listening
	broadcast_code(ACTION_EXIT, 0, 0, 0);
	//not assembled yet, do assembling
	if (assembly_round > round){
		do_assembly(round);
	}
	// if the final contig size is 0, then report the previous round
	while (round > 1) {
		logger->info("Checking the final contigs assembled in round " + int2str(round) + " ...");
		prepare_final_contig_file(round);
		if (get_file_size(final_contig_file) == 0){
			logger->info("... no contigs found in round " + int2str(round));
			round--;
		}
		else {
			break;
		}
	}
	if (get_file_size(final_contig_file) == 0) {
		return;
	}
	do_spliced_alignment();
	do_gene_finding();
	//RM HERE
	this->get_spliced_aligner()->clean_files(this->final_contig_file);
	output_summary(round);
	output_spliced_alignment();
	if (file_exists(this->gene_finding_output_file)) {
		output_gene_finding();
	}
	ofstream outFile(summary_file.c_str());
	outFile << output_content << endl;
	outFile.close();
	//RM HERE
	run_shell_command("rm -rf " + query_file + ".*");
}

void SRAssemblerMaster::clean_tmp_files(int round){
	if (round == 0) return;
	//else: remove data of previous round
	string cmd = "rm -f " + tmp_dir + "/vmatch_" + "r" + int2str(round) + "_*";
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = "rm -f " + tmp_dir + "/matched_reads_left_" + "r" + int2str(round) + "_*";
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = "rm -f " + tmp_dir + "/matched_reads_right_" + "r" + int2str(round) + "_*";
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = "rm -f " + tmp_dir + "/matched_reads_" + "r" + int2str(round) + "_*";
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = "rm -f " + tmp_dir + "/query-vs-contig_" + "r" + int2str(round) + ".*";
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = "rm -f " + tmp_dir + "/hit_contigs_" +"r" + int2str(round) + ".*";
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = "rm -f " + tmp_dir + "/long_contig_candidate_" + "r" + int2str(round) + ".*";
	logger->debug(cmd);
	run_shell_command(cmd);
}

void SRAssemblerMaster::save_query_list(){
	string fn = tmp_dir + "/query_list";
	//RM HERE
	run_shell_command("rm -rf " + fn);
	if (query_list.size() > 0){
		ofstream fs(fn.c_str());
		for (unsigned i=0;i<query_list.size();i++)
			fs << query_list[i] << endl;
		fs.close();
	}
}

void SRAssemblerMaster::load_query_list(){
	string fn = tmp_dir + "/query_list";
	if (file_exists(fn)){
		ifstream fs(fn.c_str());
		query_list.clear();
		string line;
		while (getline(fs, line))
			query_list.push_back(line);
		fs.close();
	}
}

void SRAssemblerMaster::load_saved_contigs(){
	string fn = get_saved_contig_file_name();
	this->contig_number = 1;
	if (start_round == 1)
		//RM HERE
		run_shell_command("rm -rf " + fn);
	else {
		if (file_exists(fn)){
			ifstream fs(fn.c_str());
			string line;
			while (getline(fs, line)){
				if (line.substr(0,1)==">")
					this->contig_number++;
			}
			fs.close();
		}
	}
}

string SRAssemblerMaster::get_saved_contig_file_name(){
	return tmp_dir + "/saved_contig.fasta";
}

int SRAssemblerMaster::do_assembly(int round) {
	logger->info("Doing assembly, round: " + int2str(round));
	int best_k = 0;
	unsigned int max_longest_contig = 0;
	int total_k = (end_k-start_k)/step_k + 1;
	int from;
	int i = 0;
	int completed = 0;
	int code_value;
	if (mpiSize == 1){
		for (i=1; i<=total_k; i++)
			SRAssembler::do_assembly(round, start_k + (i-1)*step_k);
	} else {
		if (total_k < mpiSize-1){
			for (i=1; i<=total_k; i++){
				sleep(0.5);
				send_code(i, ACTION_ASSEMBLY, round, start_k + (i-1)*step_k, 0);
			}
			while(completed < total_k){
				mpi_receive(code_value, from);
				completed++;
			}
		}
		else {
			for (i=1;i<mpiSize;i++){
				sleep(0.5);
				send_code(i, ACTION_ASSEMBLY, round, start_k + (i-1)*step_k, 0);
			}
			while(completed < total_k){
				mpi_receive(code_value, from);
				completed++;
				if (i <= total_k) {
					send_code(from, ACTION_ASSEMBLY, round, start_k + (i-1)*step_k, 0);
					i++;
				}
			}
		}
	}
	summary_n50 += int2str(round) + "\t";
	summary_n90 += int2str(round) + "\t";
	summary_max += int2str(round) + "\t";
	summary_total += int2str(round) + "\t";
	summary_best += int2str(round) + "\t";
	vector<Assembly_stats> stats;
	for (int k=start_k;k<=end_k;k+=step_k) {
		Assembly_stats kstats = get_assembly_stats(round, k);
		if (kstats.longest_contig > max_longest_contig){
			best_k = k;
			max_longest_contig = kstats.longest_contig;
		}
		stats.push_back(kstats);
	}
	i = 0;
	for (int k=start_k;k<=end_k;k+=step_k) {
		Assembly_stats kstats = stats[i];
		if (kstats.longest_contig > 0) {
			logger->info("The longest contig assembled in round\t" + int2str(round) + " is of length\t" + int2str(kstats.longest_contig) + " with k =\t" + int2str(k));
		}
		else {
			logger->info("No contig of the specified minimum length has been assembled by round\t" + int2str(round) + " with k =\t" + int2str(k));
		}
		string ffont = "";
		string rfont = "";
		if (k == best_k){
			ffont = "<B><font color=\"red\">";
			rfont = "</font></B>";
		}
		summary_n50 += ffont + int2str(kstats.n50) + rfont + "\t";
		summary_n90 += ffont + int2str(kstats.n90) + rfont + "\t";
		summary_max += ffont + int2str(kstats.longest_contig) + rfont + "\t";
		summary_total += ffont + int2str(kstats.total_contig) + rfont + "\t";
		i++;
	}
	summary_n50 += "\n";
	summary_n90 += "\n";
	summary_max += "\n";
	summary_total += "\n";
	if (best_k > 0) {
		logger->info("The best k-value (corresponding to the longest assembled contig) in round\t" + int2str(round) + " is k =\t" + int2str(best_k));
	}
	else {
		logger->info("No contig of the specified minimum length has been assembled by round\t" + int2str(round));
	}
	summary_best += int2str(best_k) + "\t";
	if (best_k > 0) {
		process_long_contigs(round, best_k);
	}
	//RM HERE
	get_assembler()->clean_files(tmp_dir);
	return max_longest_contig;
}

void SRAssemblerMaster::load_long_contigs() {
	string long_contig_file = tmp_dir + "/long_contig.fasta";
	string saved_contig_file_name = get_saved_contig_file_name();
	ifstream in_contig(long_contig_file.c_str());
	ofstream saved_contig_file;
	saved_contig_file.open(saved_contig_file_name.c_str(), fstream::app);
	string contig_fasta;
	string line;
	while (getline(in_contig, line)){
		contig_fasta = line + "\n";
		getline(in_contig, line);
		contig_fasta = contig_fasta + line + "\n";
		//long_contigs.push_back(contig_fasta);
		saved_contig_file << contig_fasta;
	}
	saved_contig_file.close();
}

void SRAssemblerMaster::process_long_contigs(int round, int k) {
	string long_contig_candidate_file = tmp_dir + "/long_contig_candidate_" + "r" + int2str(round-1)+ ".fasta";
	string long_contig_candidate_next_file = tmp_dir + "/long_contig_candidate_" + "r" + int2str(round)+ ".fasta";
	string long_contig_file = tmp_dir + "/long_contig_original.fasta";
	string long_contig_trimmed_file = tmp_dir + "/long_contig.fasta";
	string saved_contig_file_name = get_saved_contig_file_name();
	boost::unordered_set<string> candidate_ids;
	boost::unordered_set<string> long_contig_ids;
	ofstream out_long_contig_trimmed;
	out_long_contig_trimmed.open((long_contig_trimmed_file.c_str()), ofstream::app);
	ofstream saved_contig_file;
	saved_contig_file.open(saved_contig_file_name.c_str(), fstream::app);
	string line;
	if (file_exists(long_contig_candidate_file) && get_file_size(long_contig_candidate_file) > 0){
		this->get_aligner(round)->align_long_contigs(long_contig_candidate_file, tmp_dir, this->get_assembly_file_name(round, k), this->max_contig_lgth, candidate_ids, long_contig_ids);
		//add candidate contigs to accepted long contigs
		ifstream candidate_file(long_contig_candidate_file.c_str());
		while (getline(candidate_file, line)){
			vector<string> tokens;
			tokenize(line.substr(1), tokens, " ");
			if (candidate_ids.find(tokens[0]) != candidate_ids.end()){
				string contig_fasta = ">contig" + int2str(contig_number++) + " length " + int2str(this->max_contig_lgth) + " ";
				for (unsigned int i=3;i<tokens.size();i++){
					contig_fasta = contig_fasta + tokens[i] + " ";
				}
				contig_fasta = contig_fasta + "\n";
				getline(candidate_file, line);
				contig_fasta = contig_fasta + line + "\n";
				saved_contig_file << contig_fasta;
				out_long_contig_trimmed << contig_fasta;
			}
		}
		candidate_file.close();
	}

	ifstream in_contig(get_assembly_file_name(round, k).c_str());
	ofstream out_contig(get_query_fasta_file_name(round+1).c_str());
	ofstream out_candidate_contig(long_contig_candidate_next_file.c_str());
	ofstream out_long_contig;
	out_long_contig.open((long_contig_file.c_str()));
	string header = "";
	string seq = "";
	while (getline(in_contig, line)){
		if (line.substr(0,1) == ">"){
			if (header.length() > 0) {
				vector<string> tokens;
				tokenize(header.substr(1), tokens, " ");
				if (long_contig_ids.find(tokens[0]) != long_contig_ids.end())
					out_long_contig << header << endl << seq <<endl;
				else {
					if (seq.length() > this->ini_contig_size) {
						out_contig << header << endl << seq << endl;
					}
					if (seq.length() > this->max_contig_lgth) {
						out_candidate_contig << header << endl << seq.substr((seq.length() - max_contig_lgth) / 2,max_contig_lgth) << endl;
					}
				}
			}
			header = line;
			seq = "";
		}
		else
			seq.append(line);
	}
	vector<string> tokens;
	tokenize(header.substr(1), tokens, " ");
	if (long_contig_ids.find(tokens[0]) != long_contig_ids.end())
		out_long_contig << header << endl << seq <<endl;
	else {
		if (seq.length() > this->ini_contig_size)
			out_contig << header << endl << seq << endl;
		if (seq.length() > this->max_contig_lgth) {
			out_candidate_contig << header << endl << seq.substr((seq.length() - max_contig_lgth) / 2,max_contig_lgth) << endl;
		}
	}
	in_contig.close();
	out_contig.close();
	out_candidate_contig.close();
	out_long_contig.close();

	if (long_contig_ids.size() > 0){
		logger->info("Processing contigs longer than the specified maximal contig size " + int2str(this->max_contig_lgth) + " ...");
		string reads_on_contigs = tmp_dir + "/matched_reads.sam";
		//remove associated reads
		string cmd = "bowtie-build " + long_contig_file + " " + tmp_dir + "/long_contig >> " + logger->get_log_file();
		logger->debug(cmd);
		run_shell_command(cmd);
		for (unsigned int lib_idx=0; lib_idx < this->libraries.size(); lib_idx++){
			Library lib = this->libraries[lib_idx];
			//string type_option = (lib.get_format() == FORMAT_FASTQ)? "-q" : "-f";
			string type_option = "-f";
			if (lib.get_paired_end())
				cmd = "bowtie " + type_option + " -v 2 --suppress 2,3,4,5,6,7,8 " + tmp_dir + "/long_contig " + lib.get_matched_left_read_filename() + "," + lib.get_matched_right_read_filename() + " " + reads_on_contigs + " >> " + logger->get_log_file() + " 2>&1";
			else
				cmd = "bowtie "  + type_option + " -v 2 --suppress 2,3,4,5,6,7,8 " + tmp_dir + "/long_contig " + lib.get_matched_left_read_filename() + " " + reads_on_contigs + " >> " + logger->get_log_file() + " 2>&1";
			logger->debug(cmd);
			run_shell_command(cmd);

			//read sam file
			boost::unordered_set<string> mapped_reads;
			ifstream reads_file(reads_on_contigs.c_str());
			while (getline(reads_file, line)) {
				line.erase(line.find_last_not_of(" \n\r\t")+1);
				mapped_reads.insert(line);
			}
			reads_file.close();

			string left_matched = lib.get_matched_left_read_filename();
			string right_matched = lib.get_matched_right_read_filename();
			string left_processed = tmp_dir + "/matched_reads_processed_left." + lib.get_file_extension();
			string right_processed = tmp_dir + "/matched_reads_processed_right." + lib.get_file_extension();

			ifstream left_matched_stream(left_matched.c_str());
			ofstream left_processed_stream(left_processed.c_str());
			ifstream right_matched_stream;
			ofstream right_processed_stream;
			if (lib.get_paired_end()) {
				right_matched_stream.open(right_matched.c_str(), ios_base::in);
				right_processed_stream.open(right_processed.c_str(), ios_base::out);
			}
			string left_header = "";
			string right_header = "";
			string left_seq = "";
			string right_seq = "";
			string left_qual = "";
			string right_qual = "";
			string plus;
			while (getline(left_matched_stream, left_header)){
				string left_seq_id = "";
				string right_seq_id = "";
				string lead_chr = ">";
				if (left_header.substr(0,1) == lead_chr){
					unsigned int pos = left_header.find_first_of(" ");
					if (pos ==string::npos)
						left_seq_id = left_header;
					else
						left_seq_id = left_header.substr(1, pos-1);
					left_seq_id.erase(left_seq_id.find_last_not_of("\n\r\t")+1);
					getline(left_matched_stream, left_seq);
					if (lib.get_paired_end()){
						getline(right_matched_stream, right_header);
						getline(right_matched_stream, right_seq);
						pos = right_header.find_first_of(" ");
						if (pos ==string::npos)
							right_seq_id = right_header;
						else
							right_seq_id = right_header.substr(1, pos-1);
						right_seq_id.erase(right_seq_id.find_last_not_of(" \n\r\t")+1);
					}
					if (mapped_reads.find(left_seq_id) == mapped_reads.end() && (lib.get_paired_end() && mapped_reads.find(right_seq_id) == mapped_reads.end())){
						left_processed_stream << left_header << endl << left_seq << endl;
						if (lib.get_paired_end()){
							right_processed_stream << right_header << endl << right_seq << endl;
						}
					}
				}
			}
			left_matched_stream.close();
			right_matched_stream.close();
			left_processed_stream.close();
			right_processed_stream.close();
			cmd = "cp " + left_processed + " " + left_matched;
			logger->debug(cmd);
			run_shell_command(cmd);
			cmd = "cp " + left_processed + " " + lib.get_matched_left_read_filename(round);
			logger->debug(cmd);
			run_shell_command(cmd);
			if (lib.get_paired_end()){
				cmd = "cp " + right_processed + " " + right_matched;
				logger->debug(cmd);
				run_shell_command(cmd);
				cmd = "cp " + right_processed + " " + lib.get_matched_right_read_filename(round);
				logger->debug(cmd);
				run_shell_command(cmd);
			}
		}
	}
	// For right now masking is done by NCBI's dustmasker.
	//TODO This should be more modular
	//TODO Masking is happening before cleaning in cleaning rounds.
	string cmd;
	string contig_file;
	contig_file = get_contig_file_name(round);
	string masked_file = contig_file + ".masked";
	//cmd = "printf '\e[38;5;002mPROCESS_LONG_CONTIGS\e[0m\n'; dustmasker -in " + contig_file + " -outfmt fasta -out - | sed '/^[^>]/s/[a-z]/N/g' > " + masked_file;
	cmd = "dustmasker -in " + contig_file + " -outfmt fasta -out - | sed '/^[^>]/s/[a-z]/N/g' > " + masked_file;
	run_shell_command(cmd);
}

void SRAssemblerMaster::remove_hit_contigs(vector<string> &contig_list, int round){
	logger->debug("remove hit contigs of round " + int2str(round));
	string contig_file = get_contig_file_name(round);
	string tmp_file = tmp_dir + "/contig_tmp_" + "r" + int2str(round) + ".fasta";
	string saved_contig_file_name = get_saved_contig_file_name();
	string line;
	string header = "";
	string seq = "";
	ifstream contig_file_stream(contig_file.c_str());
	ofstream tmp_file_stream(tmp_file.c_str());
	ofstream saved_contig_file;
	saved_contig_file.open(saved_contig_file_name.c_str(), fstream::app);
	fstream filestr;
	bool hit_seq = false;
	while (getline(contig_file_stream, line)){
		if (line.substr(0,1) == ">"){
			if (seq != "") {
				if (hit_seq)
					saved_contig_file << ">contig" + int2str(contig_number++) + header.substr(header.find_first_of(" ")) + "\n" + seq + "\n";
				else
					tmp_file_stream << ">" << header << endl << seq << endl;
			}
			header = line.substr(1);
			string contig_id = header.substr(0, header.find_first_of(" "));
			hit_seq = (find(contig_list.begin(), contig_list.end(), contig_id) != contig_list.end());
			seq = "";
		}
		else
			seq.append(line);
	}
	if (seq != "") {
		if (hit_seq)
			saved_contig_file << ">contig" + int2str(contig_number++) + header.substr(header.find_first_of(" ")) + "\n" + seq + "\n";
		else
			tmp_file_stream << ">" << header << endl << seq << endl;
	}
	contig_file_stream.close();
	tmp_file_stream.close();
	saved_contig_file.close();
	string cmd = "cp " + tmp_file + " " + contig_file;
	run_shell_command(cmd);
	//RM HERE
	cmd = "rm " + tmp_file;
	run_shell_command(cmd);
}

void SRAssemblerMaster::remove_no_hit_contigs(unordered_set<string> &hit_list, int round){
	logger->debug("remove contigs without hits against query sequences in round " + int2str(round));
	string contig_file = get_contig_file_name(round);
	string tmp_file = tmp_dir + "/contig_tmp_" + "r" + int2str(round) + ".fasta";
	string line;
	string header = "";
	string seq = "";
	ifstream contig_file_stream(contig_file.c_str());
	ofstream tmp_file_stream(tmp_file.c_str());
	bool hit_seq = false;
	int contig_number = 0;
	while (getline(contig_file_stream, line)){
		if (line.substr(0,1) == ">"){
			string contig_id = int2str(contig_number);
			if (hit_seq) {
			 	tmp_file_stream << ">" << header << endl << seq << endl;
			}
			header = line.substr(1);
			hit_seq = (find(hit_list.begin(), hit_list.end(), contig_id) != hit_list.end());
			seq = "";
			contig_number++;
		}
		else
			seq.append(line);
	}
	if (hit_seq)
		tmp_file_stream << ">" << header << endl << seq << endl;
	contig_file_stream.close();
	tmp_file_stream.close();
	string cmd = "cp " + tmp_file + " " + contig_file;
	run_shell_command(cmd);
	//RM HERE
	cmd = "rm " + tmp_file;
	run_shell_command(cmd);
}

void SRAssemblerMaster::prepare_final_contig_file(int round){
	logger->debug("prepare final contig file of round " + int2str(round));
	string contig_file = get_contig_file_name(round);
	string saved_contig_file_name = get_saved_contig_file_name();
	//contig_file = this->final_scaf_file;
	//final_long_contig_file = results_dir + "/contigs_long.fasta";
	ifstream last_round_contig(contig_file.c_str());
	ifstream saved_contig(saved_contig_file_name.c_str());
	ofstream final_contig(final_contig_file.c_str());
	string line;
	while (getline(saved_contig, line))
		final_contig << line << endl;
	saved_contig.close();
	string header = "";
	string seq = "";
	while (getline(last_round_contig, line)){
		if (line.substr(0,1) == ">"){
			if (header.length() > 0) {
				vector<string> header_tokens;
				tokenize(header.substr(1), header_tokens, " ");
				header = ">contig" + int2str(contig_number++) + " ";
				for (unsigned int i=1;i<header_tokens.size();i++){
					header = header + header_tokens[i] + " ";
				}
				final_contig << header << endl << seq << endl;
			}
			header = line;
			seq = "";
		}
		else
			seq.append(line);
	}
	if (header.length() > 0) {
		vector<string> header_tokens;
		tokenize(header.substr(1), header_tokens, " ");
		string header = ">contig" + int2str(contig_number++) + " ";
		for (unsigned int i=1;i<header_tokens.size();i++){
			header = header + header_tokens[i] + " ";
		}
		final_contig << header << endl << seq << endl;
	}
}

void SRAssemblerMaster::create_folders(){
	string cmd;
	if (file_exists(tmp_dir)){
		cmd = "rm -rf " + tmp_dir;
		run_shell_command(cmd);
	}
	if (file_exists(results_dir)){
		cmd = "rm -rf " + results_dir;
		run_shell_command(cmd);
	}
	if (file_exists(intermediate_dir)){
		cmd = "rm -rf " + intermediate_dir;
		run_shell_command(cmd);
	}
	for (unsigned i=0;i<this->libraries.size();i++){
		string dir = data_dir + "/lib" + int2str(i+1);
		cmd = "mkdir -p " + dir;
		run_shell_command(cmd);
	}
	cmd = "mkdir " + results_dir;
	run_shell_command(cmd);
	cmd = "mkdir " + intermediate_dir;
	run_shell_command(cmd);
	cmd = "mkdir " + tmp_dir;
	run_shell_command(cmd);
}

void SRAssemblerMaster::remove_no_hit_contigs(int round){
	logger->info("Removing contigs without hits ...");
	string contig_file = get_contig_file_name(round);
	Aligner* aligner = get_aligner(round);
	aligner->create_index(tmp_dir + "/qindex", type, query_file);
	string program_name = aligner->get_program_name();
	program_name += "_" + get_type(1) + "_vs_contig";
	Params params = read_param_file(program_name);
	string out_file = tmp_dir + "/query_contig.aln";
	aligner->do_alignment(tmp_dir + "/qindex", type, get_match_length(1), get_mismatch_allowed(1), contig_file, params, out_file);
	unordered_set<string> hits = aligner->get_hit_list(out_file);
	remove_no_hit_contigs(hits, round);
	string cmd = "rm -f " + tmp_dir + "/qindex*";
	logger->debug(cmd);
	run_shell_command(cmd);
}

void SRAssemblerMaster::clean_unmapped_reads(int round){
//TODO improve the speed of the read finding. May be able to parse bowtie output directly.
	remove_no_hit_contigs(round);
	string contig_file = get_contig_file_name(round);
	//remove associated reads
	string index_name = tmp_dir + "/contig_index_" + int2str(round);
	string cmd = "bowtie-build " + contig_file + " " + index_name + " >> " + logger->get_log_file();
	logger->debug(cmd);
	run_shell_command(cmd);
	for (unsigned int lib_idx=0; lib_idx < this->libraries.size(); lib_idx++){
		Library lib = this->libraries[lib_idx];
		//string type_option = (lib.get_format() == FORMAT_FASTQ)? "-q" : "-f";
		string type_option = "-f";
		string reads_on_contigs = tmp_dir + "/matched_reads" + int2str(lib_idx+1) + ".sam";
		if (lib.get_paired_end())
			cmd = "bowtie " + type_option + " -v 2 --suppress 2,3,4,5,6,7,8 " + index_name + " " + lib.get_matched_left_read_filename() + "," + lib.get_matched_right_read_filename() + " " + reads_on_contigs + " >> " + logger->get_log_file() + " 2>&1";
		else
			cmd = "bowtie "  + type_option + " -v 2 --suppress 2,3,4,5,6,7,8 " + index_name + " " + lib.get_matched_left_read_filename() + " " + reads_on_contigs + " >> " + logger->get_log_file() + " 2>&1";
		logger->debug(cmd);
		run_shell_command(cmd);
		//read sam file
		boost::unordered_set<string> mapped_reads;
		ifstream reads_file(reads_on_contigs.c_str());
		string line;
		while (getline(reads_file, line)) {
			line.erase(line.find_last_not_of(" \n\r\t")+1);
			mapped_reads.insert(line);
		}
		reads_file.close();
		string left_matched = lib.get_matched_left_read_filename();
		string right_matched = lib.get_matched_right_read_filename();
		string left_processed = tmp_dir + "/matched_reads_processed_left." + "fasta";
		string right_processed = tmp_dir + "/matched_reads_processed_right." + "fasta";

		ifstream left_matched_stream(left_matched.c_str());
		ofstream left_processed_stream(left_processed.c_str());
		ifstream right_matched_stream;
		ofstream right_processed_stream;
		if (lib.get_paired_end()) {
			right_matched_stream.open(right_matched.c_str(), ios_base::in);
			right_processed_stream.open(right_processed.c_str(), ios_base::out);
		}
		string left_header = "";
		string right_header = "";
		string left_seq = "";
		string right_seq = "";
		string left_qual = "";
		string right_qual = "";
		string plus;
		while (getline(left_matched_stream, left_header)){
			string left_seq_id = "";
			string right_seq_id = "";
			if (left_header.length()==0)
				continue;
			string lead_chr = ">";
			if (left_header.substr(0,1) == lead_chr){
				size_t pos = left_header.find_first_of(" ");
				if (pos ==string::npos)
					left_seq_id = left_header.substr(1);
				else
					left_seq_id = left_header.substr(1, pos-1);
				left_seq_id.erase(left_seq_id.find_last_not_of("\n\r\t")+1);
				getline(left_matched_stream, left_seq);
				if (lib.get_paired_end()){
					getline(right_matched_stream, right_header);
					getline(right_matched_stream, right_seq);
					pos = right_header.find_first_of(" ");
					if (pos ==string::npos)
						right_seq_id = right_header.substr(1);
					else
						right_seq_id = right_header.substr(1, pos-1);
					right_seq_id.erase(right_seq_id.find_last_not_of("\n\r\t")+1);
				}
				if (mapped_reads.find(left_seq_id) != mapped_reads.end() && (lib.get_paired_end() || mapped_reads.find(right_seq_id) != mapped_reads.end())){
					left_processed_stream << left_header << endl << left_seq << endl;
					if (lib.get_paired_end()){
						right_processed_stream << right_header << endl << right_seq << endl;
					}
				}
			}
		}
		left_matched_stream.close();
		right_matched_stream.close();
		left_processed_stream.close();
		right_processed_stream.close();
		cmd = "cp " + left_processed + " " + left_matched;
		logger->debug(cmd);
		run_shell_command(cmd);
		cmd = "cp " + left_processed + " " + lib.get_matched_left_read_filename(round);
		logger->debug(cmd);
		run_shell_command(cmd);
		if (lib.get_paired_end()){
			cmd = "cp " + right_processed + " " + right_matched;
			logger->debug(cmd);
			run_shell_command(cmd);
			cmd = "cp " + right_processed + " " + lib.get_matched_right_read_filename(round);
			logger->debug(cmd);
			run_shell_command(cmd);
		}
	}
	//RM HERE
	cmd = "rm -f " + index_name + "*";
	logger->debug(cmd);
	run_shell_command(cmd);
}

SRAssemblerMaster::~SRAssemblerMaster() {
	// TODO Auto-generated destructor stub
}
