/*
 * SRAssemblerMaster.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: hchou
 */
#include "SRAssemblerMaster.h"


SRAssemblerMaster::SRAssemblerMaster() {
	// Auto-generated constructor stub
}

int SRAssemblerMaster::init(int argc, char * argv[], int rank, int mpiSize) {
	string command;
	string cmd;
	int ret = SRAssembler::init(argc, argv, rank, mpiSize);
	if (ret == -1)
		return ret;
	best_hits.insert(make_pair("coverage", std::make_tuple(0, 0.0))); // <round, coverage>
	if (!get_aligner(Aligner::PROTEIN_ALIGNER)->is_available()) return -1;
	if (!get_aligner(Aligner::DNA_ALIGNER)->is_available()) return -1;
	if (!get_assembler()->is_available()) return -1;
	if (!get_spliced_aligner()->is_available()) return -1;
	this->start_round = (over_write)?1:get_start_round();
	if (this->start_round == 1)
		create_folders();
	// Start the msg.log with the command used to run SRAsssembler.
	command = "Command: ";
	for (int i=0; i<argc;i++){
		command.append(argv[i]).append(" ");
	}
	logger->info(command);
	if (!preprocessing_only){
		// Add to the msg.log the Parameters file contents for reproducibility.
		logger->debug("Parameter file contents:");
		cmd = "cat " + param_file + " >> " + logger->get_log_file();
		logger->debug(cmd);
		run_shell_command(cmd);
	}
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
		ifstream probe_file(this->probe_file.c_str());
		string line;
		while (getline(probe_file, line)){
			if (line.substr(0,1) == ">"){
				vector<string> tokens;
				tokenize(line.substr(1), tokens, " 	");
				query_list.push_back(tokens[0]);
			}
		}
		probe_file.close();
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
	logger->info("We have " + int2str(libraries.size()) + " libraries");
	for (unsigned int i=0;i<this->libraries.size();i++){
		logger->info("library " + int2str(i+1) + ": " + libraries[i].get_library_name());
		logger->info("insert size: " + int2str(libraries[i].get_insert_size()));
		logger->info("left read: " + libraries[i].get_left_read());
		logger->info("right read: " + libraries[i].get_right_read());
		logger->info("reversed: " + int2str(libraries[i].get_reversed()));
		logger->info("Paired-end: " + int2str(libraries[i].get_paired_end()));
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
	output_content += "<PRE>\n\nThe final contigs file(<a href=\"" + get_file_name(this->final_contigs_file) + "\">" + get_file_name(this->final_contigs_file) + "</a>) are generated from round " + int2str(round) + " with k=" + int2str(best_k) + "\n";
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
	SplicedAligner* spliced_aligner = get_spliced_aligner();
	output_content += "</PRE><H3>Spliced alignment - " + spliced_aligner->get_program_name() + "</H3><HR>";
	output_content += "<PRE>\n";
	string spliced_filename = get_file_name(this->spliced_alignment_output_file);
	string hit_contigs_name = get_file_name(this->hit_contigs_file);
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
		// Test if split files have been generated.
		if (file_exists(lib->get_split_file_name(1, LEFT_READ))){
			lib->set_num_parts(get_file_count(lib->get_split_read_prefix(lib->get_left_read()) + "*.fasta"));
			// Test if split reads have been indexed by looking for the final index file that we expect to exist from the left reads.
			if (file_exists(lib->get_read_part_index_name(lib->get_num_parts(), LEFT_READ) + ".skp")){
				logger->info("Using previously split files for read library " + int2str(lib_index+1));
				broadcast_code(ACTION_TOTAL_PARTS, lib_index, lib->get_num_parts(), 0);
				continue;
			}
		}
		logger->info("Splitting read library " + int2str(lib_index+1) + " ...");
		// RM here
		// Remove any pre-existing files in case of an incomplete earlier pre-processing.
		cmd = "rm -f " + data_dir + "/lib" + int2str(lib_index+1) + "/" + get_file_base_name(lib->get_left_read()) + "* " + data_dir + "/lib" + int2str(lib_index+1) + "/" + get_file_base_name(lib->get_right_read()) + "*";
		run_shell_command(cmd);
		// 'source' variable is for identifying where an MPI message came from.
		int source;
		int part = 0;
		long long code_value;
		mpi_code code;
		// File splitting is handled by the Library.
		if (lib->get_paired_end() && mpiSize > 2){
			send_code(1, ACTION_SPLIT, lib_index, 1, 0);
			send_code(2, ACTION_SPLIT, lib_index, 2, 0);
			mpi_receive(code_value, source);
			mpi_receive(code_value, source);
		}
		else {
			lib->do_split_files(LEFT_READ, this->reads_per_file);
			if (lib->get_paired_end())
				lib->do_split_files(RIGHT_READ, this->reads_per_file);
		}
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
					send_code(part, ACTION_PRE_PROCESSING, lib_index, part, 0);
				}
				while(completed < lib->get_num_parts()){
					mpi_receive(code_value, source);
					completed++;
				}
			}
			else {
				for (part=1;part<mpiSize;part++){
					send_code(part, ACTION_PRE_PROCESSING, lib_index, part, 0);
				}
				while (completed < lib->get_num_parts()){
					mpi_receive(code_value, source);
					code = get_mpi_code(code_value);
					completed++;
					if (part <= lib->get_num_parts()){
						send_code(source, ACTION_PRE_PROCESSING, lib_index, part, 0);
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
	if (file_exists(aux_dir)){
		// This is counting backwards from the maximum round number until i=2.
		for (int i=this->num_rounds;i>1;i--) {
			bool found_previous = true;
			int mpiSize = (this->mpiSize == 0)? 1: this->mpiSize;
				// If reads file and next round index file exist (means assembled), we can continue from here.
				for (unsigned int lib_idx=0;lib_idx<this->libraries.size();lib_idx++) {
					Library lib = this->libraries[lib_idx];
					found_previous = file_exists(get_query_fasta_file_name(i+1));
				}
			if (found_previous) {
				int procID=getpid();
				broadcast_code(ACTION_MEMDIR, 0, procID, 0);
				// Make sure there is a tmp_dir.
				this->tmp_dir = this->tmp_loc + "/SRAssemblermem" + int2str(procID);
				string cmd = "\\rm -rf " + tmp_dir + "; mkdir " + tmp_dir;
				logger->debug(cmd);
				run_shell_command(cmd);
				// Make sure that the existence of the tmp_dir is obvious in case of disrupted run.
				cmd = "ln --symbolic --target-directory=" + out_dir + " " + tmp_dir;
				logger->debug(cmd);
				run_shell_command(cmd);
				// Make sure the query is indexed for cleaning rounds.
				SRAssembler::create_index(1);
				start_round = i+1;
				if (start_round > this->num_rounds)
					return start_round;
				logger->info("Previous results found. SRAssembler starts from round " + int2str(start_round));
				// Clean the temp results if it is not complete.
				run_shell_command("find " + aux_dir + " -name \"matched_reads_left_" + "r" + int2str(start_round) + "*\" -o -name \"matched_reads_right_" + "r" + int2str(start_round) + "*\" -delete");
				for (unsigned int lib_idx=0;lib_idx<this->libraries.size();lib_idx++) {
					Library lib = this->libraries[lib_idx];
					run_shell_command("cp " + lib.get_matched_left_reads_filename(i) + " " + lib.get_matched_left_reads_filename());
					if(lib.get_paired_end()) {
						run_shell_command("cp " + lib.get_matched_right_reads_filename(i) + " " + lib.get_matched_right_reads_filename());
					}
				}
				if (mpiSize > 1) {
					for (int i=1;i<mpiSize;i++)
						send_code(i, ACTION_LOAD_PREVIOUS, start_round - 1, 0, 0);
					int completed = 0;
					long long code_value = 0;
					int source = 0;
					while(completed < mpiSize - 1) {
						mpi_receive(code_value, source);
						completed++;
					}
				} else {
					load_found_reads(start_round - 1);
				}
				load_long_contigs();
				break;
			}
		}
	}
	return start_round;
}

void SRAssemblerMaster::do_walking() {
	if (preprocessing_only) {
		logger->info("Do pre-processing of reads only. The chromosome walking is skipped.");
		broadcast_code(ACTION_EXIT, 0, 0, 0);
		return;
	}

	logger->info("Start chromosome walking ...");
	logger->info("Total processors: " + int2str(mpiSize));
	int source;
	int read_part = 0;
	long long code_value;
	mpi_code code;
	int round = this->start_round;
	if (round > this->num_rounds){
		logger->info("The previous assembly has been completed");
		broadcast_code(ACTION_EXIT, 0, 0, 0);
		return;
	}
	output_summary_header();
	bool assembled;
	// If this is a continuation of a previous run, cleaning may not happen on the same schedule as it would if the round were restarted.
	int rounds_since_cleaning = 0;
	// Walking begins.
	while(true){
		logger->info("Starting round " + int2str(round) + " ...");
		int new_reads_count = 0;
		assembled = false;
		rounds_since_cleaning += 1;

		// Index the query in round 1 for the dnavsprot vmatch step.
		if (round == 1) {
			create_index(1);
		} else {
			// Mask the previous round's contigs for later searches.
			if (assembly_round < round) {
				mask_contigs(round-1);
			}
		}

		// For each library, align the reads to the queries (the query file in round 1, previously found reads if the assembly round has not yet been reached, or the assembled contigs from the previous round).
		for (unsigned lib_idx=0;lib_idx<this->libraries.size();lib_idx++){
			int completed = 0;
			Library lib = this->libraries[lib_idx];
			if (mpiSize == 1){
				for (read_part=1; read_part<=lib.get_num_parts(); read_part++){
					new_reads_count += do_alignment(round, lib_idx, read_part);
				}
			} else {
				// If there are fewer split read files than processors
				if (lib.get_num_parts() < mpiSize){
					for (read_part=1; read_part<=lib.get_num_parts(); read_part++){
						send_code(read_part, ACTION_ALIGNMENT, round, read_part, lib_idx);
					}
					while(completed < lib.get_num_parts()){
						mpi_receive(code_value, source);
						code = get_mpi_code(code_value);
						int found_new_reads = code.value2;
						new_reads_count += found_new_reads;
						completed++;
					}
				// If there are more split read files than processors
				} else {
					for (read_part=1;read_part<mpiSize;read_part++){
						send_code(read_part, ACTION_ALIGNMENT, round, read_part, lib_idx);
					}
					while (completed < lib.get_num_parts()){
						mpi_receive(code_value, source);
						code = get_mpi_code(code_value);
						int found_new_reads = code.value2;
						int file_idx = code.value1;
						new_reads_count += found_new_reads;
						completed++;
						// As files are completed, new files are sent to slaves to be aligned.
						int next_file_idx = file_idx + mpiSize - 1;
						if (next_file_idx <= lib.get_num_parts())
							send_code(source, ACTION_ALIGNMENT, round, next_file_idx, lib_idx);
					}
				}
			}
		}
		// At the end of the round, save the found reads in case you want to restart the run.
		if (mpiSize == 1){
			save_found_reads(round);
		} else {
			// Have each slave save its own found reads.
			int slave;
			for (slave=1; slave < mpiSize; slave++) {
				send_code(slave, ACTION_SAVE, round, 0, 0);
				// Wait until all the slaves have saved their found reads.
				mpi_receive(code_value, source);
			}
		}
		if (new_reads_count == 0) {
			logger->info("The walking is terminated: No new reads found.");
			break;
		}
		merge_mapped_files(round);
		long read_count = get_total_read_count(round);
		logger->info("Found new reads: " + int2str(new_reads_count) + " \tTotal matched reads: " + int2str(read_count));

		// Assemble the reads.
		if (assembly_round <= round){
			unsigned int longest_contig = do_assembly(round);
			assembled = true;
			summary_best += int2str(read_count) + "\n";
			bool cleaned = false;

			// This is a hack to avoid slowdown due to the assembly of an unreasonable number of contigs.
			// This typically only happens in the event of a repeat element in an intron.
			if (! ignore_contig_explosion) {
				string contig_line_count = run_shell_command_with_return("wc -l " + get_contig_file_name(round));
				int contig_count = str2int(contig_line_count) / 2;
				// If there are too many contigs, first try cleaning them of bad ones.
				if (contig_count > contig_limit) {
					logger->debug("Alarmingly high (" + int2str(contig_count) + ") number of contigs assembled, attempting clean.");
					remove_no_hit_contigs(round);
					contig_line_count = run_shell_command_with_return("wc -l " + get_contig_file_name(round));
					contig_count = str2int(contig_line_count) / 2;
					// If cleaning didn't work, bail on this run.
					if (contig_count > contig_limit) {
						logger->info("The walking is terminated: " + int2str(contig_count) + " contigs produced in round " + int2str(round) + ". This is too many to be a good run. Consider adjusting parameters such as Vmatch_protein_vs_contigs or increasing -i initial_contig_min. You can also rerun with the -d 0 argument to ignore contig numbers.");
						broadcast_code(ACTION_EXIT, 0, 0, 0);
						// RM HERE
						clean_tmp_files(round-1);
						return;
					}
					// If cleaning no_hit_contigs did work, finish cleaning by removing the reads that don't match the contigs that were kept.
					remove_unmapped_reads(round);
					cleaned = true;
					rounds_since_cleaning = 0;
				}
			}

			bool no_reads = true;
			// if no reads found, stop
			for (unsigned int lib_idx=0; lib_idx < this->libraries.size(); lib_idx++){
				if (get_file_size(libraries[lib_idx].get_matched_left_reads_filename()) > 0) {
					no_reads = false;
					break;
				}
			}
			if (no_reads){
				logger->info("The walking is terminated: No new reads found after removing reads associated with the assembled contigs.");
				break;
			}

			// if no contigs found, stop
			if (get_file_size(get_contig_file_name(round)) == 0) {
				logger->info("The walking is terminated: No contigs found.");
				break;
			}
			// If maximum round is reached, stop
			logger->info("Round " + int2str(round) + " is done.");
			if (round == num_rounds) {
				logger->info("The walking is terminated: The maximum round (" + int2str(num_rounds) + ") has been reached.");
				break;
			}
			// If it is time to clean out useless contigs and reads, do so.
			// This happens on a schedule based on the -b option, but can be disrupted by assembling too many contigs, or by restarting a previous round.
			if (round > 1 && !cleaned) {
				if (rounds_since_cleaning == clean_round) {
					// Assembled contigs that don't have some degree of hit to the probe are removed.
					remove_no_hit_contigs(round);
					remove_unmapped_reads(round);
					cleaned = true;
					rounds_since_cleaning = 0;
				}
			}
			// If we haven't yet found a contig that meets the required length, we try the next round.
			if (longest_contig < min_contig_lgth) {
				// RM HERE
				clean_tmp_files(round-1);
				round++;
				continue;
			}

			// do spliced alignment and remove the probe sequences already assembled
			if (round >= assembly_round && check_gene_assembled){
				string_map query_map = do_spliced_alignment(round);
				logger->debug("Best coverage so far is in round " + int2str(std::get<0>(best_hits["coverage"])) + " with coverage " + double2str(std::get<1>(best_hits["coverage"])) + ".");
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
					// Assembled contigs that don't have some degree of hit to the remaining probes are removed.
					remove_no_hit_contigs(round);
					remove_unmapped_reads(round);
					cleaned = true;
					rounds_since_cleaning = 0;
				}
			}
		// If not assembly round yet.
		} else {
			// Collect found reads into a file to be used as query in next round.
			string joined_file = aux_dir + "/matched_reads_joined.fasta";
			string cmd;
			// For each library, append the matched reads to the matched_reads_joined.fasta
			for (unsigned i=0;i<this->libraries.size();i++){
				string left_file = aux_dir + "/matched_reads_left_" + "lib" + int2str(i+1) + ".fasta";
				string right_file = aux_dir + "/matched_reads_right_" + "lib" + int2str(i+1) + ".fasta";
				if (libraries[i].get_paired_end()){
					cmd = "cat " + left_file + " " + right_file + " >> " + joined_file;
					logger->debug(cmd);
					run_shell_command(cmd);
				} else {
					cmd = "cat " + left_file + " >> " + joined_file;
					logger->debug(cmd);
					run_shell_command(cmd);
				}
			}
			logger->info("Round " + int2str(round) + " is done.");
			if (round == num_rounds) {
				logger->info("The walking is terminated: The maximum round (" + int2str(num_rounds) + ") has been reached.");
				break;
			}
		}
		// RM HERE
		clean_tmp_files(round-1);
		round++;
	}
	// Walking ends.
	// RM HERE
	clean_tmp_files(round-1);
	// If this round has not been assembled yet, do assembling.
	if ( ! assembled && assembly_round > round){
		do_assembly(round);
	}
	// Notify all slaves to stop listening.
	broadcast_code(ACTION_EXIT, 0, 0, 0);
	// if the final contig size is 0, then report the previous round
	while (round > 0) {
		logger->info("Checking the final contigs assembled in round " + int2str(round) + " ...");
		prepare_final_contigs_file(round);
		if (get_file_size(final_contigs_file) == 0){
			logger->info("... no contigs found in round " + int2str(round));
			round--;
		}
		else {
			break;
		}
	}
	if (get_file_size(final_contigs_file) == 0) {
		return;
	}
	do_spliced_alignment();
	do_gene_finding();
	// RM HERE
	this->get_spliced_aligner()->clean_files(this->final_contigs_file);
	output_summary(round);
	output_spliced_alignment();
	if (file_exists(this->gene_finding_output_file)) {
		output_gene_finding();
	}
	ofstream outFile(summary_file.c_str());
	outFile << output_content << endl;
	outFile.close();
	// RM HERE
	// Now that we're done, clean up unneccessary temporary files and the link to the tmp_dir
	string cmd = "rm -rf " + probe_file + ".* " + aux_dir + "/qindex.* " + aux_dir + "/cindex.* " + out_dir + "/" + get_file_name(tmp_dir) + " " + tmp_dir;
	logger->debug(cmd);
	run_shell_command(cmd);
}

void SRAssemblerMaster::clean_tmp_files(int round){
	if (round == 0) return;
	if (!tidy) return;
	// Remove unneccessary files from the previous round.
	string cmd;
	// TODO make into one command to run.
	logger->debug("Clean tmp files from round " + int2str(round));
	cmd = "rm -f " + aux_dir + "/matched_reads_left_" + "r" + int2str(round) + "_part* ";
	cmd += aux_dir + "/matched_reads_right_" + "r" + int2str(round) + "_part* ";
	cmd += aux_dir + "/matched_reads_" + "r" + int2str(round) + "_* ";
	cmd += aux_dir + "/query-vs-contig_" + "r" + int2str(round) + ".* ";
	cmd += aux_dir + "/query-vs-contig_k*" + "r" + int2str(round) + ".aln ";
	cmd += aux_dir + "/hit_contigs_" +"r" + int2str(round) + ".* ";
	cmd += aux_dir + "/long_contig_candidate_" + "r" + int2str(round) + ".* ";
	logger->debug(cmd);
	run_shell_command(cmd);
}

void SRAssemblerMaster::save_query_list(){
	string fn = aux_dir + "/query_list";
	// RM HERE
	run_shell_command("rm -rf " + fn);
	if (query_list.size() > 0){
		ofstream fs(fn.c_str());
		for (unsigned i=0;i<query_list.size();i++)
			fs << query_list[i] << '\n';
		fs.close();
	}
}

void SRAssemblerMaster::load_query_list(){
	string fn = aux_dir + "/query_list";
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
		// RM HERE
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

// This is used for saving contigs that reach the max length, and for contigs that have matched one query when there are multiple queries.
string SRAssemblerMaster::get_saved_contig_file_name(){
	return aux_dir + "/saved_contig.fasta";
}

int SRAssemblerMaster::do_assembly(int round) {
	logger->info("Doing assembly, round: " + int2str(round));
	int best_k = 0;
	unsigned int max_longest_contig = 0;
	unsigned int best_longest_contig = 0;
	int spliced_align_length = 0;
	int max_spliced_align = 0;
	int total_k = (end_k-start_k)/step_k + 1;
	int source;
	int i = 0;
	int completed = 0;
	long long code_value;
	if (mpiSize == 1){
		for (i=1; i<=total_k; i++)
			SRAssembler::do_assembly(round, start_k + (i-1)*step_k, 1);
	} else {
		if (total_k < mpiSize-1){
			// TODO Multithreading this doesn't work because SOAPdenovo2 can't just use the nodes that aren't in use.
			//int threads = (mpiSize - 1) / total_k;
			for (i=1; i<=total_k; i++){
				send_code(i, ACTION_ASSEMBLY, round, start_k + (i-1)*step_k, 1);
			}
			while(completed < total_k){
				mpi_receive(code_value, source);
				completed++;
			}
		} else {
			for (i=1;i<mpiSize;i++){
				send_code(i, ACTION_ASSEMBLY, round, start_k + (i-1)*step_k, 1);
			}
			while(completed < total_k){
				mpi_receive(code_value, source);
				completed++;
				if (i <= total_k) {
					send_code(source, ACTION_ASSEMBLY, round, start_k + (i-1)*step_k, 1);
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
			max_longest_contig = kstats.longest_contig;
		}
		stats.push_back(kstats);
		// Use the spliced alignment length to determine best k.
		spliced_align_length = do_spliced_alignment(round, k);
		if (spliced_align_length > max_spliced_align) {
			best_k = k;
			max_spliced_align = spliced_align_length;
			best_longest_contig = kstats.longest_contig;
		// Use the longest contig as a tie-breaker if necessary.
		} else if (spliced_align_length == max_spliced_align && kstats.longest_contig > best_longest_contig) {
			best_k = k;
			best_longest_contig = kstats.longest_contig;
		}
	}
	i = 0;
	for (int k=start_k;k<=end_k;k+=step_k) {
		Assembly_stats kstats = stats[i];
		if (kstats.longest_contig > 0) {
			logger->info("The longest contig (out of " + int2str(kstats.total_contig) + ") assembled in round\t" + int2str(round) + " is of length\t" + int2str(kstats.longest_contig) + " with k =\t" + int2str(k));
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
		logger->info("The best k-value (corresponding to the longest spliced alignment to the query) in round\t" + int2str(round) + " is k =\t" + int2str(best_k));
	}
	else {
		logger->info("No contig of the specified minimum length has been assembled by round\t" + int2str(round));
	}
	summary_best += int2str(best_k) + "\t";
	if (best_k > 0) {
		process_long_contigs(round, best_k);
	}
	// RM HERE
	if (tidy) {
		get_assembler()->clean_files(aux_dir);
	}
	return max_longest_contig;
}

void SRAssemblerMaster::load_long_contigs() {
	string long_contig_file = aux_dir + "/long_contig.fasta";
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
		saved_contig_file << contig_fasta;
	}
	saved_contig_file.close();
}

// TODO This should probably be refactored
// This function doesn't just process long contigs, it processes all contigs for length, including removing the short ones.
void SRAssemblerMaster::process_long_contigs(int round, int k) {
	string long_contig_candidate_file = aux_dir + "/long_contig_candidate_" + "r" + int2str(round-1)+ ".fasta";
	string long_contig_candidate_next_file = aux_dir + "/long_contig_candidate_" + "r" + int2str(round)+ ".fasta";
	string long_contig_file = aux_dir + "/long_contig_original.fasta";
	string long_contig_trimmed_file = aux_dir + "/long_contig.fasta";
	string saved_contig_file_name = get_saved_contig_file_name();
	boost::unordered_set<string> candidate_ids;
	boost::unordered_set<string> long_contig_ids;
	ofstream out_long_contig_trimmed;
	out_long_contig_trimmed.open((long_contig_trimmed_file.c_str()), ofstream::app);
	ofstream saved_contig_file;
	saved_contig_file.open(saved_contig_file_name.c_str(), fstream::app);
	string line;
	if (file_exists(long_contig_candidate_file) && get_file_size(long_contig_candidate_file) > 0){
		// Candidate long contigs from last round are aligned against this round's assembled contigs.
		// The ids of candidates that hit are stored in candidate_ids list.
		// The id of the matching contig from this round is stored in long_contig_ids list.
		this->get_aligner(round)->align_long_contigs(long_contig_candidate_file, aux_dir, this->get_assembly_file_name(round, k), this->max_contig_lgth, candidate_ids, long_contig_ids);
		// The matched candidate long contigs are added to the save file for contigs and to the accepted long contigs file.
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

	// Check the contigs of this round for any that exceed the max_length and make them candidate long contigs for the next round.
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
				// Contigs from this round that matched a candidate long contig from last round go into long_contig_original.fasta
				if (long_contig_ids.find(tokens[0]) != long_contig_ids.end())
					out_long_contig << header << '\n' << seq << '\n';
				else {
					// If a contig meets the length minimum for searching it goes into this round's list of assembled contigs and is a query in the next round's search
					if (seq.length() > this->query_contig_min) {
						out_contig << header << '\n' << seq << '\n';
					}
					// If a contig exceed the maximum contig length, a substring is stored in the candidate file to check and see if it is assembled again next round.
					if (seq.length() > this->max_contig_lgth) {
						// A substring of the sequence is taken from the middle.
						out_candidate_contig << header << '\n' << seq.substr((seq.length() - max_contig_lgth) / 2,max_contig_lgth) << '\n';
					}
				}
			}
			header = line;
			seq = "";
		}
		else
			seq.append(line);
	}
	// Finish processing the last contig.
	vector<string> tokens;
	tokenize(header.substr(1), tokens, " ");
	if (long_contig_ids.find(tokens[0]) != long_contig_ids.end())
		out_long_contig << header << '\n' << seq << '\n';
	else {
		if (seq.length() > this->query_contig_min)
			out_contig << header << '\n' << seq << '\n';
		if (seq.length() > this->max_contig_lgth) {
			out_candidate_contig << header << '\n' << seq.substr((seq.length() - max_contig_lgth) / 2,max_contig_lgth) << '\n';
		}
	}
	in_contig.close();
	out_contig.close();
	out_candidate_contig.close();
	out_long_contig.close();

	if (long_contig_ids.size() > 0){
		logger->info("Processing contigs longer than the specified maximal contig size " + int2str(this->max_contig_lgth) + " ...");
		// Remove associated reads.
		// TODO make this a function that remove_unmapped_reads also uses
		string cmd;
		string long_contig_index = aux_dir + "/long_cindex";
		Aligner* aligner = get_aligner(round);
		aligner->create_index(long_contig_index, "dna", long_contig_file);
		for (unsigned int lib_idx=0; lib_idx < this->libraries.size(); lib_idx++) {
			Library lib = this->libraries[lib_idx];
			// Index current matched reads for extraction purposes
			string left_matched_reads = lib.get_matched_left_reads_filename();
			string right_matched_reads;
			if (lib.get_paired_end()) {
				 right_matched_reads = lib.get_matched_right_reads_filename();
			}
			aligner->create_index(aux_dir + "/left_reads_index", "dna", left_matched_reads);
			if (lib.get_paired_end()) {
				aligner->create_index(aux_dir + "/right_reads_index", "dna", right_matched_reads);
			}

			// Use the found reads as queries against the long contigs to identify matchy reads
			string program_name = aligner->get_program_name();
			program_name += "_reads_vs_contigs";
			Params params = get_parameters(program_name);
			string vmatch_outfile = aux_dir + "/reads_vs_long_contigs.lib" + int2str(lib_idx+1) + ".round" + int2str(round) + ".vmatch";
			aligner->do_alignment(long_contig_index, "reads", 0, 2, left_matched_reads, params, vmatch_outfile);
				if (lib.get_paired_end()) {
					aligner->do_alignment(long_contig_index, "reads", 0, 2, right_matched_reads, params, vmatch_outfile);
				}

			// Use vseqselect to AVOID matchy reads

			// The .prj index file includes the total number of reads in the index
			int readcount;
			string reads_index_prj = aux_dir + "/left_reads_index.prj";
			ifstream reads_index_prj_stream(reads_index_prj.c_str());
			while (getline(reads_index_prj_stream, line)){
				if (line.substr(0,17) == "numofdbsequences="){
					readcount = str2int(line.substr(17));
					break;
				}
			}
			reads_index_prj_stream.close();

			// Complement the set of matched reads
			string vmatch_complement = aux_dir + "/reads_vs_long_contigs.lib" + int2str(lib_idx+1) + ".round" + int2str(round) + ".complement";
			ifstream vmatch_stream(vmatch_outfile.c_str());
			ofstream complement_stream(vmatch_complement.c_str());
			getline(vmatch_stream, line);
			// For each potential read number, add to the complement if not in the vmatch output
			for (int i=0; i < readcount; i++) {
				if (i < str2int(line)) {
					complement_stream << i << '\n';
					continue;
				}
				// Because the vmatch output is sorted, we only iterate the vmatch_stream when we find matching numbers
				if (i == str2int(line)) {
					getline(vmatch_stream, line);  // This will eventually run out
					continue;
				}
				// If vmatch_stream is empty, line will have no value and we keep all remaining read numbers
				complement_stream << i << '\n';
			}
			vmatch_stream.close();
			complement_stream.close();

			logger->debug("Keep reads without hits against long_contigs in round " + int2str(round));
			cmd = "vseqselect -seqnum " + vmatch_complement + " " + aux_dir + "/left_reads_index | awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0} END { printf n }' > " + left_matched_reads;
			logger->debug(cmd);
			run_shell_command(cmd);
			cmd = "cp " + left_matched_reads + " " + lib.get_matched_left_reads_filename(round);
			logger->debug(cmd);
			run_shell_command(cmd);
			if (lib.get_paired_end()) {
				cmd = "vseqselect -seqnum " + vmatch_complement + " " + aux_dir + "/right_reads_index | awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0} END { printf n }' > " + right_matched_reads;
				logger->debug(cmd);
				run_shell_command(cmd);
				cmd = "cp " + right_matched_reads + " " + lib.get_matched_right_reads_filename(round);
				logger->debug(cmd);
				run_shell_command(cmd);
			}
		}
		// RM here
		cmd = "rm -f " + aux_dir + "/left_reads_index* " + aux_dir + "/right_reads_index*";
		run_shell_command(cmd);
	}
}

void SRAssemblerMaster::remove_hit_contigs(vector<string> &contig_list, int round){
	logger->debug("Remove hit contigs of round " + int2str(round));
	string contig_file = get_contig_file_name(round);
	string tmp_file = aux_dir + "/nonhit_contigs_" + "r" + int2str(round) + ".fasta";
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
					tmp_file_stream << ">" << header << '\n' << seq << '\n';
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
			tmp_file_stream << ">" << header << '\n' << seq << '\n';
	}
	contig_file_stream.close();
	tmp_file_stream.close();
	saved_contig_file.close();
	string cmd = "cp " + tmp_file + " " + contig_file;
	run_shell_command(cmd);
	// RM HERE
	if (tidy) {
		cmd = "rm " + tmp_file;
		run_shell_command(cmd);
	}
}

void SRAssemblerMaster::prepare_final_contigs_file(int round){
	logger->debug("Prepare final contig file of round " + int2str(round));
	string contig_file = get_contig_file_name(round);
	string saved_contig_file_name = get_saved_contig_file_name();
	ifstream last_round_contig(contig_file.c_str());
	ifstream saved_contig(saved_contig_file_name.c_str());
	ofstream final_contig(final_contigs_file.c_str());
	string line;
	while (getline(saved_contig, line))
		final_contig << line << '\n';
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
				final_contig << header << '\n' << seq << '\n';
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
		final_contig << header << '\n' << seq << '\n';
	}
}

void SRAssemblerMaster::create_folders(){
	string cmd;
	if (file_exists(aux_dir)){
		cmd = "rm -rf " + aux_dir;
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
	// Remove pre-existing library symlinks.
	// RM here
	cmd = "rm -f " + data_dir + "/lib{0..9}*";
	run_shell_command(cmd);
	for (unsigned i=0;i<this->libraries.size();i++){
		string dir = data_dir + "/" + libraries[i].get_library_name();
		cmd = "mkdir -p " + dir + "; ";
		string slink = data_dir + "/lib" + int2str(i+1);
		cmd += "ln -s " + libraries[i].get_library_name() + " " + slink;
		logger->info(cmd);
		run_shell_command(cmd);
	}
	// If pre-processing only, don't bother making useless directories.
	if (preprocessing_only){
		run_shell_command("mkdir " + results_dir);
		return;
	}

	cmd = "mkdir " + results_dir + " " + intermediate_dir + " " + aux_dir;
	run_shell_command(cmd);

	// Set unique directory for temporary files, ideally stored in RAM (/dev/shm).
	// If the run is disrupted, these files will remain until a computer reboot, potentially slowing down the computer.
	int procID=getpid();
	broadcast_code(ACTION_MEMDIR, 0, procID, 0);
	this->tmp_dir = this->tmp_loc + "/SRAssemblermem" + int2str(procID);
	// If a disrupted run left a conflicting file behind, it should be removed first.
	cmd = "\\rm -rf " + tmp_dir + "; mkdir " + tmp_dir;
	logger->debug(cmd);
	run_shell_command(cmd);
	// Make sure that the existence of the tmp_dir is obvious in case of disrupted run.
	cmd = "ln --symbolic --target-directory=" + out_dir + " " + tmp_dir;
	logger->debug(cmd);
	run_shell_command(cmd);
}

void SRAssemblerMaster::remove_no_hit_contigs(int round){
	logger->info("Removing contigs without hits ...");
	string cmd;
	string alignment_type;
	string contig_file = get_contig_file_name(round);
	string contig_index = aux_dir + "/cindex";
	run_shell_command("rm -f " + contig_index + "*");
	run_shell_command("cp " + contig_file + " " + contig_file + ".beforeclean");
	Aligner* aligner = get_aligner(round);
	// Index contigs for easy extraction of hit contigs
	aligner->create_index(contig_index, "dna", contig_file);
	string program_name = aligner->get_program_name();
	program_name += "_" + get_type(1) + "_vs_contigs";
	Params params = get_parameters(program_name);
	string out_file = aux_dir + "/query_vs_contig.round" + int2str(round) + ".vmatch";
	if (this->probe_type == "dna") {
		// The "reads" type alignment ensures that we keep the hit from the query (in this case, the contigs), not the index (the dna probe).
		alignment_type = "reads";
	} else {
		alignment_type = "protein";
	}
	// Use the index of the probe_file (qindex) created in the first round.
	aligner->do_alignment(aux_dir + "/qindex", alignment_type, get_match_length(1), get_mismatch_allowed(1), contig_file, params, out_file);
	logger->debug("Removing contigs without hits against query sequences in round " + int2str(round));
	cmd = "vseqselect -seqnum " + out_file + " " + contig_index + " | awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0} END { printf n }' > " + contig_file;
	logger->debug(cmd);
	run_shell_command(cmd);
	// RM here
	// TODO Restore allowing multiple queries again? //string cmd = "rm -f " + aux_dir + "/qindex*";
	cmd = "rm -f " + out_file + " " + contig_index + "*";
	run_shell_command(cmd);
	// Create a new index of the good contigs for remove_unmapped_reads to use. Happens here because remove_unmapped_reads might be parallel.
	aligner->create_index(contig_index, "dna", contig_file);
}

void SRAssemblerMaster::remove_unmapped_reads(int round){
	logger->info("Removing found reads without matched contigs ...");
	int source;
	unsigned int completed = 0;
	long long code_value;
	mpi_code code;
	if (mpiSize == 1){
		for (unsigned int lib_idx=0; lib_idx < this->libraries.size(); lib_idx++) {
		SRAssembler::remove_unmapped_reads(lib_idx, round);
		}
	} else {
		if (int(this->libraries.size()) < mpiSize){
			for (unsigned int lib_idx = 0; lib_idx < this->libraries.size(); lib_idx++){
				send_code(lib_idx + 1, ACTION_CLEAN, lib_idx, round, 0);
			}
			while (completed < this->libraries.size()){
				mpi_receive(code_value, source);
				completed++;
			}
		// If there are more libraries than processors
		} else {
			for (int lib_idx = 0; lib_idx < mpiSize - 1; lib_idx++){ // Zero indexes libraries
				send_code(lib_idx + 1, ACTION_CLEAN, lib_idx, round, 0);
			}
			while (completed < this->libraries.size()){
				mpi_receive(code_value, source);
				code = get_mpi_code(code_value);
				int lib_idx = code.value1;
				completed++;
				// As libraries are completed, new libraries are sent to slaves to be cleaned
				int next_lib_idx = lib_idx + mpiSize - 1;
				if (next_lib_idx < int(this->libraries.size())) {
					send_code(source, ACTION_CLEAN, next_lib_idx, round, 0);
				}
			}
		}
	}
}

SRAssemblerMaster::~SRAssemblerMaster() {
	// Auto-generated destructor stub
}
