/*
 * VmatchAligner.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#include "VmatchAligner.h"
#include <boost/unordered_set.hpp>

using namespace std;

VmatchAligner::VmatchAligner(int log_level, string log_file):Aligner(log_level, log_file) {

}

/*
 * Parse Vmatch output file.
 * This function needs to :
 * Add new found reads to the found_reads_list
 * Add sequences of found reads to out_left_read and out_right_read
 *
 * Parameters:
 * output_file: the Vmatch alignment report file
 * found_reads: the reads currently found in this node
 * lib_idx: the integer number of the library for this part
 * read_part: the integer count of the split reads file for this part
 * left_read_index: the index name of the left reads of this library
 * right_read_index: the index name of the right reads of this library
 * out_left_read: the matched left-end reads file name
 * out_right_read: the matched left-end reads file name
 *
 */
int VmatchAligner::parse_output(const string& output_file, unordered_set<string>& found_reads, const int lib_idx, const int read_part, const string& left_read_index, const string& right_read_index, const string& out_left_read, const string& out_right_read) {
	// Parse the output file and get the mapped reads that had not been found yet.
	bool paired_end = (out_right_read != "");
	ifstream report_file_stream(output_file.c_str());
	int read_found = 0;
	int new_read_count = 0;
	string line;
	string cmd;
	string part_string = int2str(read_part);
	string lib_string = int2str(lib_idx);
	string tmpvseqselectfile = out_left_read + "-tmp";
	// TODO put the temporary file in tmp_dir
	ofstream tmp_file_stream(tmpvseqselectfile.c_str());

	while (getline(report_file_stream, line)) {
		read_found++;
		string seq_number = line;
		string seq_id = "lib" + lib_string + ",part" + part_string + ",read" + seq_number;
		// boost::unordered_set.find() produces past-the-end pointer if a key isn't found.
		if (found_reads.find(seq_id) == found_reads.end()) {
			new_read_count += 1;
			found_reads.insert(seq_id);
			tmp_file_stream << seq_number << '\n';
		}
	}
	logger->debug("Matched " + int2str(read_found) + " reads and " + int2str(new_read_count / 2) + " new read (pairs) in library " + int2str(lib_idx + 1) + ", part " + part_string);
	report_file_stream.close();
	tmp_file_stream.close();

	// Use awk to strip out linebreaks from within the sequence.
	cmd = "vseqselect -seqnum " + tmpvseqselectfile + " " + left_read_index + " | awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0} END { printf n }' >> " + out_left_read;
	logger->debug(cmd);
	run_shell_command(cmd);
	if (paired_end) {
		cmd = "vseqselect -seqnum " + tmpvseqselectfile + " " + right_read_index + " | awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0} END { printf n }' >> " + out_right_read;
		logger->debug(cmd);
		run_shell_command(cmd);
	}
	// TODO make this part of a cleanup function.
	// RM here
	cmd = "\\rm " + tmpvseqselectfile;
	run_shell_command(cmd);

	// We are catching two new reads if the library is paired end
	if (paired_end) {
		new_read_count *= 2;
	}
	return new_read_count;
}

void VmatchAligner::create_index(const string& index_name, const string& db_type, const string& fasta_file) {
	string cmd = "mkvtree -" + db_type + " -db " + fasta_file + " -pl -indexname " + index_name + " -allout >> " + logger->get_log_file();
	logger->debug(cmd);
	run_shell_command(cmd);
}

/* This function APPENDS the sequence numbers of hits (usually reads) to the output_file.
 * This is suitable input for vseqselect to pull those hits out of the mkvtree index.
 */
void VmatchAligner::do_alignment(const string& index_name, const string& alignment_type, int match_length, int mismatch_allowed, const string& query_file, const Params& params, const string& output_file) {
	// Are mismatches allowed by default? If not, empty string.
	string e_option = "";
	string l_option;
	if (mismatch_allowed > 0) {
		e_option = " -e " + int2str(mismatch_allowed);
	}
	// Other parameters are set up here.
	string param_list = "";
	for ( Params::const_iterator it = params.begin(); it != params.end(); ++it ){
			// Allowed number of mismatches may be handled by arguments or from params
			if (it->first == "e") {
				if (str2int(it->second) > 0){
					e_option = " -e " + it->second;
				}
				continue;
			}
			if (it->first == "l")
				match_length = str2int(it->second);
			else
				param_list += " -" + it->first + " " + it->second;
	}
	// If match_length is 0, use "-complete" for a complete length match.
	if (match_length > 0) {
		l_option = " -l " + int2str(match_length);
	} else {
		l_option = " -complete";
	}
	string cmd;
	string tmpvmfile = output_file + "-tmp";
	/* Vmatch output is appended to the output file so that left and right read searches for one part go into the same output file.
	 * In the DNA searches, the "-d -p" options search both strands.
	 * AWK selects only the column containing the sequence number. It is column 2 for "proteins" and "reads" because in those cases reads are being used as queries, not subjects.
	 * The results are sorted and each unique hit reported only once. This is not piped because the output_file may already have left reads in it.
	 */
	if (alignment_type == "protein" ) {
		cmd = "vmatch -dnavsprot 1 -q " + query_file + " -d" + l_option + e_option + " " + param_list + " -nodist -noevalue -noscore -noidentity " + index_name + " | awk '$0 !~ /^#.*/ {print $6}' >> " + output_file + "; sort -nu " + output_file + " > " + tmpvmfile + "; \\mv " + tmpvmfile + " " + output_file;
	} else if (alignment_type == "dna" ) {
		cmd = "vmatch -q " + query_file + " -d -p" + l_option + e_option + " " + param_list + " -nodist -noevalue -noscore -noidentity " + index_name + " | awk '$0 !~ /^#.*/ {print $2}' >> " + output_file + "; sort -nu " + output_file + " > " + tmpvmfile + "; \\mv " + tmpvmfile + " " + output_file;
	} else if (alignment_type == "reads") {
		cmd = "vmatch -q " + query_file + " -d -p" + l_option + e_option + " " + param_list + " -nodist -noevalue -noscore -noidentity " + index_name + " | awk '$0 !~ /^#.*/ {print $6}' >> " + output_file + "; sort -nu " + output_file + " > " + tmpvmfile + "; \\mv " + tmpvmfile + " " + output_file;
	}
	logger->debug(cmd);
	run_shell_command(cmd);
}

int VmatchAligner::get_format() {
	return FORMAT_FASTA;
}

bool VmatchAligner::is_available() {
	int ret = system("vmatch > /dev/null 2>&1");
	if (WEXITSTATUS(ret) != 1) {
		cout << "Cannot find vmatch, check your PATH variable!" << endl;
		return false;
	}
	return true;
}

string VmatchAligner::get_program_name() {
	return "Vmatch";
}


// TODO refactor
void VmatchAligner::align_long_contigs(const string& long_contig_candidate_file, const string& aux_dir, const string& contig_file, const int max_contig_size, unordered_set<string>& candidate_ids, unordered_set<string>& long_contig_ids) {
	if (!file_exists(long_contig_candidate_file)){
		return;
	}

	string indexname = aux_dir + "/contigs_index";
	string vmatch_out_file = aux_dir + "/long_contigs_candidates.vmatch";
	string cmd = "mkvtree -dna -db " + contig_file + " -pl -indexname " + indexname + " -allout >> " + logger->get_log_file();
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = "vmatch -q " + long_contig_candidate_file + " -d -p" + " -l " + int2str(max_contig_size) + " -showdesc 0 -nodist -noevalue -noscore -noidentity " + indexname + " | awk '{print $1,$2,$3,$4,$5,$6}' | uniq -f5 > " + vmatch_out_file;
	logger->debug(cmd);
	run_shell_command(cmd);
	ifstream vmatch_file_stream(vmatch_out_file.c_str());
	// parse the output file and remove the long contigs and associated reads.
	string line;
	while (getline(vmatch_file_stream, line)){
		if (line[0] != '#'){
			vector<string> tokens;
			tokenize(line, tokens, " ");
			string long_contig_id = tokens[1];
			string candidate_contig_id = tokens[5];
			long_contig_ids.insert(long_contig_id);
			candidate_ids.insert(candidate_contig_id);
		}
	}
	vmatch_file_stream.close();
	// RM here
	cmd = "rm " + vmatch_out_file + " " + indexname + "*";
	run_shell_command(cmd);
}

VmatchAligner::~VmatchAligner() {
	// TODO Auto-generated destructor stub
}
