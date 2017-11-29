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

unordered_set<string> VmatchAligner::get_hit_list(const string& output_file) {
	logger->debug("get hit list from output file " + output_file);
	ifstream report_file_stream(output_file.c_str());
	boost::unordered_set<string> current_mapped_reads;
	string line;
	while (getline(report_file_stream, line)) {
		if (line[0] != '#') {
			vector<string> tokens;
			tokenize(line, tokens, " ");
			string seq_id = tokens[0];
			current_mapped_reads.insert(seq_id);
		}
	}
	report_file_stream.close();
	return current_mapped_reads;
}
/*
 * Parse Vmatch output file.
 * First this function takes the Vmatch output files and pulls the read IDs out of the MATCH lines
 * Then it reads through ALL of the reads for that library and pulls out the reads whose IDs were found
 * TODO change the read search into something using an index.
 * Parameters:
 * output_file: the Vmatch alignment report file
 * mapped_reads: the reads currently found in this node
 * source_read: the split reads file name for this part
 * out_left_read: the matched left-end reads file name
 * out_right_read: the matched left-end reads file name
 * joined_read: the matched joined reads file name
 * fastq_format: interleaved or split
 * format: fastq or fasta
 *
 */
//int VmatchAligner::parse_output(const string& output_file, unordered_set<string>& mapped_reads, const string& read_source, const string& out_left_read, const string& out_right_read, int fastq_format, int format) {
	//logger->debug("parsing output file " + output_file);
	//ifstream report_file_stream(output_file.c_str());
	//string seq_id;
	//// current_mapped_reads are the reads that show up in this run of parse_output()
	//boost::unordered_set<string> current_mapped_reads;
	//int found_new_read = 0;
	////parse the output file and get the mapped reads have not found yet
	//string line;
	//// Find mapped read IDs
	//while (getline(report_file_stream, line)) {
		//if (line[0] != '#'){
			//// The line should be a Vmatch output line. The read ID is column 6 when doing reads as query.
			//vector<string> tokens;
			//tokenize(line, tokens, " ");
			//string seq_id = tokens[5];
			///* boost::unordered_set.find() produces past-the-end pointer if a key isn't found
			 //* so we can confirm that a read isn't previously found with this if statement
			 //*/
			//if (mapped_reads.find(seq_id) == mapped_reads.end()) {
				//found_new_read = 1;
				//current_mapped_reads.insert(seq_id);
				//mapped_reads.insert(seq_id);
			//}
		//}
	//}
	//report_file_stream.close();
	//// The mapped read IDs for this have now been found.

	////fetch sequences
	//// paired_end set to true if out_right_read not empty string
	//bool paired_end = (out_right_read != "");

	//ifstream read_source_stream(read_source.c_str());
	//ofstream out_left_read_stream(out_left_read.c_str());
	////run_shell_command("printf '\e[38;5;002m" "OUT_LEFT_READ_STREAM:" + out_left_read + "\e[0m\n'");
	//ofstream out_right_read_stream;
	//if (paired_end)
		//out_right_read_stream.open(out_right_read.c_str(), ios_base::out);
	////ofstream joined_read_stream(joined_read.c_str());
	//string left_header = "";
	//string right_header = "";
	//string left_seq = "";
	//string right_seq = "";
	//string left_qual = "";
	//string right_qual = "";
	//string plus;
	//while (getline(read_source_stream, left_header)) {
		//string left_seq_id = "";
		//string right_seq_id = "";
		//string lead_chr = (format == FORMAT_FASTQ)? "@" : ">";
		////string lead_chr = ">";
		//if (left_header.substr(0,1) == lead_chr){
			//unsigned int pos = left_header.find_first_of(" ");
			//if (pos == string::npos)
				//left_seq_id = left_header;
			//else
				//left_seq_id = left_header.substr(1, pos-1);

			////run_shell_command("printf '\e[38;5;002mLEFT_SEQ_ID:" + left_seq_id + "\e[0m\n'");
			//getline(read_source_stream, left_seq);
			//// This can be removed if we don't use FASTQ internally
			//if (format == FORMAT_FASTQ) {
				//getline(read_source_stream, plus);
				//getline(read_source_stream, left_qual);
			//}
			//// Is FASTA_interleaved a thing? If so, variable name should be changed
			////TODO How do we determine if interleaved? It would be nice if SRAssembler could do either format.
			//if (paired_end && fastq_format == FASTQ_INTERLEAVED) {
				//getline(read_source_stream, right_header);
				//pos = right_header.find_first_of(" ");
				//if (pos ==string::npos)
					//right_seq_id = right_header;
				//else
					//right_seq_id = right_header.substr(1, pos-1);
				//getline(read_source_stream, right_seq);
				//if (format == FORMAT_FASTQ) {
					//getline(read_source_stream, plus);
					//getline(read_source_stream, right_qual);
				//}
			//}
			//// If seq_id IS in current_mapped_reads, add read (or read pair) to the output reads
			////TODO speed this up by using indexing? Just a dictionary of the reads would be faster.
			//if (current_mapped_reads.find(left_seq_id) != current_mapped_reads.end() || (paired_end && current_mapped_reads.find(right_seq_id) != current_mapped_reads.end())){
				////run_shell_command("printf '\e[38;5;002m" "got a hit" "\e[0m\n'");
				//if (paired_end){
					////run_shell_command("printf '\e[38;5;002m" "got a paired-end hit" "\e[0m\n'");
					//if (fastq_format == FASTQ_INTERLEAVED){
						////run_shell_command("printf '\e[38;5;002m" "got an interleaved paired-end hit" "\e[0m\n'");
						////out_left_read_stream << "@" << left_seq_id << "_1" << endl << left_seq << endl << "+" << endl << left_qual << endl;
						////out_right_read_stream << "@" << right_seq_id << "_2" << endl << right_seq << endl << "+" << endl << right_qual << endl;
						//if (format == FORMAT_FASTA){
							////run_shell_command("printf '\e[38;5;002m" "got an interleaved paired-end FASTA hit" "\e[0m\n'");
							//out_left_read_stream << ">" << left_seq_id << endl << left_seq << endl;
							//out_right_read_stream << ">" << right_seq_id << endl << right_seq << endl;
						//} else {
							//out_left_read_stream << "@" << left_seq_id << endl << left_seq << endl << "+" << endl << left_qual << endl;
							//out_right_read_stream << "@" << right_seq_id << endl << right_seq << endl << "+" << endl << right_qual << endl;
						//}
						////joined_read_stream << ">" << left_seq_id << endl << left_seq << endl << ">" << right_seq_id << endl << right_seq << endl;
					//}
					//if (fastq_format == FASTQ_JOINED) {
						//string left_seq_part = left_seq.substr(0, left_seq.length()/2);
						//string right_seq_part = left_seq.substr(left_seq.length()/2);
						//string left_qual_part = left_qual.substr(0, left_qual.length()/2);
						//string right_qual_part = left_qual.substr(left_qual.length()/2);
						////out_left_read_stream << "@" << left_seq_id << "_1" << endl << left_seq_part << endl << "+" << endl << left_qual_part << endl;
						////out_right_read_stream << "@" << right_seq_id << "_2" << endl << right_seq_part << endl << "+" << endl << right_qual_part << endl;
						//if (format == FORMAT_FASTA){
							//out_left_read_stream << ">" << left_seq_id << endl << left_seq_part << endl;
							//out_right_read_stream << ">" << right_seq_id << endl << right_seq_part << endl;
						//} else {
							//out_left_read_stream << "@" << left_seq_id << endl << left_seq_part << endl << "+" << endl << left_qual_part << endl;
							//out_right_read_stream << "@" << right_seq_id << endl << right_seq_part << endl << "+" << endl << right_qual_part << endl;
						//}
						////joined_read_stream << ">" << left_seq_id << endl << left_seq_part << right_seq_part << endl;
					//}
				//}
				//else {
					//if (format == FORMAT_FASTA)
						//out_left_read_stream << ">" << left_seq_id << endl << left_seq << endl;
					//else
						//out_left_read_stream << "@" << left_seq_id << endl << left_seq << endl << "+" << endl << left_qual << endl;
					////joined_read_stream << ">" << left_seq_id << endl << left_seq << endl;
				//}
			//}
		//}
	//}
	//read_source_stream.close();
	//out_left_read_stream.close();
	//if (paired_end)
		//out_right_read_stream.close();
	////joined_read_stream.close();
	//return found_new_read;
//}

int VmatchAligner::parse_output(const string& output_file, unordered_set<string>& mapped_reads, const string& read_source, const string& out_left_read, const string& out_right_read, int fastq_format, int format) {
	logger->debug("parsing output file " + output_file);
	ifstream report_file_stream(output_file.c_str());
	string seq_id;
	// current_mapped_reads are the reads that show up in this run of parse_output()
	boost::unordered_set<string> current_mapped_reads;
	int found_new_read = 0;
	//parse the output file and get the mapped reads have not found yet
	string line;
	// Find mapped read IDs
	while (getline(report_file_stream, line)) {
		if (line[0] != '#'){
			// The line should be a Vmatch output line. The read ID is column 6 when doing reads as query.
			vector<string> tokens;
			tokenize(line, tokens, " ");
			string seq_id = tokens[0];
			/* boost::unordered_set.find() produces past-the-end pointer if a key isn't found
			 * so we can confirm that a read isn't previously found with this if statement
			 */
			if (mapped_reads.find(seq_id) == mapped_reads.end()) {
				found_new_read = 1;
				current_mapped_reads.insert(seq_id);
				mapped_reads.insert(seq_id);
			}
		}
	}
	report_file_stream.close();
	// The mapped read IDs for this have now been found.

	//fetch sequences
	// paired_end set to true if out_right_read not empty string
	bool paired_end = (out_right_read != "");

	ifstream read_source_stream(read_source.c_str());
	ofstream out_left_read_stream(out_left_read.c_str());
	//run_shell_command("printf '\e[38;5;002m" "OUT_LEFT_READ_STREAM:" + out_left_read + "\e[0m\n'");
	ofstream out_right_read_stream;
	if (paired_end)
		out_right_read_stream.open(out_right_read.c_str(), ios_base::out);
	//ofstream joined_read_stream(joined_read.c_str());
	string left_header = "";
	string right_header = "";
	string left_seq = "";
	string right_seq = "";
	string left_qual = "";
	string right_qual = "";
	string plus;
	while (getline(read_source_stream, left_header)) {
		string left_seq_id = "";
		string right_seq_id = "";
		string lead_chr = (format == FORMAT_FASTQ)? "@" : ">";
		//string lead_chr = ">";
		if (left_header.substr(0,1) == lead_chr){
			unsigned int pos = left_header.find_first_of(" ");
			if (pos == string::npos)
				left_seq_id = left_header;
			else
				left_seq_id = left_header.substr(1, pos-1);

			//run_shell_command("printf '\e[38;5;002mLEFT_SEQ_ID:" + left_seq_id + "\e[0m\n'");
			getline(read_source_stream, left_seq);
			// This can be removed if we don't use FASTQ internally
			if (format == FORMAT_FASTQ) {
				getline(read_source_stream, plus);
				getline(read_source_stream, left_qual);
			}
			// Is FASTA_interleaved a thing? If so, variable name should be changed
			//TODO How do we determine if interleaved? It would be nice if SRAssembler could do either format.
			if (paired_end && fastq_format == FASTQ_INTERLEAVED) {
				getline(read_source_stream, right_header);
				pos = right_header.find_first_of(" ");
				if (pos ==string::npos)
					right_seq_id = right_header;
				else
					right_seq_id = right_header.substr(1, pos-1);
				getline(read_source_stream, right_seq);
				if (format == FORMAT_FASTQ) {
					getline(read_source_stream, plus);
					getline(read_source_stream, right_qual);
				}
			}
			// If seq_id IS in current_mapped_reads, add read (or read pair) to the output reads
			//TODO speed this up by using indexing? Just a dictionary of the reads would be faster.
			if (current_mapped_reads.find(left_seq_id) != current_mapped_reads.end() || (paired_end && current_mapped_reads.find(right_seq_id) != current_mapped_reads.end())){
				//run_shell_command("printf '\e[38;5;002m" "got a hit" "\e[0m\n'");
				if (paired_end){
					//run_shell_command("printf '\e[38;5;002m" "got a paired-end hit" "\e[0m\n'");
					if (fastq_format == FASTQ_INTERLEAVED){
						//run_shell_command("printf '\e[38;5;002m" "got an interleaved paired-end hit" "\e[0m\n'");
						//out_left_read_stream << "@" << left_seq_id << "_1" << endl << left_seq << endl << "+" << endl << left_qual << endl;
						//out_right_read_stream << "@" << right_seq_id << "_2" << endl << right_seq << endl << "+" << endl << right_qual << endl;
						if (format == FORMAT_FASTA){
							//run_shell_command("printf '\e[38;5;002m" "got an interleaved paired-end FASTA hit" "\e[0m\n'");
							out_left_read_stream << ">" << left_seq_id << endl << left_seq << endl;
							out_right_read_stream << ">" << right_seq_id << endl << right_seq << endl;
						} else {
							out_left_read_stream << "@" << left_seq_id << endl << left_seq << endl << "+" << endl << left_qual << endl;
							out_right_read_stream << "@" << right_seq_id << endl << right_seq << endl << "+" << endl << right_qual << endl;
						}
						//joined_read_stream << ">" << left_seq_id << endl << left_seq << endl << ">" << right_seq_id << endl << right_seq << endl;
					}
					if (fastq_format == FASTQ_JOINED) {
						string left_seq_part = left_seq.substr(0, left_seq.length()/2);
						string right_seq_part = left_seq.substr(left_seq.length()/2);
						string left_qual_part = left_qual.substr(0, left_qual.length()/2);
						string right_qual_part = left_qual.substr(left_qual.length()/2);
						//out_left_read_stream << "@" << left_seq_id << "_1" << endl << left_seq_part << endl << "+" << endl << left_qual_part << endl;
						//out_right_read_stream << "@" << right_seq_id << "_2" << endl << right_seq_part << endl << "+" << endl << right_qual_part << endl;
						if (format == FORMAT_FASTA){
							out_left_read_stream << ">" << left_seq_id << endl << left_seq_part << endl;
							out_right_read_stream << ">" << right_seq_id << endl << right_seq_part << endl;
						} else {
							out_left_read_stream << "@" << left_seq_id << endl << left_seq_part << endl << "+" << endl << left_qual_part << endl;
							out_right_read_stream << "@" << right_seq_id << endl << right_seq_part << endl << "+" << endl << right_qual_part << endl;
						}
						//joined_read_stream << ">" << left_seq_id << endl << left_seq_part << right_seq_part << endl;
					}
				}
				else {
					if (format == FORMAT_FASTA)
						out_left_read_stream << ">" << left_seq_id << endl << left_seq << endl;
					else
						out_left_read_stream << "@" << left_seq_id << endl << left_seq << endl << "+" << endl << left_qual << endl;
					//joined_read_stream << ">" << left_seq_id << endl << left_seq << endl;
				}
			}
		}
	}
	read_source_stream.close();
	out_left_read_stream.close();
	if (paired_end)
		out_right_read_stream.close();
	//joined_read_stream.close();
	return found_new_read;
}

void VmatchAligner::create_index(const string& index_name, const string& type, const string& fasta_file) {
	string db_type = type;
	if (db_type == "cdna") db_type = "dna";
	string cmd = "mkvtree -" + db_type + " -db " + fasta_file + " -pl -indexname " + index_name + " -allout -v >> " + logger->get_log_file();
	logger->debug(cmd);
	run_shell_command(cmd);
}

void VmatchAligner::do_alignment(const string& index_name, const string& type, int match_length, int mismatch_allowed, const string& query_file, const Params& params, const string& output_file) {
	// Is this a protein query (round 1 only)? If not, empty string.
	string align_type = (type == "protein") ? "-dnavsprot 1": "";
	// Are mismatches allowed? If not, empty string.
	string e_option = "";
	if (mismatch_allowed > 0)
		e_option = " -e " + int2str(mismatch_allowed);
	// Other parameters are set up here.
	string param_list = "";
	for ( Params::const_iterator it = params.begin(); it != params.end(); ++it ){
			// Is "-e" option handled here, or above? Something seems vestigial.
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
	/*TODO
	 * change these commands to search the read index using the contigs.
	 * Check to see if the output order is still sorted and can be easily collapsed with uniq
	 */
	if (type == "protein" ) {
		string cmd = "vmatch " + align_type + " -q " + query_file + " -d" + " -l " + int2str(match_length) + " " + e_option + " " + param_list + " -showdesc 0 -nodist -noevalue -noscore -noidentity " + index_name + " | awk '$0 !~ /^#.*/ {print $6}' | uniq > " + output_file;
		logger->debug(cmd);
		run_shell_command(cmd);
	} else if (type == "cdna" ) {
		//cmd = "vmatch " + align_type + " -q " + query_file + " -d -p" + " -l " + int2str(match_length) + " " + e_option + " " + param_list + " -nodist -noevalue -noscore -noidentity " + index_name + " | awk '{print $1,$2,$3,$4,$5,$6}' | uniq -f2 >> " + output_file;
		string cmd = "vmatch " + align_type + " -q " + query_file + " -d -p" + " -l " + int2str(match_length) + " " + e_option + " " + param_list + " -showdesc 0 -nodist -noevalue -noscore -noidentity " + index_name + " | awk '$0 !~ /^#.*/ {print $2}' | sort -u >> " + output_file;
		logger->debug(cmd);
		run_shell_command(cmd);
	}
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

void VmatchAligner::align_long_contigs(const string& long_contig_candidate_file, const string& tmp_dir, const string& contig_file, const int max_contig_size, unordered_set<string>& candidate_ids, unordered_set<string>& long_contig_ids) {
	if (!file_exists(long_contig_candidate_file)){
		return;
	}

	string indexname = tmp_dir + "/contigs";
	string alignment_out_file = tmp_dir + "/output.aln";
	string cmd = "mkvtree -dna -db " + contig_file + " -pl -indexname " + indexname + " -allout -v >> " + logger->get_log_file();
	logger->debug(cmd);
	run_shell_command(cmd);
	cmd = "vmatch -q " + long_contig_candidate_file + " -l " + int2str(max_contig_size) + " -showdesc 0 -nodist -noevalue -noscore -noidentity " + indexname + " | awk '{print $1,$2,$3,$4,$5,$6}' | uniq -f5 > " + alignment_out_file;
	logger->debug(cmd);
	run_shell_command(cmd);
	ifstream out_file_stream(alignment_out_file.c_str());
	//parse the output file and remove the long contigs and associated reads.
	string line;
	while (getline(out_file_stream, line)){
		if (line[0] != '#'){
			vector<string> tokens;
			tokenize(line, tokens, " ");
			string long_contig_id = tokens[1];
			string candidate_contig_id = tokens[5];
			long_contig_ids.insert(long_contig_id);
			candidate_ids.insert(candidate_contig_id);
		}
	}
	out_file_stream.close();
}

VmatchAligner::~VmatchAligner() {
	// TODO Auto-generated destructor stub
}
