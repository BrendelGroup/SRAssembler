/*
 * Library.cpp
 *
 *  Created on: Feb 27, 2012
 *      Author: hchou
 */

#include "Library.h"

Library::Library(unsigned int lib_idx, string data_dir, string aux_dir, Logger* logger) {
	this->library_name = "";
	this->lib_idx = lib_idx;
	this->data_dir = data_dir;
	this->aux_dir = aux_dir;
	this->reversed = DIRECTION_FR;
	this->left_read = "";
	this->right_read = "";
	this->format = FORMAT_FASTQ;
	this->file_extension = "fastq";
	this->insert_size = INSERT_SIZE;
	this->num_parts = 0;
	this->paired_end = false;
	this->logger = logger;
}

string Library::get_library_name(){
	return this->library_name;
}

int Library::get_insert_size(){
	return this->insert_size;
}

int Library::get_num_parts(){
	return this->num_parts;
}

int Library::get_reversed(){
	return this->reversed;
}

int Library::get_format(){
	return this->format;
}

bool Library::get_paired_end(){
	return this->paired_end;
}

string Library::get_left_read(){
	return this->left_read;
}

string Library::get_right_read(){
	return this->right_read;
}

string Library::get_file_extension(){
	return this->file_extension;
}

void Library::set_library_name(string library_name){
	this->library_name = library_name;
}

void Library::set_insert_size(int insert_size){
	this->insert_size = insert_size;
}

void Library::set_num_parts(int num_parts){
	this->num_parts = num_parts;
}

void Library::set_reversed(int reversed){
	this->reversed = reversed;
}

void Library::set_format(int format){
	this->format = format;
	this->file_extension = (format == FORMAT_FASTQ)? "fastq" : "fasta";
}

void Library::set_paired_end(bool paired_end){
	this->paired_end = paired_end;
}

void Library::set_left_read(string left_read){
	this->left_read = left_read;
}

void Library::set_right_read(string right_read){
	this->right_read = right_read;
}

string Library::get_matched_left_reads_filename(){
	return aux_dir + "/matched_reads_left_" + "lib" + int2str(lib_idx+1) + ".fasta";
}

string Library::get_matched_left_reads_filename(int round){
	return aux_dir + "/matched_reads_left_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + ".fasta";
}

string Library::get_matched_left_reads_filename(int round, int part){
	return aux_dir + "/matched_reads_left_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + "_" + "part" + int2str(part) + ".fasta";
}

string Library::get_matched_right_reads_filename(){
	return aux_dir + "/matched_reads_right_" + "lib" + int2str(lib_idx+1) + ".fasta";
}

string Library::get_matched_right_reads_filename(int round){
	return aux_dir + "/matched_reads_right_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + ".fasta";
}

string Library::get_matched_right_reads_filename(int round, int part){
	if (paired_end)
		return aux_dir + "/matched_reads_right_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + "_" + "part" + int2str(part) + ".fasta";
	return "";
}

string Library::get_split_file_name(int file_part, int read_direction){
	string read_file = (read_direction == LEFT_READ)? this->left_read:this->right_read;
	return data_dir + "/lib" + int2str(lib_idx+1) + "/" + get_file_base_name(read_file) + "_" + "part" + int2str(file_part) + ".fasta";
}

string Library::get_read_part_index_name(int file_part, int read_direction){
	string read_file = (read_direction == LEFT_READ)? this->left_read:this->right_read;
	return data_dir + "/lib" + int2str(lib_idx+1) + "/" + get_file_base_name(read_file) + "_" + "part" + int2str(file_part);
}

string Library::get_split_read_prefix(string src_read){
	return data_dir + "/lib" + int2str(lib_idx+1) + "/" + get_file_base_name(src_read) + "_" + "part";
}

void Library::do_split_files(int read_direction, int reads_per_file){
	string read_file = (read_direction == LEFT_READ)? this->left_read:this->right_read;
	ifstream in_stream(read_file.c_str());
	string out_file;
	int part = 1;
	ofstream out_stream;
	string line;
	// Priming line
	getline(in_stream, line);
	// While not at end of input file
	while(! in_stream.eof()) {
		int linecount = 0;
		out_file = get_split_read_prefix(read_file) + int2str(part) + ".fasta";
		out_stream.open(out_file.c_str());
		while( linecount < reads_per_file * 2){
			// Write the header
			out_stream << ">" << line.substr(1) << '\n';
			linecount++;
			// Get the read sequence and write it
			getline(in_stream, line);
			out_stream << line << '\n';
			linecount++;
			if (get_format() == FORMAT_FASTQ) {
				// Get the plus line
				getline(in_stream, line);
				// Dump the plus line, get the quality score
				getline(in_stream, line);
			}
			// Dump the quality line, get the next header. End writing early of last outfile.
			if (! getline(in_stream, line)) {
				break;
			}
		}
		out_stream.close();
		part++;
	}
	out_stream.close();
	in_stream.close();
}

Library::~Library() {
	// Auto-generated destructor stub
}
