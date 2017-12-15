/*
 * Library.cpp
 *
 *  Created on: Feb 27, 2012
 *      Author: hchou
 */

#include "Library.h"

Library::Library(unsigned int lib_idx, string data_dir, string tmp_dir, Logger* logger) {
	this->lib_idx = lib_idx;
	this->data_dir = data_dir;
	this->tmp_dir = tmp_dir;
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

void Library::set_insert_size(int insert_size){
	this->insert_size = insert_size;
}

void Library::set_num_parts(int num_parts){
	cerr << "Set library " + int2str(this->lib_idx + 1) + " to " + int2str(num_parts) + " parts." << endl;
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

string Library::get_matched_left_read_filename(){
	return tmp_dir + "/matched_reads_left_" + "lib" + int2str(lib_idx+1) + ".fasta";
}

string Library::get_matched_left_read_filename(int round){
	return tmp_dir + "/matched_reads_left_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + ".fasta";
}

string Library::get_matched_left_read_filename(int round, int part){
	return tmp_dir + "/matched_reads_left_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + "_" + "part" + int2str(part) + ".fasta";
}

string Library::get_matched_right_read_filename(){
	return tmp_dir + "/matched_reads_right_" + "lib" + int2str(lib_idx+1) + ".fasta";
}

string Library::get_matched_right_read_filename(int round){
	return tmp_dir + "/matched_reads_right_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + ".fasta";
}

string Library::get_matched_right_read_filename(int round, int part){
	if (paired_end)
		return tmp_dir + "/matched_reads_right_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + "_" + "part" + int2str(part) + ".fasta";
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

//void Library::do_split_files(int read_direction, int reads_per_file){
	//string read_file = (read_direction == LEFT_READ)? this->left_read:this->right_read;
	//string cmd = "split -l " + int2str(reads_per_file * 4) + " " + read_file + " " + get_split_read_prefix(read_file);
	//logger->debug(cmd);
	//run_shell_command(cmd);
//}
void Library::do_split_files(int read_direction, int reads_per_file){
	//run_shell_command("printf '\e[38;5;002m" "do_split_files INVOKED" "\e[0m\n'");
	string read_file = (read_direction == LEFT_READ)? this->left_read:this->right_read;
	string cmd;
	// multiplier is based on the number of lines in a read. Assumes FASTQ to start.
	if (get_format() == FORMAT_FASTQ) {
		//run_shell_command("printf '\e[38;5;002m" "split FASTQ files" "\e[0m\n'");
		// Use sed magic to turn fastq into fasta. Assumes single-line reads.
		// Use awk to split in a controllable way that gives nice suffixes
		cmd = "< " + read_file + " sed -n -e '1~4s/^@/>/p;2~4p' | awk -v prefix=" + get_split_read_prefix(read_file) + " -v lines=" + int2str(reads_per_file * 2) + " 'NR%lines==1 {++i; file = prefix i \".fasta\"} {print > file}'";
	} else {
		//run_shell_command("printf '\e[38;5;002m" "split FASTA files" "\e[0m\n'");
		cmd = "< " + read_file + " awk -v prefix=" + get_split_read_prefix(read_file) + " -v lines=" + int2str(reads_per_file * 2) + " 'NR%lines==1 {++i; file = prefix i \".fasta\"} {print > file}'";
	}
	logger->debug(cmd);
	run_shell_command(cmd);
}
Library::~Library() {
	// TODO Auto-generated destructor stub
}
