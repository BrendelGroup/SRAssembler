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

string Library::get_matched_left_read_name(){
	return tmp_dir + "/matched_reads_left_" + "l" + int2str(lib_idx+1) + "." + this->file_extension;
}

string Library::get_matched_left_read_name(int round){
	return tmp_dir + "/matched_reads_left_" + "r" + int2str(round) + "_" + "l" + int2str(lib_idx+1) + "." + this->file_extension;
}

string Library::get_matched_left_read_name(int round, int idx){
	return tmp_dir + "/matched_reads_left_" + "r" + int2str(round) + "_" + "l" + int2str(lib_idx+1) + "_" + "s" + int2str(idx) + "." + this->file_extension;
}

string Library::get_matched_right_read_name(){
	return tmp_dir + "/matched_reads_right_" + "l" + int2str(lib_idx+1) + "." + this->file_extension;
}

string Library::get_matched_right_read_name(int round){
	return tmp_dir + "/matched_reads_right_" + "r" + int2str(round) + "_" + "l" + int2str(lib_idx+1) + "." + this->file_extension;
}

string Library::get_matched_right_read_name(int round, int idx){
	if (paired_end)
		return tmp_dir + "/matched_reads_right_" + "r" + int2str(round) + "_" + "l" + int2str(lib_idx+1) + "_" + "s" + int2str(idx) + "." + this->file_extension;
	return "";
}

string Library::get_joined_read_name(int round, int idx, int file_type){
//TODO: I think that this is a deprecated function. It doesn't appear to be used anywhere.
	string extension = (file_type == FORMAT_FASTQ)? "fastq" : "fasta";
	return tmp_dir + "/matched_reads_joined_" + "r" + int2str(round) + "_" + "l" + int2str(lib_idx+1) + "_" + "s" + int2str(idx) + "." + extension;
}

string Library::get_split_file_name(int idx, int file_type){
	string extension = (file_type == FORMAT_FASTQ)? "fastq" : "fasta";
	return data_dir + "/lib" + int2str(lib_idx+1) + "/" + get_file_base_name(left_read) + "_" + int2str(idx) + "." + extension;
}

string Library::get_prefix_split_src_file(string src_read){
	return tmp_dir + "/" + get_file_base_name(src_read) + "_" + int2str(lib_idx+1) + "_split_";
}

void Library::do_split_files(int read_type, int reads_per_file){
	string read_file = (read_type == LEFT_READ)? this->left_read:this->right_read;
	string cmd = "split -l " + int2str(reads_per_file * 4) + " " + read_file + " " + get_prefix_split_src_file(read_file);
	logger->debug(cmd);
	run_shell_command(cmd);
}
Library::~Library() {
	// TODO Auto-generated destructor stub
}
