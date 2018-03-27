/*
 * Utility.h
 *
 *  Created on: Jul 27, 2010
 *      Author: hchou
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include <string>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iomanip>      // std::setw
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include "Const.h"
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/foreach.hpp>

using namespace std;

typedef boost::unordered_map<std::string, string> string_map;
typedef boost::unordered_map<std::string, std::tuple<int, double>> tuple_map;

string readfile(const string& fn);
void tokenize(const string& str, vector<string>& tokens, const string& delimiters);
void run_shell_command(const string& cmd);
string run_shell_command_with_return(const string& cmd);
int str2int(const string &str);
string int2str(const int n);
string int2str(const int n, const int length);
string long2str(const int n);
double str2double(const string &str);
string double2str(const double n);
double randgen(double low, double high);
int randgen(int low, int high);
bool file_exists(const char* filename);
bool file_exists(const string &filename);
long get_file_Size(const char* filename);
long get_file_size(const string &filename);
long get_read_count(const string& filename, int format);
string get_file_base_name(const string& path);
string get_file_name(const string& path);
string trim(const string& str);
string string_format(const std::string fmt, ...);
void fastq2fasta(const string& fq, const string& fa);
void remove_duplicate_reads(const string& filename, int read_format);
unsigned int count_letters(string &str);

#endif /* UTILITY_H_ */
