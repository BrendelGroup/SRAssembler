/*
 * Utility.cpp
 *
 *  Created on: Jul 27, 2010
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#include "Utility.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <sys/stat.h> // stat, stat()

using namespace std;

string readfile(const string& fn) {
	ifstream inFile(fn.c_str());
	string filestring;
	string line;
	while(getline(inFile,line)) {
		filestring.append(line + "\n");
	}
	inFile.close();
	return filestring;
};

void fastq2fasta(const string& fq, const string& fa) {
	ifstream in_file(fq.c_str());
	ofstream out_file(fa.c_str());
	string line;
	while(getline(in_file,line)) {
		if (line.substr(0,1) == "@"){
			out_file << ">" << line.substr(1) << '\n';
			getline(in_file, line);
			out_file << line << '\n';
			getline(in_file, line);
			getline(in_file, line);
		}

	}
	in_file.close();
	out_file.close();
};

void tokenize(const string& line, vector<string>& tokens, const string& delimiters) {
	// Skip delimiters at beginning.
	string::size_type lastPos = line.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = line.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
		// Found a token, add it to the vector.
		tokens.push_back(line.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = line.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = line.find_first_of(delimiters, lastPos);
	}
}

void run_shell_command(const string& cmd) {
	char c[2048];
	strcpy(c,cmd.c_str());
	system(c);
}

string run_shell_command_with_return(const string& cmd)
{
	// setup
	string data;
	FILE *stream;
	int buf_size = 2048;
	char buffer[buf_size];

	// do it
	stream = popen(cmd.c_str(), "r");
	while ( fgets(buffer, buf_size, stream) != NULL )
	data.append(buffer);
	pclose(stream);

	// exit
	return data;
}

int str2int (const string &str) {
	int n;
	stringstream ss(str);
	ss >> n;
	return n;
};
string int2str (const int n) {
	stringstream ss;
	ss << n;
	return ss.str();
}
string int2str (const int n, const int length) {
	stringstream ss;
	ss << setw(length) << setfill('0') << n;
	return ss.str();
}
string long2str(const int n){
	stringstream ss;
	ss << n;
	return ss.str();
}

double str2double(const string &str){
	double n;
	stringstream ss(str);
	ss >> n;
	return n;
}
string double2str(const double n){
	stringstream ss;
	ss << n;
	return ss.str();
}

double randgen(double XMin, double XMax){
	return XMin + rand() * (XMax - XMin) / RAND_MAX;
}

int randgen(int XMin, int XMax){
	return XMin + rand() % (XMax - XMin + 1);
}

bool file_exists(const char* filename) {
	struct stat info;
	int ret = -1;
	// Get the file attributes.
	ret = stat(filename, &info);
	return (ret == 0);
}
bool file_exists(const string& filename) {
	return file_exists(filename.c_str());
}

long get_file_size(const char* filename) {
	struct stat info;
	// Get the file attributes.
	stat(filename, &info);
	return info.st_size;
}
long get_file_size(const string& filename) {
	return get_file_size(filename.c_str());
}

long get_read_count(const string& filename, int format) {
	ifstream inFile(filename.c_str());
	string line;
	long line_count = 0;
	while (getline(inFile,line)) {
		line_count++;
	}
	inFile.close();
	int line_per_read = (format == FORMAT_FASTQ)?4:2;
	return (line_count/line_per_read);
}

string get_file_name(const string& path) {
	int last_slash = path.find_last_of("/");
	return path.substr(last_slash + 1);
}

string get_file_base_name(const string& path) {
	int last_slash = path.find_last_of("/");
	int last_dot = path.find_last_of(".");
	return (last_dot == -1 || last_slash > last_dot)? path.substr(last_slash + 1):path.substr(last_slash + 1, last_dot - last_slash - 1);
}

string trim(const string& str) {
	size_t start = str.find_first_not_of(" \t\n\r");
	if(start == string::npos) return "";
	return str.substr(start, str.find_last_not_of(" \t\n\r") - start + 1);
}

string string_format(const string fmt, ...) {
	int size=100;
	string str;
	va_list ap;
	while (1) {
		str.resize(size);
		va_start(ap, fmt);
		int n = vsnprintf((char *)str.c_str(), size, fmt.c_str(), ap);
		va_end(ap);
		if (n > -1 && n < size) {
			str.resize(n);
			return str;
		}
		if (n > -1)
			size=n+1;
		else
			size*=2;
	}
}

unsigned int count_letters(string &str) {
	unsigned int count = 0;
	char const constexpr alpha[] = R"(abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ)";
	for(std::string::size_type i = 0; i < str.size(); ++i) {
		char test = str[i];
		if (any_of(begin(alpha), end(alpha), [test](char c){return c == test;})) {
			count++;
		}
	}
	return count;
}
