/*
 * Logger.cpp
 *
 *  Created on: Oct 22, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#include "Logger.h"
#include "Utility.h"
#include <stdio.h>
#include <time.h>

Logger* Logger::logger=NULL;

Logger::Logger(int level, string log_file):log_level(level), log_file(log_file)	 {
	// Auto-generated constructor stub
}

Logger::~Logger() {
	// Auto-generated destructor stub
}

Logger *Logger::getInstance(int level, string log_file)
{
	if (logger == NULL) {
		logger = new Logger(level, log_file);
	}
	return logger;
}

void Logger::fatal(const string& msg)
{
	print_message(msg, "[FATAL]", (log_level <= Logger::LEVEL_FATAL));
}

void Logger::info(const string& msg)
{
	print_message(msg, " [INFO]", (log_level <= Logger::LEVEL_INFO));
}

void Logger::running(const string& msg)
{
	print_message(msg, "  [RUN]", (log_level <= Logger::LEVEL_INFO));
}

void Logger::mpi(const string& msg)
{
	print_message(msg, "  [MPI]", (log_level <= Logger::LEVEL_INFO));
}

void Logger::warn(const string& msg)
{
	print_message(msg, " [WARN] ", (log_level <= Logger::LEVEL_WARN));
}

void Logger::error(const string& msg)
{
	print_message(msg, "[ERROR]", (log_level <= Logger::LEVEL_ERROR));
}

void Logger::debug(const string& msg)
{
	print_message(msg, "[DEBUG]", (log_level <= Logger::LEVEL_DEBUG));
}

string Logger::get_log_file(){
	return this->log_file;
}

int Logger::get_log_level(){
	return this->log_level;
}

void Logger::print_message(const string &msg, const string &level, bool to_std_out){

	time_t rawtime;
	struct tm * timeinfo;
	char buffer [40];
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	strftime (buffer,40,"[%F %T]",timeinfo);

	ofstream log_file_stream;
	log_file_stream.open(log_file.c_str(), ios::out | ios::app );
	log_file_stream << buffer << level + " " + msg << endl;
	log_file_stream.close();
	if (to_std_out)
		cout << buffer << level + " " + msg << endl;
}

void Logger::safe_run_shell_command(const string& cmd) {
// If the shell command does not exit properly, SRAssembler logs this and continues.
	int exitcode;
	exitcode = run_shell_command(cmd);
	if (exitcode != 0) {
		error("Command <<" + cmd + ">> errored out with status:" + int2str(exitcode));
	}
}

void Logger::fragile_run_shell_command(const string& cmd) {
// If the shell command does not exit properly, SRAssembler logs this and throws an error.
	int exitcode;
	exitcode = run_shell_command(cmd);
	if (exitcode != 0) {
		error("Command <<" + cmd + ">> errored out with status:" + int2str(exitcode));
		throw exitcode;
	}
}
