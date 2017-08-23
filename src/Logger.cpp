/*
 * Logger.cpp
 *
 *  Created on: Oct 22, 2011
 *      Author: hchou
 */

#include "Logger.h"
#include <stdio.h>
#include <time.h>

Logger* Logger::logger=NULL;

Logger::Logger(int level, string log_file):log_level(level), log_file(log_file)	 {
	// TODO Auto-generated constructor stub

}

Logger::~Logger() {
	// TODO Auto-generated destructor stub
}

Logger *Logger::getInstance(int level, string log_file)
{
	if (logger == NULL)
		logger = new Logger(level, log_file);
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

	strftime (buffer,40,"[%c]",timeinfo);

	ofstream log_file_stream;
	log_file_stream.open(log_file.c_str(), ios::out | ios::app );
	log_file_stream << buffer << level << " " << msg << endl;
	log_file_stream.close();
	if (to_std_out)
	    cout << buffer << level << " " << msg << endl;
}

