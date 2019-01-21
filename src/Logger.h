/*
 * Logger.h
 *
 *  Created on: Oct 22, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#ifndef LOGGER_H_
#define LOGGER_H_

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class Logger {
public:
	Logger(int, string);
	virtual ~Logger();
	static const int LEVEL_DEBUG = 0;
	static const int LEVEL_INFO = 1;
	static const int LEVEL_WARN = 2;
	static const int LEVEL_ERROR = 3;
	static const int LEVEL_FATAL = 4;
	static Logger* getInstance(int level, string log_file);
	void debug(const string& msg);
	void info(const string& msg);
	void running(const string& msg);
	void mpi(const string& msg);
	void warn(const string& msg);
	void error(const string& msg);
	void fatal(const string& msg);
	string get_log_file();
	int get_log_level();
	void safe_run_shell_command(const string& cmd);
	void fragile_run_shell_command(const string& cmd);
private:
	static Logger* logger;
	int log_level;
	string log_file;
	void print_message(const string& msg, const string& level, bool to_std_out);
};

#endif /* LOGGER_H_ */
