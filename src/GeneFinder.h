/*
 * GeneFinder.h
 *
 *  Created on: Nov 14, 2012
 *      Author: hchou
 */

#ifndef GENEFINDER_H_
#define GENEFINDER_H_

#include <string>
#include "Utility.h"
#include "Const.h"
#include "Logger.h"

using namespace std;

class GeneFinder {
public:
	GeneFinder(int, string);
	virtual ~GeneFinder();
	virtual void do_gene_finding(const string& genomic_file, const string& species, const Params& params, const string& output_file)=0;
	virtual string get_output_summary()=0;
	virtual string get_program_name()=0;
	virtual bool is_available()=0;
	static const int NONE = 0;
	static const int SNAP = 1;
	static GeneFinder* getInstance(int, int, string);
	void set_log_file(std::string log_file);
protected:
	Logger* logger;
	string output_string;
private:
	static GeneFinder* gene_finder;
};

#endif /* GENEFINDER_H_ */
