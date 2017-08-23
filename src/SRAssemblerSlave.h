/*
 * SRAssemblerSlave.h
 *
 *  Created on: Oct 13, 2011
 *      Author: hchou
 */

#ifndef SRASSEMBLERSLAVE_H_
#define SRASSEMBLERSLAVE_H_

#include "SRAssembler.h"

class SRAssemblerSlave: public SRAssembler {
public:
	SRAssemblerSlave();
	int init(int argc, char * argv[], int rank, int size);
	void show_usage();
	void print_message(const string&);
	void do_preprocessing();
	void do_walking();
	virtual ~SRAssemblerSlave();
private:
	void process_message();
};

#endif /* SRASSEMBLERSlave_H_ */
