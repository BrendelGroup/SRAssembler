/*
 * SRAssemblerSlave.h
 *
 *  Created on: Oct 13, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#ifndef SRASSEMBLERSLAVE_H_
#define SRASSEMBLERSLAVE_H_

#include "SRAssembler.h"

class SRAssemblerSlave: public SRAssembler {
public:
	SRAssemblerSlave();
	int init(int argc, char * argv[], int rank, int mpiSize);
	void show_usage();
	void print_message(const string&);
	void do_preprocessing();
	void remove_taboo_reads();
	void do_walking();
	virtual ~SRAssemblerSlave();
private:
	void process_message();
};

#endif /* SRASSEMBLERSlave_H_ */
