/*
 * Assembler.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#include "Assembler.h"
#include "AbyssAssembler.h"
#include "SOAPDenovoAssembler.h"

Assembler* Assembler::assembler = NULL;

Assembler::Assembler(int log_level, string log_file) {
	logger = Logger::getInstance(log_level, log_file);

}

Assembler::~Assembler() {
	// TODO Auto-generated destructor stub
}


//singleton implementation
Assembler* Assembler::getInstance(int type, int log_level, string log_file){
	if (type == ABYSS_ASSEMBLER){
		if (assembler == NULL)
			assembler = new AbyssAssembler(log_level, log_file);
		return assembler;
	}
	if (type == SOAPDENOVO_ASSEMBLER){
		if (assembler == NULL)
			assembler = new SOAPDenovoAssembler(log_level, log_file);
		return assembler;
	}
	return NULL;
}
