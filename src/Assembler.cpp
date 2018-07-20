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
	// Auto-generated destructor stub
}

// singleton implementation
Assembler* Assembler::getInstance(int assembler_type, int log_level, string log_file){
	if (assembler_type == ABYSS_ASSEMBLER){
		if (assembler == NULL)
			assembler = new AbyssAssembler(log_level, log_file);
		return assembler;
	}
	if (assembler_type == SOAPDENOVO_ASSEMBLER){
		if (assembler == NULL)
			assembler = new SOAPDenovoAssembler(log_level, log_file);
		return assembler;
	}
	return NULL;
}
