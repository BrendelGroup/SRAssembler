/*
 * Aligner.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#include "Aligner.h"
#include "VmatchBWAfastmapAligner.h"

Aligner* Aligner::protein_aligner=NULL;
Aligner* Aligner::dna_aligner=NULL;


// TODO Auto-generated constructor stub
Aligner::Aligner(int log_level, string log_file) {
	logger = Logger::getInstance(log_level, log_file);
}

Aligner::~Aligner() {
	// TODO Auto-generated destructor stub
}

//singleton implementation
Aligner* Aligner::getInstance(int type, int log_level, string log_file){
	if (type == PROTEIN_ALIGNER){
		if (protein_aligner == NULL)
			protein_aligner = new VmatchBWAfastmapAligner(log_level, log_file);
		return protein_aligner;
	}
	if (type == DNA_ALIGNER){
		if (dna_aligner == NULL)
			dna_aligner = new VmatchBWAfastmapAligner(log_level, log_file);
		return dna_aligner;
	}
	return NULL;
}
