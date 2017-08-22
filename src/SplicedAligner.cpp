/*
 * SplicedAligner.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: hchou
 */

#include "SplicedAligner.h"
#include "GTHAligner.h"
#include "GSQAligner.h"
#include "ExonerateAligner.h"

SplicedAligner* SplicedAligner::spliced_aligner = NULL;

SplicedAligner::SplicedAligner(int log_level, string log_file) {
	logger = Logger::getInstance(log_level, log_file);

}

SplicedAligner::~SplicedAligner() {
	// TODO Auto-generated destructor stub
}

string_map SplicedAligner::get_aligned_query_list(){
	return aligned_query_list;
}
int SplicedAligner::get_match_num(){
	return this->num_matches;
}
string SplicedAligner::get_output_summary(){
	return this->output_string;
}
//singleton implementation
SplicedAligner* SplicedAligner::getInstance(int type, int log_level, string log_file){
	if (type == SplicedAligner::GENOME_THREADER){
		if (spliced_aligner == NULL)
			spliced_aligner = new GTHAligner(log_level, log_file);
		return spliced_aligner;
	}
	if (type == SplicedAligner::GENE_SEQER){
		if (spliced_aligner == NULL)
			spliced_aligner = new GSQAligner(log_level, log_file);
		return spliced_aligner;
	}
	if (type == SplicedAligner::EXONERATE){
		if (spliced_aligner == NULL)
			spliced_aligner = new ExonerateAligner(log_level, log_file);
		return spliced_aligner;
	}
	return NULL;
}
