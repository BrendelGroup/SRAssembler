/*
 * SplicedAligner.cpp
 *
 *  Created on: Oct 15, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#include "SplicedAligner.h"
#include "GTHAligner.h"
#include "GSQAligner.h"

SplicedAligner* SplicedAligner::spliced_aligner = NULL;

SplicedAligner::SplicedAligner(int log_level, string log_file) {
	logger = Logger::getInstance(log_level, log_file);

}

SplicedAligner::~SplicedAligner() {
	// Auto-generated destructor stub
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
// singleton implementation
SplicedAligner* SplicedAligner::getInstance(int type, int log_level, string log_file){
	if (type == SplicedAligner::GENOMETHREADER){
		if (spliced_aligner == NULL)
			spliced_aligner = new GTHAligner(log_level, log_file);
		return spliced_aligner;
	}
	if (type == SplicedAligner::GENESEQER){
		if (spliced_aligner == NULL)
			spliced_aligner = new GSQAligner(log_level, log_file);
		return spliced_aligner;
	}
	return NULL;
}
