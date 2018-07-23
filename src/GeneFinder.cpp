/*
 * GeneFinder.cpp
 *
 *  Created on: Nov 14, 2012
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#include "GeneFinder.h"
#include "SnapGeneFinder.h"

GeneFinder* GeneFinder::gene_finder = NULL;

GeneFinder::GeneFinder(int log_level, string log_file) {
	logger = Logger::getInstance(log_level, log_file);
}

GeneFinder::~GeneFinder() {
	// Auto-generated destructor stub
}

// singleton implementation
GeneFinder* GeneFinder::getInstance(int type, int log_level, string log_file){
	if (type == GeneFinder::NONE)
		return NULL;
	if (type == GeneFinder::SNAP){
		if (gene_finder == NULL)
			gene_finder = new SnapGeneFinder(log_level, log_file);
		return gene_finder;
	}
	return NULL;
}
