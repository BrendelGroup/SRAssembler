/*
 * SRAssemblerSlave.cpp
 *
 *  Created on: Oct 13, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#include "SRAssemblerSlave.h"

SRAssemblerSlave::SRAssemblerSlave() {
	// Auto-generated constructor stub
}

int SRAssemblerSlave::init(int argc, char * argv[], int rank, int mpiSize) {
	return SRAssembler::init(argc, argv, rank, mpiSize);
}

void SRAssemblerSlave::print_message(const string& msg){
}

void SRAssemblerSlave::show_usage(){

}

void SRAssemblerSlave::do_preprocessing(){
	// SRAssembler Slaves move on to do_walking and then process messages until they get an exit code.
}

void SRAssemblerSlave::remove_taboo_reads(){
	// SRAssembler Slaves move on to do_walking and then process messages until they get an exit code.
}

void SRAssemblerSlave::do_walking(){
	process_message();
}

void SRAssemblerSlave::process_message(){
	int action;
	long long code_value;
	int value1;
	int value2;
	int value3;
	int source;
	mpi_code code;
	while(1){
		mpi_receive(code_value, source);
		code = get_mpi_code(code_value);
		action = code.action;
		value1 = code.value1;
		value2 = code.value2;
		value3 = code.value3;
		if (action == ACTION_EXIT){
			// Slave stops waiting for messages and process finishes.
			break;
		}
		if (action == ACTION_TOTAL_CHUNKS){
			// Slave tells its Libraries how many chunks they have.
			this->libraries[value1].set_num_chunks(value2);
		}
		if (action == ACTION_SPLIT){
			// Slave manages library file splitting process.
			if (value2 == 1)
				this->libraries[value1].do_split_files(LEFT_READ, this->reads_per_file);
			if (value2 == 2)
				this->libraries[value1].do_split_files(RIGHT_READ, this->reads_per_file);
			send_code(source, ACTION_RETURN, 0, 0, 0);
		}
		if (action == ACTION_PRE_PROCESSING){
			// Slave manages indexing of reads in one library chunk file.
			SRAssembler::preprocess_read_chunk(value1, value2);
			send_code(source, ACTION_RETURN, 0, value2, 0);
		}
		if (action == ACTION_ASSEMBLY){
			// Slave manages assembly with one kmer value.
			do_assembly(value1, value2, value3);
			send_code(source, ACTION_RETURN, 0, 0, 0);
		}
		if (action == ACTION_ALIGNMENT){
			// Slave manages alignment.
			int found_new_reads = do_alignment(value1, value3, value2);
			send_code(source, ACTION_RETURN, value2, found_new_reads, 0);
		}
		if (action == ACTION_LOAD_PREVIOUS){
			// If continuing a previous SRAssembler run, the list of previously found reads needs to be loaded to each Slave so they don't capture them again.
			load_found_reads(value1);
			send_code(source, ACTION_RETURN, 0, 0, 0);
		}
		if (action == ACTION_MEMDIR){
			// Slaves need to know where the shared directory for storing temporary files is.
			this->tmp_dir = this->tmp_loc + "/SRAssemblertemp" + int2str(value2);
		}
		if (action == ACTION_SAVE){
			// Save any reads found by this Slave in case this SRAssembler run is started again.
			save_found_reads(value1);
			send_code(source, ACTION_RETURN, 0, 0, 0);
		}
		if (action == ACTION_CLEAN){
			// Slave manages aligning of found reads pool to contigs kept after cleaning, and retains only the hit reads.
			SRAssembler::remove_unmapped_reads(value1, value2);
			send_code(source, ACTION_RETURN, value1, 0, 0);
		}
	}
}

SRAssemblerSlave::~SRAssemblerSlave() {
	// Auto-generated destructor stub
}
