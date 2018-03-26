/*
 * SRAssemblerSlave.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: hchou
 */

#include "SRAssemblerSlave.h"

SRAssemblerSlave::SRAssemblerSlave() {
	// TODO Auto-generated constructor stub

}

int SRAssemblerSlave::init(int argc, char * argv[], int rank, int mpiSize) {
	return SRAssembler::init(argc, argv, rank, mpiSize);
}

void SRAssemblerSlave::print_message(const string& msg){
}

void SRAssemblerSlave::show_usage(){

}

void SRAssemblerSlave::do_preprocessing(){
	//process_message();
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
	int from;
	mpi_code code;
	while(1){
		mpi_receive(code_value, from);
		code = get_mpi_code(code_value);
		action = code.action;
		value1 = code.value1;
		value2 = code.value2;
		value3 = code.value3;
		if (action == ACTION_EXIT){    //action terminated
			break;
		}
		if (action == ACTION_TOTAL_PARTS){    //action terminated
			this->libraries[value1].set_num_parts(value2);
		}
		if (action == ACTION_SPLIT){  //do_split
			if (value2 == 1)
				this->libraries[value1].do_split_files(LEFT_READ, this->reads_per_file);
			if (value2 == 2)
				this->libraries[value1].do_split_files(RIGHT_READ, this->reads_per_file);
			send_code(from, ACTION_RETURN, 0, 0, 0);
		}
		if (action == ACTION_PRE_PROCESSING){  //do_split
			SRAssembler::preprocess_read_part(value1, value2);
			send_code(from, ACTION_RETURN, 0, value2, 0);
		}
		if (action == ACTION_ASSEMBLY){  //do_assmebly
			do_assembly(value1, value2, value3);
			send_code(from, ACTION_RETURN, 0, 0, 0);
		}
		if (action == ACTION_ALIGNMENT){  //do alignment
			int found_new_reads = do_alignment(value1, value3, value2);
			send_code(from, ACTION_RETURN, value2, found_new_reads, 0);
		}
		if (action == ACTION_LOAD_PREVIOUS){  //load previous reads
			load_found_reads(value1);
			send_code(from, ACTION_RETURN, 0, 0, 0);
		}
		if (action == ACTION_MEMDIR) {
			this->mem_dir="/dev/shm/SRAssembler" + int2str(value2);
		}
		if (action == ACTION_SAVE){  //load previous reads
			save_found_reads(value1);
			send_code(from, ACTION_RETURN, 0, 0, 0);
		}
		if (action == ACTION_CLEAN){
			SRAssembler::remove_unmapped_reads(value1, value2);
			send_code(from, ACTION_RETURN, value1, 0, 0);
		}
	}
}

SRAssemblerSlave::~SRAssemblerSlave() {
	// TODO Auto-generated destructor stub
}
