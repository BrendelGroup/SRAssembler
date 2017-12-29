/*
 * MPIWrapper.h
 *
 *  Created on: Oct 16, 2011
 *      Author: hchou
 */

#ifndef MPIWRAPPER_H_
#define MPIWRAPPER_H_

#include <string>

typedef struct m_code {
		int action;
		int value1;
		int value2;
		int value3;
} mpi_code;

void mpi_send(const char* msg, const int& to );
void mpi_receive(char* msg, int& from );
void mpi_send( const int& code, const int& to );
void mpi_bcast(const int& code);
void mpi_receive( int& code, int& from );
unsigned long get_mpi_code_value(mpi_code code);
mpi_code get_mpi_code(int code_value);
void mpi_init(int argc, char * argv[] );
int mpi_get_rank();
int mpi_get_size();
void mpi_finalize();
#endif /* MPIWRAPPER_H_ */
