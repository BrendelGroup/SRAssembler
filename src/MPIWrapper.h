/*
 * MPIWrapper.h
 *
 *  Created on: Oct 16, 2011
 *     Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#ifndef MPIWRAPPER_H_
#define MPIWRAPPER_H_

#include <string>

typedef struct m_code {
		int action;
		long value1;
		long value2;
		int value3;
} mpi_code;

void mpi_send(const char* msg, const int& to );
void mpi_receive(char* msg, int& source );
void mpi_send( const long long& code_value, const int& to );
void mpi_bcast(const long long& code_value);
void mpi_receive( long long& code_value, int& source );
long long get_mpi_code_value(mpi_code code);
mpi_code get_mpi_code(long long code_value);
void mpi_init(int argc, char * argv[] );
int mpi_get_rank();
int mpi_get_size();
void mpi_finalize();
#endif /* MPIWRAPPER_H_ */
