/*
 * MPIWrapper.cpp
 *
 *  Created on: Oct 16, 2011
 *      Author: hchou
 */

#include "MPIWrapper.h"
#include <string.h>
#ifdef MPIMODE
#include <mpi.h>
#endif

unsigned long get_mpi_code_value(mpi_code code){
	return 100000000 * code.action + 1000000 * code.value1 + 100* code.value2 + code.value3;
}

mpi_code get_mpi_code(int code_value){
	mpi_code code;
	code.action = code_value / 100000000;
	code.value1 = code_value / 1000000 % 100;
	code.value2 = code_value / 100 % 10000;;
	code.value3 = code_value % 100;
	return code;
}

void mpi_send( const char* msg, const int& to ) {
#ifdef MPIMODE
	int len = strlen( ( char* )msg );
	MPI::COMM_WORLD.Send( &len, 1,   MPI::INT,  to, ( int )1 );
	MPI::COMM_WORLD.Send(  msg, len, MPI::CHAR, to, ( int )1 );
#endif
}

void mpi_receive( char* msg, int& from ) {
#ifdef MPIMODE
	int len;
	MPI::Status status;

	MPI::COMM_WORLD.Recv( &len, 1, MPI::INT, MPI::ANY_SOURCE,MPI::ANY_TAG,  status );
	from = ( ( int )status.Get_source(  ) );
	MPI::COMM_WORLD.Recv(  msg, len, MPI::CHAR, from, MPI::ANY_TAG );
	msg[len] = '\0';
#endif
}
void mpi_send( const int& code, const int& to ) {
#ifdef MPIMODE
	MPI::COMM_WORLD.Send( &code, 1,   MPI::INT,  to, ( int )1 );
#endif
}
void mpi_bcast(const int& code) {
	int mpiSize = mpi_get_size();
	int rank = mpi_get_rank();
	for (int i=0;i<mpiSize;i++)
		if (rank != i)
			mpi_send(code, i);
	//MPI::COMM_WORLD.Bcast( &code, 1,   MPI::INT, 0);
}
void mpi_receive( int& code, int& from ) {
#ifdef MPIMODE
	MPI::Status status;

	MPI::COMM_WORLD.Recv( &code, 1, MPI::INT, MPI::ANY_SOURCE,MPI::ANY_TAG,  status );
	from = ( ( int )status.Get_source(  ) );
#endif
}
void mpi_init(int argc, char * argv[] ) {
#ifdef MPIMODE
	MPI::Init(argc,argv);
#endif
}
int mpi_get_rank(){
#ifdef MPIMODE
	return MPI::COMM_WORLD.Get_rank();
#else
	return 0;
#endif
}
int mpi_get_size(){
#ifdef MPIMODE
	return MPI::COMM_WORLD.Get_size();
#else
	return 1;
#endif
}
void mpi_finalize(){
#ifdef MPIMODE
	MPI::Finalize();
#endif
}
