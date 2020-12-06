#ifndef CHECK_PARALLEL_H
#define CHECK_PARALLEL_H

// ===================================================================================================
// Class to check errors in parallel codes: compare results from 1 processor to multiple processors
// ==================================================================================================
#include "Overture.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;


class CheckParallel
{
public:

CheckParallel();
~CheckParallel();


// check the difference in an array between two parallel runs
real checkDiff( const RealArray & x, const aString & label, FILE *file = NULL );
real checkDiff( const IntegerArray & x, const aString & label, FILE *file = NULL );

// set the default file for debug info
int setDefaultFile( FILE *file );

protected:

int counter;

// The database is the new place to store parameters
mutable DataBase dbase; 

};


// Use these functions to call checkParallel from Fortran
real checkParallelArrayReal( CheckParallel *pcp, char *label_,
                         int & nd1a, int & nd1b, int & nd2a, int & nd2b, int & nd3a, int & nd3b, int & nd4a, int & nd4b, 
                         int & n1a, int & n1b, int & n2a, int & n2b, int & n3a, int & n3b, int & n4a, int & n4b,
                         real & x_, 
                         int & labelLength );

real checkParallelArrayInt( CheckParallel *pcp, char *label_,
                         int & nd1a, int & nd1b, int & nd2a, int & nd2b, int & nd3a, int & nd3b, int & nd4a, int & nd4b, 
                         int & n1a, int & n1b, int & n2a, int & n2b, int & n3a, int & n3b, int & n4a, int & n4b,
                         int & x_, 
                         int & labelLength );

#endif
  
