#ifndef COMPLEX_ARRAY
#define COMPLEX_ARRAY

//  Class for a simple 2x2 array of complex numbers

#include "Overture.h"
#include <complex>
typedef ::real LocalReal;

typedef std::complex<LocalReal> MyComplex;

// ============================================================================================
// This Complex class is used to hide std::complex from Overture header files 
// ============================================================================================
class Complex : public std::complex<LocalReal> 
{
public:  

Complex & operator=( const LocalReal & x ){ (MyComplex&)(*this) = x; return *this;}   // 
Complex & operator=( const std::complex<LocalReal> & x ){ (MyComplex&)(*this) =x; return *this;}   // 
};


// ============================================================================================
// Class for a simple 2x2 array of complex numbers
// ============================================================================================
class ComplexArray
{
public:

// constructor: 
ComplexArray( int n1, int n2 ){ nd[0]=n1; nd[1]=n2; pA = new Complex [n1*n2]; } // 

// destructor:
~ComplexArray(){ delete [] pA; }

// Copy constructor (deep copy)
ComplexArray( const ComplexArray & x ) // copy constructor
  {
    *this = x;
  }

// operator = 
ComplexArray & operator=( const ComplexArray & x )
  {
    if( nd[0]!=x.nd[0] || nd[1] != x.nd[1] )
      redim(x.nd[0],x.nd[1]);
    for( int i=0; i<nd[0]*nd[1]; i++ )
      pA[i]=x.pA[i];

    return *this;
    
    }

// return A(i1,i2) (lvalue for assignment to)
inline Complex & operator()(int i1, int i2 ){ return pA[(i1)+nd[0]*(i2)]; }    // lvalue 

// return A(i1,i2) (rvalue for use on rhs)
inline Complex operator()(int i1, int i2 ) const { return pA[(i1)+nd[0]*(i2)]; } // rvalue 


void display( FILE *file=stdout ) const;


// redimension the array
void redim( int n1, int n2 ){ delete [] pA; nd[0]=n1; nd[1]=n2;  pA = new Complex [n1*n2]; } //

protected:

// data 
int nd[2];
Complex *pA;

};

  
#endif
