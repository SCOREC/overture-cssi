#include "ComplexArray.h"



// ==============================================================
/// \brief Display the matrix
// ==============================================================
void ComplexArray::
display( FILE *file /* =stdout */ ) const 
{
  printF("Entering ComplexArray::display\n");
  const ComplexArray & A = *this;
  for( int i1=0; i1<nd[0]; i1++ )
  {
    fPrintF(file," i=%2d: ",i1);
    for( int i2=0; i2<nd[1]; i2++ )
    {
      Complex a = A(i1,i2);
      fPrintF(file,"[%10.4e %+10.4e I] ", std::real(a),imag(a));
    }
    fPrintF(file,"\n");
  }
    
    
}
