#include "Overture.h"
#include "display.h"

int
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

  aString nameOfOGFile;
  cout << "Enter the name of the overlapping grid data base file " << endl;
  cin >> nameOfOGFile;
  
  // create and read in a CompositeGrid
  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile);
  cg.update(MappedGrid::THEvertex | MappedGrid::THEvertexBoundaryNormal | MappedGrid::THEinverseVertexDerivative );   

  Index I1,I2,I3;                         // A++ Index object
  
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )  // loop over component grids
  {
    MappedGrid & mg = cg[grid];
    printF("------------------ GRID %i (%s) ------------------\n",grid,(const char*)mg.getName());
    
    getIndex(cg[grid].indexRange(),I1,I2,I3);                    // assign I1,I2,I3

    const realArray & x = mg.vertex();
    const realArray & rx = mg.inverseVertexDerivative();

    ::display(x,"vertex","%.6e ");
    ::display(rx,"rx","%.6e ");
  }    
  
  Overture::finish();          
  return 0;  
}
