#include "Cgcssi.h"

void Cgcssi::
addForcing(MappedGrid& mg, realMappedGridFunction & dvdt, int iparam[], real rparam[],
		realMappedGridFunction & dvdtImplicit = Overture::nullRealMappedGridFunction())
{
}


// int Cgcssi:: 
// applyBoundaryConditions(const real & t, realMappedGridFunction & u, 
// 			realMappedGridFunction & gridVelocity,
// 			const int & grid,
// 			const int & option=-1,
// 			realMappedGridFunction *puOld=NULL, 
// 			realMappedGridFunction *pGridVelocityOld=NULL,
// 			const real & dt=-1.)
// {
//   return 0;
// }

int Cgcssi::
getUt(const realMappedGridFunction & v, 
	  const realMappedGridFunction & gridVelocity, 
	  realMappedGridFunction & dvdt, 
	  int iparam[], real rparam[],
	  realMappedGridFunction & dvdtImplicit = Overture::nullRealMappedGridFunction(),
	  MappedGrid *pmg2=NULL,
      const realMappedGridFunction *pGridVelocity2= NULL)
{
  return 0;
}



void Cgcssi::
getTimeSteppingEigenvalue(MappedGrid & mg, 
 			       realMappedGridFunction & u, 
 			       realMappedGridFunction & gridVelocity,  
 			       real & reLambda,
 			       real & imLambda, 
			  const int & grid)
{
}
