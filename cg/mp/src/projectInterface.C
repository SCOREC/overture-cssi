#include "Cgmp.h"
#include "AdvanceOptions.h"
#include "MpParameters.h"

// -- project the temperature at interfaces
int Cgmp::
projectInterface( real t, real dt, std::vector<int> & gfIndex )
{

  const int projectInterfaceTemperature = parameters.dbase.get<bool>("projectInterfaceTemperature");
  if( !projectInterfaceTemperature )
    return 0;

  ForDomain(d)
  {
    // printF("  domain %d: %s\n",d,(const char*)domainSolver[d]->getClassName());
    if( domainSolver[d] )
      domainSolver[d]->projectInterface(  t, dt, domainSolver[d]->gf[gfIndex[d]] );

  }

  return 0;
}