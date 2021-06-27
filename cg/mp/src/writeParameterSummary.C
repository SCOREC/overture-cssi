#include "Cgmp.h"

// ===================================================================================================================
/// \brief Output run-time parameters for the header.
/// \param file (input) : write output to this file.
///
// ===================================================================================================================
void Cgmp::
writeParameterSummary( FILE * file )
{

  const Parameters::TimeSteppingMethod & timeSteppingMethod= parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod");

  MpParameters::MultiDomainAlgorithmEnum multiDomainAlgorithm=
                 parameters.dbase.get<MpParameters::MultiDomainAlgorithmEnum>("multiDomainAlgorithm");
  
  const Parameters::InterfaceCommunicationModeEnum & interfaceCommunicationMode= 
    parameters.dbase.get<Parameters::InterfaceCommunicationModeEnum>("interfaceCommunicationMode");
  const int & interfaceProjectionGhostOption = parameters.dbase.get<int>("interfaceProjectionGhostOption");
  const bool & hasHeatFluxInterfaces = parameters.dbase.get<bool>("hasHeatFluxInterfaces");
  const bool & hasTractionInterfaces = parameters.dbase.get<bool>("hasTractionInterfaces");
  
    
  fPrintF(file,"\n"
          "******************************************************************\n");
  fPrintF(file,
          "             %s Version 0.2                                 \n"
          "             -----------------                              \n",
          (const char*)getClassName()   );
  

  fPrintF(file,"\n"
          " cfl = %f, tFinal=%e, tPrint = %e \n"
          " time stepping method: %s.\n"
          " hasHeatFluxInterfaces=%d, hasTractionInterfaces=%d\n"
          " number of PC corrections = %i.\n"
          " solve coupled interface equations = %i. ( =1: for explicit time-stepping, solve coupled jump conditions for T).\n"
          " use %s interface transfer. (useNewInterfaceTransfer=%d)\n"
          " relax correction steps = %i.\n"
          " multi-domain algorithm = %s.\n"
          " interface communication mode= %s,\n"
          " project interface = %s. (interfaceProjectionOption=%i, interface-ghost=%s)\n"
          " project initial conditions = %i\n"
          ,
          parameters.dbase.get<real >("cfl"),
          parameters.dbase.get<real >("tFinal"),
          parameters.dbase.get<real >("tPrint"),
          (const char*)Parameters::timeSteppingName[timeSteppingMethod],
          (int)hasHeatFluxInterfaces,(int)hasTractionInterfaces,
          (int)parameters.dbase.get<int>("numberOfPCcorrections"),
          (int)parameters.dbase.get<bool>("solveCoupledInterfaceEquations"),
          (parameters.dbase.get<bool>("useNewInterfaceTransfer") ? "new" : "old"),
          (int)parameters.dbase.get<bool>("useNewInterfaceTransfer"),
          (int)parameters.dbase.get<bool>("relaxCorrectionSteps"),
          (multiDomainAlgorithm==MpParameters::defaultMultiDomainAlgorithm          ? "default" :
           multiDomainAlgorithm==MpParameters::stepAllThenMatchMultiDomainAlgorithm ? "step all then match" : 
           multiDomainAlgorithm==MpParameters::multiStageAlgorithm                  ? "multi-stage" : 
                                                                                      "unknown"),
          (interfaceCommunicationMode==Parameters::autoRequestInterfaceData       ? "autoRequestInterfaceData" :
           interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded ? "requestInterfaceDataWhenNeeded" :
                                                                                    "unknown"),
          (parameters.dbase.get<bool>("projectInterface")? "true" : "false"),
          parameters.dbase.get<int>("interfaceProjectionOption"),
          (interfaceProjectionGhostOption==0 ? "extrapolate" : 
           interfaceProjectionGhostOption==1 ? "compatibility" : 
           interfaceProjectionGhostOption==2 ? "exact" : "domain BC" ),
           parameters.dbase.get<bool>("projectMultiDomainInitialConditions")
           );


  if( parameters.dbase.get<bool >("twilightZoneFlow") )
  {
    fPrintF(file," Twilight zone flow\n");
  }

  InterfaceList & interfaceList = parameters.dbase.get<InterfaceList>("interfaceList");
  for( int inter=0; inter < interfaceList.size(); inter++ )
  {
    InterfaceDescriptor & interfaceDescriptor = interfaceList[inter]; 

  
    fPrintF(file,"\n  ------ Interface %i is an interface between domain1=%i (%s,%s) and domain2=%i (%s,%s) ------- \n",
        inter,interfaceDescriptor.domain1,
        (const char*)domainSolver[interfaceDescriptor.domain1]->getName(),
        (const char*)domainSolver[interfaceDescriptor.domain1]->getClassName(),
        interfaceDescriptor.domain2,
        (const char*)domainSolver[interfaceDescriptor.domain2]->getClassName(),
        (const char*)domainSolver[interfaceDescriptor.domain2]->getName());

    // there may be multiple grid faces that lie on the interface:     
    for( int face=0; face<interfaceDescriptor.gridListSide1.size(); face++ )
    {
      GridFaceDescriptor & gridDescriptor1 = interfaceDescriptor.gridListSide1[face];
      GridFaceDescriptor & gridDescriptor2 = interfaceDescriptor.gridListSide2[face];
      
      const int d1=gridDescriptor1.domain, grid1=gridDescriptor1.grid, side1=gridDescriptor1.side, dir1=gridDescriptor1.axis;
      const int d2=gridDescriptor2.domain, grid2=gridDescriptor2.grid, side2=gridDescriptor2.side, dir2=gridDescriptor2.axis;
      

      if( face==0 )
      {
        const IntegerArray & interfaceType1 = domainSolver[d1]->parameters.dbase.get<IntegerArray >("interfaceType");
        const IntegerArray & interfaceType2 = domainSolver[d2]->parameters.dbase.get<IntegerArray >("interfaceType");
        assert( interfaceType1(side1,dir1,grid1) == interfaceType2(side2,dir2,grid2) );

        if( interfaceType1(side1,dir1,grid1)==Parameters::heatFluxInterface )
        {
          fPrintF(file,"  InterfaceType = heatFluxInterface\n");
        }
        else if( interfaceType1(side1,dir1,grid1)==Parameters::tractionInterface )
        {
          fPrintF(file,"  InterfaceType = tractionInterface\n");
        }

        // -- face info ---
        fPrintF(file,"  face=%d: [d1,grid1,side1,dir1]=[%d,%d,%d,%d] [d2,grid2,side2,dir2]=[%d,%d,%d,%d] \n",face,d1,grid1,side1,dir1,d2,grid2,side2,dir2);
      }
    } // end for face

    fPrintF(file,
        "  interface solve: tol=%9.3e, omega=%9.3e, max iterations=%i, relaxCorrectionSteps=%i\n",
        interfaceDescriptor.interfaceTolerance,
        interfaceDescriptor.interfaceOmega,
        interfaceDescriptor.maximumNumberOfIntefaceIterations,
        (int)parameters.dbase.get<bool>("relaxCorrectionSteps")
       );

  } // end for interface
  fPrintF(file,"\n See interfaceFile.log for more interface information.\n");


  fPrintF(file,"******************************************************************\n\n");
  

}