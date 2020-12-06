#include "Maxwell.h"
#include "Ogshow.h"
#include "HDF_DataBase.h"
#include "display.h"
#include "BodyForce.h"
#include "DispersiveMaterialParameters.h"

// static int restartNumber=-1;

//\begin{>>MaxwellInclude.tex}{\subsection{saveShow}} 
void Maxwell::
saveShow( int current, real t, real dt )
//=========================================================================================
// /Description:
//    Save a solution in the show file. This routine will also save restart files.
//
// /current (input) : save this grid function.
//
//\end{MaxwellInclude.tex}  
//=========================================================================================
{

  if( show==NULL )
    return;

  real cpu0=getCPU();

  const int myid = Communication_Manager::My_Process_Number;

/* ---
  if( parameters.saveRestartFile )
  {
    // keep two restart files, just in case we crash while writing one of them
    restartNumber = (restartNumber+1) % 2;
    saveRestartFile(gf0, restartNumber==0 ? "ob1.restart" : "ob2.restart" );
  }
  --- */
  
  if( debug & 1 )
    printF("saving a solution in the show file...\n");
  
  int & numberSavedToShowFile = parameters.dbase.get<int>("numberSavedToShowFile");

  Ogshow & showFile = *show;
  if( numberSavedToShowFile==-1 )
  {
    saveParametersToShowFile();
  }
  
  assert( cgp!=NULL );
  CompositeGrid & cg= *cgp;

  bool isDispersiveSingleDomain=false;
  bool isDispersiveMultiDomain=false;
  if( cg.numberOfDomains()>1 )
  {
    for( int domain=0; domain<cg.numberOfDomains(); domain++ )
    {
      const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
      // const bool isDispersive = dmp.isDispersiveMaterial() ||  (method==bamx && totalNumberOfPolarizationComponents(grid)>0);
      isDispersiveMultiDomain = isDispersiveMultiDomain || dmp.isDispersiveMaterial();
      if( isDispersiveMultiDomain ) break;
    }
    
    if( isDispersiveMultiDomain )
      printF("\n +++++ saveShow: MULTI-DOMAIN PROBLEM AND DISPERSIVE +++++++\n\n");

  }
  else
  {
    int domain=0;
    const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
    isDispersiveSingleDomain = dmp.isDispersiveMaterial();
    
    if( isDispersiveSingleDomain )
      printF("\n +++++ saveShow: SINGLE-DOMAIN PROBLEM AND DISPERSIVE +++++++\n\n");
  }
  


  showFile.startFrame();


  const int appendToOldShowFile = parameters.dbase.get<int >("appendToOldShowFile");
  if( numberSavedToShowFile==-1 && appendToOldShowFile )
  {
    // -- do not save first solution if we are appending to an existing show file.
    // This solution is already there.
    numberSavedToShowFile=0;
    // There is no need to save the grid if we are appending:
    parameters.dbase.get<int >("saveGridInShowFile")=false;
    return;
  }
  if( numberSavedToShowFile==-1 )
  {
    // first call -- save general parameters
    numberSavedToShowFile=0;
    // parameters.saveParametersToShowFile();
  }
  numberSavedToShowFile++;

  
  HDF_DataBase *dbp=NULL;
  #ifdef OV_USE_HDF5
    bool putToDataBase=true;    // hdf5  -- put on all processors
  #else
    bool putToDataBase= parameters.dbase.get<int >("myid")==0; // hdf4 - only put on processor 0
  #endif

  // printf(" ***** Maxwell:: saveShow: myid=%i, putToDataBase=%i\n",myid,(int)putToDataBase);
    

  if( putToDataBase )
  {
    dbp = showFile.getFrame();
    assert( dbp!=NULL );
    // save parameters that go in this frame
    HDF_DataBase & db = *dbp;
    db.put(t,"time");
    db.put(dt,"dt");

    // --- Save body/boundary force regions in the first frame ---
    if( numberSavedToShowFile==1 )
    {
      if( parameters.dbase.get<bool >("turnOnBoundaryForcing") )
      {
	printF("\n @@@ saveShow: Save material region boundaries to the show file @@@@ \n\n");

	// Here is the array of boundary forcings:
	std::vector<BodyForce*> & boundaryForcings =  parameters.dbase.get<std::vector<BodyForce*> >("boundaryForcings");
	// -- save the number of regions:
	db.put((int)boundaryForcings.size(),"numberOfBoundaryForceRegions");

	for( int bf=0; bf<boundaryForcings.size(); bf++ )
	{
	  const BodyForce & boundaryForce = *boundaryForcings[bf];
	  boundaryForce.put(db,sPrintF("BoundaryForce%i",bf));
	}
      }
      else
      {
	db.put(0,"numberOfBoundaryForceRegions"); // there are no boundary force regions
      }

    }
  }
  
  
  char buffer[80]; 
  aString showFileTitle[5];

  showFileTitle[0]=sPrintF(buffer,"Maxwell %s",(const char *)methodName);
  showFileTitle[1]=sPrintF(buffer,"t=%4.3f, dt=%8.2e",t,dt);
  showFileTitle[2]="";  // marks end of titles

  for( int i=0; showFileTitle[i]!=""; i++ )
    showFile.saveComment(i,showFileTitle[i]);


  if( mgp==NULL )
  {
    // save a CompositeGridFunction...
    if( t<=0. )
    {
      // ** CHECK ME **
      getErrors( current, t, dt ); // *wdh* Jan 21, 2019 -- make sure errors are computed at t=0
    }
    
    
    realCompositeGridFunction v;
    realCompositeGridFunction & u = getAugmentedSolution(current,v,t);

    // realCompositeGridFunction & u = getCGField(HField,current);

    if( saveGridInShowFile )
    {
      if( debug & 4 ) 
         printF("***Save grid in the show file: numberOfComponentGrids=%i\n",cg.numberOfComponentGrids());
      showFile.saveSolution( u,"u",Ogshow::useCurrentFrame );  // save the grid and the grid function
      saveGridInShowFile=false;
      showFileFrameForGrid=showFile.getNumberOfFrames();
    }
    else
    {
      if( debug & 4 )
        printF("***Save solution in the show file: showFileFrameForGrid=%i (-2=default)\n",
	     showFileFrameForGrid);
      showFile.saveSolution( u,"u",showFileFrameForGrid );  // save the grid function
    }
  }
  else
  {
    // Save a MappedGridFunction...
    realMappedGridFunction & u = fields[current];
    if( saveGridInShowFile )
    {
      if( debug & 4 )
        printF("***Save grid in the show file: numberOfComponentGrids=%i\n",cg.numberOfComponentGrids());
      showFile.saveSolution( u,"u",Ogshow::useCurrentFrame );  // save the grid and the grid function
      saveGridInShowFile=false;
      showFileFrameForGrid=showFile.getNumberOfFrames();
    }
    else
    {
      if( debug & 4 )
        printF("***Save solution in the show file: showFileFrameForGrid=%i (-2=default)\n",
	     showFileFrameForGrid);
      showFile.saveSolution( u,"u",showFileFrameForGrid );  // save the grid function
    }
  }
  
  // Here we save time sequences to the show file
  // Only save if this is the last frame in a subFile
  if( dbase.get<bool >("saveSequencesEveryTime") && showFile.isLastFrameInSubFile() )
  {  
//    #ifndef USE_PPP
    // fix me for parallel -- adding this causes the program to hang when closing the show file.
    saveSequencesToShowFile();
//    #endif
  }

  // *wdh* 2019/12/07 Comment this out to avoid error at end when saving sequences: (not found in common/src/saveShow.bC)
  // showFile.endFrame();  

  timing(timeForShowFile)+=getCPU()-cpu0;
}


int Maxwell::
saveParametersToShowFile()
// =================================================================================================
// /Description:
//     Save PDE specific parameters in the show file.
//     These parameters can be used for a restart. They can also be used, for example,
//     by the user defined derived functions (when viewing the show file with plotStuff).
// 
//\end{OB_ParametersInclude.tex}  
// =================================================================================================
{
  assert( show!=NULL );

  ListOfShowFileParameters showFileParams;

  // save parameters
  showFileParams.push_back(ShowFileParameter("Maxwell's Equations","pde"));
    
  // printF("\n ****** saveParametersToShowFile: ez=%d\n\n",ez);
  // OV_ABORT("stop here");
  
  showFileParams.push_back(ShowFileParameter("exFieldComponent",ex));
  showFileParams.push_back(ShowFileParameter("eyFieldComponent",ey));
  showFileParams.push_back(ShowFileParameter("ezFieldComponent",ez));

  showFileParams.push_back(ShowFileParameter("hxFieldComponent",hx));
  showFileParams.push_back(ShowFileParameter("hyFieldComponent",hy));
  showFileParams.push_back(ShowFileParameter("hzFieldComponent",hz));

  showFileParams.push_back(ShowFileParameter("eps",eps));
  showFileParams.push_back(ShowFileParameter("mu",mu));

  show->saveGeneralParameters(showFileParams);
    

  return 0;
}


//\begin{>>MaxwellInclude.tex}{\subsection{saveSequenceInfo}} 
int Maxwell::
saveSequenceInfo( real t0, RealArray & sequenceData )
//=========================================================================================
// /Description:
//    Save info into the time history arrays.
// 
//\end{MaxwellInclude.tex}  
//=========================================================================================
{
  if( show==NULL )
    return 0;

  if( sequenceCount >= timeSequence.getLength(0) )
  {
    int num=timeSequence.getLength(0);
    Range R(0,num-1),all;
    RealArray seq;  seq=sequence;
    num=int(num*1.5+100);
    timeSequence.resize(num);
    sequence.redim(num,numberOfSequences);
    sequence(R,all)=seq;
  }

  timeSequence(sequenceCount)=t0;
  for( int n=sequenceData.getBase(0); n<=sequenceData.getBound(0); n++ )
    sequence(sequenceCount,n)=sequenceData(n); 

  sequenceCount++;
  return 0;
}



//\begin{>>MaxwellInclude.tex}{\subsection{saveSequencesToShowFile}} 
int Maxwell::
saveSequencesToShowFile()
//=========================================================================================
// /Description:
//
//\end{CompositeGridSolverInclude.tex}  
//=========================================================================================
{
  if( show==NULL || sequenceCount<=0 )
    return 0;
  
  assert( cgp!=NULL );
  CompositeGrid & cg= *cgp;
  const int numberOfDimensions = cg.numberOfDimensions();
 
  Range I(0,sequenceCount-1);
  Range N=sequence.dimension(1);
  
  // Is this next line correct?
  // int numberOfComponents = method==nfdtd ? hz-ex+1 : max(ey,ez) + hz + 1 -ex +1;
  int numberOfComponents=0;
  if( method==nfdtd )
  {
    numberOfComponents = 3;
  }
  else if( method==sosup )
  {
    // 2D: ex,ey,hz, ext,eyt,hzt  
    // 3D: ex,ey,ez, ext,eyt,ezt  
    numberOfComponents = 6; 
  }
  else if( method==bamx )
  {
    const int & solveForAllFields = dbase.get<int>("solveForAllFields");
    numberOfComponents = (numberOfDimensions ==3 || solveForAllFields==1) ? 6 : 3; 
  }
  else if( method==yee )
  {
    numberOfComponents = numberOfDimensions==2 ? 3 : 6;
  }
  else
  {
    numberOfComponents = max(ey,ez) + hz + 1 -ex +1;
  }
  
  aString *name = new aString [numberOfSequences]; 
  for( int n=0; n<numberOfComponents; n++ )
  {
    if( cgerrp!=NULL )
    {
      name[n]=cgerrp[0].getName(n);
    }
    else
    {
      name[n]=sPrintF("error%i",n);
    }
    for( int i=0; i<name[n].length(); i++ )
    { // change blanks to underscores (for matlab)
      if( name[n][i]==' ' )
      {
	name[n][i]='_';
      }
    }
  }


  if( computeEnergy && numberOfSequences>numberOfComponents )
  {
    name[numberOfComponents  ]="Energy";
    name[numberOfComponents+1]="Delta_E";
  }
  
  // display(sequence(I,N),"saveSequencesToShowFile: sequence(I,N)");
  
  if( false )
  {
    printf("saveSequencesToShowFile() myid=%i sequenceCount=%i\n",myid,sequenceCount);
    fflush(0);
  }
  

  // NOTE: This function must be called by ALL processors in parallel
  show->saveSequence("errors",timeSequence(I),sequence(I,N),name);
  
  delete [] name;
  
  return 0;
}
