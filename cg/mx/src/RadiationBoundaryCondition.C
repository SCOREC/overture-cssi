#include "RadiationBoundaryCondition.h"
#include "RadiationKernel.h"
#include "OGFunction.h"
#include "display.h"
#include "ParallelUtility.h"

#define exmax EXTERN_C_NAME(exmax)
#define wpulse EXTERN_C_NAME(wpulse)
#define wdpulse EXTERN_C_NAME(wdpulse)
#define radEval EXTERN_C_NAME(radeval)

extern "C"
{

void exmax(double&Ez,double&Bx,double&By,const int &nsources,const double&xs,const double&ys,
           const double&tau,const double&var,const double&amp, const double&a,
           const double&x,const double&y,const double&time);

void wpulse(double&w,const int &nsources,const double&xs,const double&ys,
           const double&tau,const double&var,const double&amp, const double&a,
           const double&x,const double&y,const double&time);

void wdpulse(double&w,double&wt,double&wx,const int &nsources,const double&xs,const double&ys,
           const double&tau,const double&var,const double&amp, const double&a,
           const double&x,const double&y,const double&time);

void radEval(const int &nd, const int &nd1a,const int &nd1b,
             const int &nd2a,const int &nd2b,const int &nd3a,const int &nd3b,
             const int & gridIndexRange,real& u1,real&  u, real & xy, real & rsxy,
             const int &boundaryCondition, 
             const int &md1a,const int &md1b, const int &md2a,const int &md2b, real& huv,
             const int &ld1a,const int &ld1b, const int &ld2a,const int &ld2b, const int &ld3a,const int &ld3b, real& huv3d,

             const int &sd1a,const int &sd1b, const int &sd2a,const int &sd2b, 
             const int &sd3a,const int &sd3b, const int &sd4a,const int &sd4b, const int &sd5a,const int &sd5b,
	     real& uSave,
	     
             const int &bd1a,const int &bd1b, const int &bd2a,const int &bd2b, 
             const int &bd3a,const int &bd3b, const int &bd4a,const int &bd4b, const int &bd5a,const int &bd5b,
	     real& bd,
	     
             const int &ipar,const real& rpar, const int &ierr );

}

// ****** NOTE: these values should be the same as those in forcing.h
static const int nsources=1;
static double xs[nsources]={0.}, ys[nsources]={1.e-8*1./3.}, tau[nsources]={-.95}, var[nsources]={30.}, 
              amp[nsources]={1.};
static double period= 1.;  // period in y

#define getExactSolution EXTERN_C_NAME(getexactsolution)
extern "C"
{
void getExactSolution( real&x,real&y,real&t,int&ndx,int&ndy,int&ndt, real&wd )
// evaluate a derivative of the exact solution
//  ndx,ndy,ndt : number of x,y,t derivatives, e.g. ndx=1,ndy=0,ndt=1 returns w.xt
{
  real w,wt,wx,wy;
  
  if( ndx==0 && ndy==0 && ndt==0 )
  {
    wpulse(wd,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t);
  }
  else if( ndx==0 && ndy==0 && ndt==1 )
  {
    wdpulse(w,wt,wx,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t);
    wd=wt;
  }
  else if( ndx==1 && ndy==0 && ndt==0 )
  {
    wdpulse(w,wt,wx,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t);
    wd=wx;
  }
  else if( ndx==0 && ndy==0 && ndt==3 )
  {
    const real eps=1.e-3;
    real w1,w2,w3,w4;
    
    wpulse(w1,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t-2*eps);
    wpulse(w2,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t-  eps);
    wpulse(w3,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t+  eps);
    wpulse(w4,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t+2*eps);

    wd = (-w1 +2.*w2-2.*w3+w4)/(2.*eps*eps*eps);
    
  }
  else if( ndx==0 && ndy==0 && ndt==4 )
  {
    const real eps=1.e-2;

    real w1,w2,w3,w4,w5;
    
    wpulse(w1,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t-2*eps);
    wpulse(w2,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t-  eps);
    wpulse(w3,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t      );
    wpulse(w4,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t+  eps);
    wpulse(w5,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,t+2*eps);

    wd = (w1 -4.*w2 +6.*w3 -4.*w4 +w5)/(eps*eps*eps*eps);
    
  }
  else if( ndx==1 && ndy==0 && ndt==1 )
  {
    const real eps=1.e-2;

    real w1,w2,w3,w4;
    
    wpulse(w1,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x-  eps,y,t-  eps);
    wpulse(w2,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x+  eps,y,t-  eps);
    wpulse(w3,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x-  eps,y,t+  eps);
    wpulse(w4,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x+  eps,y,t+  eps);

    wd = (w4-w3-w2+w1)/(4.*eps*eps);
    
  }
  else
  {
    printF("getExactSolution:ERROR: ndx,ndy,ndt=%i,%i,%i\n",ndx,ndy,ndt);
  }
}
}



int RadiationBoundaryCondition::debug =0;
real RadiationBoundaryCondition::cpuTime=0.;


RadiationBoundaryCondition::RadiationBoundaryCondition(int orderOfAccuracy_ /* =4 */ )
// =========================================================================================
// /Description:
//     This object can be used to apply a radition boundary condition to one face of a grid.
// 
// =========================================================================================
{
  orderOfAccuracy=orderOfAccuracy_;
  numberOfDerivatives=orderOfAccuracy;
  currentTimeLevel=-1;
  currentTime=0.;
  numberOfTimeLevels=0;
  rside=-1;
  raxis=-1;
  numberOfPoles=-1;
  
  tz=NULL;
  radiationKernel=NULL;
  
  dbase.put<int>("numberOfDimensions")=2;
  dbase.put<bool>("useParallelVersion")=false;

  dbase.put<FILE*>("pDebugFile")=stdout;  // debug file (for parallel)
  
}



RadiationBoundaryCondition::~RadiationBoundaryCondition()
{
  delete radiationKernel;
}


// =========================================================================================
// \brief Set the order of accuracy.
// 
// =========================================================================================
int RadiationBoundaryCondition::
setOrderOfAccuracy(int orderOfAccuracy_ )
{
  orderOfAccuracy=orderOfAccuracy_;
  return 0;
}

// =========================================================================================
// \brief Set the number of poles in the RadiationKernel 
// 
// =========================================================================================
int RadiationBoundaryCondition::
setNumberOfPoles( int numPoles )
{
  numberOfPoles=numPoles;
  return 0;
}

// =========================================================================================
// \brief Set the debug bit flag.
// =========================================================================================
int RadiationBoundaryCondition::setDebug( int debugFlag, FILE *pDebugFile /* =stdout */  )
{
  debug=debugFlag;
  dbase.get<FILE*>("pDebugFile")=pDebugFile;
  return 0;
}


// =========================================================================================
// \brief Uuse new parallel version.
// =========================================================================================
int RadiationBoundaryCondition::useParallelVersion( bool trueOrFalse )
{
  dbase.get<bool>("useParallelVersion")=trueOrFalse;
}



int RadiationBoundaryCondition::
initialize( realMappedGridFunction & u, 
	    int side, int axis,
	    int nc1_/* =0 */, int nc2_/* =0 */, 
	    real c_ /* =1. */, real period_ /* = -1. */, 
	    int numberOfModes_ /* =-1 */, 
	    int orderOfTimeStepping_ /* =-1 */, int numberOfPoles_ /* =-1 */ )
// ===============================================================================================
// /Description:
//     Initialize the boundary conditions.
//
// /mg (input) : apply BC to a face of this grid. (Fastest results if the number of grid points in 
//   the tangential direction is 2**M+1 )
// /side,axis : apply the radiation BC to this face of the grid.
// /nc1,nc2 : apply BC to these components n=nc1,nc1+1,...,nc1
// /numberOfGridPoints : the number of distinct grid points in the periodic direction. 
//                       Fastest results if this is a power of 2.
// /period (input) : solution is periodic on this interval ( set to <0 to have the period computed automatically)
// /c (input) : wave speed.
// /numberOfModes (input) : specify the number of Fourier modes to use in the BC. 
// 
// ===============================================================================================
{
  MappedGrid & mg = *u.getMappedGrid();

  int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
  numberOfDimensions = mg.numberOfDimensions();


  rside=side;
  raxis=axis;

  nc1=nc1_;
  nc2=nc2_;
  assert( nc2-nc1 >=0 );

  const int axisp1=(axis+1) % mg.numberOfDimensions();
  if( !mg.isPeriodic(axisp1) )
  {
    printF("RadiationBoundaryCondition::initialize:ERROR: raxis=%i but mg.isPeriodic(axisp1=%i)=%i\n",
           raxis,axisp1,mg.isPeriodic(axisp1));
    OV_ABORT("error");
  }
  const int axisp2=(axis+2) % mg.numberOfDimensions();
  if( numberOfDimensions==3 && !mg.isPeriodic(axisp2) )
  {
    printF("RadiationBoundaryCondition::initialize:ERROR: raxis=%i but mg.isPeriodic(axisp2=%i)=%i\n",
           raxis,axisp2,mg.isPeriodic(axisp2));
    OV_ABORT("error");
  }  

  if( radiationKernel==NULL )
    radiationKernel = new RadiationKernel; 
  

  // if( numberOfDimensions != 2 )
  // {
  //   OV_ABORT("RadiationBoundaryCondition:: FINISH ME : numberOfDimensions");
  // }
  
  const bool & useParallelVersion = dbase.get<bool>("useParallelVersion");
  if( useParallelVersion )
    printF("**** RadiationBC:: initialize USE PARALLEL VERSION *****\n");

  // This is really the number of points in the Fourier transform:
  numberOfGridPoints1=mg.gridIndexRange(1,axisp1)-mg.gridIndexRange(0,axisp1);
  numberOfModes1=numberOfModes_;
  const int maxModes1=numberOfGridPoints1/2-1;  // this is the max allowed
  if( numberOfModes1<=0 )
  {
    // guess the number of modes to use
    //  Interior scheme converges like  C1*h**M   where M=orderOfAccuracy
    //  Boundary scheme is spectral:    C2*h2**N   where N=1/h2 and h2=spacing for boundary: h2=L/numberOfModes
    //    Set:  C1*h**M = C2*h2**N  -> Nlog(h2) + log(C2) = M*log(h) + log(C1)
    //             Nlog(N) = M*log(1/h) + log(C1/C2)
    if( false )
      numberOfModes1=min(maxModes1,max(8,numberOfGridPoints1/4));  // *test* July 5, 2019
    else
      numberOfModes1=min(maxModes1,max(8,int(sqrt(numberOfGridPoints1))));  


    printF(" **** RadiationBC: numberOfModes1=%i (numberOfGridPoints1=%i) ****\n",numberOfModes1,numberOfGridPoints1);
    
  }
  
  numberOfGridPoints2=mg.gridIndexRange(1,axisp2)-mg.gridIndexRange(0,axisp2);
  numberOfModes2=numberOfModes_;
  const int maxModes2=numberOfGridPoints2/2-1;  // this is the max allowed
  if( numberOfModes2<=0 )
  {
    numberOfModes2=min(maxModes2,max(8,int(sqrt(numberOfGridPoints2))));  
    printF(" **** RadiationBC: numberOfModes2=%i (numberOfGridPoints2=%i) ****\n",numberOfModes2,numberOfGridPoints2);
    
  }

// // **** TEST
//   numberOfModes1=maxModes1;
//   numberOfModes2=maxModes2;


  period1=period_;
  period2=period_; 
  c=c_;
  radius=-1.;
  
  const bool isRectangular = mg.isRectangular();

  if( radiationKernel->getKernelType() == RadiationKernel::planar )
  {
    if( !isRectangular )
    {
      printF("RaditionBoundaryCondition::ERROR: rectangular grid expected for RadiationKernel::planar\n");
      OV_ABORT("ERROR");
    }
      
    if( orderOfAccuracy==2 )
    {
      numberOfDerivatives = 1 + 1;  // we save u,ux
    }
    else
    {
      numberOfDerivatives = 1 + 1 + 1 + 1;  // we save u,ux,uxx,uxxx
    }

    if( period1<0 )
    {
      // determine the period (periodic interval)
      real dx[3];
      mg.getDeltaX(dx);

      period1=dx[axisp1]*numberOfGridPoints1; 
      printF("RadiationBoundaryCondition:INFO the periodic interval was computed to be period1=%9.3e\n",period1);

      period2=dx[axisp2]*numberOfGridPoints2; 
      if( mg.numberOfDimensions()==3 )
	printF("RadiationBoundaryCondition:INFO the periodic interval was computed to be period2=%9.3e\n",period2);

    }
  }
  else if( radiationKernel->getKernelType() == RadiationKernel::cylindrical )
  {
    // radiationKernel->setKernelType(RadiationKernel::cylindrical);
    if( isRectangular )
    {
      printF("RaditionBoundaryCondition::ERROR: curvilinear grid expected for RadiationKernel::cylindrical.\n");
      OV_ABORT("ERROR");
    }
    

    if( orderOfAccuracy==2 )
    {
      numberOfDerivatives = 1 + 2;  // we save u,ux,uy  
    }
    else
    {
      numberOfDerivatives = 1 + 2 + 3 + 4;  // we save u, ux,uy, uxx,uxy,uyy, uxxx,uxxy,uxyy,uyyy 
    }

    // determine the radius 
    Mapping & map = mg.mapping().getMapping();
    realArray r(1,3), x(1,3);
    r=0.;
    r(0,axis)=(real)side;
      
    map.map(r,x);
      
    radius = sqrt( x(0,0)*x(0,0) + x(0,1)*x(0,1) );
      
    printF("RaditionBoundaryCondition:INFO kernelType is cylindrical and the radius was computed to be %10.4e\n",
	   radius);


  }
  else
  {
    printF("RaditionBoundaryCondition::ERROR: un-expected kernelType=%d\n",(int)radiationKernel->getKernelType());
    OV_ABORT("ERROR");
  }
  
  
  orderOfTimeStepping=orderOfTimeStepping_;
  if( orderOfTimeStepping==-1 )
    orderOfTimeStepping=orderOfAccuracy+1;
  
  numberOfPoles=numberOfPoles_;
  if( numberOfPoles!=21 && numberOfPoles!=31 )
  {
    if( numberOfPoles>0 )
    {
      printF("RadiationBoundaryCondition::initialize:ERROR: numberOfPoles=%d is not available. Choosing numberOfPoles=21.\n",numberOfPoles);
    }
    
    numberOfPoles=21;
  }

  if( !(numberOfPoles==21 || numberOfPoles==31 ) )
  {
    printF("RadiationBoundaryCondition::initialize: ERROR: numberOfPoles=%d not supported\n",numberOfPoles);
    OV_ABORT("error");
  }
  
  printF("\n ++++  RadiationBoundaryCondition::initialize: SETTING numberOfPoles=%d +++++\n",numberOfPoles);
  
  // **TEST**
  // numberOfPoles=31;
  


  const int numberOfComponents=nc2-nc1+1;
  const int numberOfFields=numberOfDerivatives*numberOfComponents;
  numberOfTimeLevels=orderOfAccuracy+1;
  
  if( useParallelVersion )
  {
    // *new way* using parallel version -- Oct/Nov 2020
    radiationKernel->setUseFourierTransformClass(true);
    radiationKernel->initialize( u,side,axis,numberOfFields,numberOfModes1, numberOfModes2,c,orderOfTimeStepping, numberOfPoles );

    // -- We store a number of time levels of u and some of it's normal derivatives on the boundary --
    if( !dbase.has_key("boundaryData") )
      dbase.put<RealArray>("boundaryData");

    // printF("\n **** RADBC: INIT boundaryData: numberOfDerivatives=%d numberOfComponents=%d numberOfFields=%d\n",
    // 	   numberOfDerivatives,numberOfComponents,numberOfFields);
    

    RealArray & boundaryData = dbase.get<RealArray>("boundaryData");
    Range Nf=numberOfFields, Nt=numberOfTimeLevels;
    Index Rv[3], &R1=Rv[0], &R2=Rv[1], &R3=Rv[2];
    R1=u.dimension(0); R2=u.dimension(1); R3=u.dimension(2); 
    Rv[axis]=Range(mg.gridIndexRange(side,axis),mg.gridIndexRange(side,axis));

    OV_GET_SERIAL_ARRAY(real,u,uLocal);
    int includeGhost=1;  // what should this be? 
    bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,R1,R2,R3,includeGhost);

    if( ok )
    {
      boundaryData.redim(R1,R2,R3,Nf,Nt);
      boundaryData=0.;
    }
    
      
  }
  else
  {
    // *old serial code*
    radiationKernel->initialize( numberOfDimensions,
				 numberOfGridPoints1, numberOfGridPoints2,
				 numberOfFields, 
				 numberOfModes1, numberOfModes2,
				 period1, period2, c, 
				 orderOfTimeStepping, numberOfPoles, radius );


    if( numberOfDimensions==2 )
      uSave.redim(numberOfGridPoints1,numberOfTimeLevels,numberOfDerivatives,numberOfComponents);
    else
      uSave.redim(numberOfGridPoints1,numberOfGridPoints2,numberOfTimeLevels,numberOfDerivatives,numberOfComponents);

    uSave=0.;
  }
  
  currentTimeLevel=0;
  currentTime=0.;
  
}



int RadiationBoundaryCondition::
assignBoundaryConditions( realMappedGridFunction & u, real t, real dt,
			  realMappedGridFunction & u2 )
// ==============================================================================================================
// /Description:
//     Assign the radiation boundary conditions.
//
//  u (input/output) : u at time t at interior and boundary points. This routine will assign ghost points of u
//  u2 (input): u at the previous time, u(t-dt)
// ==============================================================================================================
{
  // debug=1;  // *** TEMP **
  

  
  if( t<0 ) return 0;  // Do this for now -- assume the solution is zero near the radiation boundary
  if( t==0 && t==currentTime ) return 0;

  const int myid = max(0,Communication_Manager::My_Process_Number);
  const int np   = max(1,Communication_Manager::Number_Of_Processors);

  real time0=getCPU();
  if( t<= currentTime )
  {
    printF("RBC:assignBC: Error: t=%9.3e <= currentTime=%9.3e. dt=%8.3e\n",t,currentTime,dt);
    return 0;
  }
  assert( t>currentTime );  // 
  assert( numberOfTimeLevels>0 );

  currentTime=t;
  currentTimeLevel = (currentTimeLevel +1 ) % numberOfTimeLevels;

  const bool & useParallelVersion = dbase.get<bool>("useParallelVersion");
  if( debug>0 && t<=5*dt )
  {
    printF("RadiationBoundaryCondition::assignBoundaryConditions t=%9.3e dt=%8.2e rside=%i raxis=%i useParallelVersion=%d\n",
	   t,dt,rside,raxis,useParallelVersion);
  }
  
 FILE *pDebugFile= dbase.get<FILE*>("pDebugFile");


  const RadiationKernel::KernelTypeEnum kernelType = radiationKernel->getKernelType();
  
  MappedGrid & mg = *( u.getMappedGrid() );
  const int numberOfDimensions = mg.numberOfDimensions();
  
  const bool isRectangular = mg.isRectangular();
  if( kernelType==RadiationKernel::planar ) 
    assert( isRectangular );

  real dx[3];
  mg.getDeltaX(dx);

  // Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
  int i1,i2,i3;


  OV_GET_SERIAL_ARRAY(real,u,uLocal);
  OV_GET_SERIAL_ARRAY(real,u2,u2Local);


  const int side=rside, axis=raxis;
  const int numberOfComponents=nc2-nc1+1;
  const int numberOfFields=numberOfDerivatives*numberOfComponents;
  
  
  RealArray & ub = uSave; 
  
  Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
  getBoundaryIndex(mg.indexRange(),side,axis,I1,I2,I3);
  int includeGhost=0;  // No ghost *wdh* Nov 27, 2020
  bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3,includeGhost);

  int n1a,n1b,n2a,n2b,n3a,n3b;
  if( true )
  {
    n1a=I1.getBase(); n1b=I1.getBound();
    n2a=I2.getBase(); n2b=I2.getBound();
    n3a=I3.getBase(); n3b=I3.getBound();
  }
  else
  {
     
    if( axis==0 )
    {
      n1a=mg.gridIndexRange(side,0); n1b=n1a;
      n2a=mg.gridIndexRange(0,1);    n2b=mg.gridIndexRange(1,1)-1;
      n3a=mg.gridIndexRange(0,2);    n3b= numberOfDimensions==2 ? mg.gridIndexRange(1,2) : mg.gridIndexRange(1,2)-1;
    }
    else
    {
      n1a=mg.gridIndexRange(0,0);    n1b=mg.gridIndexRange(1,0)-1;
      n2a=mg.gridIndexRange(side,1); n2b=n2a;
      n3a=mg.gridIndexRange(0,2);    n3b= numberOfDimensions==2 ? mg.gridIndexRange(1,2) : mg.gridIndexRange(1,2)-1;
    }
  }
  if( false )
  {
    printF("RadBC: [n1a,n1b]=[%d,%d] [n2a,n2b]=[%d,%d] [n3a,n3b]=[%d,%d]  I1=[%d,%d] I2=[%d,%d] I3=[%d,%d]\n",
           n1a,n1b,n2a,n2b,n3a,n3b, I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound());
  }
  
  const int axisp1=(axis+1) % mg.numberOfDimensions();
  const int axisp2=(axis+2) % mg.numberOfDimensions();

  Range R1(mg.dimension(0,axisp1),mg.dimension(1,axisp1));
  Range R2(mg.dimension(0,axisp2),mg.dimension(1,axisp2));

	  
  const int m=currentTimeLevel;
  const int mm1= (m-1+numberOfTimeLevels) % numberOfTimeLevels;
  const int mm2= (m-2+numberOfTimeLevels) % numberOfTimeLevels;
  const int mm3= (m-3+numberOfTimeLevels) % numberOfTimeLevels;
  // const int mm4= (m-4+numberOfTimeLevels) % numberOfTimeLevels;

  RealArray vv, huv;
  if(  useParallelVersion )
  {
    // *new* parallel version -- Oct 28, 2020
    // Works in serial or parallel but generally slower

    Range N=numberOfFields;
    RealArray & boundaryData = dbase.get<RealArray>("boundaryData");
    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    I1=u.dimension(0); I2=u.dimension(1); I3=u.dimension(2); 
    Iv[axis]=Range(mg.gridIndexRange(side,axis),mg.gridIndexRange(side,axis));

    int includeGhost=1;  // what should this be? 
    bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3,includeGhost);

    // Set this to true to compare the new and old versions of the RadiationKernel
    bool useOldRadEval= false; // true && np==1;
    RealArray vv2,huv2; // ***** clean this up ****

    if( ok )
    {
      // Index C(0,numberOfFields,numberOfDerivatives);
      // --- Save u boundary data at current time level m (derivative data is saved below in radEval)
      for( int n=nc1; n<=nc2; n++ )
      {
	int mc = 0 + numberOfDerivatives*(n);
	boundaryData(I1,I2,I3,mc,m) = uLocal(I1,I2,I3,n); 
      }
    
      vv.redim(I1,I2,I3,N);
      huv.redim(I1,I2,I3,N);

      vv(I1,I2,I3,N) = boundaryData(I1,I2,I3,N,mm1);  //   all (u,ux,...) boundary data at previous time level mm1

      if( useOldRadEval )
      {
	vv2.redim(vv); huv2.redim(huv);
	vv2=vv;  huv2=huv;
	// reshape to support old way : 
	if( numberOfDimensions==2 )
	{
	  if( axis==0 )
	  {
	    vv2.reshape(I2,N); huv2.reshape(I2,N);
	  }
	  else
	  {
	    vv2.reshape(I1,N); huv2.reshape(I1,N);
	  }
	}
	else
	{
	  if( axis==0 )
	  {
	    vv2.reshape(I2,I3,N); huv2.reshape(I2,I3,N);
	  }
	  else if( axis==1 )
	  {
	    vv2.reshape(I1,I3,N); huv2.reshape(I1,I3,N);
	  }
	  else
	  {
	    vv2.reshape(I1,I2,N); huv2.reshape(I1,I2,N);
	  }
	}
      }
    }
    // end if ok 

    if( tz==NULL )
    {
      // ----------------------------------------
      // -------- EVALUATE RADIATION KERNEL -----
      // ----------------------------------------
      // this next call will assign periodic images on hu also.
      if( useOldRadEval )
	radiationKernel->evaluateKernel( dt, vv2, huv2 ); 

      radiationKernel->evaluateKernelParallel( dt, vv, huv ); 

    }
    else
    {
        huv=0; // do this for now for TZ 
    }

    if( ok )
    {
      // reshape to support old way :   (**fix me ?) 
      if( numberOfDimensions==2 )
      {
	if( axis==0 )
	  huv.reshape(I2,N);
	else
	  huv.reshape(I1,N);
      }
      else
      {
	if( axis==0 )
	  huv.reshape(I2,I3,N);
	else if( axis==1 )
	  huv.reshape(I1,I3,N);
	else
	  huv.reshape(I1,I2,N);
      }
      if( useOldRadEval )
      {
	if( false )
	{
	  RealArray err(huv);
	  err=huv-huv2;
	  ::display(huv2,"huv2","%9.3e ");
	  ::display(huv,"huv","%9.3e ");
	  ::display(err,"err","%9.3e ");
	}
      
	real maxDiff = max(fabs(huv-huv2));
	printF("RBC: t=%9.3e: |huv-huv2| = %9.3e\n",t,maxDiff);
      }
    }
    

  }
  else
  {

    // ----- OLD SERIAL VERSION (faster) ----

    
    // --------- SAVE THE CURRENT SOLUTION ON THE BOUNDARY ARRAY  -----
    i3=n3a;
    if( axis==0 )
    {
      i1=side==0 ? n1a : n1b;
      for( int n=nc1; n<=nc2; n++ )
      {
	if( numberOfDimensions==2 )
	{
	  for( int i2=n2a; i2<=n2b; i2++ )
	    ub(i2,m,0,n)=uLocal(i1,i2,i3,n);              // save solution at time t, m=currentLevel 
	}
	else
	{
	  for( int i3=n3a; i3<=n3b; i3++ )
	  {
	    for( int i2=n2a; i2<=n2b; i2++ )
	    {
	      ub(i2,i3,m,0,n)=uLocal(i1,i2,i3,n);         // save solution at time t
	    }
	  }
	}
      }// end for n 

      if( false )
      {
	printF("RBC: t=%9.3e: [n2a,n2b][n3a,n3b]=[%d,%d][%d,%d]\n",t,n2a,n2b,n3a,n3b);
	int ey=1;
	display(uLocal(i1,R1,R2,ey),sPrintF("Ey on side=%d at t=%9.3e",side,t),"%9.2e ");
      }
    
    }
    else if( axis==1 )
    {
      i2=side==0 ? n2a : n2b;
      for( int n=nc1; n<=nc2; n++ )
      {
	if( numberOfDimensions==2 )
	{
	  for( int i1=n1a; i1<=n1b; i1++ )
	    ub(i1,m,0,n)=uLocal(i1,i2,i3,n);         // save solution at time t
	}
	else
	{
	  for( int i3=n3a; i3<=n3b; i3++ )
	    for( int i1=n1a; i1<=n1b; i1++ )
	      ub(i1,i3,m,0,n)=uLocal(i1,i2,i3,n);         // save solution at time t
	}
      }
    
    }
    else
    {
      OV_ABORT("FIX ME");
    }
  
  

    // ***** Evaluate the kernel ***

    Range N=numberOfFields; 
    // RealArray vv, huv;  // do all kernel evaluations at once
    if( numberOfDimensions==2 )
    {
      vv.redim(R1,N); huv.redim(R1,N);
      // vv=-99.;  // for testing 
      vv=-9999999.;  huv=0.; // for testing 
    }
    else
    {
      vv.redim(R1,R2,N); huv.redim(R1,R2,N);
      // do this for now: 
      vv=0.;  huv=0.;  // for testing 
    }
  

    #define vvn(i,m,n) vv(i,m+numberOfDerivatives*(n-nc1))
    #define vvn3(i,j,m,n) vv(i,j,m+numberOfDerivatives*(n-nc1))

    if( false )
      display(ub,"ub","%9.2e ");


    // --- set vv = u(t-dt) on the boundary ---
    

    // Assign periodic images too seems to be necessary  *wdh* **Nov 7, 2020**
    Index I1(mg.gridIndexRange(0,axisp1),mg.gridIndexRange(1,axisp1)); 
    Index I2(mg.gridIndexRange(0,axisp2),mg.gridIndexRange(1,axisp2)); 

    //Index I1(mg.gridIndexRange(0,axisp1),mg.gridIndexRange(1,axisp1)-1); 
    //Index I2(mg.gridIndexRange(0,axisp2),mg.gridIndexRange(1,axisp2)-1); 

    for( int n=nc1; n<=nc2; n++ )  // --- loop over components ---
    {
      for( int i=0; i<numberOfDerivatives; i++ )  // loop over derivatives 
      {
	if( numberOfDimensions==2 )
	{
	  vvn(I1,i,n)=ub(I1,mm1,i,n); 
	  // vvn(R1,i,n)=ub(R1,mm1,i,n);   // assign over R1 just so there are some values in ghost when debugging 
	}
	else
	{
	  vvn3(I1,I2,i,n)=ub(I1,I2,mm1,i,n); 
	}
	
      }
    }
    
    if( false )
    {
      printF("\n +++++ RBC: after assigning vv, currentTime=%d, m=%d, numberOfTimeLevels=%d\n",currentTimeLevel,m,numberOfTimeLevels);
      display(vv, sPrintF("vv  BEFORE evaluateKernel side=%d, t=%9.2e",side,t),"%9.2e ");
    }
    

    // ----------------------------------------
    // -------- EVALUATE RADIATION KERNEL -----
    // ----------------------------------------
    if( tz==NULL )
    {

      // this next call will assign periodic images on hu also.
      radiationKernel->evaluateKernel( dt, vv, huv ); 

    }
    else
    {
      huv=0; // do this for now

      // we should set
      // huv(i,0) = -[ u.t + u.r + u/(2r) ]/c
    }

  
    #undef vvn 

  }  // eval old serial way 

  
  if( debug & 8 || debug & 8  )
  {
    display(vv, sPrintF("vv  after evaluateKernel side=%d, t=%9.2e",side,t),pDebugFile,"%9.2e ");
    display(huv,sPrintF("huv after evaluateKernel side=%d, t=%9.2e",side,t),pDebugFile,"%9.2e ");
  }
  
  // -------------------------------------------------------------------------
  // ------------ FILL GHOST POINTS ------------------------------------------
  // -------------------------------------------------------------------------

  real *uptr = uLocal.getDataPointer();
  real *u2ptr = u2Local.getDataPointer();
  real *puSave = uSave.getDataPointer();

  real *pxy = uptr;  // set default if not used.
  real *prsxy=uptr;  

  if( !isRectangular )
  {
    mg.update(MappedGrid::THEvertex | MappedGrid::THEinverseVertexDerivative); 

    OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
    OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal);
    
    // pxy=mg.vertex().getDataPointer(); 
    // prsxy=mg.inverseVertexDerivative().getDataPointer(); 
    pxy=xLocal.getDataPointer(); 
    prsxy=rxLocal.getDataPointer(); 

  }
  else if( debug>0 )
  {
    mg.update(MappedGrid::THEvertex); 
    OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);

    // pxy=mg.vertex().getDataPointer(); 
    pxy=xLocal.getDataPointer(); 
  }
  
  const int base=huv.getBase(0);
  real *phuv= huv.getDataPointer(); 
    

  const IntegerArray & dim = mg.dimension();
  int ierr=0;
  int ipar[30];
  real rpar[30];
			  
  int grid=0;
  int gridType=0;  // rectangular

  real eps=1., mu=1.; // not used

  real alpha=orderOfAccuracy==2 ? 1. : 1.;  // dissipation term -- 10 is probably too large

  assert( n1a==I1.getBase() && n1b==I1.getBound() );
  assert( n2a==I2.getBase() && n2b==I2.getBound() );
  assert( n3a==I3.getBase() && n3b==I3.getBound() );

  ipar[0] = side;
  ipar[1] = axis;
  ipar[2] = grid;
  ipar[3] = n1a;
  ipar[4] = n1b;
  ipar[5] = n2a;
  ipar[6] = n2b;
  ipar[7] = n3a;
  ipar[8] = n3b;
  ipar[9] = nc1;
  ipar[10]= nc2; 
  ipar[11]= currentTimeLevel;
  ipar[12]= numberOfTimeLevels;

  ipar[13]= gridType;
  ipar[14]= orderOfAccuracy;
  ipar[15]= debug;
  ipar[16]= (int)kernelType;
  ipar[17]= useParallelVersion;
  
     	      			   
  rpar[0] = dx[0];
  rpar[1] = dx[1];
  rpar[2] = dx[2];
  rpar[3] = mg.gridSpacing(0);
  rpar[4] = mg.gridSpacing(1);
  rpar[5] = mg.gridSpacing(2);
  rpar[6] = alpha;
  // assert( tz!=NULL );
  rpar[7]=(real &)tz;  // twilight zone pointer
     	      			   
  rpar[10]= t;
  rpar[11]= dt;
  rpar[12]= eps;
  rpar[13]= mu;
  rpar[14]= c;

  
  // boundaryData: 
  RealArray & bd = useParallelVersion ? dbase.get<RealArray>("boundaryData") : uSave;
  if( ok )
  {
    radEval( mg.numberOfDimensions(), 
	     // dim(0,0),dim(1,0),dim(0,1),dim(1,1),dim(0,2),dim(1,2),
	     uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),uLocal.getBase(2),uLocal.getBound(2),
	     mg.gridIndexRange(0,0), *uptr, *u2ptr, *pxy, *prsxy,
	     mg.boundaryCondition(0,0),
	     huv.getBase(0),huv.getBound(0),                               0,numberOfDerivatives-1,*phuv,
	     huv.getBase(0),huv.getBound(0),huv.getBase(1),huv.getBound(1),0,numberOfDerivatives-1,*phuv,
	     uSave.getBase(0),uSave.getBound(0),
	     uSave.getBase(1),uSave.getBound(1),
	     uSave.getBase(2),uSave.getBound(2),
	     uSave.getBase(3),uSave.getBound(3),
	     uSave.getBase(4),uSave.getBound(4),
	     *puSave,
	     bd.getBase(0),bd.getBound(0),
	     bd.getBase(1),bd.getBound(1),
	     bd.getBase(2),bd.getBound(2),
	     bd.getBase(3),bd.getBound(3),
	     bd.getBase(4),bd.getBound(4),
	     *bd.getDataPointer(),
	     ipar[0], rpar[0], ierr );
  }
  

  if( false )
  {
    // Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];
    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);

    // ::display(uLocal(Ig1,Ig2,Ig3,0),"u-ghost After radEval","%6.3f ");
    ::display(uLocal,"u After radEval","%6.3f ");
  }
  


  cpuTime+=getCPU()-time0;
  return 0;
  
}

