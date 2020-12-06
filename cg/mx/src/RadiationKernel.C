#include "RadiationKernel.h"
#include "display.h"
#include "ParallelUtility.h"

#include "FourierTransform.h"


#define bcperq21d EXTERN_C_NAME(bcperq21d)
#define bcperq31d EXTERN_C_NAME(bcperq31d)
#define bcper3dq21 EXTERN_C_NAME(bcper3dq21)
#define bccyld EXTERN_C_NAME(bccyld)
#define AMCOF EXTERN_C_NAME(amcof)
#define RFFTI EXTERN_C_NAME(rffti)
#define RFFTF EXTERN_C_NAME(rfftf)
#define RFFTB EXTERN_C_NAME(rfftb)
#define ZFFTI EXTERN_C_NAME(zffti)
#define ZFFTF EXTERN_C_NAME(zfftf)
#define ZFFTB EXTERN_C_NAME(zfftb)

extern "C"
{
  void bcperq21d( const int&nda,const int&ndb,double & p,double &f,double &ploc,double &c,double &len,double &dt,
                int&n,int&m,int&md,int&ns,int&ord,double &fold,double &phi,double &amc,
		double &fftsave,int&bcinit);

  void bcperq31d( const int&nda,const int&ndb,double & p,double &f,double &ploc,double &c,double &len,double &dt,
                int&n,int&m,int&md,int&ns,int&ord,double &fold,double &phi,double &amc,
		double &fftsave,int&bcinit);

  void bccyld(const int&nda,const int&ndb,
             double &p,double &f,double &ploc,double &c,double &r,double &dt,int&n,int&m,int&md,int&ns,int&ord,
             double &fold, double &phi,
             int&npoles,double &alpha,double &beta,double &amc,double &fftsave,int&bcinit);

  void bcper3dq21( const int&nda1,const int&ndb1,const int&nda2,const int&ndb2,
		   double &p,double &f,double &pl,double &zl,double &pl1,double &zl2,double &c,
                   const double&len1,const double&len2,double &dt,
		   const int&n1,const int&n2,const int&m,const int&md1,const int&md2,const int&ns1,const int&ns2,
                   const int&ord,double &fold,double &phi,double &amc,double &fftsave1,double &fftsave2,
                   int& bcinit);

  void AMCOF( real & amc, const int & orderInTime);
  void RFFTI( const int & n, real & fftsave );
  void RFFTF( const int & n1, real & f, real & fftsave );
  void RFFTB( const int & n1, real & f, real & fftsave );

  void ZFFTI( const int & n, real & fftsave );
  void ZFFTF( const int & n1, real & f, real & fftsave );
  void ZFFTB( const int & n1, real & f, real & fftsave );

}

real RadiationKernel::cpuTime=0.;

RadiationKernel::RadiationKernel()
// ======================================================================================
// This class is used to evaluate the Kernel of a Radiation boundary condition.
// ======================================================================================
{
  kernelType=planar;
  
  ploc=NULL;
  fold=NULL;
  phi=NULL;
  amc=NULL;
  fftsave1=NULL;
  fftsave2=NULL;

  zl=NULL;
  pl1=NULL;
  zl2=NULL;  

  radius=1.;
  alpha=NULL;
  beta=NULL;
  npoles=NULL;
  
  bcinit=-1;
  
  dbase.put<int>("numberOfDimensions")=2;
  dbase.put<int>("bcInitComplex")=0;  

  dbase.put<int>("currentTimeLevel")=-1;    // current time level for time integral by quadrature 

  dbase.put<int>("rside")=0;  
  dbase.put<int>("raxis")=0;

  dbase.put<bool>("useFourierTransformClass")=false; // if true, use new FourierTransform class
  
  dbase.put<int>("debug")=0;
  dbase.put<aString>("debugFileName")="radKernel";

}




// =====================================================================================================
/// \brief Set the debug flag, non-zero to output info
/// \note Set to 1 for some output and 1+2=3 for more
// =====================================================================================================
int RadiationKernel::setDebug( int value )
{
  int & debug = dbase.get<int>("debug");
  debug=value;
  return 0;
}

// =====================================================================================================
/// \brief Set the debug file name 
// =====================================================================================================
int RadiationKernel::setDebugFileName( const aString & fileName )
{
  aString & debugFileName = dbase.get<aString>("debugFileName");
  debugFileName = fileName;
  
  return 0;
}

// ============================================================================
// /Description:
//   Set the kernel type. Call this before initialize.
// ============================================================================
int RadiationKernel::setKernelType( KernelTypeEnum type )
{
  kernelType=type;
  return 0;
}

RadiationKernel::KernelTypeEnum RadiationKernel::
getKernelType() const
{
  return kernelType;
}


// ===========================================================================================
/// \brief: Specify whether to use the new FourierTransform class (Supports parallel FFTs)
// ===========================================================================================
int RadiationKernel::setUseFourierTransformClass( bool trueOrFalse /* = true */ )
{
  bool & useFourierTransformClass = dbase.get<bool>("useFourierTransformClass");
  useFourierTransformClass=trueOrFalse;

  return 0;
}



// ===========================================================================================
/// \brief:  New initialize routine for parallel
/// \param mg, rside, raxis  (input) : radiation BC is on face (rside,raxis) of this mapped grid.
/// \param numberOfFields (input) : number of fields.
/// \param orderOfTimeStepping (input): order of accuracy of time-stepping of the Kernel equations.
/// \param numberOfPoles (input) : number of poles in Kernel expansion.
// ===========================================================================================
int RadiationKernel::initialize( realMappedGridFunction & u, 
				 int rside, int raxis,
				 int numberOfFields_, 
				 int numberOfModes1_,   int numberOfModes2_,
				 real c_, 
				 int orderOfTimeStepping_, int numberOfPoles_,
                                 real radius_ /* =1. */ )				 
{
  const int myid = max(0,Communication_Manager::My_Process_Number);
  const int np=max(1,Communication_Manager::Number_Of_Processors);


  MappedGrid & mg = *u.getMappedGrid();

  int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
  numberOfDimensions = mg.numberOfDimensions();


  const int & debug = dbase.get<int>("debug");
  const aString & debugFileName = dbase.get<aString>("debugFileName");
  if( !dbase.has_key("debugFile") )
  {
    dbase.put<FILE*>("debugFile")=NULL;

    FILE *& debugFile=dbase.get<FILE*>("debugFile");
    if( debug !=0 && debugFile==NULL )
    {
      char fileName[40];
      sprintf(fileName,"%s%i.debug",(const char*)debugFileName,myid);
      debugFile= fopen(fileName,"w");
    }
  }
  FILE *debugFile=dbase.get<FILE*>("debugFile");



  if( debug & 1 )
    printF("RadiationKernel::initialize... numberOfDimensions=%d, np=%d\n",numberOfDimensions,np);
  
  OV_GET_SERIAL_ARRAY(real,u,uLocal);
  if( debug & 1 )
  {
    fprintf(debugFile,"RadKer: myid=%d: uLocal=[%d,%d][%d,%d][%d,%d]\n",myid,
	   uLocal.getBase(0),uLocal.getBound(0),
	   uLocal.getBase(1),uLocal.getBound(1),
	   uLocal.getBase(2),uLocal.getBound(2));
  }
  

  IndexBox uLocalArrayBox;
  CopyArray::getLocalArrayBox( myid, u, uLocalArrayBox );

  if( debug & 1 )
    fprintf(debugFile,"RadKer:init: myid=%i uLocal bounds=[%i,%i][%i,%i][%i,%i][%i,%i] (no ghost)\n",
	    myid,
	    uLocalArrayBox.base(0),uLocalArrayBox.bound(0),
	    uLocalArrayBox.base(1),uLocalArrayBox.bound(1),
	    uLocalArrayBox.base(2),uLocalArrayBox.bound(2),
	    uLocalArrayBox.base(3),uLocalArrayBox.bound(3));
    

  Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
  getIndex(mg.indexRange(),I1,I2,I3);
  int includeGhost=0;
  bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3,includeGhost);
  if( debug & 1 )
    fprintf(debugFile,"RadiationKernel:init: myid=%i indexRange localBounds=[%i,%i][%i,%i][%i,%i]\n",
	    myid,I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound());

  dbase.get<int>("rside")=rside;  
  dbase.get<int>("raxis")=raxis;  

  // const int axisp1 = (raxis+1) % numberOfDimensions;
  // const int axisp2 = (raxis+2) % numberOfDimensions;
  int axisp1, axisp2;
  
  const IntegerArray & gid = mg.gridIndexRange();
  if( raxis==0 )
  {
    axisp1=1; axisp2=2;  // tangential directions
  }
  else if( raxis==1 )
  {
    axisp1=0; axisp2=2;
  }
  else 
  {
    axisp1=0; axisp2=1;
  }

  numberOfGridPoints1 = gid(1,axisp1)-gid(0,axisp1);
  numberOfGridPoints2 = gid(1,axisp2)-gid(0,axisp2);

  if( debug & 1 )
    fprintf(debugFile,"RadiationKernel: axisp1=%d: numberOfGridPoints1=%d, axisp2=%d: numberOfGridPoints2=%d"
	   " (power of 2 is best)\n",
	   axisp1,numberOfGridPoints1,axisp2,numberOfGridPoints2);

  // int numberOfModes1=numberOfModes1_;
  // int numberOfModes2=numberOfModes2_;
  

  // --- save the local box of grid points used in the transform ----
  IndexBox & fbox = dbase.put<IndexBox>("fbox");
  fbox.setBounds( uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
		  uLocal.getBase(2),uLocal.getBound(2),0,0 );
  if( raxis>=0 )
  {
    // restrict box to a face
    fbox.setBaseBound( raxis,gid(rside,raxis),gid(rside,raxis) );
  }
  if( debug & 1 )
  {
    fprintf(debugFile,"RadKer: fBox =[%3i,%3i][%3i,%3i][%3i,%3i][%3i,%3i]\n",
	   fbox.base(0),fbox.bound(0),
	   fbox.base(1),fbox.bound(1),
	   fbox.base(2),fbox.bound(2),
	   fbox.base(3),fbox.bound(3));
  }


  const bool isRectangular = mg.isRectangular();

  if( getKernelType() == planar )
  {
    if( !isRectangular )
    {
      printF("RadiationKernel::ERROR: rectangular grid expected for RadiationKernel::planar\n");
      OV_ABORT("ERROR");
    }

    // determine the period (periodic interval)
    real dx[3];
    mg.getDeltaX(dx);

     // -- check periodicity ---
    if( mg.isPeriodic(axisp1)!=Mapping::derivativePeriodic )
    {
      printF("RaditionKernel::ERROR: grid is not `derivativePeriodic' along axisp1=%d\n",axisp1);
      OV_ABORT("ERROR");
    }
    if( numberOfDimensions==3 && mg.isPeriodic(axisp2)!=Mapping::derivativePeriodic )
    {
      printF("RaditionKernel::ERROR: grid is not `derivativePeriodic' along axisp2=%d\n",axisp2);
      OV_ABORT("ERROR");
    }
    
    period1=dx[axisp1]*numberOfGridPoints1; 
    period2=dx[axisp2]*numberOfGridPoints2; 
      
    if( debug & 1 )
    {
      printF("RadiationKernel:INFO the periodic interval was computed to be period1=%9.3e\n",period1);
      if( mg.numberOfDimensions()==3 )
	printF("RadiationKernel:INFO the periodic interval was computed to be period2=%9.3e\n",period2);
    }
    
  }
  else if( getKernelType() == cylindrical )
  {
    if( isRectangular )
    {
      printF("RaditionKernel::ERROR: curvilinear grid expected for RadiationKernel::cylindrical.\n");
      OV_ABORT("ERROR");
    }
    
    // ---- determine the radius ----
    // **FIX ME FOR PARALLEL**
    Mapping & map = mg.mapping().getMapping();
    realArray r(1,3), x(1,3);
    r=0.;
    r(0,raxis)=(real)rside;
      
    map.map(r,x);
      
    radius = sqrt( x(0,0)*x(0,0) + x(0,1)*x(0,1) );
      
    printF("RaditionKernel:INFO kernelType is cylindrical and the radius was computed to be %10.4e\n",
	   radius);


  }
  else
  {
    printF("RaditionKernel::ERROR: un-expected kernelType=%d\n",(int)getKernelType());
    OV_ABORT("ERROR");
  }

  

  // if( np>1 )
  // {
  //   OV_ABORT("RadiationKernel::initialize: finish me");
  // }


  const bool & useFourierTransformClass = dbase.get<bool>("useFourierTransformClass");
  //  printF("******  useFourierTransformClass=%d ******\n",(int) useFourierTransformClass);
  
  if( useFourierTransformClass )
  {
    // ----- class to compute FFT's in parallel ------

    dbase.put<FourierTransform>("fourierTransform");
    FourierTransform & fourierTransform = dbase.get<FourierTransform>("fourierTransform");

    // Turn off physical transforms (reverse of default forward/backward FFTs)
    fourierTransform.usePhysicalTransforms(false);
    int fftDebug=debug; // 
    fourierTransform.setDebug(fftDebug);
   
    const int ndfft = mg.numberOfDimensions()-1;  // 1D FFTs in 2D,  2D FFTs in 3D 
    fourierTransform.initialize( u, ndfft, rside,raxis );


    // --- Here is the size of the LOCAL FFT ---
    dbase.put<IndexBox>("fftBox");
    IndexBox & fftBox = dbase.get<IndexBox>("fftBox");
    fourierTransform.getLocalIndexBox( fftBox );


    if( debug & 1 )
    {
      fprintf(debugFile,"RadKer:INIT: fftBox =[%3i,%3i][%3i,%3i][%3i,%3i][%3i,%3i]\n",
	      fftBox.base(0),fftBox.bound(0),
	      fftBox.base(1),fftBox.bound(1),
	      fftBox.base(2),fftBox.bound(2),
	      fftBox.base(3),fftBox.bound(3));
    }

  }

  // do this for now 
  return initialize(mg.numberOfDimensions(),
		    numberOfGridPoints1, numberOfGridPoints2,
		    numberOfFields_, 
		    numberOfModes1_,   numberOfModes2_,
		    period1,  period2,
		    c_, 
		    orderOfTimeStepping_, numberOfPoles_,
		    radius_  );

  
  
  return 0;
}




// ============================================================================
/// \brief: initialize the RadiationKernel (version for 2d).
// ============================================================================
int RadiationKernel::initialize( int numberOfGridPoints1_, 
				 int numberOfFields_, 
				 int numberOfModes1_, 
				 real period1_,  
				 real c_, 
				 int orderOfTimeStepping_, int numberOfPoles_,
                                 real radius_ /* =1. */ )
{
  int numberOfDimensions_=2;
  int numberOfGridPoints2_=1;
  int numberOfModes2_=1;
  real period2_=1.;
  
  return initialize(numberOfDimensions_,
		    numberOfGridPoints1_, numberOfGridPoints2_,
		    numberOfFields_, 
		    numberOfModes1_,   numberOfModes2_,
		    period1_,  period2_,
		    c_, 
		    orderOfTimeStepping_, numberOfPoles_,
		    radius_  );
  
}


// ============================================================================
/// \brief: initialize the RadiationKernel (version for 3d).
// ============================================================================
int RadiationKernel::initialize( int numberOfDimensions_,
                                 int numberOfGridPoints1_, int numberOfGridPoints2_,
				 int numberOfFields_, 
				 int numberOfModes1_,   int numberOfModes2_,
				 real period1_,  real period2_,
				 real c_, 
				 int orderOfTimeStepping_, int numberOfPoles_,
                                 real radius_ /* =1. */ )
{
  const int & debug = dbase.get<int>("debug");
  
  int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
  numberOfDimensions=numberOfDimensions_;

  numberOfGridPoints1=numberOfGridPoints1_;
  numberOfGridPoints2=numberOfGridPoints2_;

  numberOfFields=numberOfFields_;

  numberOfModes1=numberOfModes1_;
  numberOfModes2=numberOfModes2_;

  period1=period1_;
  period2=period2_;

  c=c_;
  orderOfTimeStepping=orderOfTimeStepping_;
  if( debug & 1 )
    printF("\n ############ RadiationKernel::initialize numberOfDimensions=%d, orderOfTimeStepping=%d ############\n"
	   "       numberOfGridPoints1=%d, numberOfGridPoints2=%d, numberOfModes1=%d, numberOfModes2=%d\n",
	   numberOfDimensions,orderOfTimeStepping,numberOfGridPoints1,numberOfGridPoints2,numberOfModes1,numberOfModes2);
  
  if( numberOfModes1>numberOfGridPoints1/2 )
  {
    printF("RadiationKernel::initialize: error: numberOfModes1=%d > numberOfGridPoints1/2=%d\n",
	   numberOfModes1,numberOfGridPoints1/2);
    OV_ABORT("ERROR");
  }
  if( numberOfDimensions==3 && numberOfModes2>numberOfGridPoints2/2 )
  {
    printF("RadiationKernel::initialize: error: numberOfModes2=%d > numberOfGridPoints2/2=%d\n",
	   numberOfModes2,numberOfGridPoints2/2);
    OV_ABORT("ERROR");
  }
  


  numberOfPoles=numberOfPoles_;
  
  int totalGridPoints= numberOfDimensions==2 ? numberOfGridPoints1 : numberOfGridPoints1*numberOfGridPoints2;

  delete [] ploc;
  ploc = new double [totalGridPoints];
  
  ns1=2*numberOfGridPoints1+15;
  ns2=4*numberOfGridPoints2+15;
  
  if( kernelType==planar )
  {

    if( numberOfDimensions==2 )
    {

      delete [] fold;
      fold= new double [(orderOfTimeStepping-1)*numberOfModes1*numberOfPoles*numberOfFields*2]; // complex
  
      delete [] phi;
      phi = new double [numberOfModes1*numberOfPoles*numberOfFields*2];  // complex

      delete [] amc;
      amc = new double[orderOfTimeStepping+1];
    
      delete [] fftsave1;
      fftsave1 = new double [ns1];
    }
    else
    {
      // --- bcper3dq21
      // complex*16 alpha(21),beta(21),phi(0:md1,-md2:md2,21,m)
      // complex*16 fold(0:ord-2,0:md1,-md2:md2,21,m),adon,phat,xfact,ii 
      // integer nda1,ndb1,nda2,ndb2 
      // double precision p(nda1:ndb1,nda2:ndb2,m),f(nda1:ndb1,nda2:ndb2,m)
      // double precision pl(n1,n2),pl1(n1)
      // complex*16 zl(n1,n2),zl2(n2) 
      // double precision fftsave1(ns1),fftsave2(ns2)
      // double precision amc(-1:ord-2),scl1,scl2,dt,w,len1,len2,c
	
      delete [] fold;
      fold= new double [(orderOfTimeStepping-1)*(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields*2]; // complex
  
      delete [] phi;
      phi = new double [(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields*2];  // complex

      delete [] amc;
      amc = new double[orderOfTimeStepping+1];
    
      delete [] fftsave1;
      delete [] fftsave2;
      fftsave1 = new double [ns1];
      fftsave2 = new double [ns2];

      delete [] pl1;
      pl1 = new double [numberOfGridPoints1];
      
      delete [] zl;
      zl = new double [totalGridPoints*2];  // complex

      delete [] zl2;
      zl2 = new double [numberOfGridPoints2*2];  // complex
      
      
    }
    
  }

//  double precision: p(n,m) - m fields to which the operator should be applied
//
//  double precision: f(n,m) - the m results
//
//  double precision: ploc(n) - workspace 
//
//  double precision: c - the wave speed
//
//  double precision: r - the radius
//
//  double precision: dt - the time step
//
//  integer: n the number of grid points - most efficient if it has small
//           prime factors, preferably even
//
//  integer: m the number of fields
//
//  integer: md the maximum mode used in the bc md < n/2
//
//  integer: ns>=2*n+15 
//
//  integer: ord - time-stepping order - note that the stability domain for
//                 Adams-Moulton methods gets small if this is too big
// 
//  complex*16: fold(0:ord-2,0:md,44,m) - stored values for time-stepping
//
//  complex*16: phi(0:md,44,m) - the auxiliary functions computed here
//
//  integer: npoles(0:md) - #poles - computed here
//
//  complex*16: alpha(0:md,44) - the amplitudes computed here
//
//  complex*16: beta(0:md,44) - the poles computed here
//
//  double precision: amc(-1:ord-2) - Adams-Moulton coefficients (computed here)
//                              use amcof.f
//
//  double precision: fftsave(ns) - used by fftpack - link to rffti,rfftf,rfftb

  else if( kernelType==cylindrical )
  {
    const int mdp1=numberOfModes1+1;

    numberOfPoles=44;  // for now all modes get 44 poles

    delete [] fold;
    fold= new double [(orderOfTimeStepping-1)*mdp1*numberOfPoles*numberOfFields*2]; // complex
  
    delete [] phi;
    phi = new double [mdp1*numberOfPoles*numberOfFields*2];  // complex

    delete [] amc;
    amc = new double [orderOfTimeStepping+1];
    
    delete [] fftsave1;
    fftsave1 = new double [ns1];
    
    delete [] npoles;
    npoles = new int [mdp1];
    
    delete [] alpha;
    alpha = new double [mdp1*numberOfPoles*2];  // complex

    delete [] beta;
    beta = new double [mdp1*numberOfPoles*2];  // complex

    
  }
  
  bcinit=0;
  
}


int RadiationKernel::evaluateKernel( double dt, RealArray & u, RealArray & hu )
// ========================================================================================
//
//  Assign the radiation kernel. This version used Ton Hagstrom Fortran routines
//            bcperq21   -- 2D, 21 pole approximation
//            bcperq31   -- 2D, 31 pole approximation 
//            bcper3dq21 -- 3D, 21 pole approximation
//
// This routine always assumes the arrays u and hu have the 
// appropriate 1D or 2D dimensions:
//   raxis=0 :  u(0:ny-1) or u(0:ny-1,0:nz-1)
//   raxis=1 :  u(0:nx-1) or u(0:nx-1,0:nz-1)
//   raxis=2 :  u(0:nx-1) or u(0:nx-1,0:ny-1)
//
//  This routine assumes that the input values run from u(0),...,u(numberOfGridPoints-1)
//
// 
//   The output hu(i) is defined for all i (with periodic images assigned too)
// ========================================================================================
{
  real time=getCPU();
  
  const int & debug = dbase.get<int>("debug");
  // FILE *debugFile=dbase.get<FILE*>("debugFile");
  const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");

  assert( ploc!=NULL );
  
  double *pu = u.getDataPointer();
  double *pf = hu.getDataPointer();
  
  // const int rside = dbase.get<int>("rside");
  // const int raxis = dbase.get<int>("raxis");

  // This routine always assumes the arrays u and hu have the 
  // appropriate 1D or 2D dimensions:
  //   raxis=0 :  u(0:ny-1) or u(0:ny-1,0:nz-1)
  //   raxis=1 :  u(0:nx-1) or u(0:nx-1,0:nz-1)
  //   raxis=2 :  u(0:nx-1) or u(0:nx-1,0:ny-1)
  int axis=0;  // normal axis 
  int axisp1 = 1;
  int axisp2 = 2;
  
  // make sure u and hu are the same size: 
  assert( u.getBase(0)==hu.getBase(0) && u.getBound(0)==hu.getBound(0) );
  assert( u.getBase(axisp1)==hu.getBase(axisp1) && u.getBound(axisp1)==hu.getBound(axisp1) );
  if( numberOfDimensions==3 )
  {
    assert( u.getBase(axisp2)==hu.getBase(axisp2) && u.getBound(axisp2)==hu.getBound(axisp2) );
  }
  
  

  // Leading dimensions of the u and hu arrays: 
  const int nda1=u.getBase(0)+1, ndb1=u.getBound(0)+1;  // add one (base one assumed in bcperq21d)
  const int nda2=u.getBase(1)+1, ndb2=u.getBound(1)+1;  // add one (base one assumed in bcperq21d)
  
  // printF("RadiationKernel::evaluateKernel: nda1=%d, ndb1=%d (base 1), numberOfGridPoints1=%d numberOfModes1=%d, u=[%d,%d]\n",
  // 	 nda1,ndb1,numberOfGridPoints1,numberOfModes1,
  // 	 u.getBase(0),u.getBound(0));
	 

  assert( bcinit>=0 );
  
  if( kernelType==planar )
  {
    if( numberOfDimensions==2 )
    {
      if( numberOfPoles==21 )
      {
	// 21 -pole approximation: advance the "numberOfModes" auxillary variables in time and
	// evaluate the Kernel 
	bcperq21d( nda1,ndb1,
		   *pu, *pf, *ploc, c,period1,dt,numberOfGridPoints1,numberOfFields,numberOfModes1,ns1,
		   orderOfTimeStepping,*fold,*phi,*amc, *fftsave1,bcinit);
      }
      else if( numberOfPoles==31 )
      {
	// here is the 31 pole approximation
	bcperq31d( nda1,ndb1,
		   *pu, *pf, *ploc, c,period1,dt,numberOfGridPoints1,numberOfFields,numberOfModes1,ns1,
		   orderOfTimeStepping,*fold,*phi,*amc, *fftsave1,bcinit);
      }
      else
      {
	OV_ABORT("RadiationKernel::evaluateKernel:ERROR: numberOfPoles");
      }

      if( false )
      {
	printF("RadiationKernel::AFTER bcperq21d, numberOfGridPoints1=%d, numberOfFields=%d, bcinit=%d\n",numberOfGridPoints1,numberOfFields,bcinit);
	display(u,"RadiationKernel:: u","%9.2e ");
	display(hu,"RadiationKernel:: hu","%9.2e ");
	// OV_ABORT("stop here for now");
      }
	

    }
    else
    {
      // ---- 3D ----
      if( numberOfPoles==21 )
      {
	// 21 -pole approximation: advance the "numberOfModes" auxillary variables in time and
	// evaluate the Kernel 
	// bcperq21d( nda,ndb,
	// 	   *pu, *pf, *ploc, c,period1,dt,numberOfGridPoints1,numberOfFields,numberOfModes1,ns1,
	// 	   orderOfTimeStepping,*fold,*phi,*amc, *fftsave,bcinit);

        // int nda1=nda, ndb1=nda;
        // int nda2=nda, ndb2=nda;  // fix me 
	if( debug & 1 )
	  printF("RadiationKernel:evaluateKernel **3D** period1=%g, period2=%g\n",period1,period2);

        bcper3dq21(nda1,ndb1,nda2,ndb2, *pu, *pf, *ploc, *zl, *pl1, *zl2,
		   c,
		   period1,period2, dt,
		   numberOfGridPoints1,numberOfGridPoints2,
		   numberOfFields, numberOfModes1,numberOfModes2,
		   ns1,ns2,
		   orderOfTimeStepping,*fold,*phi,*amc, *fftsave1, *fftsave2,
		   bcinit);
	if( false )
	{
	  printF("RadiationKernel::AFTER bcper3dq21\n");
	  display(u,"RadiationKernel:: u","%9.2e ");
	  display(hu,"RadiationKernel:: hu","%9.2e ");
	  // OV_ABORT("stop here for now");
	}
	
  
      }
      else
      {
	OV_ABORT("RadiationKernel::evaluateKernel:ERROR: numberOfPoles");
      }

    }
    
  }
  else if( kernelType==cylindrical )
  {
    bccyld(nda1,ndb1,
           *pu, *pf, *ploc, c,radius, dt,numberOfGridPoints1,numberOfFields,numberOfModes1,ns1,
           orderOfTimeStepping, *fold,*phi,
	   *npoles, *alpha, *beta,*amc,*fftsave1,bcinit);
  }
  else
  {
    printF("ERROR: un-implemented kernelType=%i\n",kernelType);
    Overture::abort("error");
  }
  

  if( false )
  {
    printF("RadKer:evaluateKernel:bcperq: numberOfFields=%d, numberOfGridPoints1=%d\n",
           numberOfFields,numberOfGridPoints1);
  }
  
  if( debug & 2 )
  {
    // fprintf(debugFile,"RadKer:evaluateKernel:bcperq\n");
    // fprintf(debugFile,"fold=");
    // // fold= new double [(orderOfTimeStepping-1)*mdp1*numberOfPoles*numberOfFields*2]; // complex
    // int orderInTimeMinus1=orderOfTimeStepping-1;
    // #define foldv(l,k,pole,mc) fold2c[(l)+orderInTimeMinus1*( (k)+ numberOfModes1*( (pole) + numberOfPoles*(mc)))]
    // for( int l=0; l<
    // for( int k=0; k<numberOfModes1; k++ )
      
  }
  


  // --- assign periodic images ---
  if( numberOfDimensions==2 )
  {
    for( int i=hu.getBase(0); i<0; i++ )
    {
      for( int j=0; j<numberOfFields; j++ )
	hu(i,j)=hu(i+numberOfGridPoints1,j);
    }
    for( int i=numberOfGridPoints1; i<=hu.getBound(0); i++ )
    {
      for( int j=0; j<numberOfFields; j++ )
	hu(i,j)=hu(i-numberOfGridPoints1,j);
    }
  }
  else
  {
    // --- assign periodic images ----
    // ** CHECK ME**

    Range M = numberOfFields;

    Range I3= Range(hu.getBase(1),hu.getBound(1));
    Range I2= Range(hu.getBase(0),-1);
    hu(I2,I3,M) = hu(I2+numberOfGridPoints1,I3,M);    // left-ghost = right
    I2 = Range(numberOfGridPoints1,hu.getBound(0));
    hu(I2,I3,M) = hu(I2-numberOfGridPoints1,I3,M);    // right-ghost = left 
    
    I2= Range(hu.getBase(0),hu.getBound(0));
    I3= Range(hu.getBase(1),-1);
    hu(I2,I3,M) = hu(I2,I3+numberOfGridPoints2,M);    // bottom-ghost = top 
    I3 = Range(numberOfGridPoints2,hu.getBound(1));
    hu(I2,I3,M) = hu(I2,I3-numberOfGridPoints2,M);    // top-ghost = bottom
    
    
  }
  
  cpuTime+=getCPU()-time;
  return 0;
}


// ------------- Keep routines that use complex at the bottom of the file to avoid name conflicts (real) ---------

#include <complex>
typedef ::real LocalReal;
typedef ::real OV_real;
typedef std::complex<LocalReal> Complex;


// ======================================================================================================
/// \brief Destructor
// ======================================================================================================
RadiationKernel::~RadiationKernel()
{

  if( dbase.has_key("debugFile") )
  {
    FILE * debugFile=dbase.get<FILE*>("debugFile");
    if( debugFile!=NULL )
      fclose(debugFile);
  }


  delete [] ploc;
  delete [] fold;
  delete [] phi;
  delete [] amc;
  delete [] fftsave1;
  delete [] fftsave2;
  delete [] alpha;
  delete [] beta;
  delete [] npoles;

  delete [] zl;
  delete [] pl1;
  delete [] zl2;

  if( dbase.has_key("phic") )
    delete [] dbase.get<Complex*>("phic");

  if( dbase.has_key("foldc") )
    delete [] dbase.get<Complex*>("foldc");

  if( dbase.has_key("alphac") )
    delete [] dbase.get<Complex*>("alphac");

  if( dbase.has_key("betac") )
    delete [] dbase.get<Complex*>("betac");

  if( dbase.has_key("phi2c") )
    delete [] dbase.get<Complex*>("phi2c");

  if( dbase.has_key("fold2c") )
    delete [] dbase.get<Complex*>("fold2c");

}


// =====================================================================================================
/// \brief initialize alpha and beta arrays defining the poles in the non-local kernel approximation
// ======================================================================================================
int RadiationKernel::initializePoles()
{
  // ---- alpha and beta define the fit of the convolution kernel in Laplace transform space ----
  dbase.put<Complex*>("alphac")=NULL;
  Complex *& alphac = dbase.get<Complex*>("alphac");
  alphac = new Complex [numberOfPoles];                    

  dbase.put<Complex*>("betac")=NULL;
  Complex *& betac = dbase.get<Complex*>("betac");
  betac = new Complex [numberOfPoles];                     

#define alpha(i) alphac[(i-1)]
#define beta(i)  betac[(i-1)]

 Complex I(0.,1.);

  if( numberOfPoles == 21 )
  {
    // values from bcperq21.f (base 1)
    alpha(1) = -.2410467618025768E-06 + (  -.2431987763837349E-06)*I;
    alpha(2) = -.2410467618025768E-06 + (   .2431987763837349E-06)*I; 
    alpha(3) = -.1617695923999794E-05 + (  -.1638622585172068E-05)*I;
    alpha(4) = -.1617695923999794E-05 + (   .1638622585172068E-05)*I;
    alpha(5) = -.7723476507531262E-05 + (  -.7878743138182415E-05)*I;
    alpha(6) = -.7723476507531262E-05 + (   .7878743138182415E-05)*I;
    alpha(7) = -.3400304516975200E-04 + (  -.3510673092397324E-04)*I;
    alpha(8) = -.3400304516975200E-04 + (   .3510673092397324E-04)*I; 
    alpha(9) = -.1454893381589074E-03 + (  -.1535469093409158E-03)*I;
    alpha(10)= -.1454893381589074E-03 + (   .1535469093409158E-03)*I;
    alpha(11)= -.6104572904148162E-03 + (  -.6733883694898616E-03)*I;
    alpha(12)= -.6104572904148162E-03 + (   .6733883694898616E-03)*I;
    alpha(13)= -.2473202929583869E-02 + (  -.3011442350813045E-02)*I;
    alpha(14)= -.2473202929583869E-02 + (   .3011442350813045E-02)*I;
    alpha(15)= -.8964957513027030E-02 + (  -.1398751873403249E-01)*I;
    alpha(16)= -.8964957513027030E-02 + (   .1398751873403249E-01)*I;
    alpha(17)= -.1846252520037211E-01 + (  -.6565858806543060E-01)*I;
    alpha(18)= -.1846252520037211E-01 + (   .6565858806543060E-01)*I;
    alpha(19)=  .9181095934161065E-01 + (  -.2076825633238755E+00)*I;
    alpha(20)=  .9181095934161065E-01 + (   .2076825633238755E+00)*I;
    alpha(21)=  .3787484004895032E+00 + (   .0000000000000000E+00)*I;

    beta(1) = -.4998142304334231E-04 + (   .9999998607359947E+00)*I;
    beta(2) = -.4998142304334231E-04 + (  -.9999998607359947E+00)*I;
    beta(3) = -.2501648855535112E-03 + (   .9999990907954994E+00)*I;
    beta(4) = -.2501648855535112E-03 + (  -.9999990907954994E+00)*I;
    beta(5) = -.8021925048752190E-03 + (   .9999958082358295E+00)*I;
    beta(6) = -.8021925048752190E-03 + (  -.9999958082358295E+00)*I;
    beta(7) = -.2263515963206483E-02 + (   .9999820162287431E+00)*I;
    beta(8) = -.2263515963206483E-02 + (  -.9999820162287431E+00)*I;
    beta(9) = -.6112737916031916E-02 + (   .9999224860282032E+00)*I;
    beta(10)= -.6112737916031916E-02 + (  -.9999224860282032E+00)*I;
    beta(11)= -.1625071664643320E-01 + (   .9996497460330479E+00)*I;
    beta(12)= -.1625071664643320E-01 + (  -.9996497460330479E+00)*I;
    beta(13)= -.4295328074381198E-01 + (   .9982864080633248E+00)*I;
    beta(14)= -.4295328074381198E-01 + (  -.9982864080633248E+00)*I;
    beta(15)= -.1129636068874967E+00 + (   .9907617913485537E+00)*I;
    beta(16)= -.1129636068874967E+00 + (  -.9907617913485537E+00)*I;
    beta(17)= -.2902222956062986E+00 + (   .9462036470847180E+00)*I;
    beta(18)= -.2902222956062986E+00 + (  -.9462036470847180E+00)*I;
    beta(19)= -.6548034445533449E+00 + (   .7077228221122372E+00)*I;
    beta(20)= -.6548034445533449E+00 + (  -.7077228221122372E+00)*I;
    beta(21)= -.9345542777004186E+00 + (   .0000000000000000E+00)*I;

  }
  else if( numberOfPoles==31 )
  {
    // values from bcperq31.f (base 1)
    alpha(1) = -.7248698190879261E-07 + (-.7266856459887401E-07)*I;
    alpha(2) = -.7248698190879261E-07 + ( .7266856459887401E-07)*I;
    alpha(3) = -.3761817129252043E-06 + (-.3769314065715886E-06)*I;
    alpha(4) = -.3761817129252043E-06 + ( .3769314065715886E-06)*I;
    alpha(5) = -.1264056865588544E-05 + (-.1266668072310468E-05)*I;
    alpha(6) = -.1264056865588544E-05 + ( .1266668072310468E-05)*I;
    alpha(7) = -.3738466040644283E-05 + (-.3758052767416740E-05)*I;
    alpha(8) = -.3738466040644283E-05 + ( .3758052767416740E-05)*I;
    alpha(9) = -.1053054276380811E-04 + (-.1064172701253104E-04)*I;
    alpha(10)= -.1053054276380811E-04 + ( .1064172701253104E-04)*I;
    alpha(11)= -.2902787565325144E-04 + (-.2948662740491706E-04)*I;
    alpha(12)= -.2902787565325144E-04 + ( .2948662740491706E-04)*I;
    alpha(13)= -.7907555696587643E-04 + (-.8101521788775704E-04)*I;
    alpha(14)= -.7907555696587643E-04 + ( .8101521788775704E-04)*I;
    alpha(15)= -.2136714788000509E-03 + (-.2223118737627421E-03)*I;
    alpha(16)= -.2136714788000509E-03 + ( .2223118737627421E-03)*I;
    alpha(17)= -.5718490556660803E-03 + (-.6120497162181238E-03)*I;
    alpha(18)= -.5718490556660803E-03 + ( .6120497162181238E-03)*I;
    alpha(19)= -.1500648374605461E-02 + (-.1702427459760661E-02)*I;
    alpha(20)= -.1500648374605461E-02 + ( .1702427459760661E-02)*I;
    alpha(21)= -.3777846647633256E-02 + (-.4822679646038999E-02)*I;
    alpha(22)= -.3777846647633256E-02 + ( .4822679646038999E-02)*I;
    alpha(23)= -.8612337122259847E-02 + (-.1397170663788553E-01)*I;
    alpha(24)= -.8612337122259847E-02 + ( .1397170663788553E-01)*I;
    alpha(25)= -.1436234309979386E-01 + (-.4057035360912333E-01)*I;
    alpha(26)= -.1436234309979386E-01 + ( .4057035360912333E-01)*I;
    alpha(27)=  .5115479570066042E-02 + (-.1052205533094523E+00)*I;
    alpha(28)=  .5115479570066042E-02 + ( .1052205533094523E+00)*I;
    alpha(29)=  .1353590959088723E+00 + (-.1619011637672531E+00)*I;
    alpha(30)=  .1353590959088723E+00 + ( .1619011637672531E+00)*I;
    alpha(31)=  .2773935622389255E+00 + ( .0000000000000000E+00)*I;

    beta(1) = -.2295682215845089E-04 + ( .9999999803841022E+00)*I;
    beta(2) = -.2295682215845089E-04 + (-.9999999803841022E+00)*I;
    beta(3) = -.1023267035480556E-03 + ( .9999999247031288E+00)*I;
    beta(4) = -.1023267035480556E-03 + (-.9999999247031288E+00)*I;
    beta(5) = -.2744661422580783E-03 + ( .9999998277644731E+00)*I;
    beta(6) = -.2744661422580783E-03 + (-.9999998277644731E+00)*I;
    beta(7) = -.6185331183100160E-03 + ( .9999993934564934E+00)*I;
    beta(8) = -.6185331183100160E-03 + (-.9999993934564934E+00)*I;
    beta(9) = -.1293014263135570E-02 + ( .9999973246515538E+00)*I;
    beta(10)= -.1293014263135570E-02 + (-.9999973246515538E+00)*I;
    beta(11)= -.2607493030160015E-02 + ( .9999910719461364E+00)*I;
    beta(12)= -.2607493030160015E-02 + (-.9999910719461364E+00)*I;
    beta(13)= -.5163854131302370E-02 + ( .9999735586993596E+00)*I;
    beta(14)= -.5163854131302370E-02 + (-.9999735586993596E+00)*I;
    beta(15)= -.1013381010454508E-01 + ( .9999195856921692E+00)*I;
    beta(16)= -.1013381010454508E-01 + (-.9999195856921692E+00)*I;
    beta(17)= -.1979469547015627E-01 + ( .9997465694790964E+00)*I;
    beta(18)= -.1979469547015627E-01 + (-.9997465694790964E+00)*I;
    beta(19)= -.3855416231134041E-01 + ( .9991412761425118E+00)*I;
    beta(20)= -.3855416231134041E-01 + (-.9991412761425118E+00)*I;
    beta(21)= -.7487966608350244E-01 + ( .9968954030259501E+00)*I;
    beta(22)= -.7487966608350244E-01 + (-.9968954030259501E+00)*I;
    beta(23)= -.1446954921242388E+00 + ( .9885229324572602E+00)*I;
    beta(24)= -.1446954921242388E+00 + (-.9885229324572602E+00)*I;
    beta(25)= -.2757301824230285E+00 + ( .9580367554955076E+00)*I;
    beta(26)= -.2757301824230285E+00 + (-.9580367554955076E+00)*I;
    beta(27)= -.5026369164103634E+00 + ( .8538756767749224E+00)*I;
    beta(28)= -.5026369164103634E+00 + (-.8538756767749224E+00)*I;
    beta(29)= -.8041374504940593E+00 + ( .5573022536695248E+00)*I;
    beta(30)= -.8041374504940593E+00 + (-.5573022536695248E+00)*I;
    beta(31)= -.9694745780110367E+00 + ( .0000000000000000E+00)*I;

  }
  else
  {
    printF("RadiationKernel::initializePoles: ERROR: numberOfPoles=%d not supported\n",numberOfPoles);
    OV_ABORT("error");
  }
	


      
  return 0;
}



// ========================================================================================
//
/// \brief  Evaluate the radiation kernel. *New* serial version, gett ready for for parallel 
///       THIS VERSION NOT USED
/// \param u (input) :  
/// \param hu (output) : the output hu(i) is defined for all i (with periodic images assigned too)
//  
//
//  This routine assumes that the input values run from u(0),...,u(numberOfGridPoints-1)
//
/// \Note: See paper
//      author="B. Alpert and L. Greengard and T. Hagstrom",
//      title="Nonreflecting Boundary Conditions for the Time-Dependent Wave Equation",
//      journal=JCP,
//      volume="180",
//      year="2002",
//      pages="270--296")
// 
// ========================================================================================
int RadiationKernel::evaluateKernelNew( double dt, RealArray & u, RealArray & hu )
{
  LocalReal time=getCPU();
  
  const int & debug = dbase.get<int>("debug");
  const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");

  int & bcInitComplex = dbase.get<int>("bcInitComplex");

  assert( ploc!=NULL );
  
  double *pu = u.getDataPointer();
  double *pf = hu.getDataPointer();
  
  // We assume normal direction is axis=0 for now:
  int axis=0;  // normal axis 
  int axisp1 = 1;
  int axisp2 = 2;
  
  // make sure u and hu are the same size: 
  assert( u.getBase(0)==hu.getBase(0) && u.getBound(0)==hu.getBound(0) );
  assert( u.getBase(axisp1)==hu.getBase(axisp1) && u.getBound(axisp1)==hu.getBound(axisp1) );
  if( numberOfDimensions==3 )
  {
    assert( u.getBase(axisp2)==hu.getBase(axisp2) && u.getBound(axisp2)==hu.getBound(axisp2) );
  }
  

  // Leading dimensions of the u and hu arrays: 
  const int nda1=u.getBase(0)+1, ndb1=u.getBound(0)+1;  // add one (base one assumed in bcperq21d)
  const int nda2=u.getBase(1)+1, ndb2=u.getBound(1)+1;  // add one (base one assumed in bcperq21d)
  
  assert( bcInitComplex>=0 );
  
  if( kernelType==planar )
  {
    // ------------- Planar Kernel : flat boundary -------------

    // **NEW WAY**
    // Start from code in bcperq21.f
    // We solve for auxilary vaiables phi_j 
    //    (d/dt - c beta_j w) phi_j = c alpha_j w^2 phat 
    //    phat(t) = 
    //    fha(t)t = sum phi_j   :    FourierTranform of solution 

    Complex I(0.,1.);

    const int orderInTime=orderOfTimeStepping; 
    const int orderInTimeMinus1=orderInTime-1;

    int md = numberOfModes1;    // number of Fourier modes to keep, max = n1/2
    if( numberOfDimensions==3 )
      md *= numberOfModes2;

    int numComp=numberOfFields;  
    int n1=numberOfGridPoints1;  // number of grid points
    int ns1 =2*n1+15;            // for real FFT
    int ns1c=4*n1+15;            // for complex FFT
       
    int n2=numberOfGridPoints2;  // number of grid points
    int ns2 =2*n2+15;            // for real FFT
    int ns2c=4*n2+15;            // for complex FFT
    


    if( !dbase.has_key("phic") )
    {
      dbase.put<Complex*>("phic")=NULL;
      Complex *& phic = dbase.get<Complex*>("phic");
      if( numberOfDimensions==2 )
        phic = new Complex [md*numberOfPoles*numComp];            
      else 
        phic = new Complex [(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields];            

      dbase.put<Complex*>("foldc")=NULL;
      Complex *&foldc = dbase.get<Complex*>("foldc");
      if( numberOfDimensions==2 )
        foldc = new Complex[orderInTimeMinus1*md*numberOfPoles*numComp];
      else
	foldc = new Complex[orderInTimeMinus1*(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields];


      // phi = new double [(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields*2];  // complex
      // fold= new double [(orderOfTimeStepping-1)*(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields*2]; // complex



      dbase.put<RealArray>("amc");
      RealArray & amc = dbase.get<RealArray>("amc");
      amc.redim(Range(-1,orderInTime-2));

      dbase.put<RealArray>("fftsave");
      RealArray & fftsave = dbase.get<RealArray>("fftsave");
      fftsave.redim(ns1);

      dbase.put<RealArray>("fftsavec");
      RealArray & fftsavec = dbase.get<RealArray>("fftsavec");
      fftsavec.redim(ns1c);

      if( numberOfDimensions==3 )
      {
	dbase.put<RealArray>("fftsave2");
	RealArray & fftsave2 = dbase.get<RealArray>("fftsave2");
	fftsave2.redim(ns2);

	dbase.put<RealArray>("fftsave2c");
	RealArray & fftsave2c = dbase.get<RealArray>("fftsave2c");
	fftsave2c.redim(ns2c);

      }
      

      // dbase.put<Complex*>("fftsavec")=NULL;
      // Complex *&fftsavec = dbase.get<Complex*>("fftsavec");
      // fftsavec = new Complex[4*n1+15]; // *** delete me ***

      // --- initialize alpha and beta arrays defining the poles in the non-local kernel approximation
      initializePoles();

    }
      
    // ---- alpha and beta define the fit of the convolution kernel in Laplace transform space ----
    Complex *& alphac = dbase.get<Complex*>("alphac");
    Complex *& betac = dbase.get<Complex*>("betac");
#define alpha(i) alphac[(i-1)]
#define beta(i)  betac[(i-1)]



    if( numberOfDimensions==2 )  // *new way* Finish me for 3D
    {


      Complex *& phic = dbase.get<Complex*>("phic");
#define phi(k,pole,mc) phic[(k)+md*( (pole) + numberOfPoles*(mc))]
       
      // ComplexArray fold(0:orderInTime-1,md,numberOfPoles,numComp) : 
      Complex *&foldc = dbase.get<Complex*>("foldc");
#define fold(l,k,pole,mc) foldc[(l)+orderInTimeMinus1*( (k)+md*( (pole) + numberOfPoles*(mc)))]

       
      RealArray & fftsave = dbase.get<RealArray>("fftsave"); 
      RealArray &fftsavec = dbase.get<RealArray>("fftsavec");     

      RealArray & amc = dbase.get<RealArray>("amc");
       
       
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int k=0; k<numberOfModes1; k++ )
	    {
	      phi(k,pole,mc)=0.;
	      for( int l=0; l<=orderInTime-2; l++ )
	      {
		fold(l,k,pole,mc)=0.;
	      }
	    }
	  }
	}
	 
	RFFTI(n1,fftsave(0));

        // complex FFT: works-space: nsc > =4*n1+15 
        ZFFTI(n1,fftsavec(0));

	bcInitComplex=1;
	 
      } // end initialize

      RealArray ploc(n1); // , p(n1,numComp);
      // RealArray f(n1);  // put answer here for now 
       
      // For complex FFT's
      //  zlocr(i) = real-part
      //  zloci(i) = imag-part
      const int n1c=2*n1;
      LocalReal *pzloc = new LocalReal [n1*2];  // holds complex transforms 
#define zlocr(i) pzloc[2*(i)]
#define zloci(i) pzloc[1+2*(i)]

      LocalReal scl=twoPi/period1;
      // printF("scl=%20.14e\n",scl);
      Range I1=n1;
       
      // printF("u.getBase(0)=%d\n",u.getBase(0));
      

      // loop over components: 
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {
	// ploc(I1)=p(I1,mc);    // solution for this component 
	ploc(I1)=u(I1,mc);    // solution for this component 
	RFFTF( n1,ploc(0),fftsave(0) );   // ploc = uHat 

        for( int i=0; i<n1; i++ )
	{
	  zlocr(i) = u(i,mc); // real part of data
	  zloci(i)=0.;        // imag part of data 
	}
	
	  
	ZFFTF( n1,pzloc[0],fftsavec(0) );   
	if( false )
	{
	  printf("After Forward: n1=%d (zloc holds results for uHat(k) and uHat(-k)=conj(uHat(k))\n ploc=[",n1);
	  for( int i=0; i<n1; i++ ){ printF("%13.6e,",ploc(i)); }   // 
          printF("]\n");
	  printf(" zloc=[");
	  for( int i=0; i<n1; i++ ){ printF("%13.6e,%13.6e,",zlocr(i),zloci(i)); }   // 
          printF("]\n");
	  
	}
	
        // convert complex FFT to real:
	//   pHat( k) = pzloc[k] ,    k=0,1,2,...,n1/2-1
	//   pHat(-k) = pzloc[n-k]    k=0,1,2,...,n1/2


	for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	{
	  for( int k=0; k<numberOfModes1; k++ )   // "k=1" loop over actual Fourier modes k that we compute 
	  {
	    const int kp = k+1;
	    LocalReal w = scl*kp;
	    Complex expBeta = exp(c*beta(pole+1)*dt*w);  // beta's are complex 
	    Complex xfact = expBeta;
	    Complex adon = xfact*phi(k,pole,mc);         // add-on 
	    // printF("* adon=[%12.4e,%12.4e]\n",std::real(adon),std::imag(adon));
	    for( int l=0; l<=orderInTime-2; l++ )
	    {
	      adon += xfact*amc(l)*fold(l,k,pole,mc);
	      xfact *= expBeta;
	      // printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
	    }
	    Complex phat = ploc(2*kp-1) + ploc(2*kp)*I;  // use real FFT

            phat = zlocr(kp) + zloci(kp)*I;              // use complex FFT
	     

	    adon += c*dt*amc(-1)*alpha(pole+1)*w*w*phat;
	    phi(k,pole,mc)=adon;

	    // printF("expBeta=[%12.4e,%12.4e]\n",std::real(expBeta),std::imag(expBeta));
	    //printF("phi(%d,%d,%d)=[%12.4e,%12.4e]\n",k,pole,mc,std::real(phi(k,pole,mc)),std::imag(phi(k,pole,mc)));
	    
	    for( int l=orderInTime-2; l>=1; l-- )
	    {
	      fold(l,k,pole,mc)=fold(l-1,k,pole,mc);
	    }
	    fold(0,k,pole,mc)=c*dt*w*w*alpha(pole+1)*phat;   // alpha's are complex 
	  }  // end for k
	}  // end for pole

	ploc=0.;  // fill in hHat 
	for( int pole=0; pole<numberOfPoles; pole++ )
	{
	  for( int k=0; k<numberOfModes1; k++ )   // loop over Fourier modes we compute 
	  {
	    const int kp = k+1;
	    ploc(2*kp-1) += std::real(phi(k,pole,mc));    
	    ploc(2*kp  ) += std::imag(phi(k,pole,mc)); 

	  }
	}

        // real FFT:      cos(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  
        // complex FFT:   cos(0*x) sin(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  

        // convert complex FFT to real:  (for n1 even)
	//   pHat( k) = pzloc[k] ,    k=0,1,2,...,n1/2-1
	//   pHat(-k) = pzloc[n-k]    k=  1,2,...,n1/2
        
	for( int i=0; i<2*n1; i++ )
	{
	  pzloc[i]=0.;
	}
	
        zlocr(0) = ploc(0);
        zloci(0) = 0.;
        int j=1;
        for( int i=1; i<n1/2; i++ )
	{
          // k=1,2,...,n1/2-1
	  zlocr(i) = ploc(j  );
	  zloci(i) = ploc(j+1);
          // negative k are the complex conjugate:
	  zlocr(n1-i) =  ploc(j  );
	  zloci(n1-i) = -ploc(j+1);

	  j+=2;
	  
	}
        // if( (n1 % 2) != 0 )
	// {
	//   // printF("WARNING: n1=%d is NOT EVEN\n",n1);
        //   zlocr(n1)=0.;
        //   zloci(n1)=0.;
	if( false )
	{
	  printF("Before backward:\n ploc=[");
	  for( int i=0; i<n1; i++ ){ printF("%16.9e,",ploc(i)); }   // 
          printF("]\n");
	  printF(" zloc=[");
	  for( int i=0; i<n1; i++ ){ printF("%16.9e,%16.9e,",zlocr(i),zloci(i)); }   // 
          printF("]\n");
	}


	RFFTB( n1, ploc(0),fftsave(0) );
        // use real FFT result 
	hu(I1,mc) = ploc(I1)/n1;
	     
	  
	// }

        ZFFTB( n1,pzloc[0],fftsavec(0) );  	       
	if( false )
	{
	  printf("AFter backward: (print only real(zloc))\n ploc=[");
	  for( int i=0; i<n1; i++ ){ printF("%16.9e,",ploc(i)); }   // 
          printF("]\n");
	  printf("zlocr=[");
	  for( int i=0; i<n1; i++ ){ printF("%16.9e,",zlocr(i)); }   // 
          printF("]\n");
	  
	}
	// use complex FFT result 
	for( int i=0; i<n1; i++ )
	{
          hu(i,mc) = zlocr(i)/n1;
	}
	

      }  // end for component mc
        
     
      delete [] pzloc;

    }
    else if( numberOfDimensions==3 )
    {
      // ---------- *NEW* 3D version for parallel ------------

      if( debug & 1 )
	printF("$$$$$$$$$$$$ ++ RADIATION KERNEL -- NEW 3D VERSION bcInitComplex=%d $$$$$$$$$$$$$\n",bcInitComplex);

      #ifdef USE_PPP
        printF(" ----> RUNNING IN PARALLEL <-----\n");
      #endif

      // phi(0:numberOfModes1,-numberOfModes2:numberOfModes2,numberOfPoles,numFields)
      const int md2=2*numberOfModes2+1;
      const int numberOfModes1p1 = numberOfModes1+1;
      
      Complex *& phic = dbase.get<Complex*>("phic");
#define phi(k1,k2,pole,mc) phic[(k1)+numberOfModes1p1*((k2+numberOfModes2)+md2*( (pole) + numberOfPoles*(mc)))]
       
      // ComplexArray fold(0:orderInTime-1,md,numberOfPoles,numComp) : 
      Complex *&foldc = dbase.get<Complex*>("foldc");
#define fold(l,k1,k2,pole,mc) foldc[(l)+orderInTimeMinus1*( (k1)+numberOfModes1p1*( (k2+numberOfModes2)+md2*( (pole) + numberOfPoles*(mc))))]

       
      RealArray & fftsave = dbase.get<RealArray>("fftsave"); 
      RealArray &fftsavec = dbase.get<RealArray>("fftsavec");     
      RealArray &fftsave2c = dbase.get<RealArray>("fftsave2c");     

      RealArray & amc = dbase.get<RealArray>("amc");
       
       
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int k2=-numberOfModes2; k2<=numberOfModes2; k2++ )
	    {
	      for( int k1=0; k1<=numberOfModes1; k1++ )
	      {
		phi(k1,k2,pole,mc)=0.;
		for( int l=0; l<=orderInTime-2; l++ )
		{
		  fold(l,k1,k2,pole,mc)=0.;
		}
	      }
	    }
	  }
	}
	 
	RFFTI(n1,fftsave(0));
        ZFFTI(n2,fftsave2c(0));

	bcInitComplex=1;
	 
      } // end initialize

      

      RealArray pl1(n1);   // hold temp values 

      // complex*16 zl(n1,n2),zl2(n2)        
      Complex * pzl = new Complex [n1*n2];  
#define zl(i1,i2) pzl[(i1)+n1*(i2)]

      // Could not pas this to fortran (?)
//       Complex * pzl2 = new Complex [n2]; 
// #define zl2(i2) pzl2[(i2)]

      LocalReal *pzl2 = new LocalReal [n2*2];  // holds complex transforms 
#define zl2r(i) pzl2[2*(i)]
#define zl2i(i) pzl2[1+2*(i)]      

//       // For complex FFT's
//       //  zlocr(i) = real-part
//       //  zloci(i) = imag-part
//       const int n1c=2*n1;
//       LocalReal *pzloc = new LocalReal [n1*2];  // holds complex transforms 
// #define zlocr(i) pzloc[2*(i)]
// #define zloci(i) pzloc[1+2*(i)]

      LocalReal scl1=twoPi/period1;
      LocalReal scl2=twoPi/period2;

      // printF("scl=%20.14e\n",scl);
      Range I1=n1, I2=n2;
       

      // loop over components: 
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {
        // ---- FFT's in direction 1------
	for( int i2=0; i2<n2; i2++ )
	{
	  pl1(I1)=u(I1,i2,mc);    // solution for this component 

	  RFFTF( n1,pl1(0),fftsave(0) );   // ploc = uHat 

          zl(0,i2)=pl1(0);
          if( (n1 % 2) == 0 )
	  {
            for( int i1=1; i1<n1/2; i1++ )
	    {
              zl(i1,i2)         =pl1(2*i1-1) + I*pl1(2*i1);  // check me 
              zl((n1/2)+i1-1,i2)=pl1(2*i1-1) - I*pl1(2*i1);
	    }
            zl(n1-1,i2)=pl1(n1-1);
	  }
          else
	  {
	    for( int i1=1; i1<(n1/2)+1; i1++ )
	    {
              zl(i1,i2)       =pl1(2*i1-1) + I*pl1(2*i1);   // check me 
              zl((n1/2)+i1,i2)=pl1(2*i1-1) - I*pl1(2*i1);
	    }
	  }

	}
	if( false )
	{
	  printF("Stage 1: after initial FFTs, zl:\n");
	  for( int i2=0; i2<n2; i2++ )
	  {
	    for( int i1=0; i1<n1; i1++ )
	    {
	      printF("[%9.3e,%9.3e] ",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	    }
	    printF("\n");
	  }
	}
	

        // ---- FFT's in direction 2------
	for( int i1=0; i1<n1; i1++ )
	{
	  for( int i2=0; i2<n2; i2++ )
	  {
	    // zl2(i2)=zl(i1,i2);
	    zl2r(i2) = std::real(zl(i1,i2));
	    zl2i(i2) = std::imag(zl(i1,i2));
	    
	  }
	  
	  // ZFFTF( n2,zl2[0],fftsave2c(0) );
	  ZFFTF( n2,pzl2[0],fftsave2c(0) );

          for( int i2=0; i2<n2; i2++ )
	  {
            // zl(i1,i2)=zl2(i2);
            zl(i1,i2)=zl2r(i2) + zl2i(i2)*I;
	  }
	  
	}
	if( false )
	{
	  printF("Stage 1: after 2nd FFTs, zl:\n");
	  for( int i2=0; i2<n2; i2++ )
	  {
	    for( int i1=0; i1<n1; i1++ )
	    {
	      printF("[%9.3e,%9.3e] ",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	    }
	    printF("\n");
	  }
	}
	
	for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	{
	  for( int k2=0; k2<n2; k2++ )   // "k=1" loop over actual Fourier modes k that we compute 
	  {
            int k2p=k2+1;
            int k2t;
	    if( k2p<=(numberOfModes2+1) )
	    {
	      k2t=k2p-1;
	    }
	    else if( (k2p-1-n2) >= (-numberOfModes2) )
	    {
	      k2t=k2p-1-n2;
	    }
	    else
	    {
	      continue;
	    }
	    for( int k1t=0; k1t<=numberOfModes1; k1t++ )
	    {
              int k1=k1t;  // check ??
	      int k1p=k1t+1;

	      LocalReal w= sqrt( SQR(scl1*k1t) + SQR(scl2*k2t) );

	      Complex expBeta = exp(c*beta(pole+1)*dt*w);  // beta's are complex 
	      Complex xfact = expBeta;
	      Complex adon = xfact*phi(k1t,k2t,pole,mc);         // add-on
	      for( int l=0; l<=orderInTime-2; l++ )
	      {
		adon += xfact*amc(l)*fold(l,k1t,k2t,pole,mc);
		xfact *= expBeta;
		// printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
	      }
	      adon += c*dt*amc(-1)*alpha(pole+1)*w*w*zl(k1,k2);  // check k1,k2
              phi(k1t,k2t,pole,mc)=adon;

	      for( int l=orderInTime-2; l>=1; l-- )
	      {
		fold(l,k1t,k2t,pole,mc)=fold(l-1,k1t,k2t,pole,mc);
	      }
	      fold(0,k1t,k2t,pole,mc)=c*dt*w*w*alpha(pole+1)*zl(k1,k2);   // alpha's are complex 
	    }
	  } // end for k2 
	} // end for pole 

	for( int i1=0; i1<n1; i1++ )
	{
          for( int i2=0; i2<n2; i2++ )
            zl(i1,i2)=0.;
	}

        // ---- zl(k1,k2) = SUM_pole  phi(k1t,k2t,pole,mc); 
	for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	{
	  for( int k1t=0; k1t<=numberOfModes1; k1t++ )
	  {
	    int k1=k1t;
	    for( int k2t=0; k2t<=numberOfModes2; k2t++ )
	    {
	      int k2=k2t;
	      zl(k1,k2) += phi(k1t,k2t,pole,mc);
	    }
	    for( int k2t=-numberOfModes2; k2t<=-1; k2t++ )
	    {
              int k2 = n2 + k2t;  // check 
	      zl(k1,k2) += phi(k1t,k2t,pole,mc);
	    }
	  }
	}
	// end for pole 
	
        // ------------- INVERSE FFTs in DIRECTION 2  -----------------
        for( int k1t=0; k1t<=numberOfModes1; k1t++ )
	{
	  int k1=k1t;
          for( int k2=0; k2<n2; k2++ )
	  {
            // zl2(k2)=zl(k1,k2);
	    zl2r(k2) = std::real(zl(k1,k2));
	    zl2i(k2) = std::imag(zl(k1,k2));

	  }

          // ZFFTB(n2,zl2[0],fftsave2c(0));
          ZFFTB(n2,pzl2[0],fftsave2c(0));

	  for( int k2=0; k2<n2; k2++ )
	  {
	    // zl(k1,k2)=zl2(k2);
	    zl(k1,k2)=zl2r(k2) + zl2i(k2)*I;
	  }
        }
      
        // ------------- INVERSE FFTs in DIRECTION 1  -----------------
	for( int i2=0; i2<n2; i2++ )
	{
          pl1(0)=std::real(zl(0,i2));
	  for( int k1=1; k1<=numberOfModes1; k1++ ) // check 
	  {
            pl1(2*k1-1)=std::real(zl(k1,i2));  // check 
            pl1(2*k1  )=std::imag(zl(k1,i2));
	  }
	  for( int k1=2*numberOfModes1+1; k1<n1; k1++ )
	  {
	    pl1(k1)=0.;
	  }

          RFFTB( n1,pl1(0),fftsave(0));

	  for( int i1=0; i1<n1; i1++ )
	  {
	    hu(i1,i2,mc)=pl1(i1)/(n1*n2);
	  }
	} // end for i2


      }  // end for mc 
	
      delete [] pzl;
      delete [] pzl2;



    }
    else
    {
      OV_ABORT("ERROR: unexpected numberOfDimensions");
    }
    
    
  }
  else if( kernelType==cylindrical )
  {
    bccyld(nda1,ndb1,
           *pu, *pf, *ploc, c,radius, dt,numberOfGridPoints1,numberOfFields,numberOfModes1,ns1,
           orderOfTimeStepping, *fold,*phi,
	   *npoles, *alpha, *beta,*amc,*fftsave1,bcinit);
  }
  else
  {
    printF("ERROR: un-implemented kernelType=%i\n",kernelType);
    Overture::abort("error");
  }
  

  // --- assign periodic images ---
  if( numberOfDimensions==2 )
  {
    for( int i=hu.getBase(0); i<0; i++ )
    {
      for( int j=0; j<numberOfFields; j++ )
	hu(i,j)=hu(i+numberOfGridPoints1,j);
    }
    for( int i=numberOfGridPoints1; i<=hu.getBound(0); i++ )
    {
      for( int j=0; j<numberOfFields; j++ )
	hu(i,j)=hu(i-numberOfGridPoints1,j);
    }
  }
  else
  {
    // --- assign periodic images ----
    // ** CHECK ME**

    Range M = numberOfFields;

    Range I3= Range(hu.getBase(1),hu.getBound(1));
    Range I2= Range(hu.getBase(0),-1);
    hu(I2,I3,M) = hu(I2+numberOfGridPoints1,I3,M);    // left-ghost = right
    I2 = Range(numberOfGridPoints1,hu.getBound(0));
    hu(I2,I3,M) = hu(I2-numberOfGridPoints1,I3,M);    // right-ghost = left 
    
    I2= Range(hu.getBase(0),hu.getBound(0));
    I3= Range(hu.getBase(1),-1);
    hu(I2,I3,M) = hu(I2,I3+numberOfGridPoints2,M);    // bottom-ghost = top 
    I3 = Range(numberOfGridPoints2,hu.getBound(1));
    hu(I2,I3,M) = hu(I2,I3-numberOfGridPoints2,M);    // top-ghost = bottom
    
    
  }
  
  cpuTime+=getCPU()-time;
  return 0;
}



#define evalRadKernel EXTERN_C_NAME(evalradkernel)

extern "C"
{

// -- optimized routine for computing H(u) and updating auxilary functions --
void evalRadKernel( const int & ipar, const LocalReal & rpar,
		    const LocalReal & amc, const Complex &  alpha, const Complex & beta,
		    const int&nd1a, const int&nd1b, const int&nd2a, const int&nd2b,
		    const int&nd3a, const int&nd3b, const int&nd4a, const int&nd4b, const int&nd5a, const int&nd5b,
		    Complex & fold2d, Complex & fold3d, Complex & phi2d, Complex & phi3d,
		    const int&ndzl1a, const int&ndzl1b, const int&ndzl2a, const int&ndzl2b,
		    Complex & zl2d, Complex & zl3d );
}


// ========================================================================================
//
/// \brief  Evaluate the radiation kernel. This version works in serial or parallel.
/// 
/// \param upLocal (input) :  input data on the boundary.
/// \param hupLocal (output) : the output hu(i) is defined for all i (with periodic images assigned too)
//  
//
/// \Note: See the notes in CG/DMX/radbc/radbc.pdf 
/// \Note: See paper
//      author="B. Alpert and L. Greengard and T. Hagstrom",
//      title="Nonreflecting Boundary Conditions for the Time-Dependent Wave Equation",
//      journal=JCP,
//      volume="180",
//      year="2002",
//      pages="270--296")
// 
// ========================================================================================
int RadiationKernel::evaluateKernelParallel( double dt, RealArray & upLocal, RealArray & hupLocal )
{

  
  // printF(">>RadiationKernel::evaluateKernelParallel **NEWEST**\n");
  

  LocalReal time=getCPU();
  
  const int myid = max(0,Communication_Manager::My_Process_Number);
  const int np=max(1,Communication_Manager::Number_Of_Processors);


  bool useComplexTransforms= true;   // true = use complex transforms where real transforms could be used instead
  // if( useComplexTransforms )
  //   printF(" RadiationKernel::evaluateKernelParallel: useComplexTransforms\n");
      
  const bool & useFourierTransformClass = dbase.get<bool>("useFourierTransformClass");
  assert( useFourierTransformClass );
  
  const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
  const int rside = dbase.get<int>("rside");
  const int raxis = dbase.get<int>("raxis");

  const int orderInTime=orderOfTimeStepping; 
  const int orderInTimeMinus1=orderInTime-1;
  const int numberOfTimeLevels = orderInTime-1; // we store this many time levels 

  int & currentTimeLevel = dbase.get<int>("currentTimeLevel");    // current time level for time integral by quadrature 
  currentTimeLevel = (currentTimeLevel + 1) % numberOfTimeLevels;

  const int & debug = dbase.get<int>("debug");
  FILE *debugFile=dbase.get<FILE*>("debugFile");
  if( debug & 1 )
  {
    fprintf(debugFile,"evaluateKernelParallel: numberOfDimensions=%d\n",numberOfDimensions);
  }
  

  int & bcInitComplex = dbase.get<int>("bcInitComplex");

  // assert( ploc!=NULL );
  

  
  assert( bcInitComplex>=0 );
  
  FourierTransform & fourierTransform = dbase.get<FourierTransform>("fourierTransform");

  // -- Here is the computational box ---
  Range I1,I2,I3;
  IndexBox & fbox = dbase.get<IndexBox>("fbox");
  I1=Range(fbox.base(0),fbox.bound(0));
  I2=Range(fbox.base(1),fbox.bound(1));
  I3=Range(fbox.base(2),fbox.bound(2));
	  
  if( np>1 && debug >0 )
  {
    fprintf(debugFile,"evaluateKernelParallel fbox=[%i,%i][%i,%i][%i,%i]\n",
	    fbox.base(0),fbox.bound(0),fbox.base(1),fbox.bound(1),fbox.base(2),fbox.bound(2));
    fflush(debugFile);
      
    // fclose(debugFile);
      
    // OV_ABORT("evaluateKernelParallel: stop here for now");
      
  }

  const int n1=numberOfGridPoints1;  // number of grid points
  const int n2=numberOfGridPoints2;  // number of grid points
  const int n1by2 = n1/2;
  const int n2by2 = n2/2;

  // --- Local array sizes for the FFT: 
  IndexBox & fftBox = dbase.get<IndexBox>("fftBox");
  int axisp1, axisp2;
  if( raxis==0 )
  {
    axisp1=1; axisp2=2;
  }
  else if( raxis==1 )
  {
    axisp1=0; axisp2=2;
  }
  else
  {
    axisp1=0; axisp2=1;
  }

  // Local index bounds on loops
  const int i1a=fftBox.base(axisp1), i1b=fftBox.bound(axisp1);  // 0:n1 -> i1a:i1b
  const int i2a=fftBox.base(axisp2), i2b=fftBox.bound(axisp2);  // 0:n2 -> i2a:i2b.

  const int n1f = fftBox.bound(axisp1)-fftBox.base(axisp1)+1;
  const int n2f = fftBox.bound(axisp2)-fftBox.base(axisp2)+1;


  bool isEmpty = fftBox.isEmpty();

  if( !isEmpty )
  {
    if( debug & 1 )
    {
      fprintf(debugFile,"evalKernel: fftBox =[%3i,%3i][%3i,%3i][%3i,%3i][%3i,%3i]\n",
	      fftBox.base(0),fftBox.bound(0),
	      fftBox.base(1),fftBox.bound(1),
	      fftBox.base(2),fftBox.bound(2),
	      fftBox.base(3),fftBox.bound(3));
   
      fprintf(debugFile,"n1=%d, n2=%d\n",n1,n2);
      fprintf(debugFile,"local array sizes: [i1a,i1b][i2a,i2b]=[%3d,%3d][%3d,%3d] \n",i1a,i1b,i2a,i2b);
      fprintf(debugFile,"local array sizes: n1f=%d, n2f=%d, isEmpty=%d \n",n1f,n2f,(int)isEmpty);
      fflush(debugFile);
    
    }

    if( upLocal.getBase(axisp1)>i1a || upLocal.getBound(axisp1)<i1b ||
	( numberOfDimensions==3 && upLocal.getBase(axisp2)>i2a || upLocal.getBound(axisp2)<i2b ) )
    {
      if( debug>0 )
      {
	fprintf(debugFile,"RadiationKernel::evaluateKernelParallel:ERROR: array upLocal is not of the expected size!\n");
	fprintf(debugFile," axisp1=%d: [i1a,i1b]=[%d,%d] upLocal=[%d,%d]\n",axisp1,i1a,i1b,
		upLocal.getBase(axisp1),upLocal.getBound(axisp1));
	fprintf(debugFile," axisp2=%d: [i2a,i2b]=[%d,%d] upLocal=[%d,%d]\n",axisp2,i2a,i2b,
		upLocal.getBase(axisp2),upLocal.getBound(axisp2));
	// ::display(upLocal,"upLocal");
	fclose(debugFile);
      }
    
      OV_ABORT("ERROR- upLocal not of the expected size -- see debug file for more info");
    }
  }
  else
  {
    if( debug & 1 )
    {
      fprintf(debugFile,"RadEval::evaluateKernelParallel: Index box is empty.\n");
    }
    
  }
  

  if( kernelType==planar )
  {
    // ------------- Planar Kernel : flat boundary -------------

    // The code below is based on the Fortran routine bcperq21.f
    // 
    // We solve for auxilary vaiables phi_j 
    //    (d/dt - c beta_j w) phi_j = c alpha_j w^2 phat 
    //    phat(t) = 
    //    fhat(t) = sum phi_j   :    FourierTranform of solution 

    Complex I(0.,1.);

    int md = 2*numberOfModes1+1;    // total number of Fourier modes to keep, max = n1/2
    if( numberOfDimensions==3 )
      md *= 2*numberOfModes2+1;

    int numComp=numberOfFields;  


    if( !dbase.has_key("phi2c") )
    {
      // ---- initialize work spaces and data ----

      // Here are the new versions with storage for "k" based on the FFT format 
      dbase.put<Complex*>("phi2c")=NULL;
      Complex *& phi2c = dbase.get<Complex*>("phi2c");
      if( numberOfDimensions==2 )
        phi2c = new Complex [n1f*numberOfPoles*numComp];            
      else 
        phi2c = new Complex [(n1f*n2f)*numberOfPoles*numberOfFields];            

      dbase.put<Complex*>("fold2c")=NULL;
      Complex *&fold2c = dbase.get<Complex*>("fold2c");
      if( numberOfDimensions==2 )
        fold2c = new Complex[orderInTimeMinus1*(n1f)*numberOfPoles*numComp];
      else
	fold2c = new Complex[orderInTimeMinus1*(n1f*n2f)*numberOfPoles*numberOfFields];


      // array to hold Adams-Moulton coefficients 
      dbase.put<RealArray>("amc");
      RealArray & amc = dbase.get<RealArray>("amc");
      amc.redim(Range(-1,orderInTime-2));


      // --- initialize alpha and beta arrays defining the poles in the non-local kernel approximation
      initializePoles();

    }

    Complex *& alphac = dbase.get<Complex*>("alphac");
    Complex *& betac = dbase.get<Complex*>("betac");
    #define alpha(i) alphac[(i-1)]
    #define beta(i)  betac[(i-1)]


    const LocalReal scl1=twoPi/period1;
    const LocalReal scl2= numberOfDimensions==3 ? twoPi/period2 : 1.0;

    const bool useOptEval=true; // false; // true; // true   // *******

    // ipar and rpar are for evalRadKernel 
    int ipar[] = {
      numberOfDimensions,
      numberOfPoles,
      orderInTime,
      n1,n2,
      numberOfModes1,numberOfModes2,
      i1a,i1b,i2a,i2b,
      0, // mc : filled in below 
      currentTimeLevel,
      numberOfTimeLevels
    }; // 
    LocalReal rpar[]={ scl1,scl2,c,dt};  // 


    if( numberOfDimensions==2 )  
    {
      // ================ PLANAR TWO DIMENSIONS ======================

      // Here are the new versions with storage for "k" based on the FFT format 
      Complex *& phi2c = dbase.get<Complex*>("phi2c");
      #define phi2(kfft,pole,mc) phi2c[((kfft)-i1a)+ n1f*( (pole) + numberOfPoles*(mc))]
       
      // ComplexArray fold(0:orderInTime-1,md,numberOfPoles,numComp) : 
      Complex *&fold2c = dbase.get<Complex*>("fold2c");
      #define fold2(l,kfft,pole,mc) fold2c[(l)+orderInTimeMinus1*( ((kfft)-i1a)+ n1f*( (pole) + numberOfPoles*(mc)))]

      RealArray & amc = dbase.get<RealArray>("amc");
      
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int kfft=i1a; kfft<=i1b; kfft++ )  // loop over modes on this processor 
	    {
	      phi2(kfft,pole,mc)=0.;
	      for( int l=0; l<=orderInTime-2; l++ )
	      {
		fold2(l,kfft,pole,mc)=0.;
	      }
	    }
	  }
	}
	 
	 
	bcInitComplex=1;
	 
      } // end initialize



      Complex *pzl = new Complex[n1f];
      #define zl(i1) pzl[i1-i1a]
      
      LocalReal scl=twoPi/period1;

      // Array dimensions for evalOpt: 
      int nd1a=0, nd1b=orderInTimeMinus1-1;
      int nd2a=i1a, nd2b=i1b;
      int nd3a=0, nd3b=numberOfPoles-1;
      int nd4a=0, nd4b=numComp-1;
      int nd5a=0, nd5b=0;
      int ndzl1a=i1a, ndzl1b=i1b;
      int ndzl2a=i2a, ndzl2b=i2b;

      // ------- loop over components -------------
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {

	fourierTransform.forwardTransform( upLocal,mc, pzl  );

	if( !isEmpty )
	{
	  if( debug & 1 )
	  {
	    fprintf(debugFile,"After fourierTransform.forwardTransform, mc=%d  **NEWER** \n",mc);
	    for( int i=0; i<n1; i++ ){ fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i)),std::imag(zl(i))); }   // 
	    fprintf(debugFile,"]\n");
	  }

	  if( useOptEval )
	  {
            // ---- optimized evaluation of H(u) and auxilary functions ----
	    ipar[11]=mc;
            evalRadKernel( ipar[0],rpar[0],amc(-1),alphac[0],betac[0],
			   nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,nd5a,nd5b,
			   fold2c[0], fold2c[0], phi2c[0], phi2c[0],
			   ndzl1a,ndzl1b,ndzl2a,ndzl2b, pzl[0], pzl[0] );
	    
	  }
	  else if( true )
          {
            // --- non-optimized evaluate of H(u) -----
            
            // --- loop over kfft ---
            for( int kfft=i1a; kfft<=i1b; kfft++ )
            {
              const int k = kfft<n1by2 ? kfft : kfft-n1;  // true Fourier mode k
              if( abs(k) <= numberOfModes1 )
              {
                // LocalReal w = scl*k;
                LocalReal w = scl*abs(k); // *wdh* Oct 30, 2020 : sqrt( k^2 )
                Complex phat;
                phat = zl(kfft);   // use kfft here 

                for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
                {

                  Complex expBeta = exp(c*beta(pole+1)*dt*w);  // beta's are complex 
                  Complex xfact = expBeta;
                  Complex adon = xfact*phi2(kfft,pole,mc);         // add-on 
                  // printF("* adon=[%12.4e,%12.4e]\n",std::real(adon),std::imag(adon));
                  for( int l=0; l<=orderInTime-2; l++ )
                  {
                    adon += xfact*amc(l)*fold2(l,kfft,pole,mc);
                    xfact *= expBeta;
                    // printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
                  }


                  // phat = zlocr(kfft) + zloci(kfft) *I;   // use kfft here 
                  // if( k<0 )
                  //   phat=std::conj(phat);
		
                  adon += c*dt*amc(-1)*alpha(pole+1)*w*w*phat;
		
                  phi2(kfft,pole,mc)=adon;
		  
                  // printF("expBeta=[%12.4e,%12.4e]\n",std::real(expBeta),std::imag(expBeta));
                  // printF("kfft=%d, phi=[%18.10e,%18.10e]\n",kfft,std::real(phi2(kfft,pole,mc)),std::imag(phi2(kfft,pole,mc)));
                  // printF("phat=[%18.10e,%18.10e] alpha=[%18.10e,%18.10e] beta=[%18.10e,%18.10e] \n",
                  //        std::real(phat),std::imag(phat),
                  //        std::real(alpha(pole+1)),std::imag(alpha(pole+1)),
                  //        std::real(beta(pole+1)),std::imag(beta(pole+1)));
	    
                  for( int l=orderInTime-2; l>=1; l-- )
                  {
                    fold2(l,kfft,pole,mc)=fold2(l-1,kfft,pole,mc);
                  }
                  fold2(0,kfft,pole,mc)=c*dt*w*w*alpha(pole+1)*phat;   // alpha's are complex 
		
                  // printF("radKer: k=%2d, phat=[%12.5e,%12.5e] phi=[%12.5e,%12.5e] fold(0,%d,%d,%d)=(%12.5e,%12.5e)\n",k,
                  //        std::real(phat),std::imag(phat),
                  //        std::real(phi2(kfft,pole,mc)),std::imag(phi2(kfft,pole,mc)),
                  //        kfft,pole,mc,std::real(fold2(0,kfft,pole,mc)),
                  //        std::imag(fold2(0,kfft,pole,mc)));

                
                }  // end for pole
              }
              
            }  // end for k

            // --- loop over kfft ---
            for( int kfft=i1a; kfft<=i1b; kfft++ )
            {
              zl(kfft)=0.;
              int k = kfft<n1by2 ? kfft : kfft-n1;  // true k
              if( abs(k)<=numberOfModes1 )
              {
                for( int pole=0; pole<numberOfPoles; pole++ ) 
                {
                  zl(kfft) += phi2(kfft,pole,mc);
                }
              }
            }

          }
          else 
	  {
            // --- loop over kfft ---
	    for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	    {
	      for( int kfft=i1a; kfft<=i1b; kfft++ )
	      {
		const int k = kfft<n1by2 ? kfft : kfft-n1;  // true Fourier mode k
		if( abs(k) > numberOfModes1 )
		  continue;   // do not compute this mode
	      
		// LocalReal w = scl*k;
		LocalReal w = scl*abs(k); // *wdh* Oct 30, 2020 : sqrt( k^2 )

		Complex expBeta = exp(c*beta(pole+1)*dt*w);  // beta's are complex 
		Complex xfact = expBeta;
		Complex adon = xfact*phi2(kfft,pole,mc);         // add-on 
		// printF("* adon=[%12.4e,%12.4e]\n",std::real(adon),std::imag(adon));
		for( int l=0; l<=orderInTime-2; l++ )
		{
		  adon += xfact*amc(l)*fold2(l,kfft,pole,mc);
		  xfact *= expBeta;
		  // printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
		}


		Complex phat;
		phat = zl(kfft);   // use kfft here 
		// phat = zlocr(kfft) + zloci(kfft) *I;   // use kfft here 
		// if( k<0 )
		//   phat=std::conj(phat);
		
		adon += c*dt*amc(-1)*alpha(pole+1)*w*w*phat;
		
		phi2(kfft,pole,mc)=adon;
		  
		// printF("expBeta=[%12.4e,%12.4e]\n",std::real(expBeta),std::imag(expBeta));
		// printF("kfft=%d, phi=[%18.10e,%18.10e]\n",kfft,std::real(phi2(kfft,pole,mc)),std::imag(phi2(kfft,pole,mc)));
		// printF("phat=[%18.10e,%18.10e] alpha=[%18.10e,%18.10e] beta=[%18.10e,%18.10e] \n",
                //        std::real(phat),std::imag(phat),
                //        std::real(alpha(pole+1)),std::imag(alpha(pole+1)),
                //        std::real(beta(pole+1)),std::imag(beta(pole+1)));
	    
		for( int l=orderInTime-2; l>=1; l-- )
		{
		  fold2(l,kfft,pole,mc)=fold2(l-1,kfft,pole,mc);
		}
		fold2(0,kfft,pole,mc)=c*dt*w*w*alpha(pole+1)*phat;   // alpha's are complex 
		
		// printF("radKer: k=%2d, phat=[%12.5e,%12.5e] phi=[%12.5e,%12.5e] fold(0,%d,%d,%d)=(%12.5e,%12.5e)\n",k,
                //        std::real(phat),std::imag(phat),
                //        std::real(phi2(kfft,pole,mc)),std::imag(phi2(kfft,pole,mc)),
		//        kfft,pole,mc,std::real(fold2(0,kfft,pole,mc)),
		//        std::imag(fold2(0,kfft,pole,mc)));

	      }  // end for k
	    }  // end for pole

	 
	  
            // --- loop over kfft ---
            for( int kfft=i1a; kfft<=i1b; kfft++ )
            {
              zl(kfft)=0.;
              int k = kfft<n1by2 ? kfft : kfft-n1;  // true k
              if( abs(k)<=numberOfModes1 )
              {
                for( int pole=0; pole<numberOfPoles; pole++ ) 
                {
                  zl(kfft) += phi2(kfft,pole,mc);
                }
              }
            }

	  } // end if !useOptEval 
	  
	}
	// end if !isEmpty

	fourierTransform.backwardTransform( pzl, hupLocal,mc  );

	if( debug & 1 )
	{
	  ::display(hupLocal,sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWEST**",mc),
		    debugFile,"%13.6e ");
	  

	  // ::display(hupLocal(I1,I2,I3),sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),
	  // 	    debugFile,"%13.6e ");
	}
	  
      }  // end for component mc
        
     
      // delete [] pzloc;
      delete [] pzl;

    }
    else if( numberOfDimensions==3 )
    {
      // ---------- PLANAR THREE-DIMENSIONS PARALLEL ------------

      if( debug & 1 )
      {
	printF("$$$$$$$$$$$$ RADIATION KERNEL -- NEWEST 3D PARALLEL VERSION bcInitComplex=%d $$$$$$$$$$$$$\n",bcInitComplex);
	fprintf(debugFile," n1=%d, numberOfModes1=%d, n2=%d, numberOfModes2=%d\n",n1,numberOfModes1,n2,numberOfModes2);
      }
      
      // #ifdef USE_PPP
      //   printF(" ----> RUNNING IN PARALLEL <-----\n");
      // #endif

      
      Complex *& phi2c = dbase.get<Complex*>("phi2c");
      #define phi2(k1fft,k2fft,pole,mc) phi2c[(k1fft-i1a) + n1f*( (k2fft-i2a) + n2f*( (pole) + numberOfPoles*(mc)))]
       
      Complex *&fold2c = dbase.get<Complex*>("fold2c");
      #define fold2(l,k1fft,k2fft,pole,mc) fold2c[(l)+orderInTimeMinus1*( (k1fft-i1a)+n1f*( (k2fft-i2a) + n2f*( (pole) + numberOfPoles*(mc))))]

       
      RealArray & amc = dbase.get<RealArray>("amc");
       
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int k2fft=i2a; k2fft<=i2b; k2fft++ )
	    {
	      for( int k1fft=i1a; k1fft<=i1b; k1fft++ )
	      {
		phi2(k1fft,k2fft,pole,mc)=0.;
		for( int l=0; l<=orderInTime-2; l++ )
		{
		  fold2(l,k1fft,k2fft,pole,mc)=0.;
		}
	      }
	    }
	  }
	}

	bcInitComplex=1;
	 
      } // end initialize

      
      // Array dimensions for evalOpt: 
      int nd1a=0, nd1b=orderInTimeMinus1-1;
      int nd2a=i1a, nd2b=i1b;
      int nd3a=i2a, nd3b=i2b;
      int nd4a=0, nd4b=numberOfPoles-1;
      int nd5a=0, nd5b=numComp-1;
      int ndzl1a=i1a, ndzl1b=i1b;
      int ndzl2a=i2a, ndzl2b=i2b;

      // complex*16 zl(n1,n2),zl2(n2)        
      Complex *pzl = new Complex [n1f*n2f];  
      #define zl(i1,i2) pzl[(i1-i1a)+n1f*(i2-i2a)]

      LocalReal scl1=twoPi/period1;
      LocalReal scl2=twoPi/period2;

      // printF("scl=%20.14e\n",scl);
      // Range I1=n1, I2=n2;
       
      // loop over components: 
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {

	fourierTransform.forwardTransform( upLocal,mc, pzl  );
	  
	if( !isEmpty )
	{
	  if( debug & 1 )
	  {
	    fprintf(debugFile,"After fourierTransform.forwardTransform, mc=%d  **NEWER** \n",mc);
	    for( int i2=i2a; i2<=i2b; i2++ )
	    {
	      fprintf(debugFile,"i2=%d: ",i2);
	      for( int i1=i1a; i1<=i1b; i1++ )
	      {
		fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	      }
	      fprintf(debugFile,"]\n");
	    }
	  }

	  if( useOptEval )
	  {
            // ---- optimized evaluation of H(u) and auxilary functions ----
	    ipar[11]=mc;
            evalRadKernel( ipar[0],rpar[0],*amc.getDataPointer(),alphac[0],betac[0],
			   nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,nd5a,nd5b,
			   fold2c[0], fold2c[0], phi2c[0], phi2c[0],
			   ndzl1a,ndzl1b,ndzl2a,ndzl2b, pzl[0], pzl[0] );
	    
	  }
	  else 
	  {
            // ---- Evaluation of H(u) and auxilary functions ----

            // --- loop over k1fft and k2fft ---
	    for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	    {
	      for( int k2fft=i2a; k2fft<=i2b; k2fft++ )
	      {
		const int k2 = k2fft<n2by2 ? k2fft : k2fft-n2;  // true Fourier mode k2 
		if( abs(k2)>numberOfModes2 )
		  continue;   // no need to compute this mode

		for( int k1fft=i1a; k1fft<=i1b; k1fft++ )
		{
		  const int k1 = k1fft<n1by2 ? k1fft : k1fft-n1;  // true Fourier mode k1
		  if( abs(k1)>numberOfModes1 )
		    continue;    // no need to compute this mode

		  LocalReal w= sqrt( SQR(scl1*k1) + SQR(scl2*k2) );

		  Complex expBeta = exp(c*beta(pole+1)*dt*w);              // beta's are complex 
		  Complex xfact = expBeta;
		  Complex adon = xfact*phi2(k1fft,k2fft,pole,mc);                 // add-on
		  for( int l=0; l<=orderInTime-2; l++ )
		  {
		    adon += xfact*amc(l)*fold2(l,k1fft,k2fft,pole,mc);
		    xfact *= expBeta;
		    // printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
		  }
		  adon += c*dt*amc(-1)*alpha(pole+1)*w*w*zl(k1fft,k2fft);  
		  phi2(k1fft,k2fft,pole,mc)=adon;

		  for( int l=orderInTime-2; l>=1; l-- )
		  {
		    fold2(l,k1fft,k2fft,pole,mc)=fold2(l-1,k1fft,k2fft,pole,mc);
		  }
		  fold2(0,k1fft,k2fft,pole,mc)=c*dt*w*w*alpha(pole+1)*zl(k1fft,k2fft);   // alpha's are complex 
		}
	      } // end for k2 
	    } // end for pole 


            // --- loop over k1fft and k2fft ---
	    for( int k1fft=i1a; k1fft<=i1b; k1fft++ )
	    {
	      for( int k2fft=i2a; k2fft<=i2b; k2fft++ )
	      {
		zl(k1fft,k2fft)=0.;
		int k1 = k1fft<n1by2 ? k1fft : k1fft-n1;  // true k1
		int k2 = k2fft<n2by2 ? k2fft : k2fft-n2;  // true k2 
		if( abs(k1)<=numberOfModes1 && abs(k2)<=numberOfModes2 )
		{
		  for( int pole=0; pole<numberOfPoles; pole++ ) 
		  {
		    zl(k1fft,k2fft) += phi2(k1fft,k2fft,pole,mc);
		  }
		}
	      }
	    }
	  }
          

	  if( debug & 1 )
	  {
	    fprintf(debugFile,"Before backwardTransform, mc=%d  zl \n",mc);
	    for( int i2=i2a; i2<=i2b; i2++ )
	    {
	      fprintf(debugFile,"i2=%d: ",i2);
	      for( int i1=i1a; i1<=i1b; i1++ )
	      {
		fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	      }
	      fprintf(debugFile,"]\n");
	    }
	  }


	  
	} // end !isEmpty

	fourierTransform.backwardTransform( pzl, hupLocal,mc  );
	  
	if( debug  & 1 && !isEmpty )
	{
          // fprintf(debugFile,"After fourierTransform.backwardTransform, mc=%d  **NEWER**\n",mc);
	  // ::display(hupLocal,sPrintF("After fourierTransform.backwardTransform, mc=%d  hpuLocal:",mc),debugFile,"%13.6e ");
	  
          Index Iv[4], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2], &I4=Iv[3];

          I1=hupLocal.dimension(0); I2=hupLocal.dimension(1); I3=hupLocal.dimension(2); I4=hupLocal.dimension(3);
          Iv[raxis]=Range(fbox.base(raxis),fbox.base(raxis));
	  
	  if( true )
	  {
	    RealArray temp;
	    temp.redim(I1,I2,I3,I4);
	    temp=hupLocal(I1,I2,I3,I4);
	    if( raxis==0 )
	      temp.reshape(I2,I3,I4);
	    else if( raxis==1 )
	      temp.reshape(I1,I3,I4);
	    else
	      temp.reshape(I1,I2,I4);
	  
	    ::display(temp,sPrintF("After fourierTransform.backwardTransform, before periodic update, mc=%d  hpuLocal:",mc),debugFile,"%13.6e ");

	  }
	  
	  // // ** FIX ME FOR PARALLEL 

	  // for( int i3=i3a; i3<=i3b, i3++ )
	  // {
	  //   if( raxis!=2 ) fprintf("i3=%d: ");
	  //   for( int i2=i2a; i2<=i2b, i2++ )
	  //   {
	  //     if( raxis!=1 ) fprintf(debugFile,"i2=%d: ");
	  //     for( int i1=i1a; i1<=i1b, i1++ )
	  //     {
	  // 	fprintf(debugFile,"%13.6e ",hupLocal(i1,i2,i3));
	  //     }
	  //     if( raxis!=0 ) fprintf(debugFile,"\n");
	  //   }
	  //   if( raxis!=1 ) fprintf(debugFile,"\n");
	  // }
	  // if( raxis!=2 ) fprintf(debugFile,"\n");
	
	  // // ** FIX ME FOR PARALLEL 
	  // RealArray temp;
	  // temp.redim(I1,I2,I3);
	  // temp=hupLocal(I1,I2,I3);
	  // if( raxis==0 )
	  //   temp.reshape(I2,I3);
	  // else if( raxis==1 )
	  //   temp.reshape(I1,I3);
	  // else
	  //   temp.reshape(I1,I2);
	    
	  // ::display(temp,sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),debugFile,"%13.6e ");
	  //  ::display(hupLocal(I1,I2,I3),sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),"%13.6e ");
	}

      }  // end for mc 
	
      delete [] pzl;


    }
    else
    {
      OV_ABORT("ERROR: unexpected numberOfDimensions");
    }
    
    


  }
  else if( kernelType==cylindrical )
  {
    OV_ABORT("RadiationKernel:parallel: finish me for kernelType==cylindrical");
    
    // bccyld(nda1,ndb1,
    //        *pu, *pf, *ploc, c,radius, dt,numberOfGridPoints1,numberOfFields,numberOfModes1,ns1,
    //        orderOfTimeStepping, *fold,*phi,
    // 	   *npoles, *alpha, *beta,*amc,*fftsave1,bcinit);
  }
  else
  {
    printF("ERROR: un-implemented kernelType=%i\n",kernelType);
    Overture::abort("error");
  }

  // --- assign periodic images ---
  Range N=numberOfFields;
  fourierTransform.periodicUpdate( hupLocal,N ); 
  
  cpuTime+=getCPU()-time;
  return 0;
}

// ========================================================================================
//
/// \brief  Evaluate the radiation kernel. Version 3 *New* parallel version
/// \param u (input) :  
/// \param hu (output) : the output hu(i) is defined for all i (with periodic images assigned too)
//  
//
//  This routine assumes that the input values run from u(0),...,u(numberOfGridPoints-1)
//
/// \Note: See paper
//      author="B. Alpert and L. Greengard and T. Hagstrom",
//      title="Nonreflecting Boundary Conditions for the Time-Dependent Wave Equation",
//      journal=JCP,
//      volume="180",
//      year="2002",
//      pages="270--296")
// 
// ========================================================================================
int RadiationKernel::evaluateKernelParallelV3( double dt, RealArray & upLocal, RealArray & hupLocal )
{

  printF(">>RadiationKernel::evaluateKernelParallel **NEWEST**\n");
  

  LocalReal time=getCPU();
  
  const int myid = max(0,Communication_Manager::My_Process_Number);
  const int np=max(1,Communication_Manager::Number_Of_Processors);


  bool useComplexTransforms= true;   // true = use complex transforms where real transforms could be used instead
  // if( useComplexTransforms )
  //   printF(" RadiationKernel::evaluateKernelParallel: useComplexTransforms\n");
      
  const bool & useFourierTransformClass = dbase.get<bool>("useFourierTransformClass");
  assert( useFourierTransformClass );
  
  const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
  const int rside = dbase.get<int>("rside");
  const int raxis = dbase.get<int>("raxis");

  const int & debug = dbase.get<int>("debug");
  FILE *debugFile=dbase.get<FILE*>("debugFile");
  if( debug & 1 )
  {
    fprintf(debugFile,"evaluateKernelParallel: numberOfDimensions=%d\n",numberOfDimensions);
  }
  

  int & bcInitComplex = dbase.get<int>("bcInitComplex");

  assert( ploc!=NULL );
  
  // double *pu = u.getDataPointer();
  // double *pf = hu.getDataPointer();
  
  // // We assume normal direction is axis=0 for now:
  // int axis=0;  // normal axis 
  // int axisp1 = 1;
  // int axisp2 = 2;
  
  // // make sure u and hu are the same size: 
  // assert( u.getBase(0)==hu.getBase(0) && u.getBound(0)==hu.getBound(0) );
  // assert( u.getBase(axisp1)==hu.getBase(axisp1) && u.getBound(axisp1)==hu.getBound(axisp1) );
  // if( numberOfDimensions==3 )
  // {
  //   assert( u.getBase(axisp2)==hu.getBase(axisp2) && u.getBound(axisp2)==hu.getBound(axisp2) );
  // }
  

  // // Leading dimensions of the u and hu arrays: 
  // const int nda1=u.getBase(0)+1, ndb1=u.getBound(0)+1;  // add one (base one assumed in bcperq21d)
  // const int nda2=u.getBase(1)+1, ndb2=u.getBound(1)+1;  // add one (base one assumed in bcperq21d)
  
  assert( bcInitComplex>=0 );
  
  FourierTransform & fourierTransform = dbase.get<FourierTransform>("fourierTransform");

  // -- Here is the computational box ---
  Range I1,I2,I3;
  IndexBox & fbox = dbase.get<IndexBox>("fbox");
  I1=Range(fbox.base(0),fbox.bound(0));
  I2=Range(fbox.base(1),fbox.bound(1));
  I3=Range(fbox.base(2),fbox.bound(2));
	  
  if( np>1 && debug >0 )
  {
    fprintf(debugFile,"evaluateKernelParallel fbox=[%i,%i][%i,%i][%i,%i]\n",
	    fbox.base(0),fbox.bound(0),fbox.base(1),fbox.bound(1),fbox.base(2),fbox.bound(2));
    fflush(debugFile);
      
    // fclose(debugFile);
      
    // OV_ABORT("evaluateKernelParallel: stop here for now");
      
  }

  const int n1=numberOfGridPoints1;  // number of grid points
  const int n2=numberOfGridPoints2;  // number of grid points

  IndexBox & fftBox = dbase.get<IndexBox>("fftBox");
  // --- Local array sizes for the FFT: 
  // int n1f = fftBox.bound(0)-fftBox.base(0)+1;
  // int n2f = fftBox.bound(1)-fftBox.base(1)+1;
  // int n3f = fftBox.bound(2)-fftBox.base(2)+1;    
  int axisp1, axisp2;
  if( raxis==0 )
  {
    axisp1=1; axisp2=2;
  }
  else if( raxis==1 )
  {
    axisp1=0; axisp2=2;
  }
  else
  {
    axisp1=0; axisp2=1;
  }

  // Local index bounds on loops
  const int i1a=fftBox.base(axisp1), i1b=fftBox.bound(axisp1);  // 0:n1 -> i1a:i1b
  const int i2a=fftBox.base(axisp2), i2b=fftBox.bound(axisp2);  // 0:n2 -> i2a:i2b.

  const int n1f = fftBox.bound(axisp1)-fftBox.base(axisp1)+1;
  const int n2f = fftBox.bound(axisp2)-fftBox.base(axisp2)+1;

  // if( np==1 )
  // {
  //   assert( n1==n1f && n2==n2f );
  // }
  

  // const int n1by2 = min(n1/2,i1b);   // 0:n1/2 --> i1a:n1by2 
  // const int n2by2 = min(n2/2,i2b);   // 0:n1/2 --> i1a:n2by2 

  // const int numModes1 = min(numberOfModes1,i1b+1);    // 0..numberOfModes1 -> i1a:numModes1
  // const int numModes2 = min(numberOfModes2,i2b+1);    // 0..numberOfModes2 -> i1a:numModes2
  
  

  bool isEmpty = fftBox.isEmpty();

  if( debug & 1 )
  {
    fprintf(debugFile,"evalKernel: fftBox =[%3i,%3i][%3i,%3i][%3i,%3i][%3i,%3i]\n",
	    fftBox.base(0),fftBox.bound(0),
	    fftBox.base(1),fftBox.bound(1),
	    fftBox.base(2),fftBox.bound(2),
	    fftBox.base(3),fftBox.bound(3));
   
    fprintf(debugFile,"n1=%d, n2=%d\n",n1,n2);
    fprintf(debugFile,"local array sizes: [i1a,i1b][i2a,i2b]=[%3d,%3d][%3d,%3d] \n",i1a,i1b,i2a,i2b);
    fprintf(debugFile,"local array sizes: n1f=%d, n2f=%d, isEmpty=%d \n",n1f,n2f,(int)isEmpty);
    fflush(debugFile);
    
  }


  if( kernelType==planar )
  {
    // ------------- Planar Kernel : flat boundary -------------

    // **NEW WAY**
    // Start from code in bcperq21.f
    // We solve for auxilary vaiables phi_j 
    //    (d/dt - c beta_j w) phi_j = c alpha_j w^2 phat 
    //    phat(t) = 
    //    fha(t)t = sum phi_j   :    FourierTranform of solution 

    Complex I(0.,1.);

    const int orderInTime=orderOfTimeStepping; 
    const int orderInTimeMinus1=orderInTime-1;

    int md = 2*numberOfModes1+1;    // total number of Fourier modes to keep, max = n1/2
    if( numberOfDimensions==3 )
      md *= 2*numberOfModes2+1;

    int numComp=numberOfFields;  


    if( !dbase.has_key("phic") )
    {
      // ---- initialize work spaces and data ----

      // *NOTE* for now allocate space for all modes in parallel, even though not needed

      dbase.put<Complex*>("phic")=NULL;
      Complex *& phic = dbase.get<Complex*>("phic");
      // allocate space for plus and minus modes   -numberOfModes1:numberOfModes1
      if( numberOfDimensions==2 )
        phic = new Complex [2*(numberOfModes1+1)*numberOfPoles*numComp];            
      else 
        phic = new Complex [(2*numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields];            

      dbase.put<Complex*>("foldc")=NULL;
      Complex *&foldc = dbase.get<Complex*>("foldc");
      if( numberOfDimensions==2 )
        foldc = new Complex[orderInTimeMinus1*2*(numberOfModes1+1)*numberOfPoles*numComp];
      else
	foldc = new Complex[orderInTimeMinus1*(2*numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields];


      // array to hold Adams-Moulton coefficients 
      dbase.put<RealArray>("amc");
      RealArray & amc = dbase.get<RealArray>("amc");
      amc.redim(Range(-1,orderInTime-2));


      // --- initialize alpha and beta arrays defining the poles in the non-local kernel approximation
      initializePoles();

    }

    Complex *& alphac = dbase.get<Complex*>("alphac");
    Complex *& betac = dbase.get<Complex*>("betac");
    #define alpha(i) alphac[(i-1)]
    #define beta(i)  betac[(i-1)]


    


    if( numberOfDimensions==2 )  
    {
      // ================ PLANAR TWO DIMENSIONS ======================


      Complex *& phic = dbase.get<Complex*>("phic");
      #define phi(k,pole,mc) phic[(k+numberOfModes1)+md*( (pole) + numberOfPoles*(mc))]
       
      // ComplexArray fold(0:orderInTime-1,md,numberOfPoles,numComp) : 
      Complex *&foldc = dbase.get<Complex*>("foldc");
      #define fold(l,k,pole,mc) foldc[(l)+orderInTimeMinus1*( (k+numberOfModes1)+md*( (pole) + numberOfPoles*(mc)))]

      RealArray & amc = dbase.get<RealArray>("amc");
      
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    // for( int k=0; k<numberOfModes1; k++ )
	    for( int k=-numberOfModes1; k<=numberOfModes1; k++ )
	    {
	      phi(k,pole,mc)=0.;
	      for( int l=0; l<=orderInTime-2; l++ )
	      {
		fold(l,k,pole,mc)=0.;
	      }
	    }
	  }
	}
	 
	bcInitComplex=1;
	 
      } // end initialize

      Complex *pzl = new Complex[n1f];
      #define zl(i1) pzl[i1-i1a]
      
//       RealArray ploc(n1); // , p(n1,numComp);
       
//       // For complex FFT's
//       //  zlocr(i) = real-part
//       //  zloci(i) = imag-part
//       const int n1c=2*n1;
//       LocalReal *pzloc = new LocalReal [n1*2];  // holds complex transforms 
// #define zlocr(i) pzloc[2*(i-i1a)]
// #define zloci(i) pzloc[1+2*(i-i1a)]


      LocalReal scl=twoPi/period1;
      // printF("scl=%20.14e\n",scl);
      // I1=n1;
       

      // loop over components: 
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {

	fourierTransform.forwardTransform( upLocal,mc, pzl  );

	if( !isEmpty )
	{
	  if( debug & 1 )
	  {
	    fprintf(debugFile,"After fourierTransform.forwardTransform, mc=%d  **NEWER** \n",mc);
	    for( int i=0; i<n1; i++ ){ fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i)),std::imag(zl(i))); }   // 
	    fprintf(debugFile,"]\n");
	  }

	  for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	  {
	    // for( int k=0; k<numberOfModes1; k++ )   // "k=1" loop over actual Fourier modes k that we compute 
	    for( int k=-numberOfModes1; k<=numberOfModes1; k++ ) // actual Fourier modes k 
	    {
	      // const int kp = k+1;
              const int kfft = k>=0 ? k : n1+k;  // location of mode k in fft 
	      assert( kfft>=0 && kfft<n1 );
	      if( kfft<i1a || kfft>i1b )
		continue;   // no data on this processor
	      
	      LocalReal w = scl*k;
	      Complex expBeta = exp(c*beta(pole+1)*dt*w);  // beta's are complex 
	      Complex xfact = expBeta;
	      Complex adon = xfact*phi(k,pole,mc);         // add-on 
	      // printF("* adon=[%12.4e,%12.4e]\n",std::real(adon),std::imag(adon));
	      for( int l=0; l<=orderInTime-2; l++ )
	      {
		adon += xfact*amc(l)*fold(l,k,pole,mc);
		xfact *= expBeta;
		// printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
	      }


	      Complex phat;
	      phat = zl(kfft);   // use kfft here 
	      // phat = zlocr(kfft) + zloci(kfft) *I;   // use kfft here 

	      adon += c*dt*amc(-1)*alpha(pole+1)*w*w*phat;
	      phi(k,pole,mc)=adon;

	      // printF("expBeta=[%12.4e,%12.4e]\n",std::real(expBeta),std::imag(expBeta));
	      //printF("phi(%d,%d,%d)=[%12.4e,%12.4e]\n",k,pole,mc,std::real(phi(k,pole,mc)),std::imag(phi(k,pole,mc)));
	    
	      for( int l=orderInTime-2; l>=1; l-- )
	      {
		fold(l,k,pole,mc)=fold(l-1,k,pole,mc);
	      }
	      fold(0,k,pole,mc)=c*dt*w*w*alpha(pole+1)*phat;   // alpha's are complex 
	    }  // end for k
	  }  // end for pole


	  // real FFT:      cos(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  
	  // complex FFT:   cos(0*x) sin(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  

	  // convert complex FFT to real:  (for n1 even)
	  //   pHat( k) = pzloc[k] ,    k=0,1,2,...,n1/2-1
	  //   pHat(-k) = pzloc[n-k]    k=  1,2,...,n1/2
        
	  for( int kfft=i1a; kfft<=i1b; kfft++ )
	  {
	    zl(kfft)=0.;
	  }
	  
	  // ploc=0.;  // fill in hHat 
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int k=-numberOfModes1; k<=numberOfModes1; k++ )   // loop over Fourier modes we compute 
	    {
              const int kfft = k>=0 ? k : n1+k;  // location of mode k in fft 
	      assert( kfft>=0 && kfft<n1 );
	      if( kfft>=i1a &&  kfft<=i1b )
	      {
		zl(kfft) += phi(k,pole,mc);
	      }
	      
	      // const int kp = k+1;
	      // ploc(2*kp-1) += std::real(phi(k,pole,mc));    
	      // ploc(2*kp  ) += std::imag(phi(k,pole,mc)); 

	    }
	  }
	

	  // // do this for now: *fix me*
	  // // for( int i=0; i<n1; i++ ) 
	  // for( int i=i1a; i<=i1b; i++ ) 
	  // {
	  //   zlocr(i)=0.;
	  //   zloci(i)=0.;
	  // }

	  // zlocr(0) = ploc(0);
	  // zloci(0) = 0.;
	  // int j=1;
	  // // for( int i=1; i<n1/2; i++ )
	  // for( int i=max(1,i1a); i<n1by2; i++ )
	  // {
	  //   // real FFT:      cos(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  
	  //   // complex FFT:   cos(0*x) sin(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  


	  //   // convert complex FFT to real:
	  //   //   pHat( k) = pzloc[k] ,    k=0,1,2,...,n1/2-1
	  //   //   pHat(-k) = pzloc[n-k]    k=0,1,2,...,n1/2

	  //   // k=1,2,...,n1/2-1
	  //   zlocr(i) = ploc(j  );
	  //   zloci(i) = ploc(j+1);

	  //   // negative k are the complex conjugate:
          //   // *** FIX ME FOR PARALLEL ***
	  //   zlocr(n1-i) =  ploc(j  );
	  //   zloci(n1-i) = -ploc(j+1);

	  //   j+=2;
	  
	  // }
	}
	// end if !isEmpty

	fourierTransform.backwardTransform( pzl, hupLocal,mc  );

	if( debug & 1 )
	{
	  ::display(hupLocal,sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWEST**",mc),
		    debugFile,"%13.6e ");
	  // ::display(hupLocal(I1,I2,I3),sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),
	  // 	    debugFile,"%13.6e ");
	}
	  
      }  // end for component mc
        
     
      // delete [] pzloc;
      delete [] pzl;

    }
    else if( numberOfDimensions==3 )
    {
      // ---------- PLANAR THREE-DIMENSIONS PARALLEL ------------

      if( debug & 1 )
      {
	printF("$$$$$$$$$$$$ RADIATION KERNEL -- NEWEST 3D PARALLEL VERSION bcInitComplex=%d $$$$$$$$$$$$$\n",bcInitComplex);
	fprintf(debugFile," n1=%d, numberOfModes1=%d, n2=%d, numberOfModes2=%d\n",n1,numberOfModes1,n2,numberOfModes2);
      }
      
      #ifdef USE_PPP
        printF(" ----> RUNNING IN PARALLEL <-----\n");
      #endif

      
      // phi(-numberOfModes1:numberOfModes1,-numberOfModes2:numberOfModes2,numberOfPoles,numFields)
      const int md1 = 2*numberOfModes1+1;
      const int md2 = 2*numberOfModes2+1;
      
      Complex *& phic = dbase.get<Complex*>("phic");
      #define phi(k1,k2,pole,mc) phic[(k1+numberOfModes1) + md1*( (k2+numberOfModes2) + md2*( (pole) + numberOfPoles*(mc)))]
       
      // ComplexArray fold(0:orderInTime-1,-numberOfModes1:numberOfModes1,-numberOfModes2:numberOfModes2,numberOfPoles,numComp) : 
      Complex *&foldc = dbase.get<Complex*>("foldc");
      #define fold(l,k1,k2,pole,mc) foldc[(l)+orderInTimeMinus1*( (k1+numberOfModes1)+md1*( (k2+numberOfModes2)+md2*( (pole) + numberOfPoles*(mc))))]

       
      RealArray & amc = dbase.get<RealArray>("amc");
       
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int k2=-numberOfModes2; k2<=numberOfModes2; k2++ )
	    {
	      for( int k1=-numberOfModes1; k1<=numberOfModes1; k1++ )
	      {
		phi(k1,k2,pole,mc)=0.;
		for( int l=0; l<=orderInTime-2; l++ )
		{
		  fold(l,k1,k2,pole,mc)=0.;
		}
	      }
	    }
	  }
	}

	bcInitComplex=1;
	 
      } // end initialize

      

      // complex*16 zl(n1,n2),zl2(n2)        
      Complex *pzl = new Complex [n1f*n2f];  
      #define zl(i1,i2) pzl[(i1-i1a)+n1f*(i2-i2a)]

      LocalReal scl1=twoPi/period1;
      LocalReal scl2=twoPi/period2;

      // printF("scl=%20.14e\n",scl);
      // Range I1=n1, I2=n2;
       
      // loop over components: 
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {

	fourierTransform.forwardTransform( upLocal,mc, pzl  );
	  
	if( !isEmpty )
	{
	  if( debug & 1 )
	  {
	    fprintf(debugFile,"After fourierTransform.forwardTransform, mc=%d  **NEWER** \n",mc);
	    for( int i2=i2a; i2<=i2b; i2++ )
	    {
	      fprintf(debugFile,"i2=%d: ",i2);
	      for( int i1=i1a; i1<=i1b; i1++ )
	      {
		fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	      }
	      fprintf(debugFile,"]\n");
	    }
	  }

	
	  for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	  {
	    for( int k2=-numberOfModes2; k2<=numberOfModes2; k2++ )     // actual fourier modes
	    {
              const int k2fft = k2>=0 ? k2 : n2+k2;                     // location of mode k2 in fft 
	      assert( k2fft>=0 && k2fft<n2 );

	      if( k2fft<i2a || k2fft>i2b )
		continue;   // no data on this processor

	      for( int k1=-numberOfModes1; k1<=numberOfModes1; k1++ )   // actual fourier modes
	      {
		const int k1fft = k1>=0 ? k1 : n1+k1;                   // location of mode k1 in fft 
		assert( k1fft>=0 && k1fft<n1 );
		if( k1fft<i1a || k1fft>i1b )
		  continue;   // no data on this processor

		LocalReal w= sqrt( SQR(scl1*k1) + SQR(scl2*k2) );

		Complex expBeta = exp(c*beta(pole+1)*dt*w);              // beta's are complex 
		Complex xfact = expBeta;
		Complex adon = xfact*phi(k1,k2,pole,mc);                 // add-on
		for( int l=0; l<=orderInTime-2; l++ )
		{
		  adon += xfact*amc(l)*fold(l,k1,k2,pole,mc);
		  xfact *= expBeta;
		  // printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
		}
		adon += c*dt*amc(-1)*alpha(pole+1)*w*w*zl(k1fft,k2fft);  
		phi(k1,k2,pole,mc)=adon;

		for( int l=orderInTime-2; l>=1; l-- )
		{
		  fold(l,k1,k2,pole,mc)=fold(l-1,k1,k2,pole,mc);
		}
		fold(0,k1,k2,pole,mc)=c*dt*w*w*alpha(pole+1)*zl(k1fft,k2fft);   // alpha's are complex 
	      }
	    } // end for k2 
	  } // end for pole 

	  // ** do this for now: *fix me*
	  for( int i1=i1a; i1<=i1b; i1++ )
	  {
	    for( int i2=i2a; i2<=i2b; i2++ )
	    {
	      zl(i1,i2)=0.;
	    }
	  }

	  // ---- zl(k1,k2) = SUM_pole  phi(k1t,k2t,pole,mc); 
	  for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	  {
	    for( int k2=-numberOfModes2; k2<=numberOfModes2; k2++ )
	    {
              const int k2fft = k2>=0 ? k2 : n2+k2;  // location of mode k2 in fft 
	      assert( k2fft>=0 && k2fft<n2 );
	      if( k2fft<i2a || k2fft>i2b )
		continue;   // no data on this processor
	      
	      for( int k1=-numberOfModes1; k1<=numberOfModes1; k1++ )
	      {
		const int k1fft = k1>=0 ? k1 : n1+k1;  // location of mode k1 in fft 
		assert( k1fft>=0 && k1fft<n1 );
		if( k1fft>=i1a && k1fft<=i1b )
		{
		  zl(k1fft,k2fft) += phi(k1,k2,pole,mc);
		}
		
	      }
	    }
	  } // end for pole 
	
	  if( debug & 1 )
	  {
	    fprintf(debugFile,"Before backwardTransform, mc=%d  zl \n",mc);
	    for( int i2=i2a; i2<=i2b; i2++ )
	    {
	      fprintf(debugFile,"i2=%d: ",i2);
	      for( int i1=i1a; i1<=i1b; i1++ )
	      {
		fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	      }
	      fprintf(debugFile,"]\n");
	    }
	  }

	  // For real data we have the relationship: 
	  //    uHat(-k1,-k2) = conj( uHat(k1,k2) )
	  //    NOTE:   "zl(-k1,-k2)" = zl(n1-k1,n2-k2) for k1>0, k2>0  (negative freq stored here)
	  //
	  //       +------------+------------+
	  //       |            |            |
	  //       |   conj(B)  |     A      |
	  //       |            |   given    |
	  //       |            |            |
	  // k2=0  +------------+------------+
	  //       |            |            |
	  //       |            |     B      |
	  //       |  conj(A)   |   given    |
	  //       |            |            |
	  //       |            |            |
	  //       +------------+------------+
	  //   k1=-n1/2        k1=0        k1=n1/2 

	  // int i2=0;
	  // for( int i1=1; i1<=numberOfModes1; i1++ )
	  // {
	  //   zl(n1-i1,i2) = std::conj(zl(i1,i2));   // x-axis 
	  // }
	    
	  // // now fill in conj(A) and conj(B) regions
	  // for( int i2=1; i2<=numberOfModes2; i2++ )
	  // {
	  //   for( int i1=1; i1<=numberOfModes1; i1++ )
	  //   {
	  //     zl(n1-i1   ,i2) = std::conj(zl(i1,n2-i2)); 
	  //     zl(n1-i1,n2-i2) = std::conj(zl(i1,i2   )); 
	  //   }
	  // }
	  // if( debug & 1 )
	  // {
	  //   fprintf(debugFile,"Before backwardTransform, mc=%d  zl AFTER FILLING IN CONJ \n",mc);
	  //   for( int i2=0; i2<n2; i2++ )
	  //   {
	  //     fprintf(debugFile,"i2=%d: ",i2);
	  //     for( int i1=0; i1<n1; i1++ )
	  //     {
	  // 	fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	  //     }
	  //     fprintf(debugFile,"]\n");
	  //   }
	  // }
	  
	} // end !isEmpty

	fourierTransform.backwardTransform( pzl, hupLocal,mc  );
	  
	if( debug  & 1 && !isEmpty && np==1 )
	{
	  // ** FIX ME FOR PARALLEL 
	  RealArray temp;
	  temp.redim(I1,I2,I3);
	  temp=hupLocal(I1,I2,I3);
	  if( raxis==0 )
	    temp.reshape(I2,I3);
	  else if( raxis==1 )
	    temp.reshape(I1,I3);
	  else
	    temp.reshape(I1,I2);
	    
	  ::display(temp,sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),debugFile,"%13.6e ");
	  //  ::display(hupLocal(I1,I2,I3),sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),"%13.6e ");
	}

      }  // end for mc 
	
      delete [] pzl;


    }
    else
    {
      OV_ABORT("ERROR: unexpected numberOfDimensions");
    }
    
    


  }
  else if( kernelType==cylindrical )
  {
    OV_ABORT("RadiationKernel:parallel: finish me for kernelType==cylindrical");
    
    // bccyld(nda1,ndb1,
    //        *pu, *pf, *ploc, c,radius, dt,numberOfGridPoints1,numberOfFields,numberOfModes1,ns1,
    //        orderOfTimeStepping, *fold,*phi,
    // 	   *npoles, *alpha, *beta,*amc,*fftsave1,bcinit);
  }
  else
  {
    printF("ERROR: un-implemented kernelType=%i\n",kernelType);
    Overture::abort("error");
  }

  // --- assign periodic images ---
  fourierTransform.periodicUpdate( hupLocal ); 
  
  cpuTime+=getCPU()-time;
  return 0;
}

// ========================================================================================
//
/// \brief  Evaluate the radiation kernel. *New* parallel version
/// \param u (input) :  
/// \param hu (output) : the output hu(i) is defined for all i (with periodic images assigned too)
//  
//
//  This routine assumes that the input values run from u(0),...,u(numberOfGridPoints-1)
//
/// \Note: See paper
//      author="B. Alpert and L. Greengard and T. Hagstrom",
//      title="Nonreflecting Boundary Conditions for the Time-Dependent Wave Equation",
//      journal=JCP,
//      volume="180",
//      year="2002",
//      pages="270--296")
// 
// ========================================================================================
int RadiationKernel::evaluateKernelParallelV2( double dt, RealArray & upLocal, RealArray & hupLocal )
{

  LocalReal time=getCPU();
  
  const int myid = max(0,Communication_Manager::My_Process_Number);
  const int np=max(1,Communication_Manager::Number_Of_Processors);


  bool useComplexTransforms= true;   // true = use complex transforms where real transforms could be used instead
  // if( useComplexTransforms )
  //   printF(" RadiationKernel::evaluateKernelParallel: useComplexTransforms\n");
      
  const bool & useFourierTransformClass = dbase.get<bool>("useFourierTransformClass");
  assert( useFourierTransformClass );
  
  const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
  const int rside = dbase.get<int>("rside");
  const int raxis = dbase.get<int>("raxis");

  const int & debug = dbase.get<int>("debug");
  FILE *debugFile=dbase.get<FILE*>("debugFile");
  if( debug & 1 )
  {
    fprintf(debugFile,"evaluateKernelParallel: numberOfDimensions=%d\n",numberOfDimensions);
  }
  

  int & bcInitComplex = dbase.get<int>("bcInitComplex");

  assert( ploc!=NULL );
  
  // double *pu = u.getDataPointer();
  // double *pf = hu.getDataPointer();
  
  // // We assume normal direction is axis=0 for now:
  // int axis=0;  // normal axis 
  // int axisp1 = 1;
  // int axisp2 = 2;
  
  // // make sure u and hu are the same size: 
  // assert( u.getBase(0)==hu.getBase(0) && u.getBound(0)==hu.getBound(0) );
  // assert( u.getBase(axisp1)==hu.getBase(axisp1) && u.getBound(axisp1)==hu.getBound(axisp1) );
  // if( numberOfDimensions==3 )
  // {
  //   assert( u.getBase(axisp2)==hu.getBase(axisp2) && u.getBound(axisp2)==hu.getBound(axisp2) );
  // }
  

  // // Leading dimensions of the u and hu arrays: 
  // const int nda1=u.getBase(0)+1, ndb1=u.getBound(0)+1;  // add one (base one assumed in bcperq21d)
  // const int nda2=u.getBase(1)+1, ndb2=u.getBound(1)+1;  // add one (base one assumed in bcperq21d)
  
  assert( bcInitComplex>=0 );
  
  FourierTransform & fourierTransform = dbase.get<FourierTransform>("fourierTransform");

  // -- Here is the computational box ---
  Range I1,I2,I3;
  IndexBox & fbox = dbase.get<IndexBox>("fbox");
  I1=Range(fbox.base(0),fbox.bound(0));
  I2=Range(fbox.base(1),fbox.bound(1));
  I3=Range(fbox.base(2),fbox.bound(2));
	  
  if( np>1 && debug >0 )
  {
    fprintf(debugFile,"evaluateKernelParallel fbox=[%i,%i][%i,%i][%i,%i]\n",
	    fbox.base(0),fbox.bound(0),fbox.base(1),fbox.bound(1),fbox.base(2),fbox.bound(2));
    fflush(debugFile);
      
    // fclose(debugFile);
      
    // OV_ABORT("evaluateKernelParallel: stop here for now");
      
  }

  const int n1=numberOfGridPoints1;  // number of grid points
  const int n2=numberOfGridPoints2;  // number of grid points

  IndexBox & fftBox = dbase.get<IndexBox>("fftBox");
  // --- Local array sizes for the FFT: 
  // int n1f = fftBox.bound(0)-fftBox.base(0)+1;
  // int n2f = fftBox.bound(1)-fftBox.base(1)+1;
  // int n3f = fftBox.bound(2)-fftBox.base(2)+1;    
  int axisp1, axisp2;
  if( raxis==0 )
  {
    axisp1=1; axisp2=2;
  }
  else if( raxis==1 )
  {
    axisp1=0; axisp2=2;
  }
  else
  {
    axisp1=0; axisp2=1;
  }
  const int n1f = fftBox.bound(axisp1)-fftBox.base(axisp1)+1;
  const int n2f = fftBox.bound(axisp2)-fftBox.base(axisp2)+1;

  // local bounds on loops
  const int i1a=fftBox.base(axisp1), i1b=fftBox.bound(axisp1);  // 0:n1 -> i1a:i1b
  const int i2a=fftBox.base(axisp2), i2b=fftBox.bound(axisp2);  // 0:n2 -> i2a:i2b.

  const int n1by2 = min(n1/2,i1b);   // 0:n1/2 --> i1a:n1by2 
  const int n2by2 = min(n2/2,i2b);   // 0:n1/2 --> i1a:n2by2 

  const int numModes1 = min(numberOfModes1,i1b+1);    // 0..numberOfModes1 -> i1a:numModes1
  const int numModes2 = min(numberOfModes2,i2b+1);    // 0..numberOfModes2 -> i1a:numModes2
  
  

  bool isEmpty = fftBox.isEmpty();

  if( debug & 1 )
  {
    fprintf(debugFile,"evalKernel: fftBox =[%3i,%3i][%3i,%3i][%3i,%3i][%3i,%3i]\n",
	    fftBox.base(0),fftBox.bound(0),
	    fftBox.base(1),fftBox.bound(1),
	    fftBox.base(2),fftBox.bound(2),
	    fftBox.base(3),fftBox.bound(3));
   
    fprintf(debugFile,"local array sizes: [i1a,i1b][i2a,i2b]=[%3d,%3d][%3d,%3d] \n",i1a,i1b,i2a,i2b);
    fprintf(debugFile,"local array sizes: n1f=%d, n2f=%d, isEmpty=%d \n",n1f,n2f,(int)isEmpty);
    fflush(debugFile);
    
  }


  if( kernelType==planar )
  {
    // ------------- Planar Kernel : flat boundary -------------

    // **NEW WAY**
    // Start from code in bcperq21.f
    // We solve for auxilary vaiables phi_j 
    //    (d/dt - c beta_j w) phi_j = c alpha_j w^2 phat 
    //    phat(t) = 
    //    fha(t)t = sum phi_j   :    FourierTranform of solution 

    Complex I(0.,1.);

    const int orderInTime=orderOfTimeStepping; 
    const int orderInTimeMinus1=orderInTime-1;

    int md = numberOfModes1;    // number of Fourier modes to keep, max = n1/2
    if( numberOfDimensions==3 )
      md *= numberOfModes2;

    int numComp=numberOfFields;  


    if( !dbase.has_key("phic") )
    {
      // ---- initialize work spaces and data ----

      dbase.put<Complex*>("phic")=NULL;
      Complex *& phic = dbase.get<Complex*>("phic");
      if( numberOfDimensions==2 )
        phic = new Complex [md*numberOfPoles*numComp];            
      else 
        phic = new Complex [(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields];            

      dbase.put<Complex*>("foldc")=NULL;
      Complex *&foldc = dbase.get<Complex*>("foldc");
      if( numberOfDimensions==2 )
        foldc = new Complex[orderInTimeMinus1*md*numberOfPoles*numComp];
      else
	foldc = new Complex[orderInTimeMinus1*(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields];


      // array to hold Adams-Moulton coefficients 
      dbase.put<RealArray>("amc");
      RealArray & amc = dbase.get<RealArray>("amc");
      amc.redim(Range(-1,orderInTime-2));


      // --- initialize alpha and beta arrays defining the poles in the non-local kernel approximation
      initializePoles();

    }

    Complex *& alphac = dbase.get<Complex*>("alphac");
    Complex *& betac = dbase.get<Complex*>("betac");
    #define alpha(i) alphac[(i-1)]
    #define beta(i)  betac[(i-1)]


    


    if( numberOfDimensions==2 )  
    {
      // ================ PLANAR TWO DIMENSIONS ======================


      Complex *& phic = dbase.get<Complex*>("phic");
      #define phi(k,pole,mc) phic[(k)+md*( (pole) + numberOfPoles*(mc))]
       
      // ComplexArray fold(0:orderInTime-1,md,numberOfPoles,numComp) : 
      Complex *&foldc = dbase.get<Complex*>("foldc");
      #define fold(l,k,pole,mc) foldc[(l)+orderInTimeMinus1*( (k)+md*( (pole) + numberOfPoles*(mc)))]

      RealArray & amc = dbase.get<RealArray>("amc");
      
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int k=0; k<numberOfModes1; k++ )
	    {
	      phi(k,pole,mc)=0.;
	      for( int l=0; l<=orderInTime-2; l++ )
	      {
		fold(l,k,pole,mc)=0.;
	      }
	    }
	  }
	}
	 
	bcInitComplex=1;
	 
      } // end initialize

      RealArray ploc(n1); // , p(n1,numComp);
       
      // For complex FFT's
      //  zlocr(i) = real-part
      //  zloci(i) = imag-part
      const int n1c=2*n1;
      LocalReal *pzloc = new LocalReal [n1*2];  // holds complex transforms 
#define zlocr(i) pzloc[2*(i-i1a)]
#define zloci(i) pzloc[1+2*(i-i1a)]

      LocalReal scl=twoPi/period1;
      // printF("scl=%20.14e\n",scl);
      // I1=n1;
       

      // loop over components: 
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {

	fourierTransform.forwardTransform( upLocal,mc, pzloc  );

	if( !isEmpty )
	{
	  if( debug & 1 )
	  {
	    fprintf(debugFile,"After fourierTransform.forwardTransform, mc=%d  **NEWER** \n",mc);
	    for( int i=0; i<n1; i++ ){ fprintf(debugFile,"[%13.6e,%13.6e]",zlocr(i),zloci(i)); }   // 
	    fprintf(debugFile,"]\n");
	  }

	  for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	  {
	    // for( int k=0; k<numberOfModes1; k++ )   // "k=1" loop over actual Fourier modes k that we compute 
	    for( int k=i1a; k<numModes1; k++ )   // "k=1" loop over actual Fourier modes k that we compute 
	    {
	      const int kp = k+1;
	      LocalReal w = scl*kp;
	      Complex expBeta = exp(c*beta(pole+1)*dt*w);  // beta's are complex 
	      Complex xfact = expBeta;
	      Complex adon = xfact*phi(k,pole,mc);         // add-on 
	      // printF("* adon=[%12.4e,%12.4e]\n",std::real(adon),std::imag(adon));
	      for( int l=0; l<=orderInTime-2; l++ )
	      {
		adon += xfact*amc(l)*fold(l,k,pole,mc);
		xfact *= expBeta;
		// printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
	      }


	      Complex phat;
	      phat = zlocr(kp) + zloci(kp) *I;

	      adon += c*dt*amc(-1)*alpha(pole+1)*w*w*phat;
	      phi(k,pole,mc)=adon;

	      // printF("expBeta=[%12.4e,%12.4e]\n",std::real(expBeta),std::imag(expBeta));
	      //printF("phi(%d,%d,%d)=[%12.4e,%12.4e]\n",k,pole,mc,std::real(phi(k,pole,mc)),std::imag(phi(k,pole,mc)));
	    
	      for( int l=orderInTime-2; l>=1; l-- )
	      {
		fold(l,k,pole,mc)=fold(l-1,k,pole,mc);
	      }
	      fold(0,k,pole,mc)=c*dt*w*w*alpha(pole+1)*phat;   // alpha's are complex 
	    }  // end for k
	  }  // end for pole


	  // real FFT:      cos(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  
	  // complex FFT:   cos(0*x) sin(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  

	  // convert complex FFT to real:  (for n1 even)
	  //   pHat( k) = pzloc[k] ,    k=0,1,2,...,n1/2-1
	  //   pHat(-k) = pzloc[n-k]    k=  1,2,...,n1/2
        
	  ploc=0.;  // fill in hHat 
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    // for( int k=0; k<numberOfModes1; k++ )   // loop over Fourier modes we compute 
	    for( int k=i1a; k<numModes1; k++ )   // loop over Fourier modes we compute 
	    {
	      const int kp = k+1;
	      ploc(2*kp-1) += std::real(phi(k,pole,mc));    
	      ploc(2*kp  ) += std::imag(phi(k,pole,mc)); 

	    }
	  }
	

	  // do this for now: *fix me*
	  // for( int i=0; i<n1; i++ ) 
	  for( int i=i1a; i<=i1b; i++ ) 
	  {
	    zlocr(i)=0.;
	    zloci(i)=0.;
	  }

	  zlocr(0) = ploc(0);
	  zloci(0) = 0.;
	  int j=1;
	  // for( int i=1; i<n1/2; i++ )
	  for( int i=max(1,i1a); i<n1by2; i++ )
	  {
	    // real FFT:      cos(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  
	    // complex FFT:   cos(0*x) sin(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  


	    // convert complex FFT to real:
	    //   pHat( k) = pzloc[k] ,    k=0,1,2,...,n1/2-1
	    //   pHat(-k) = pzloc[n-k]    k=0,1,2,...,n1/2

	    // k=1,2,...,n1/2-1
	    zlocr(i) = ploc(j  );
	    zloci(i) = ploc(j+1);

	    // negative k are the complex conjugate:
            // *** FIX ME FOR PARALLEL ***
	    zlocr(n1-i) =  ploc(j  );
	    zloci(n1-i) = -ploc(j+1);

	    j+=2;
	  
	  }
	}
	// end if !isEmpty

	fourierTransform.backwardTransform( pzloc, hupLocal,mc  );

	if( debug & 1 )
	{
	  ::display(hupLocal,sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),
		    debugFile,"%13.6e ");
	  // ::display(hupLocal(I1,I2,I3),sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),
	  // 	    debugFile,"%13.6e ");
	}
	  
      }  // end for component mc
        
     
      delete [] pzloc;

    }
    else if( numberOfDimensions==3 )
    {
      // ---------- PLANAR THREE-DIMENSIONS PARALLEL ------------

      if( debug & 1 )
      {
	printF("$$$$$$$$$$$$ RADIATION KERNEL -- NEW 3D PARALLEL VERSION bcInitComplex=%d $$$$$$$$$$$$$\n",bcInitComplex);
	fprintf(debugFile," n1=%d, numberOfModes1=%d, n2=%d, numberOfModes2=%d\n",n1,numberOfModes1,n2,numberOfModes2);
      }
      
      #ifdef USE_PPP
        printF(" ----> RUNNING IN PARALLEL <-----\n");
      #endif

      
      // phi(0:numberOfModes1,-numberOfModes2:numberOfModes2,numberOfPoles,numFields)
      const int md2=2*numberOfModes2+1;
      const int numberOfModes1p1 = numberOfModes1+1;
      
      Complex *& phic = dbase.get<Complex*>("phic");
#define phi(k1,k2,pole,mc) phic[(k1)+numberOfModes1p1*((k2+numberOfModes2)+md2*( (pole) + numberOfPoles*(mc)))]
       
      // ComplexArray fold(0:orderInTime-1,md,numberOfPoles,numComp) : 
      Complex *&foldc = dbase.get<Complex*>("foldc");
#define fold(l,k1,k2,pole,mc) foldc[(l)+orderInTimeMinus1*( (k1)+numberOfModes1p1*( (k2+numberOfModes2)+md2*( (pole) + numberOfPoles*(mc))))]

       
      RealArray & amc = dbase.get<RealArray>("amc");
       
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int k2=-numberOfModes2; k2<=numberOfModes2; k2++ )
	    {
	      for( int k1=0; k1<=numberOfModes1; k1++ )
	      {
		phi(k1,k2,pole,mc)=0.;
		for( int l=0; l<=orderInTime-2; l++ )
		{
		  fold(l,k1,k2,pole,mc)=0.;
		}
	      }
	    }
	  }
	}

	bcInitComplex=1;
	 
      } // end initialize

      

      // complex*16 zl(n1,n2),zl2(n2)        
      Complex * pzl = new Complex [n1*n2];  
      #define zl(i1,i2) pzl[(i1)+n1*(i2)]

      LocalReal scl1=twoPi/period1;
      LocalReal scl2=twoPi/period2;

      // printF("scl=%20.14e\n",scl);
      // Range I1=n1, I2=n2;
       
      // loop over components: 
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {

	fourierTransform.forwardTransform( upLocal,mc, pzl  );
	  
	if( !isEmpty )
	{
	  if( debug & 1 )
	  {
	    fprintf(debugFile,"After fourierTransform.forwardTransform, mc=%d  **NEWER** \n",mc);
	    for( int i2=0; i2<n2; i2++ )
	    {
	      fprintf(debugFile,"i2=%d: ",i2);
	      for( int i1=0; i1<n1; i1++ )
	      {
		fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	      }
	      fprintf(debugFile,"]\n");
	    }
	  }

	
	  for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	  {
	    for( int k2=0; k2<n2; k2++ )   // "k=1" loop over actual Fourier modes k that we compute 
	    {
	      int k2p=k2+1;
	      int k2t;
	      if( k2p<=(numberOfModes2+1) )
	      {
		k2t=k2p-1;
	      }
	      else if( (k2p-1-n2) >= (-numberOfModes2) )
	      {
		k2t=k2p-1-n2;
	      }
	      else
	      {
		continue;
	      }
	      for( int k1t=0; k1t<=numberOfModes1; k1t++ )
	      {
		int k1=k1t;  // check ??
		int k1p=k1t+1;

		LocalReal w= sqrt( SQR(scl1*k1t) + SQR(scl2*k2t) );

		Complex expBeta = exp(c*beta(pole+1)*dt*w);  // beta's are complex 
		Complex xfact = expBeta;
		Complex adon = xfact*phi(k1t,k2t,pole,mc);         // add-on
		for( int l=0; l<=orderInTime-2; l++ )
		{
		  adon += xfact*amc(l)*fold(l,k1t,k2t,pole,mc);
		  xfact *= expBeta;
		  // printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
		}
		adon += c*dt*amc(-1)*alpha(pole+1)*w*w*zl(k1,k2);  // check k1,k2
		phi(k1t,k2t,pole,mc)=adon;

		for( int l=orderInTime-2; l>=1; l-- )
		{
		  fold(l,k1t,k2t,pole,mc)=fold(l-1,k1t,k2t,pole,mc);
		}
		fold(0,k1t,k2t,pole,mc)=c*dt*w*w*alpha(pole+1)*zl(k1,k2);   // alpha's are complex 
	      }
	    } // end for k2 
	  } // end for pole 

	  // ** do this for now: *fix me*
	  for( int i1=0; i1<n1; i1++ )
	  {
	    for( int i2=0; i2<n2; i2++ )
	      zl(i1,i2)=0.;
	  }

	  // ---- zl(k1,k2) = SUM_pole  phi(k1t,k2t,pole,mc); 
	  for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	  {
	    for( int k1t=0; k1t<=numberOfModes1; k1t++ )
	    {
	      int k1=k1t;
	      for( int k2t=0; k2t<=numberOfModes2; k2t++ )
	      {
		int k2=k2t;
		zl(k1,k2) += phi(k1t,k2t,pole,mc);
	      }
	      for( int k2t=-numberOfModes2; k2t<=-1; k2t++ )
	      {
		int k2 = n2 + k2t;  // check 
		zl(k1,k2) += phi(k1t,k2t,pole,mc);
	      }
	    }
	  }
	  // end for pole 
	
	  if( debug & 1 )
	  {
	    fprintf(debugFile,"Before backwardTransform, mc=%d  zl before zlNew \n",mc);
	    for( int i2=0; i2<n2; i2++ )
	    {
	      fprintf(debugFile,"i2=%d: ",i2);
	      for( int i1=0; i1<n1; i1++ )
	      {
		fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	      }
	      fprintf(debugFile,"]\n");
	    }
	  }

	  // For real data we have the relationship: 
	  //    uHat(-k1,-k2) = conj( uHat(k1,k2) )
	  //    NOTE:   "zl(-k1,-k2)" = zl(n1-k1,n2-k2) for k1>0, k2>0  (negative freq stored here)
	  //
	  //       +------------+------------+
	  //       |            |            |
	  //       |   conj(B)  |     A      |
	  //       |            |   given    |
	  //       |            |            |
	  // k2=0  +------------+------------+
	  //       |            |            |
	  //       |            |     B      |
	  //       |  conj(A)   |   given    |
	  //       |            |            |
	  //       |            |            |
	  //       +------------+------------+
	  //   k1=-n1/2        k1=0        k1=n1/2 

	  int i2=0;
	  for( int i1=1; i1<=numberOfModes1; i1++ )
	  {
	    zl(n1-i1,i2) = std::conj(zl(i1,i2));   // x-axis 
	  }
	    
	  // now fill in conj(A) and conj(B) regions
	  for( int i2=1; i2<=numberOfModes2; i2++ )
	  {
	    for( int i1=1; i1<=numberOfModes1; i1++ )
	    {
	      zl(n1-i1   ,i2) = std::conj(zl(i1,n2-i2)); 
	      zl(n1-i1,n2-i2) = std::conj(zl(i1,i2   )); 
	    }
	  }
	  if( debug & 1 )
	  {
	    fprintf(debugFile,"Before backwardTransform, mc=%d  zl AFTER FILLING IN CONJ \n",mc);
	    for( int i2=0; i2<n2; i2++ )
	    {
	      fprintf(debugFile,"i2=%d: ",i2);
	      for( int i1=0; i1<n1; i1++ )
	      {
		fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	      }
	      fprintf(debugFile,"]\n");
	    }
	  }
	  
	} // end !isEmpty

	fourierTransform.backwardTransform( pzl, hupLocal,mc  );
	  
	if( debug  & 1 && !isEmpty )
	{
	  // ** FIX ME FOR PARALLEL 
	  RealArray temp;
	  temp.redim(I1,I2,I3);
	  temp=hupLocal(I1,I2,I3);
	  if( raxis==0 )
	    temp.reshape(I2,I3);
	  else if( raxis==1 )
	    temp.reshape(I1,I3);
	  else
	    temp.reshape(I1,I2);
	    
	  ::display(temp,sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),debugFile,"%13.6e ");
	  //  ::display(hupLocal(I1,I2,I3),sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),"%13.6e ");
	}

      }  // end for mc 
	
      delete [] pzl;


    }
    else
    {
      OV_ABORT("ERROR: unexpected numberOfDimensions");
    }
    
    


  }
  else if( kernelType==cylindrical )
  {
    OV_ABORT("RadiationKernel:parallel: finish me for kernelType==cylindrical");
    
    // bccyld(nda1,ndb1,
    //        *pu, *pf, *ploc, c,radius, dt,numberOfGridPoints1,numberOfFields,numberOfModes1,ns1,
    //        orderOfTimeStepping, *fold,*phi,
    // 	   *npoles, *alpha, *beta,*amc,*fftsave1,bcinit);
  }
  else
  {
    printF("ERROR: un-implemented kernelType=%i\n",kernelType);
    Overture::abort("error");
  }

  // --- assign periodic images ---
  fourierTransform.periodicUpdate( hupLocal ); 
  
  cpuTime+=getCPU()-time;
  return 0;
}

// ========================================================================================
//   *** OLD ***
///
/// \brief  Evaluate the radiation kernel. *New* parallel version
/// \param u (input) :  
/// \param hu (output) : the output hu(i) is defined for all i (with periodic images assigned too)
//  
//
//  This routine assumes that the input values run from u(0),...,u(numberOfGridPoints-1)
//
/// \Note: See paper
//      author="B. Alpert and L. Greengard and T. Hagstrom",
//      title="Nonreflecting Boundary Conditions for the Time-Dependent Wave Equation",
//      journal=JCP,
//      volume="180",
//      year="2002",
//      pages="270--296")
// 
// ========================================================================================
int RadiationKernel::evaluateKernelParallelV1( double dt, RealArray & u, RealArray & hu,
						RealArray & upLocal, RealArray & hupLocal )
{
  LocalReal time=getCPU();
  
  bool useComplexTransforms= true;   // true = use complex transforms where real transforms could be used instead
  if( useComplexTransforms )
    printF(" RadiationKernel::evaluateKernelParallel: useComplexTransforms\n");
      
  const bool & useFourierTransformClass = dbase.get<bool>("useFourierTransformClass");
  assert( useFourierTransformClass );
  
  const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
  const int rside = dbase.get<int>("rside");
  const int raxis = dbase.get<int>("raxis");

  const int & debug = dbase.get<int>("debug");
  FILE *debugFile=dbase.get<FILE*>("debugFile");
  if( debug & 1 )
  {
    fprintf(debugFile,"evaluateKernelParallel: numberOfDimensions=%d\n",numberOfDimensions);
  }
  

  int & bcInitComplex = dbase.get<int>("bcInitComplex");

  assert( ploc!=NULL );
  
  double *pu = u.getDataPointer();
  double *pf = hu.getDataPointer();
  
  // We assume normal direction is axis=0 for now:
  int axis=0;  // normal axis 
  int axisp1 = 1;
  int axisp2 = 2;
  
  // make sure u and hu are the same size: 
  assert( u.getBase(0)==hu.getBase(0) && u.getBound(0)==hu.getBound(0) );
  assert( u.getBase(axisp1)==hu.getBase(axisp1) && u.getBound(axisp1)==hu.getBound(axisp1) );
  if( numberOfDimensions==3 )
  {
    assert( u.getBase(axisp2)==hu.getBase(axisp2) && u.getBound(axisp2)==hu.getBound(axisp2) );
  }
  

  // Leading dimensions of the u and hu arrays: 
  const int nda1=u.getBase(0)+1, ndb1=u.getBound(0)+1;  // add one (base one assumed in bcperq21d)
  const int nda2=u.getBase(1)+1, ndb2=u.getBound(1)+1;  // add one (base one assumed in bcperq21d)
  
  assert( bcInitComplex>=0 );
  
  if( kernelType==planar )
  {
    // ------------- Planar Kernel : flat boundary -------------

    // **NEW WAY**
    // Start from code in bcperq21.f
    // We solve for auxilary vaiables phi_j 
    //    (d/dt - c beta_j w) phi_j = c alpha_j w^2 phat 
    //    phat(t) = 
    //    fha(t)t = sum phi_j   :    FourierTranform of solution 

    Complex I(0.,1.);

    const int orderInTime=orderOfTimeStepping; 
    const int orderInTimeMinus1=orderInTime-1;

    int md = numberOfModes1;    // number of Fourier modes to keep, max = n1/2
    if( numberOfDimensions==3 )
      md *= numberOfModes2;

    int numComp=numberOfFields;  
    int n1=numberOfGridPoints1;  // number of grid points
    int ns1 =2*n1+15;            // for real FFT
    int ns1c=4*n1+15;            // for complex FFT
       
    int n2=numberOfGridPoints2;  // number of grid points
    int ns2 =2*n2+15;            // for real FFT
    int ns2c=4*n2+15;            // for complex FFT
    


    if( !dbase.has_key("phic") )
    {

      dbase.put<Complex*>("phic")=NULL;
      Complex *& phic = dbase.get<Complex*>("phic");
      if( numberOfDimensions==2 )
        phic = new Complex [md*numberOfPoles*numComp];            
      else 
        phic = new Complex [(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields];            

      dbase.put<Complex*>("foldc")=NULL;
      Complex *&foldc = dbase.get<Complex*>("foldc");
      if( numberOfDimensions==2 )
        foldc = new Complex[orderInTimeMinus1*md*numberOfPoles*numComp];
      else
	foldc = new Complex[orderInTimeMinus1*(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields];


      // phi = new double [(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields*2];  // complex
      // fold= new double [(orderOfTimeStepping-1)*(numberOfModes1+1)*(2*numberOfModes2+1)*numberOfPoles*numberOfFields*2]; // complex


      // array to hold Adams-Moulton coefficients 
      dbase.put<RealArray>("amc");
      RealArray & amc = dbase.get<RealArray>("amc");
      amc.redim(Range(-1,orderInTime-2));


      // --- initialize alpha and beta arrays defining the poles in the non-local kernel approximation
      initializePoles();

    }

    FourierTransform & fourierTransform = dbase.get<FourierTransform>("fourierTransform");

    Complex *& alphac = dbase.get<Complex*>("alphac");
    Complex *& betac = dbase.get<Complex*>("betac");
    #define alpha(i) alphac[(i-1)]
    #define beta(i)  betac[(i-1)]


    // --- transition to new way ---
    RealArray uLocal;
    Range I1,I2,I3;
    {
      IndexBox & fbox = dbase.get<IndexBox>("fbox");
      I1=Range(fbox.base(0),fbox.bound(0));
      I2=Range(fbox.base(1),fbox.bound(1));
      I3=Range(fbox.base(2),fbox.bound(2));
	  
      uLocal.redim(I1,I2,I3);
      uLocal=0.;
    }


    if( numberOfDimensions==2 )  // *new way* Finish me for 3D
    {


      Complex *& phic = dbase.get<Complex*>("phic");
#define phi(k,pole,mc) phic[(k)+md*( (pole) + numberOfPoles*(mc))]
       
      // ComplexArray fold(0:orderInTime-1,md,numberOfPoles,numComp) : 
      Complex *&foldc = dbase.get<Complex*>("foldc");
#define fold(l,k,pole,mc) foldc[(l)+orderInTimeMinus1*( (k)+md*( (pole) + numberOfPoles*(mc)))]

       
      // RealArray & fftsave = dbase.get<RealArray>("fftsave"); 
      // RealArray &fftsavec = dbase.get<RealArray>("fftsavec");     

      RealArray & amc = dbase.get<RealArray>("amc");
       
       
      
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int k=0; k<numberOfModes1; k++ )
	    {
	      phi(k,pole,mc)=0.;
	      for( int l=0; l<=orderInTime-2; l++ )
	      {
		fold(l,k,pole,mc)=0.;
	      }
	    }
	  }
	}
	 
	bcInitComplex=1;
	 
      } // end initialize

      RealArray ploc(n1); // , p(n1,numComp);
      // RealArray f(n1);  // put answer here for now 
       
      // For complex FFT's
      //  zlocr(i) = real-part
      //  zloci(i) = imag-part
      const int n1c=2*n1;
      LocalReal *pzloc = new LocalReal [n1*2];  // holds complex transforms 
#define zlocr(i) pzloc[2*(i)]
#define zloci(i) pzloc[1+2*(i)]

      LocalReal scl=twoPi/period1;
      // printF("scl=%20.14e\n",scl);
      // I1=n1;
       
      // printF("u.getBase(0)=%d\n",u.getBase(0));
      

      // loop over components: 
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {

	fourierTransform.forwardTransform( upLocal,mc, pzloc  );

	if( debug & 1 )
	{
	  fprintf(debugFile,"After fourierTransform.forwardTransform, mc=%d  **NEWER** \n",mc);
	  for( int i=0; i<n1; i++ ){ fprintf(debugFile,"[%13.6e,%13.6e]",zlocr(i),zloci(i)); }   // 
	  fprintf(debugFile,"]\n");
	}

	for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	{
	  for( int k=0; k<numberOfModes1; k++ )   // "k=1" loop over actual Fourier modes k that we compute 
	  {
	    const int kp = k+1;
	    LocalReal w = scl*kp;
	    Complex expBeta = exp(c*beta(pole+1)*dt*w);  // beta's are complex 
	    Complex xfact = expBeta;
	    Complex adon = xfact*phi(k,pole,mc);         // add-on 
	    // printF("* adon=[%12.4e,%12.4e]\n",std::real(adon),std::imag(adon));
	    for( int l=0; l<=orderInTime-2; l++ )
	    {
	      adon += xfact*amc(l)*fold(l,k,pole,mc);
	      xfact *= expBeta;
	      // printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
	    }
	    Complex phat;
            phat = zlocr(kp) + zloci(kp) *I;

	    adon += c*dt*amc(-1)*alpha(pole+1)*w*w*phat;
	    phi(k,pole,mc)=adon;

	    // printF("expBeta=[%12.4e,%12.4e]\n",std::real(expBeta),std::imag(expBeta));
	    //printF("phi(%d,%d,%d)=[%12.4e,%12.4e]\n",k,pole,mc,std::real(phi(k,pole,mc)),std::imag(phi(k,pole,mc)));
	    
	    for( int l=orderInTime-2; l>=1; l-- )
	    {
	      fold(l,k,pole,mc)=fold(l-1,k,pole,mc);
	    }
	    fold(0,k,pole,mc)=c*dt*w*w*alpha(pole+1)*phat;   // alpha's are complex 
	  }  // end for k
	}  // end for pole


        // real FFT:      cos(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  
        // complex FFT:   cos(0*x) sin(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  

        // convert complex FFT to real:  (for n1 even)
	//   pHat( k) = pzloc[k] ,    k=0,1,2,...,n1/2-1
	//   pHat(-k) = pzloc[n-k]    k=  1,2,...,n1/2
        
	ploc=0.;  // fill in hHat 
	for( int pole=0; pole<numberOfPoles; pole++ )
	{
	  for( int k=0; k<numberOfModes1; k++ )   // loop over Fourier modes we compute 
	  {
	    const int kp = k+1;
	    ploc(2*kp-1) += std::real(phi(k,pole,mc));    
	    ploc(2*kp  ) += std::imag(phi(k,pole,mc)); 

	  }
	}
	

	// do this for now: 
	for( int i=0; i<n1; i++ ) 
	{
	  zlocr(i)=0.;
	  zloci(i)=0.;
	}

	zlocr(0) = ploc(0);
	zloci(0) = 0.;
	int j=1;
	for( int i=1; i<n1/2; i++ )
	{
	  // real FFT:      cos(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  
	  // complex FFT:   cos(0*x) sin(0*x) cos(1*x) sin(1*x) cos(2*x) sin(2*x)  


	  // convert complex FFT to real:
	  //   pHat( k) = pzloc[k] ,    k=0,1,2,...,n1/2-1
	  //   pHat(-k) = pzloc[n-k]    k=0,1,2,...,n1/2

	  // k=1,2,...,n1/2-1
	  zlocr(i) = ploc(j  );
	  zloci(i) = ploc(j+1);
	  // negative k are the complex conjugate:
	  zlocr(n1-i) =  ploc(j  );
	  zloci(n1-i) = -ploc(j+1);

	  j+=2;
	  
	}


	fourierTransform.backwardTransform( pzloc, hupLocal,mc  );

	// ::display(uLocal,"uLocal after backward");

	// assert( rside==0 ); // *** FIX ME **
	// if( raxis==0 )
	// {
	//   for( int i=0; i<n1; i++ )
	//     hu(i,mc)=hupLocal(0,i,0);
	// }
	// else
	// {
	//   for( int i=0; i<n1; i++ )
	//     hu(i,mc)=hupLocal(i,0,0);
	// }

	// if( true )
	// {
	//   printF("After fourierTransform.backwardTransform, mc=%d  **NEWER** \n hu=[",mc);
	//   for( int i=0; i<n1; i++ ){ printF("%13.6e,",hu(i,mc)); }   // 
	//   printF("]\n");
	// }
	if( true )
	{
	  ::display(hupLocal(I1,I2,I3),sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),"%13.6e ");
	}


      }  // end for component mc
        
     
      delete [] pzloc;

    }
    else if( numberOfDimensions==3 )
    {
      // ---------- *NEW* 3D version for parallel ------------

      printF("$$$$$$$$$$$$ RADIATION KERNEL -- NEW 3D PARALLEL VERSION bcInitComplex=%d $$$$$$$$$$$$$\n",bcInitComplex);
      printF(" n1=%d, numberOfModes1=%d, n2=%d, numberOfModes2=%d\n",n1,numberOfModes1,n2,numberOfModes2);
      
      #ifdef USE_PPP
        printF(" ----> RUNNING IN PARALLEL <-----\n");
      #endif

      
      // phi(0:numberOfModes1,-numberOfModes2:numberOfModes2,numberOfPoles,numFields)
      const int md2=2*numberOfModes2+1;
      const int numberOfModes1p1 = numberOfModes1+1;
      
      Complex *& phic = dbase.get<Complex*>("phic");
#define phi(k1,k2,pole,mc) phic[(k1)+numberOfModes1p1*((k2+numberOfModes2)+md2*( (pole) + numberOfPoles*(mc)))]
       
      // ComplexArray fold(0:orderInTime-1,md,numberOfPoles,numComp) : 
      Complex *&foldc = dbase.get<Complex*>("foldc");
#define fold(l,k1,k2,pole,mc) foldc[(l)+orderInTimeMinus1*( (k1)+numberOfModes1p1*( (k2+numberOfModes2)+md2*( (pole) + numberOfPoles*(mc))))]

       
      RealArray & amc = dbase.get<RealArray>("amc");
       
      if( bcInitComplex==0 )
      {
	// --- initialize  ----    

	// get coefficients for Adams-Moulton
	AMCOF(amc(-1),orderInTime);

	for( int mc=0; mc<numComp; mc++ )
	{
	  for( int pole=0; pole<numberOfPoles; pole++ )
	  {
	    for( int k2=-numberOfModes2; k2<=numberOfModes2; k2++ )
	    {
	      for( int k1=0; k1<=numberOfModes1; k1++ )
	      {
		phi(k1,k2,pole,mc)=0.;
		for( int l=0; l<=orderInTime-2; l++ )
		{
		  fold(l,k1,k2,pole,mc)=0.;
		}
	      }
	    }
	  }
	}

	bcInitComplex=1;
	 
      } // end initialize

      

      // RealArray pl1(n1);   // hold temp values 

      // complex*16 zl(n1,n2),zl2(n2)        
      Complex * pzl = new Complex [n1*n2];  
#define zl(i1,i2) pzl[(i1)+n1*(i2)]

      Complex * pzlNew = new Complex [n1*n2];  
#define zlNew(i1,i2) pzlNew[(i1)+n1*(i2)]

      // Could not pas this to fortran (?)
//       Complex * pzl2 = new Complex [n2]; 
// #define zl2(i2) pzl2[(i2)]

      LocalReal scl1=twoPi/period1;
      LocalReal scl2=twoPi/period2;

      // printF("scl=%20.14e\n",scl);
      // Range I1=n1, I2=n2;
       
      // loop over components: 
      for( int mc=0; mc<numComp; mc++ )  // "i" 
      {

	fourierTransform.forwardTransform( upLocal,mc, pzlNew  );
	  
	if( debug & 1 )
	{
	  fprintf(debugFile,"After fourierTransform.forwardTransform, mc=%d  **NEWER** \n",mc);
	  for( int i2=0; i2<n2; i2++ )
	  {
	    fprintf(debugFile,"i2=%d: ",i2);
	    for( int i1=0; i1<n1; i1++ )
	    {
	      fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zlNew(i1,i2)),std::imag(zlNew(i1,i2)));
	    }
	    fprintf(debugFile,"]\n");
	  }
	}
	
	for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	{
	  for( int k2=0; k2<n2; k2++ )   // "k=1" loop over actual Fourier modes k that we compute 
	  {
            int k2p=k2+1;
            int k2t;
	    if( k2p<=(numberOfModes2+1) )
	    {
	      k2t=k2p-1;
	    }
	    else if( (k2p-1-n2) >= (-numberOfModes2) )
	    {
	      k2t=k2p-1-n2;
	    }
	    else
	    {
	      continue;
	    }
	    for( int k1t=0; k1t<=numberOfModes1; k1t++ )
	    {
              int k1=k1t;  // check ??
	      int k1p=k1t+1;

	      LocalReal w= sqrt( SQR(scl1*k1t) + SQR(scl2*k2t) );

	      Complex expBeta = exp(c*beta(pole+1)*dt*w);  // beta's are complex 
	      Complex xfact = expBeta;
	      Complex adon = xfact*phi(k1t,k2t,pole,mc);         // add-on
	      for( int l=0; l<=orderInTime-2; l++ )
	      {
		adon += xfact*amc(l)*fold(l,k1t,k2t,pole,mc);
		xfact *= expBeta;
		// printF("l=%d, adon=[%12.4e,%12.4e]\n",l,std::real(adon),std::imag(adon));
	      }
	      adon += c*dt*amc(-1)*alpha(pole+1)*w*w*zl(k1,k2);  // check k1,k2
              phi(k1t,k2t,pole,mc)=adon;

	      for( int l=orderInTime-2; l>=1; l-- )
	      {
		fold(l,k1t,k2t,pole,mc)=fold(l-1,k1t,k2t,pole,mc);
	      }
	      fold(0,k1t,k2t,pole,mc)=c*dt*w*w*alpha(pole+1)*zl(k1,k2);   // alpha's are complex 
	    }
	  } // end for k2 
	} // end for pole 

	for( int i1=0; i1<n1; i1++ )
	{
          for( int i2=0; i2<n2; i2++ )
            zl(i1,i2)=0.;
	}

        // ---- zl(k1,k2) = SUM_pole  phi(k1t,k2t,pole,mc); 
	for( int pole=0; pole<numberOfPoles; pole++ )   // "j" loop over number of poles in the approx. of the Kernel
	{
	  for( int k1t=0; k1t<=numberOfModes1; k1t++ )
	  {
	    int k1=k1t;
	    for( int k2t=0; k2t<=numberOfModes2; k2t++ )
	    {
	      int k2=k2t;
	      zl(k1,k2) += phi(k1t,k2t,pole,mc);
	    }
	    for( int k2t=-numberOfModes2; k2t<=-1; k2t++ )
	    {
              int k2 = n2 + k2t;  // check 
	      zl(k1,k2) += phi(k1t,k2t,pole,mc);
	    }
	  }
	}
	// end for pole 
	
	// for( int i2=0; i2<n2; i2++ )
	// {
	//   for( int i1=0; i1<n1; i1++ )
	//   {
	//     zlNew(i1,i2) = 0.;
	//   }
	// }

	// For real data:
	//    uHat(-k1,-k2) = conj( uHat(k1,k2) )
	//
	//       +------------+------------+
	//       |            |            |
	//       |   conj(B)  |     A      |
	//       |            |   given    |
	//       |            |            |
	// k2=0  +------------+------------+
	//       |            |            |
	//       |            |     B      |
	//       |  conj(A)   |   given    |
	//       |            |            |
	//       |            |            |
	//       +------------+------------+
	//   k1=-n1/2        k1=0        k1=n1/2 
	for( int i2=0; i2<n2; i2++ )
	{
	  for( int i1=numberOfModes1+1; i1<n1; i1++ )
	  {
	    zlNew(i1,i2)=0.;  // zero out (mostly unused) -- some of these will be assigned below
	  }
	  for( int i1=0; i1<=numberOfModes1; i1++ )
	  {
	    zlNew(i1,i2) = zl(i1,i2);   
	    if( i2==0 )
	      zlNew(n1-i1,i2) = std::conj(zl(i1,i2));   // x-axis 
	  }
	    
	}
	// now fill in conj(A) and conj(B) regions
	for( int i2=1; i2<=numberOfModes2; i2++ )
	{
	  for( int i1=1; i1<=numberOfModes1; i1++ )
	  {
	    zlNew(n1-i1   ,i2) = std::conj(zl(i1,n2-i2)); 
	    zlNew(n1-i1,n2-i2) = std::conj(zl(i1,i2   )); 
	  }
	}
	if( debug & 1 )
	{
	  fprintf(debugFile,"Before backwardTransform, mc=%d  **NEWER** zlNew \n",mc);
	  for( int i2=0; i2<n2; i2++ )
	  {
	    fprintf(debugFile,"i2=%d: ",i2);
	    for( int i1=0; i1<n1; i1++ )
	    {
	      fprintf(debugFile,"[%13.6e,%13.6e]",std::real(zlNew(i1,i2)),std::imag(zlNew(i1,i2)));
	    }
	    fprintf(debugFile,"]\n");
	  }
	}
	  
	fourierTransform.backwardTransform( pzlNew, hupLocal,mc  );

	if( debug  & 1 )
	{
	  RealArray temp(I1,I2,I3);
	  temp=hupLocal(I1,I2,I3);
	  if( raxis==0 )
	    temp.reshape(I2,I3);
	  else if( raxis==1 )
	    temp.reshape(I1,I3);
	  else
	    temp.reshape(I1,I2);
	    
	  ::display(temp,sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),debugFile,"%13.6e ");
	  //  ::display(hupLocal(I1,I2,I3),sPrintF("After fourierTransform.backwardTransform, mc=%d  **NEWER**",mc),"%13.6e ");
	}
	// assert( rside==0 ); // *** FIX ME **
	// assert( raxis==0 ); // *** FIX ME **

	// printF("After fourierTransform.backwardTransform, mc=%d  **NEWER** hupLocal:\n",mc);
	// for( int i2=0; i2<n2; i2++ )
	// {
	//   printF("i2=%d: ",i2);
	//   for( int i1=0; i1<n1; i1++ )
	//   {
	//     printF("%13.6e,",hupLocal(0,i1,i2,mc));
	//   }
	//   printF("\n");
	// }
	
	// end ! useFourierTransformClass

      }  // end for mc 
	
      delete [] pzl;


    }
    else
    {
      OV_ABORT("ERROR: unexpected numberOfDimensions");
    }
    
    
    // --- assign periodic images ---
    fourierTransform.periodicUpdate( hupLocal ); 


  }
  else if( kernelType==cylindrical )
  {
    bccyld(nda1,ndb1,
           *pu, *pf, *ploc, c,radius, dt,numberOfGridPoints1,numberOfFields,numberOfModes1,ns1,
           orderOfTimeStepping, *fold,*phi,
	   *npoles, *alpha, *beta,*amc,*fftsave1,bcinit);
  }
  else
  {
    printF("ERROR: un-implemented kernelType=%i\n",kernelType);
    Overture::abort("error");
  }
  


  if( numberOfDimensions==2 )
  {
    for( int i=hu.getBase(0); i<0; i++ )
    {
      for( int j=0; j<numberOfFields; j++ )
	hu(i,j)=hu(i+numberOfGridPoints1,j);
    }
    for( int i=numberOfGridPoints1; i<=hu.getBound(0); i++ )
    {
      for( int j=0; j<numberOfFields; j++ )
	hu(i,j)=hu(i-numberOfGridPoints1,j);
    }
  }
  else
  {
    // --- assign periodic images ----
    // ** CHECK ME**

    Range M = numberOfFields;

    Range I3= Range(hu.getBase(1),hu.getBound(1));
    Range I2= Range(hu.getBase(0),-1);
    hu(I2,I3,M) = hu(I2+numberOfGridPoints1,I3,M);    // left-ghost = right
    I2 = Range(numberOfGridPoints1,hu.getBound(0));
    hu(I2,I3,M) = hu(I2-numberOfGridPoints1,I3,M);    // right-ghost = left 
    
    I2= Range(hu.getBase(0),hu.getBound(0));
    I3= Range(hu.getBase(1),-1);
    hu(I2,I3,M) = hu(I2,I3+numberOfGridPoints2,M);    // bottom-ghost = top 
    I3 = Range(numberOfGridPoints2,hu.getBound(1));
    hu(I2,I3,M) = hu(I2,I3-numberOfGridPoints2,M);    // top-ghost = bottom
    
    
  }
  
  cpuTime+=getCPU()-time;
  return 0;
}





