// ============================================================================================
//
//   Class to perform Serial and Parallel FFT's in 2D and 3D
//
// ============================================================================================

#include "FourierTransform.h"
#include "ParallelUtility.h"
#include "GhostBoundaryUpdate.h"

#ifdef USE_PPP
  // in parallel you should "setenv USE_FFTMPI useFFTMPI"
  // #define USE_FFTMPI
  // Do this for now 
  // #define USE_FFTPACK
#else
   // in serial use FFTPACK routines 
  #define USE_FFTPACK
#endif

#ifdef USE_FFTMPI
// --- FFTMPI ---
#include "fft2d.h"

#include "fft3d.h"

using namespace FFTMPI_NS;

#endif 

#ifdef USE_FFTPACK
// ---- START FFTPACK ----
#define RFFTI EXTERN_C_NAME(rffti)
#define RFFTF EXTERN_C_NAME(rfftf)
#define RFFTB EXTERN_C_NAME(rfftb)
#define ZFFTI EXTERN_C_NAME(zffti)
#define ZFFTF EXTERN_C_NAME(zfftf)
#define ZFFTB EXTERN_C_NAME(zfftb)

extern "C"
{
  void RFFTI( const int & n, real & fftsave );
  void RFFTF( const int & n1, real & f, real & fftsave );
  void RFFTB( const int & n1, real & f, real & fftsave );

  void ZFFTI( const int & n, real & fftsave );
  void ZFFTF( const int & n1, real & f, real & fftsave );
  void ZFFTB( const int & n1, real & f, real & fftsave );
}

// ---- END FFTPACK ----
#endif


real FourierTransform::cpuTime=0.;



// =====================================================================================================
/// \brief FourierTransform constructor.
// =====================================================================================================
FourierTransform::FourierTransform()
{

  dbase.put<int>("debug")=0;
  dbase.put<aString>("debugFileName")="ft";

  dbase.put<int>("numberOfDimensions");
  dbase.put<int[3]>("dir");
  dbase.put<int[2]>("face");    // hold face on which an FFT will be applied.

  // Physical transforms are the reverse from the default 
  dbase.put<bool>("usePhysicalTransforms")=true;

#ifdef USE_FFTMPI
// --- START FFTMPI ---

  dbase.put<FFT2d*>("pfft2d")=NULL;
  dbase.put<FFT3d*>("pfft3d")=NULL;

  dbase.put<int>("dataSize")=0; // size of data for this processor
  dbase.put<int>("fftsize")=0;  // workspace size needed for this processor
// ---- END FFTMPI ---
#endif
  
}

// =====================================================================================================
/// \brief FourierTransform destructor.
// =====================================================================================================
FourierTransform::~FourierTransform()
{

#ifdef USE_FFTMPI
// --- START FFTMPI ---
  delete dbase.get<FFT2d*>("pfft2d");
  delete dbase.get<FFT3d*>("pfft3d");
// ---- END FFTMPI ---
#endif

}

// =====================================================================================================
/// \brief Return the Index Box that holds the local bounds on the FFT domain
// =====================================================================================================
int FourierTransform::getLocalIndexBox( IndexBox & box ) const
{
  if( !dbase.has_key("fftBox") )
  {
    printF("FourierTransform:getLocalIndexBox: box has not been created yet. Call initialize first.\n");
    return 1;
  }
  else
  {
    IndexBox & fftBox = dbase.get<IndexBox>("fftBox");

    box.setBounds( fftBox.base(0),fftBox.bound(0),
		   fftBox.base(1),fftBox.bound(1),
		   fftBox.base(2),fftBox.bound(2),
		   fftBox.base(3),fftBox.bound(3) );
  }
  return 0;
}


// =====================================================================================================
/// \brief Set the debug flag, non-zero to output info
/// \note Set to 1 for some output and 1+2=3 for more
// =====================================================================================================
int FourierTransform::setDebug( int value )
{
  int & debug = dbase.get<int>("debug");
  debug=value;
  return 0;
}

// =====================================================================================================
/// \brief Set the debug file name 
// =====================================================================================================
int FourierTransform::setDebugFileName( const aString & fileName )
{
  aString & debugFileName = dbase.get<aString>("debugFileName");
  debugFileName = fileName;
  
  return 0;
}

// =====================================================================================================
/// \brief Set the transform option, physical transforms are the reverse of the default DFTs.
// =====================================================================================================
int FourierTransform::usePhysicalTransforms( bool usePhysicalTransforms )
{
  dbase.get<bool>("usePhysicalTransforms")=usePhysicalTransforms;
  return 0;
}


// =====================================================================================================
/// \brief Initialize the Fourier Transform
// =====================================================================================================
int FourierTransform::initialize( realMappedGridFunction & u,
				  const int ndfft /* =1 */ ,
				  const int sidefft /* = -1 */ , const int axisfft /* = -1 */  )
{
  MappedGrid & mg = *u.getMappedGrid();
  
  initialize( u,mg.indexRange(),mg.gridIndexRange(),ndfft,sidefft,axisfft );

  // The GhostBoundaryUpdate class is used to update parallel ghost and periodic ghost 
  if( !dbase.has_key("ghostBoundaryUpdate") )
  {
    dbase.put<GhostBoundaryUpdate>("ghostBoundaryUpdate");
    GhostBoundaryUpdate & ghostBoundaryUpdate = dbase.get<GhostBoundaryUpdate>("ghostBoundaryUpdate");

    const int & debug = dbase.get<int>("debug");
    ghostBoundaryUpdate.setDebug(debug);

    const aString & debugFileName = dbase.get<aString>("debugFileName");
    aString name = "gbu" + debugFileName;
    ghostBoundaryUpdate.setDebugFileName( name );

    int *face = dbase.get<int[2]>("face");
    int sidefft = face[0];
    int axisfft = face[1];

    MappedGrid & mg = *u.getMappedGrid();
    ghostBoundaryUpdate.initialize( u, sidefft,axisfft );
  }

  return 0;
}


  
// =====================================================================================================
/// \brief Initialize the Fourier Transform
/// \param x (input) : P++ array that provides the parallel distribution
/// \param indexRange(0:1,0:*) : index range for active points, excluding periodic images
/// \param ndfft (input) : 1,2 or 3 for 1D, 2D or 3D FFT's 
/// \param sidefft,axisfft (input) : Specify a face to perform the FFT on.
///   For example, if the array x is 2D, and ndfft=1 then the 1D FFT will performed
///       on the face (sidefft,axisfft). E.g. (sidefft,axisfft)=(0,0) will be the left face.
///   
///       (used to determine length of transform, and parallel distribution).
// ============================================================================================
int FourierTransform::initialize( realArray & x,
                                  const IntegerArray & indexRange,
                                  const IntegerArray & gridIndexRange,
				  const int ndfft /* =1 */ ,
				  const int sidefft /* = -1 */ , const int axisfft /* = -1 */  )
{
  const real time=getCPU();

  const int myid = max(0,Communication_Manager::My_Process_Number);
  const int np   = max(1,Communication_Manager::Number_Of_Processors);

  const int & debug = dbase.get<int>("debug");
  int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
  numberOfDimensions=ndfft;

  int *dir = dbase.get<int[3]>("dir");
  int dir1=0, dir2=1, dir3=2;
  dir[0]=dir1;
  dir[1]=dir2;
  dir[2]=dir3;

  int *face = dbase.get<int[2]>("face");
  face[0]=sidefft;
  face[1]=axisfft;
  
  if( debug & 1 )
    printF(">>>> FourierTransform::initialize... np=%d, numberOfDimensions=%d, sidefft=%d, axisfft=%d\n",
	   np,numberOfDimensions,sidefft,axisfft);

  OV_GET_SERIAL_ARRAY(real,x,xLocal);

  if( debug & 2 )
    printf("FT: myid=%d: xLocal=[%d,%d][%d,%d][%d,%d]\n",myid,
	   xLocal.getBase(0),xLocal.getBound(0),
	   xLocal.getBase(1),xLocal.getBound(1),
	   xLocal.getBase(2),xLocal.getBound(2));
   
  
  // ---- Here is the Index Box for the local array of x -----
  IndexBox xBox;
  CopyArray::getLocalArrayBox( myid, x, xBox );

  if( debug & 2 )
    fprintf(stdout,"FT: myid=%i xLocal box: xbox=[%i,%i][%i,%i][%i,%i][%i,%i] (no ghost)\n",
	    myid,
	    xBox.base(0),xBox.bound(0),
	    xBox.base(1),xBox.bound(1),
	    xBox.base(2),xBox.bound(2),
	    xBox.base(3),xBox.bound(3));

  // Here is the Index box for the full array x, or the requested face of x 
  IndexBox fBox;
  IntegerArray box(2,3);
  for( int axis=0; axis<3; axis++ )
  {
    // box(0,axis)=x.getBase(axis);
    // box(1,axis)=max(x.getBase(axis),x.getBound(axis)-1);  // exclude periodic image
    box(0,axis)=indexRange(0,axis);
    box(1,axis)=indexRange(1,axis);  // exclude periodic image on right
    if( axis==axisfft )
    {
      // box(0,axis)=box(sidefft,axis);
      // box(1,axis)=box(sidefft,axis);
      // box(0,axis)=box(sidefft,axis) + sidefft;  // FFT is performed on the boundary 
      box(0,axis)=gridIndexRange(sidefft,axis);  // FFT is performed on the boundary 
      box(1,axis)=box(0,axis);
    }
  }
  fBox.setBounds( box(0,0),box(1,0), box(0,1),box(1,1), box(0,2),box(1,2) );

  if( debug & 2 )
    fprintf(stdout,"FT: myid=%i face-box: fBox=[%i,%i][%i,%i][%i,%i][%i,%i]\n",
	    myid,
	    fBox.base(0),fBox.bound(0),
	    fBox.base(1),fBox.bound(1),
	    fBox.base(2),fBox.bound(2),
	    fBox.base(3),fBox.bound(3));

  if( !dbase.has_key("fftBox") )
    dbase.put<IndexBox>("fftBox");

  IndexBox & fftBox = dbase.get<IndexBox>("fftBox");
  IndexBox::intersect( xBox, fBox, fftBox );

  if( debug & 2 )
    fprintf(stdout,"FT: myid=%i (local) fftBox=[%i,%i][%i,%i][%i,%i][%i,%i]\n",
	    myid,
	    fftBox.base(0),fftBox.bound(0),
	    fftBox.base(1),fftBox.bound(1),
	    fftBox.base(2),fftBox.bound(2),
	    fftBox.base(3),fftBox.bound(3));



  #ifdef USE_FFTPACK
  // **********************************************************************8
  // ---- START FFTPACK ----
  // **********************************************************************8

  bool & useComplexTransforms = dbase.put<bool>("useComplexTransforms");
  useComplexTransforms=true;

  int nx = fBox.bound(0)-fBox.base(0)+1;
  int ny = fBox.bound(1)-fBox.base(1)+1;
  int nz = fBox.bound(2)-fBox.base(2)+1;

  int & n1dim = dbase.put<int>("n1dim");  // length of 1D FFTs      (whether for a 1D, 2D or 3D array)
  int & n2dim = dbase.put<int>("n2dim");  // 2nd length of 2D FFTs  (whether for a 2D or 3D array)
  int & n3dim = dbase.put<int>("n3dim");  // 3rd length of 3D FFTs
  
  n1dim=nx;
  n2dim=ny;
  n3dim=nz;

  // ** CHECK ME***
  if( axisfft>=0  )
  {
    if( numberOfDimensions==1 )
    {
      // --- we are performing a 1D FFT on a 2D/3D array -----
      // NOTE: In 1D, n1dim = nx, ny or nz
      //       In 2D/3D n1dim=nx, n2dim=ny, n3dim=nz
      n1dim = max(nx,ny,nz);
      n2dim=1;
      n3dim=1;
    }
  }

  int ns1 =2*n1dim+15;            // for real FFT
  int ns1c=4*n1dim+15;            // for complex FFT
       
  int ns2 =2*n2dim+15;            // for real FFT
  int ns2c=4*n2dim+15;            // for complex FFT

  int ns3 =2*n3dim+15;            // for real FFT
  int ns3c=4*n3dim+15;            // for complex FFT


  // Workspace for real transforms
  dbase.put<RealArray>("fftsave");
  RealArray & fftsave = dbase.get<RealArray>("fftsave");

      
  // workspace for complex transforms: 
  dbase.put<RealArray>("fftsave1c");
  RealArray & fftsave1c = dbase.get<RealArray>("fftsave1c");
	
  dbase.put<RealArray>("fftsave2c");
  RealArray & fftsave2c = dbase.get<RealArray>("fftsave2c");

  dbase.put<RealArray>("fftsave3c");
  RealArray & fftsave3c = dbase.get<RealArray>("fftsave3c");

  // --- Allocate work space for transforms ---
  if( n1dim>1 )
  {
    fftsave1c.redim(ns1c);
    ZFFTI(n1dim,fftsave1c(0));
  }
  if( n2dim>1 )
  {
    fftsave2c.redim(ns2c);
    ZFFTI(n2dim,fftsave2c(0));
  }
  if( n3dim>1 )
  {
    fftsave3c.redim(ns3c);
    ZFFTI(n3dim,fftsave3c(0));
  }


  // ---- END FFTPACK ----
  #endif 



  #ifdef USE_FFTMPI
  // ***************************************************************************
  // --- START FFTMPI ---
  // ***************************************************************************

  #ifdef FFT_SINGLE
  int precision = 1;
  OV_ABORT("fix me");
  #else
  int precision = 2;  // double precision
  #endif

  // Since fftmpi needs at least 2 point in each direction we sometimes need
  // to perform a 2D transform for a 1D FFT

  int nx = fBox.bound(0)-fBox.base(0)+1;
  int ny = fBox.bound(1)-fBox.base(1)+1;
  int nz = fBox.bound(2)-fBox.base(2)+1;

  // Adjust the FFT box for fftmpi and 1D FFT's on 2D arrays
  int & adjustForOneDimension = dbase.put<int>("adjustForOneDimension");
  adjustForOneDimension=0;  // 0 = no adjustment needed

  IndexBox fftBoxa;
  
  for( int axis=0; axis<3; axis++ )
  {
    box(0,axis)=fftBox.base(axis);
    box(1,axis)=fftBox.bound(axis);
  }
  if( ndfft==1 )
  {
    // --- fftmpi has no 1D FFT ---
    // we use the 2D fft -- and there needs to be at least 2 points in each direction
    if( ny==1 )
    {
      adjustForOneDimension=2; // adjust in y direction
      
      ny=2;
      box(0,1)=0;
      box(1,1)=1;
      
      // if( sidefft==0 )
      //   box(1,1)=box(0,1) + 1;
      // else
      //   box(0,1)=box(1,1) - 1;
    }
    else if( nx==1 )
    {
      adjustForOneDimension=1; // adjust in x direction

      nx=2;
      box(0,0)=0;
      box(1,0)=1;
      // if( sidefft==0 )
      //   box(1,0)=box(0,0) + 1;
      // else
      //   box(0,0)=box(1,0) - 1;
    }
    
  }

  if( ndfft==2 && nz>1 )
  {
    // --- 2D FFT on a 3D array ---
    if( nx==1 )
    {
      // shift box bounds x <- y, y <- z 
      nx=ny; ny=nz; 
      box(0,0)=box(0,1);  box(1,0)=box(1,1);  
      box(0,1)=box(0,2);  box(1,1)=box(1,2);  
      box(0,2)=0; box(1,2)=0;
      
    }
    else if( ny==1 )
    {
      // shift box bounds y <- z 
      ny=nz;
      box(0,1)=box(0,2);  box(1,1)=box(1,2);  
      box(0,2)=0; box(1,2)=0;
    }
    
  }
  

  // if( ndfft==2 && nz>1 )
  // {
  //   if( debug & 1 )
  //   {
  //     printF("FourierTransform: Changing 2D transform to numberOfDimensions=3 since nz>1 and fftmpi needs at"
  // 	     " least 2 points in each direction\n");
  //   }
  //   numberOfDimensions=3;
  //   // dbase.put<int>("nx")=nx;
  //   // dbase.put<int>("ny")=ny;
  //   // dbase.put<int>("nz")=nz;
    
  // }
  

  fftBoxa.setBounds( box(0,0),box(1,0), box(0,1),box(1,1), box(0,2),box(1,2) );

  if( (debug & 2) &&   (ndfft==1 || (ndfft==2 && nz>1) ) )
  {
    fprintf(stdout,"FT: myid=%i (local) ADJUSTED fftBoxa=[%i,%i][%i,%i][%i,%i][%i,%i] (for 1D FFTs)\n",
	    myid,
	    fftBoxa.base(0),fftBoxa.bound(0),
	    fftBoxa.base(1),fftBoxa.bound(1),
	    fftBoxa.base(2),fftBoxa.bound(2),
	    fftBoxa.base(3),fftBoxa.bound(3));
  }
  


  // if( ndfft==1 )
  // {
  //   // fttmpi does not have 1D FFT's -- for now do a 2D FFT
  //   ny=2;   // fftmpi needs at least 2 points
  // }
  

  int & dataSize =  dbase.get<int>("dataSize"); // size of data for this processor
  
  // dataSize = 2*nxLocal*nyLocal*nzLocal;
  int numPointsLocal = (  (fftBoxa.bound(0)-fftBoxa.base(0)+1)
			 *(fftBoxa.bound(1)-fftBoxa.base(1)+1)
			 *(fftBoxa.bound(2)-fftBoxa.base(2)+1) );
  dataSize = 2*numPointsLocal;


  //printF("\n ++++++ BEFORE fft.setup\n");

  int permute=0;

  int & fftsize = dbase.get<int>("fftsize"); // workspace size needed for this processor
  int sendsize,recvsize;
  if( numberOfDimensions<=2 )
  {
    FFT2d *& pfft2d = dbase.get<FFT2d*>("pfft2d");
    pfft2d = new FFT2d(MPI_COMM_WORLD,precision);
    assert( pfft2d != NULL );

    FFT2d & fft = *pfft2d;

    fft.setup(nx,ny,
	      fftBoxa.base(0), fftBoxa.bound(0), fftBoxa.base(1),fftBoxa.bound(1),  // local dims for input
	      fftBoxa.base(0), fftBoxa.bound(0), fftBoxa.base(1),fftBoxa.bound(1),  // local dims for output
	      permute,fftsize,sendsize,recvsize);  

    if( debug & 1 )
      printf(" After fftMPI setup: Using 1D FFT library = %s\n",fft.fft1d);

    // fft.setup(nx,ny,
    // 	      inxlo,inxhi,inylo,inyhi,
    // 	      outxlo,outxhi,outylo,outyhi,
    // 	      permute,fftsize,sendsize,recvsize);  
  }
  else if( numberOfDimensions==3 )
  {
    FFT3d *& pfft3d = dbase.get<FFT3d*>("pfft3d");
    pfft3d = new FFT3d(MPI_COMM_WORLD,precision);
    assert( pfft3d != NULL );

    FFT3d & fft = *pfft3d;

    fft.setup(nx,ny,nz,
	      fftBoxa.base(0), fftBoxa.bound(0), fftBoxa.base(1),fftBoxa.bound(1), fftBoxa.base(2),fftBoxa.bound(2),
	      fftBoxa.base(0), fftBoxa.bound(0), fftBoxa.base(1),fftBoxa.bound(1), fftBoxa.base(2),fftBoxa.bound(2),
	      permute,fftsize,sendsize,recvsize);  

    if( debug & 1 )
      printf(" After fftMPI setup: Using 1D FFT library = %s\n",fft.fft1d);

    // fft.setup(nx,ny,nz,
    // 	      inxlo,inxhi,inylo,inyhi,inzlo,inzhi,
    // 	      outxlo,outxhi,outylo,outyhi,outzlo,outzhi,
    // 	      permute,fftsize,sendsize,recvsize);  
    
  }
  else
  {
    OV_ABORT("FourierTransform::initialize:ERROR: unexpected numberOfDimensions");
  }
  

  
  if( debug & 2 )
    printf("\n ++++++ After fft.setup: myid=%d, fftsize=%d, sendsize=%d, recvsize=%d \n",myid,fftsize,sendsize,recvsize);
  
  //if( useNew )
  //  OV_ABORT("FourierTranform: AFTER stup stop here for now");

/* ---
  fft->remaponly = rflag;

  fft->collective = cflag;
  fft->exchange = eflag;
  fft->packflag = pflag;

  int permute;
  if (mode == 0 || mode == 2) permute = 0;
  else permute = 2;

  // will use fftsize to allocate work buffer
  // ignore sendsize, recvsize b/c let FFT allocate remap buffers internally
  // set timesetup and timetune
  // reset nloop if tuning and user nloop = 0

  int sendsize,recvsize;

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  if (!tuneflag) {
    fft->setup(nx,ny,
               inxlo,inxhi,inylo,inyhi,outxlo,outxhi,outylo,outyhi,
               permute,fftsize,sendsize,recvsize);
  } else {
    int flag = 0;
    if (mode >= 2) flag = 1;
    fft->tune(nx,ny,
              inxlo,inxhi,inylo,inyhi,outxlo,outxhi,outylo,outyhi,
              permute,fftsize,sendsize,recvsize,
              flag,tuneper,tunemax,tuneextra);
    if (nloop == 0) nloop = fft->npertrial;
  }

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  ---- */

  // ---- END FFTMPI ---
  #endif
  
  cpuTime+=getCPU()-time;

  return 0;

} // end initialize



#include <complex>
typedef ::real LocalReal;
typedef ::real OV_real;
typedef std::complex<LocalReal> Complex;


// =========================================================================================
///
/// \brief Perform a complex forward transform
/// \note The forward transform class the 'backward' ZFFTB routine and divides by "n" since
///    this is what is normally expected in Applied Math in order to give the expected Fourier coefficients
///
/// \param u(0:n1-1,0:n2-1,mc) (input)
/// \param mc (input) : component number "mc" in the array u  
/// \param z(0:n1-1,0:n2-1) (output) : complex Fourier transform coefficients
///
// =========================================================================================
int FourierTransform::forwardTransform( const RealArray & u, const int mc, void * zl  )
{
  const LocalReal time = getCPU();

  const int myid = max(0,Communication_Manager::My_Process_Number);
  const int np   = max(1,Communication_Manager::Number_Of_Processors);

  const int & debug = dbase.get<int>("debug");
  const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");

  // The physical forward FFT uses FFTB NOT FFTF
  const bool usePhysicalTransforms = dbase.get<bool>("usePhysicalTransforms");
  const bool scaleForward=usePhysicalTransforms;  // if true divide result by "n", total number of points


  Complex I(0.,1.);

  #ifdef USE_FFTPACK
  // *************************************************************************
  // ---- START FFTPACK -----
  // *************************************************************************

  const bool & useComplexTransforms = dbase.get<bool>("useComplexTransforms");

  if( debug & 1 )
    printF("FourierTransform::forwardTransform: numberOfDimensions=%d, mc=%d, useComplexTransforms=%d\n",
	   numberOfDimensions,mc,int(useComplexTransforms));

  assert( useComplexTransforms );
  
  IndexBox & fftBox = dbase.get<IndexBox>("fftBox");

  const int i1a=fftBox.base(0), i1b=fftBox.bound(0);
  const int i2a=fftBox.base(1), i2b=fftBox.bound(1);
  const int i3a=fftBox.base(2), i3b=fftBox.bound(2);
  
  const int n1=i1b-i1a+1;  
  const int n2=i2b-i2a+1;
  const int n3=i3b-i3a+1;


  if( debug & 2 )
    printf("forwardTransform: myid=%d: [i1a,i1b][i2a,i2b][i3a,i3b]=[%d,%d][%d,%d][%d,%d]\n",myid,i1a,i1b,i2a,i2b,i3a,i3b);
  
  const int & n1dim = dbase.get<int>("n1dim");
  const int & n2dim = dbase.get<int>("n2dim");
  const int & n3dim = dbase.get<int>("n3dim");
  const LocalReal scaleFactor = scaleForward ? 1./(n1dim*n2dim*n3dim) : 1.;

  if( numberOfDimensions==1 )
  {
    // --------- ONE DIMENSIONAL FORWARD TRANSFORM -----------

    RealArray & fftsave = dbase.get<RealArray>("fftsave");
    RealArray & fftsave1c = dbase.get<RealArray>("fftsave1c");

    LocalReal *pzloc = (LocalReal*)zl;       
    //  zlocr(i) = real-part
    //  zloci(i) = imag-part
    #define zlocr(i) pzloc[2*(i)]
    #define zloci(i) pzloc[1+2*(i)]

       
    // printF("u.getBase(0)=%d\n",u.getBase(0));
      
    if( debug & 2 )
      printF("FourierTransform::forwardTransform: FFTPACK: n1dim=%d\n",n1dim);
    
    // ----- real-space data is stored in a local array with local array base and bounds -----
    int i=0;
    for( int i3=i3a; i3<=i3b; i3++ )
    {
      for( int i2=i2a; i2<=i2b; i2++ )
      {
	for( int i1=i1a; i1<=i1b; i1++ )
	{
	  zlocr(i) = u(i1,i2,i3,mc)*scaleFactor;   // real part 
	  zloci(i) = 0.;                           // imag part
          i++;
	}
      }
    }

    if( usePhysicalTransforms )
      ZFFTB( n1dim,pzloc[0],fftsave1c(0) );     // NOTE use "backward" for physical FFT
    else
      ZFFTF( n1dim,pzloc[0],fftsave1c(0) );     // NOTE use "forward" 

	
    if( false )
    {
      printf("FT:fftpack: After Forward: n1dim=%d (zloc holds results for uHat(k) and uHat(-k)=conj(uHat(k))\n",n1dim);

      printf(" zloc=");
      for( int i=0; i<n1dim; i++ ){ printF("[%10.3e,%10.3e],",zlocr(i),zloci(i)); }   // 
      printF("]\n");
	  
    }
    
  }
  else if( numberOfDimensions==2 )
  {
    // --------- TWO DIMENSIONAL FORWARD TRANSFORM -----------
    Complex *pzl = (Complex*)zl;


    RealArray & fftsave = dbase.get<RealArray>("fftsave");
    RealArray & fftsave1c = dbase.get<RealArray>("fftsave1c");
    RealArray & fftsave2c = dbase.get<RealArray>("fftsave2c");
    RealArray & fftsave3c = dbase.get<RealArray>("fftsave3c");

    if( debug & 2 )
      printf("forwardTransform: n1dim=%d, n2dim=%d, n3dim=%d\n",n1dim,n2dim,n3dim);

    #define zl(i1,i2) pzl[(i1)+n1*(i2)]

    LocalReal *pzl1 = new LocalReal [n1dim*2];  // holds complex transforms 
    #define zl1r(i) pzl1[  2*(i)]
    #define zl1i(i) pzl1[1+2*(i)]      

    LocalReal *pzl2 = new LocalReal [n2dim*2];  // holds complex transforms 
    #define zl2r(i) pzl2[  2*(i)]
    #define zl2i(i) pzl2[1+2*(i)]      

    LocalReal *pzl3 = new LocalReal [n3dim*2];  // holds complex transforms 
    #define zl3r(i) pzl3[  2*(i)]
    #define zl3i(i) pzl3[1+2*(i)]      


    // ---- FFT's in direction 1------
    #define zl3(i1,i2,i3) pzl[(i1)+n1*( (i2) + n2*(i3)) ]

    if( i3a==i3b )
    {
      // --- i1,i2 vary ---
      int i3=i3a;
      for( int i2=i2a; i2<=i2b; i2++ )
      {
	int i=0;
	for( int i1=i1a; i1<=i1b; i1++ )
	{
	  zl1r(i)=u(i1,i2,i3,mc)*scaleFactor;  // real part 
	  zl1i(i)=0.;                          // imag part
	  i++;
	}

	if( usePhysicalTransforms )
	  ZFFTB( n1dim,pzl1[0],fftsave1c(0) );     // NOTE use "backward" 
	else
	  ZFFTF( n1dim,pzl1[0],fftsave1c(0) );     // NOTE use "forward" 

	int j2=i2-i2a, j3=i3-i3a;
	for( int j1=0; j1<n1; j1++ )
	{
	  zl3(j1,j2,j3) = zl1r(j1) + zl1i(j1)*I;
	}
      }

    }
    else if( i2a==i2b )
    {
      // --- i1,i3 vary ---
      int i2=i2a;
      for( int i3=i3a; i3<=i3b; i3++ )
      {
	int i=0;
	for( int i1=i1a; i1<=i1b; i1++ )
	{
	  zl1r(i)=u(i1,i2,i3,mc)*scaleFactor;  // real part 
	  zl1i(i)=0.;                          // imag part
	  i++;
	}

	if( usePhysicalTransforms )
	  ZFFTB( n1dim,pzl1[0],fftsave1c(0) );     // NOTE use "backward" 
	else
	  ZFFTF( n1dim,pzl1[0],fftsave1c(0) );   

	int j2=i2-i2a, j3=i3-i3a;
	for( int j1=0; j1<n1; j1++ )
	{
	  zl3(j1,j2,j3) = zl1r(j1) + zl1i(j1)*I;
	}
      }
    }
    else
    {
      // --- i2,i3 vary ---
      int i1=i1a;
      for( int i3=i3a; i3<=i3b; i3++ )
      {
	int i=0;
	for( int i2=i2a; i2<=i2b; i2++ )
	{
	  zl2r(i)=u(i1,i2,i3,mc)*scaleFactor;  // real part 
	  zl2i(i)=0.;                          // imag part
	  i++;
	}

	if( usePhysicalTransforms )
	  ZFFTB( n2dim,pzl2[0],fftsave2c(0) );    // NOTE use "backward" 
	else
	  ZFFTF( n2dim,pzl2[0],fftsave2c(0) );  

	int j1=i1-i1a, j3=i3-i3a;
	for( int j2=0; j2<n2; j2++ )
	{
	  zl3(j1,j2,j3) = zl2r(j2) + zl2i(j2)*I;
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
    // ----- Fill in data -----
    if( i3a==i3b )
    {
      // --- i1,i2 vary ---
      int i3=0;
      for( int i1=0; i1<n1; i1++ )
      {
	for( int i2=0; i2<n2; i2++ )
	{
	  // zl2(i2)=zl(i1,i2);
	  zl2r(i2) = std::real(zl3(i1,i2,i3));
	  zl2i(i2) = std::imag(zl3(i1,i2,i3));
	    
	}

	if( usePhysicalTransforms )
	  ZFFTB( n2dim,pzl2[0],fftsave2c(0) );     // NOTE use "backward" 
	else
	  ZFFTF( n2dim,pzl2[0],fftsave2c(0) );    

	for( int i2=0; i2<n2; i2++ )
	{
	  zl3(i1,i2,i3)=zl2r(i2) + zl2i(i2)*I;
	}
	  
      }

    }
    else if( i2a==i2b )
    {
      // --- i1,i3 vary ---
      int i2=0;
      for( int i1=0; i1<n1; i1++ )
      {
	for( int i3=0; i3<n3; i3++ )
	{
	  zl3r(i3) = std::real(zl3(i1,i2,i3));
	  zl3i(i3) = std::imag(zl3(i1,i2,i3));
	}

	if( usePhysicalTransforms )
	  ZFFTB( n3dim,pzl3[0],fftsave3c(0) );      // NOTE use "backward" 
	else
	  ZFFTF( n3dim,pzl3[0],fftsave3c(0) );    

	for( int i3=0; i3<n3; i3++ )
	{
	  zl3(i1,i2,i3)=zl3r(i3) + zl3i(i3)*I;
	}
      }
    }
    else
    {
      // --- i2,i3 vary ---
      int i1=0;
      for( int i2=0; i2<n2; i2++ )
      {
	for( int i3=0; i3<n3; i3++ )
	{
	  zl3r(i3) = std::real(zl3(i1,i2,i3));
	  zl3i(i3) = std::imag(zl3(i1,i2,i3));
	}

	if( usePhysicalTransforms )
	  ZFFTB( n3dim,pzl3[0],fftsave3c(0) );     // NOTE use "backward" 
	else
	  ZFFTF( n3dim,pzl3[0],fftsave3c(0) );     

	for( int i3=0; i3<n3; i3++ )
	{
	  zl3(i1,i2,i3)=zl3r(i3) + zl3i(i3)*I;
	}
	  
      }
    }
    
    if( false )
    {
      printF("FourierTransform::forwardTransform: zl:\n");
      for( int i2=0; i2<n2; i2++ )
      {
	for( int i1=0; i1<n1; i1++ )
	{
	  printF("[%9.2e,%9.2e] ",std::real(zl(i1,i2)),std::imag(zl(i1,i2)));
	}
	printF("\n");
      }
    }


    delete [] pzl1;
    delete [] pzl2;
    delete [] pzl3;
    
  }  // end 2D transforms
  else if( numberOfDimensions==3 )
  {
    // --------- THREE DIMENSIONAL FORWARD TRANSFORM -----------

    assert( useComplexTransforms );
    
    RealArray & fftsave   = dbase.get<RealArray>("fftsave");
    RealArray & fftsave1c = dbase.get<RealArray>("fftsave1c");
    RealArray & fftsave2c = dbase.get<RealArray>("fftsave2c");
    RealArray & fftsave3c = dbase.get<RealArray>("fftsave3c");

    // -- Complex transform is stored here ---
    Complex *pzl = (Complex*)zl;
    #define zl(i1,i2,i3) pzl[(i1)+n1*( (i2)+n2*(i3) )]

    LocalReal *pzl1 = NULL;
    pzl1 = new LocalReal [n1*2];  // holds complex transforms 
    #define zl1r(i) pzl1[2*(i)]
    #define zl1i(i) pzl1[1+2*(i)]      

    LocalReal *pzl2 = new LocalReal [n2*2];  // holds complex transforms 
    #define zl2r(i) pzl2[2*(i)]
    #define zl2i(i) pzl2[1+2*(i)]      

    LocalReal *pzl3 = new LocalReal [n3*2];  // holds complex transforms 
    #define zl3r(i) pzl3[2*(i)]
    #define zl3i(i) pzl3[1+2*(i)]      

    // ---- FFT's in direction 1------
    for( int i3=0; i3<n3; i3++ )
    {
      for( int i2=0; i2<n2; i2++ )
      {
	// 1D complex transform 
	for( int i1=0; i1<n1; i1++ )
	{
	  // set real and imaginary parts 
	  zl1r(i1)=u(i1,i2,i3,mc)*scaleFactor;    // solution for this component ***** ASSUMES BASE 0 ****
	  zl1i(i1)=0.;
	}

	if( usePhysicalTransforms )
	  ZFFTB( n1,pzl1[0],fftsave1c(0) );       // NOTE use "backward" 
	else
	  ZFFTF( n1,pzl1[0],fftsave1c(0) );    

	for( int i1=0; i1<n1; i1++ )
	{
	  zl(i1,i2,i3) = zl1r(i1) + zl1i(i1)*I;
	}
      } // end for i2 
	
      if( false )
      {
	printF("Stage 1: after initial FFTs, zl:\n");
	for( int i3=0; i3<n3; i3++ )
	{
          printF("----- i3=%d -----\n",i3);
	  for( int i2=0; i2<n2; i2++ )
	  {
	    for( int i1=0; i1<n1; i1++ )
	    {
	      printF("[%9.3e,%9.3e] ",std::real(zl(i1,i2,i3)),std::imag(zl(i1,i2,i3)));
	    }
	    printF("\n");
	  }
	}
      }
    }
    

    // ------ FFT's in direction 2 ------
    for( int i3=0; i3<n3; i3++ )
    {
      for( int i1=0; i1<n1; i1++ )
      {
	for( int i2=0; i2<n2; i2++ )
	{
	  zl2r(i2) = std::real(zl(i1,i2,i3));
	  zl2i(i2) = std::imag(zl(i1,i2,i3));
	}

	if( usePhysicalTransforms )
	  ZFFTB( n2,pzl2[0],fftsave2c(0) );      // NOTE use "backward" 
	else
	  ZFFTF( n2,pzl2[0],fftsave2c(0) );      // NOTE use "backward" 

	for( int i2=0; i2<n2; i2++ )
	{
	  zl(i1,i2,i3)=zl2r(i2) + zl2i(i2)*I;
	}
      }
    }

    // ------ FFT's in direction 3 ------
    for( int i1=0; i1<n1; i1++ )
    {
      for( int i2=0; i2<n2; i2++ )
      {
	for( int i3=0; i3<n3; i3++ )
	{
	  zl3r(i3) = std::real(zl(i1,i2,i3));
	  zl3i(i3) = std::imag(zl(i1,i2,i3));
	}

	if( usePhysicalTransforms )
	  ZFFTB( n3,pzl3[0],fftsave3c(0) );     // NOTE use "backward" 
	else
	  ZFFTF( n3,pzl3[0],fftsave3c(0) );     // NOTE use "backward" 

	for( int i3=0; i3<n3; i3++ )
	{
	  zl(i1,i2,i3)=zl3r(i3) + zl3i(i3)*I;
	}
      }
    }
    
    if( false )
    {
      printF("FourierTransform::forwardTransform: zl:\n");
      for( int i3=0; i3<n3; i3++ )
      {
	printF("----- i3=%d -----\n",i3);
	for( int i2=0; i2<n2; i2++ )
	{
	  for( int i1=0; i1<n1; i1++ )
	  {
	    printF("[%9.2e,%9.2e] ",std::real(zl(i1,i2,i3)),std::imag(zl(i1,i2,i3)));
	  }
	  printF("\n");
	}
      }
    }
    

    delete [] pzl1;
    delete [] pzl2;
    delete [] pzl3;

  }   // end 3D transforms
  else
  {

    OV_ABORT("FourierTransform:ERROR: numberOfDimensions");
    
  }
  
  cpuTime+=getCPU()-time;

  return 0;
  // *********************************************************************************
  // ------ END FFTPACK -----
  // *********************************************************************************
  #endif




  #ifdef USE_FFTMPI
  // *********************************************************************************
  // --- START FFTMPI ---
  // *********************************************************************************

  int *dir = dbase.get<int[3]>("dir");
  const int dir1=dir[0], dir2=dir[1], dir3=dir[2];

  if( debug & 1 )
    printF("FourierTransform::forwardTransform: FFTMPI: numberOfDimensions=%d, mc=%d\n",
	   numberOfDimensions,mc);


  
  // fftsize may be different from dataSize/2 (may include storage for workspace)
  const int & fftsize = dbase.get<int>("fftsize"); // workspace size needed for this processor
  LocalReal *work = new LocalReal [2*fftsize];

  const int & dataSize =  dbase.get<int>("dataSize"); // size of data for this processor

  if( debug & 2)
    printf("forwardTransform: myid=%d: dataSize/2=%d fftsize=%d\n",myid,dataSize/2,fftsize);
  
  if( numberOfDimensions<1 || numberOfDimensions>3 )
  {
    printF("FourierTransform::forwardTransform:ERROR -- invalid numberOfDimensions=%d",numberOfDimensions);
    OV_ABORT("FourierTransform::forwardTransform:ERROR -- invalid numberOfDimensions");
  
  }

  const int & adjustForOneDimension = dbase.get<int>("adjustForOneDimension");
  IndexBox & fftBox = dbase.get<IndexBox>("fftBox");

  int i1a=fftBox.base(dir1), i1b=fftBox.bound(dir1);
  int i2a=fftBox.base(dir2), i2b=fftBox.bound(dir2);
  int i3a=fftBox.base(dir3), i3b=fftBox.bound(dir3);
  
  // if( numberOfDimensions==1 )
  //   i2b=i2a;

  if( debug & 2 )
    printf("forwardTransform: myid=%d: [i1a,i1b][i2a,i2b][i3a,i3b]=[%d,%d][%d,%d][%d,%d]\n",myid,i1a,i1b,i2a,i2b,i3a,i3b);
  

  // ----- Fill in data -----
  int i=0;
  for( int i3=i3a; i3<=i3b; i3++ )
  {
    for( int i2=i2a; i2<=i2b; i2++ )
    {
      for( int i1=i1a; i1<=i1b; i1++ )
      {
	work[i]=u(i1,i2,i3,mc); i++;  // real part 
	work[i]=0.;             i++;  // imag part
	if( adjustForOneDimension==1 )
	{
          // duplicate data for a 1D transform on a 2D array -- x-direction
	  work[i]=u(i1,i2,i3,mc); i++;  // real part 
	  work[i]=0.;             i++;  // imag part
	}
	
      }
    }
  }
  
  // For 1D FFTs we have two points in the y-direction since fftmpi does not support 1D FFTs
  // Just duplicate the data
  // printF(" 1D: i=%d =? dataSize/2=%d, dataSize=%d\n",i,dataSize/2,dataSize);
  if( adjustForOneDimension==2 )
  {
    // duplicate data for a 1D transform on a 2D array -- y-direction
    assert( i== dataSize/2 );
    const int i0=i;
    for( int j=0; j<dataSize/2; j++ ){  work[i]=work[i-i0]; i++; } // 
  }

  // printF(" 1D: i=%d =? dataSize=%d\n",i,dataSize);
  assert( i==dataSize );

  if( numberOfDimensions<=2 )
  {
    FFT2d *& pfft2d = dbase.get<FFT2d*>("pfft2d");
    if( pfft2d==NULL )
    {
      printF("FourierTransform::forward:ERROR: The FFT has not been initialized yet!\n");
      OV_ABORT("ERROR");
    }
    FFT2d & fft = *pfft2d;

    fft.compute(work,work,1);  // forward FFT
  }
  else if( numberOfDimensions==3 )
  {
    FFT3d *& pfft3d = dbase.get<FFT3d*>("pfft3d");
    if( pfft3d==NULL )
    {
      printF("FourierTransform::forward:ERROR: The FFT has not been initialized yet!\n");
      OV_ABORT("ERROR");
    }
    FFT3d & fft = *pfft3d;

    if( usePhysicalTransforms )
      fft.compute(work,work,1);   // forward FFT
    else
      fft.compute(work,work,-1);  // backward

  }
  else
  {
    OV_ABORT("FourierTransform::forwardTransform:ERROR: numberOfDimensions!");
  }
  
  // LocalReal *pzloc = (LocalReal*)zl;
  // for( int i=0; i<dataSize; i++ ){ pzloc[i]=work[i]; } //  copy work to fHat

  Complex *pzl = (Complex*)zl;
  #define zl(i1,i2,i3) pzl[ (i1) + n1*( (i2) + n2*(i3) ) ]

  int n1=i1b-i1a+1;  
  int n2=i2b-i2a+1;
  int n3=i3b-i3a+1;
  // if( numberOfDimensions==1 )  // -- for checking "1D" result ---
  // {
  //   if( n2==1) n2=2;
  // }
  
  i=0;
  for( int i3=0; i3<n3; i3++ )
  {
    for( int i2=0; i2<n2; i2++ )
    {
      for( int i1=0; i1<n1; i1++ )
      {
	zl(i1,i2,i3) = work[i] + work[i+1]*I; i+=2;
	if( adjustForOneDimension==1 )
	  i+=2;
      }
    }
  }
  
      
  if( debug & 2 )
  {
    if( numberOfDimensions==1 && np==1 )
    {
      printF("FourierTransform:AFTER forward complex: zl:\n");
      for( int i3=0; i3<n3; i3++ )
      {
	if( numberOfDimensions==3 )
	  printF("---------------- i3=%d ---------------\n",i3);
	for( int i2=0; i2<n2; i2++ )
	{
	  for( int i1=0; i1<n1; i1++ )
	  {
	    printF("[%9.2e,%9.2e] ",std::real(zl(i1,i2,i3)),std::imag(zl(i1,i2,i3)));
	  }
	  printF("\n");
	}
      }
    }
  }
  


  delete [] work;

  cpuTime+=getCPU()-time;

  return 0;
  

  // *********************************************************************************
  // ---- END FFTMPI ---
  // *********************************************************************************
  #endif


  OV_ABORT("FourierTransform::forwardTransform: UNKNOWN FFT OPTION");

  return 0;
}



// =========================================================================================
///
/// \brief Perform a complex inverse (backward) transform
// 
/// \param z(0:n1-1,0:n2-1) (input) : complex Fourier transform coefficients (NOTE: changed on output)
/// \param u(0:n1-1,0:n2-1,mc) (output) : 
/// \param mc (input) : component number "mc" in the array u  
///
// =========================================================================================
int FourierTransform::backwardTransform( void *zl, RealArray & u, const int mc  )
{
  const LocalReal time = getCPU();

  const int myid = max(0,Communication_Manager::My_Process_Number);
  const int np   = max(1,Communication_Manager::Number_Of_Processors);

  const int & debug = dbase.get<int>("debug");
  const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");

  Complex *pzl = (Complex*)zl;
  
  const bool usePhysicalTransforms = dbase.get<bool>("usePhysicalTransforms");
  const bool scaleBackward=!usePhysicalTransforms;  // if true scale backward transform by "1/n"
  
  Complex I(0.,1.);

  #ifdef USE_FFTPACK
  // *************************************************************************
  // ---- START FFTPACK -----
  // *************************************************************************

  const bool & useComplexTransforms = dbase.get<bool>("useComplexTransforms");
  if( debug & 1 )
    printF("FourierTransform::backwardTransform: numberOfDimensions=%d, mc=%d, useComplexTransforms=%d\n",
	   numberOfDimensions,mc,int(useComplexTransforms));

  IndexBox & fftBox = dbase.get<IndexBox>("fftBox");

  int i1a=fftBox.base(0), i1b=fftBox.bound(0);
  int i2a=fftBox.base(1), i2b=fftBox.bound(1);
  int i3a=fftBox.base(2), i3b=fftBox.bound(2);

  const int n1=i1b-i1a+1;  
  const int n2=i2b-i2a+1;
  const int n3=i3b-i3a+1;

  const int & n1dim = dbase.get<int>("n1dim");
  const int & n2dim = dbase.get<int>("n2dim");
  const int & n3dim = dbase.get<int>("n3dim");

  const LocalReal scaleFactor = scaleBackward ? 1./(n1dim*n2dim*n3dim) : 1.;

  if( numberOfDimensions==1 )
  {
    // ----------------- ONE DIMENSIONAL BACKWARD TRANSFORM -----------------

    RealArray & fftsave = dbase.get<RealArray>("fftsave");
    RealArray & fftsave1c = dbase.get<RealArray>("fftsave1c");


    const int & n1dim = dbase.get<int>("n1dim");
    // Range I1=n1;

    // -- RealArray & ploc = *( (RealArray*)zl); //  p(n1,numComp);
       

    // For complex FFT's
    //  zlocr(i) = real-part
    //  zloci(i) = imag-part
    // -- LocalReal *pzloc = new LocalReal [n1*2];  // holds complex transforms 

    LocalReal *pzloc = (LocalReal*)zl;
    #define zlocr(i) pzloc[2*(i)]
    #define zloci(i) pzloc[1+2*(i)]

    if( usePhysicalTransforms )
      ZFFTF( n1dim,pzloc[0],fftsave1c(0) );        // NOTE: use forward
    else 
      ZFFTB( n1dim,pzloc[0],fftsave1c(0) );        // backward

    if( false )
    {
      //printf("AFter backward: (print only real(zloc))\n ploc=[");
      //for( int i=0; i<n1; i++ ){ printF("%16.9e,",ploc(i)); }   // 
      // printF("]\n");
      printf("FT:AFTER backward:\n");
      printf("zlocr=[");
      for( int i=0; i<n1dim; i++ ){ printF("%16.9e,",zlocr(i)); }   // 
      printF("]\n");
	  
    }

    int i=0;
    for( int i3=i3a; i3<=i3b; i3++ )
    {
      for( int i2=i2a; i2<=i2b; i2++ )
      {
	for( int i1=i1a; i1<=i1b; i1++ )
	{
	  u(i1,i2,i3,mc) = zlocr(i)*scaleFactor; i++;  // real part 
	}
      }
    }

    // delete [] pzloc;

  }
  else if( numberOfDimensions==2 )
  {
    // ----------------- TWO DIMENSIONAL BACKWARD TRANSFORM -----------------

    RealArray & fftsave = dbase.get<RealArray>("fftsave");
    RealArray & fftsave1c = dbase.get<RealArray>("fftsave1c");
    RealArray & fftsave2c = dbase.get<RealArray>("fftsave2c");
    RealArray & fftsave3c = dbase.get<RealArray>("fftsave3c");

    #define zl(i1,i2) pzl[(i1)+n1*(i2)]

    LocalReal *pzl1 = new LocalReal [n1dim*2];  // holds complex transforms 
    #define zl1r(i) pzl1[2*(i)]
    #define zl1i(i) pzl1[1+2*(i)]      

    LocalReal *pzl2 = new LocalReal [n2dim*2];  // holds complex transforms 
    #define zl2r(i) pzl2[2*(i)]
    #define zl2i(i) pzl2[1+2*(i)]      

    LocalReal *pzl3 = new LocalReal [n3dim*2];  // holds complex transforms 
    #define zl3r(i) pzl3[2*(i)]
    #define zl3i(i) pzl3[1+2*(i)]      

    // ------------- INVERSE FFTs in DIRECTION 2  -----------------
    #define zl3(i1,i2,i3) pzl[(i1)+n1*( (i2) + n2*(i3)) ]

    if( i3a==i3b )
    {
      // --- i1,i2 vary ---
      int i3=0;
      for( int i1=0; i1<n1; i1++ )
      {
	for( int i2=0; i2<n2; i2++ )
	{
	  zl2r(i2) = std::real(zl3(i1,i2,i3));
	  zl2i(i2) = std::imag(zl3(i1,i2,i3));
	}

	if( usePhysicalTransforms )
	  ZFFTF(n2dim,pzl2[0],fftsave2c(0));          // NOTE: use forward
	else
	  ZFFTB(n2dim,pzl2[0],fftsave2c(0));      

	for( int i2=0; i2<n2; i2++ )
	{
	  zl3(i1,i2,i3)=zl2r(i2) + zl2i(i2)*I;
	}
      }	

    }
    else if( i2a==i2b )
    {
      // --- i1,i3 vary ---
      int i2=0;
      for( int i1=0; i1<n1; i1++ )
      {
	for( int i3=0; i3<n3; i3++ )
	{
	  zl3r(i3) = std::real(zl3(i1,i2,i3));
	  zl3i(i3) = std::imag(zl3(i1,i2,i3));
	}

	if( usePhysicalTransforms )
	  ZFFTF(n3dim,pzl3[0],fftsave3c(0));          // NOTE: use forward
	else
	  ZFFTB(n3dim,pzl3[0],fftsave3c(0));      

	for( int i3=0; i3<n3; i3++ )
	{
	  zl3(i1,i2,i3)=zl3r(i3) + zl3i(i3)*I;
	}
      }	
    }
    else
    {
      // --- i2,i3 vary ---
      int i1=0;
      for( int i2=0; i2<n2; i2++ )
      {
	for( int i3=0; i3<n3; i3++ )
	{
	  zl3r(i3) = std::real(zl3(i1,i2,i3));
	  zl3i(i3) = std::imag(zl3(i1,i2,i3));
	}

	if( usePhysicalTransforms )
	  ZFFTF(n3dim,pzl3[0],fftsave3c(0));          // NOTE: use forward
	else
	  ZFFTB(n3dim,pzl3[0],fftsave3c(0));          // NOTE: use forward

	for( int i3=0; i3<n3; i3++ )
	{
	  zl3(i1,i2,i3)=zl3r(i3) + zl3i(i3)*I;
	}
      }	
    }
    
    // ------------- INVERSE FFTs in DIRECTION 1  -----------------

    if( i3a==i3b )
    {
      // --- i1,i2 vary ---
      // LocalReal scaleFactor = 1./(n1*n2);
      for( int i2=i2a; i2<=i2b; i2++ )
      {
	int i3=0;
	for( int i1=0; i1<n1; i1++ )
	{
	  zl1r(i1)=std::real(zl3(i1,i2,i3)); 
	  zl1i(i1)=std::imag(zl3(i1,i2,i3));
	}

	// inverse transform: 
	if( usePhysicalTransforms )
	  ZFFTF( n1dim,pzl1[0],fftsave1c(0) );
	else
	  ZFFTB( n1dim,pzl1[0],fftsave1c(0) );

	i3=i3a;
	int i=0;
	for( int i1=i1a; i1<=i1b; i1++ )
	{
	  u(i1,i2,i3,mc)=zl1r(i)*scaleFactor;  // save real part 
	  i++;
	}
      } // end for i2
	
    }
    else if( i2a==i2b )
    {
      // --- i1,i3 vary ---
      // LocalReal scaleFactor = 1./(n1*n3);
      for( int i3=i3a; i3<=i3b; i3++ )
      {
	int i2=0;
	for( int i1=0; i1<n1; i1++ )
	{
	  zl1r(i1)=std::real(zl3(i1,i2,i3)); 
	  zl1i(i1)=std::imag(zl3(i1,i2,i3));
	}

	// inverse transform: 
	if( usePhysicalTransforms )
	  ZFFTF( n1dim,pzl1[0],fftsave1c(0) );
	else
	  ZFFTB( n1dim,pzl1[0],fftsave1c(0) );

	i2=i2a;
	int i=0;
	for( int i1=i1a; i1<=i1b; i1++ )
	{
	  u(i1,i2,i3,mc)=zl1r(i)*scaleFactor;  // save real part 
	  i++;
	}
      } // end for i3
    }
    else
    {
      // --- i2,i3 vary ---
      // LocalReal scaleFactor = 1./(n2*n3);
      for( int i3=i3a; i3<=i3b; i3++ )
      {
	int i1=0;
	for( int i2=0; i2<n2; i2++ )
	{
	  zl2r(i2)=std::real(zl3(i1,i2,i3)); 
	  zl2i(i2)=std::imag(zl3(i1,i2,i3));
	}

	// inverse transform: 
	if( usePhysicalTransforms )
	  ZFFTF( n2dim,pzl2[0],fftsave2c(0) );
	else
	  ZFFTB( n2dim,pzl2[0],fftsave2c(0) );

	i1=i1a;
	int i=0;
	for( int i2=i2a; i2<=i2b; i2++ )
	{
	  u(i1,i2,i3,mc)=zl2r(i)*scaleFactor;  // save real part 
	  i++;
	}
      } // end for i3
    }


    delete [] pzl1;
    delete [] pzl2;
    delete [] pzl3;

  }  // end 2D transform

  else if( numberOfDimensions==3 )
  {
    // ----------------- THREE DIMENSIONAL BACKWARD TRANSFORM -----------------

    RealArray & fftsave   = dbase.get<RealArray>("fftsave");
    RealArray & fftsave1c = dbase.get<RealArray>("fftsave1c");
    RealArray & fftsave2c = dbase.get<RealArray>("fftsave2c");
    RealArray & fftsave3c = dbase.get<RealArray>("fftsave3c");

    #define zl(i1,i2,i3) pzl[(i1)+n1*( (i2) + n2*(i3) )]

    LocalReal *pzl1 = new LocalReal [n1*2];  // holds complex transforms 
    #define zl1r(i) pzl1[2*(i)]
    #define zl1i(i) pzl1[1+2*(i)]      

    LocalReal *pzl2 = new LocalReal [n2*2];  // holds complex transforms 
    #define zl2r(i) pzl2[2*(i)]
    #define zl2i(i) pzl2[1+2*(i)]      

    LocalReal *pzl3 = new LocalReal [n3*2];  // holds complex transforms 
    #define zl3r(i) pzl3[2*(i)]
    #define zl3i(i) pzl3[1+2*(i)]      


    // ------------- INVERSE FFTs in DIRECTION 3  -----------------
    for( int i2=0; i2<n2; i2++ )
    {
      for( int i1=0; i1<n1; i1++ )
      {
	for( int i3=0; i3<n3; i3++ )
	{
	  zl3r(i3) = std::real(zl(i1,i2,i3));
	  zl3i(i3) = std::imag(zl(i1,i2,i3));
	}

        if( usePhysicalTransforms )
	  ZFFTF(n3,pzl3[0],fftsave3c(0));
	else
	  ZFFTB(n3,pzl3[0],fftsave3c(0));

	for( int i3=0; i3<n3; i3++ )
	{
	  zl(i1,i2,i3)=zl3r(i3) + zl3i(i3)*I;
	}
      }
    }
    
    // ------------- INVERSE FFTs in DIRECTION 2  -----------------
    for( int i3=0; i3<n3; i3++ )
    {
      for( int i1=0; i1<n1; i1++ )
      {
	for( int i2=0; i2<n2; i2++ )
	{
	  zl2r(i2) = std::real(zl(i1,i2,i3));
	  zl2i(i2) = std::imag(zl(i1,i2,i3));
	}

        if( usePhysicalTransforms )
  	  ZFFTF(n2,pzl2[0],fftsave2c(0));
        else
  	  ZFFTB(n2,pzl2[0],fftsave2c(0));

	for( int i2=0; i2<n2; i2++ )
	{
	  zl(i1,i2,i3)=zl2r(i2) + zl2i(i2)*I;
	}
      }
    }
    
    // ------------- INVERSE FFTs in DIRECTION 1  -----------------
    // LocalReal scaleFactor = 1./(n1*n2*n3);
    for( int i3=0; i3<n3; i3++ )
    {
      for( int i2=0; i2<n2; i2++ )
      {
	for( int i1=0; i1<n1; i1++ )
	{
	  zl1r(i1) = std::real(zl(i1,i2,i3));
	  zl1i(i1) = std::imag(zl(i1,i2,i3));
	}

        if( usePhysicalTransforms )
  	  ZFFTF(n1,pzl1[0],fftsave1c(0));
        else
  	  ZFFTB(n1,pzl1[0],fftsave1c(0));

	for( int i1=0; i1<n1; i1++ )
	{
	  // zl(i1,i2,i3)=zl1r(i1) + zl1i(i1)*I;
          u(i1,i2,i3,mc) = zl1r(i1)*scaleFactor;  // return real part 
	}
      }
    }
 

    delete [] pzl1;
    delete [] pzl2;
    delete [] pzl3;

  
  }
  else 
  {
    OV_ABORT("ERROR");
  }
  
  cpuTime+=getCPU()-time;

  return 0;
  // *********************************************************************************
  // ------ END FFTPACK -----
  // *********************************************************************************
  #endif


  #ifdef USE_FFTMPI
  // *********************************************************************************
  // --- START FFTMPI --- BACKWARD TRANSFORM ---
  // *********************************************************************************

  const int & adjustForOneDimension = dbase.get<int>("adjustForOneDimension");
  if( debug & 1 )
    printF("FourierTransform::backwardTransform: FFTMPI: numberOfDimensions=%d, mc=%d adjustForOneDimension=%d\n",
	   numberOfDimensions,mc,adjustForOneDimension);

  int *dir = dbase.get<int[3]>("dir");
  const int dir1=dir[0], dir2=dir[1], dir3=dir[2];

  
  // fftsize may be different from dataSize/2 (may include storage for workspace)
  const int & fftsize = dbase.get<int>("fftsize"); // workspace size needed for this processor
  LocalReal *work = new LocalReal [2*fftsize];

  const int & dataSize =  dbase.get<int>("dataSize"); // size of data for this processor

  if( debug & 2 )
    printf("backwardTransform: myid=%d: dataSize/2=%d fftsize=%d\n",myid,dataSize/2,fftsize);
  
  if( numberOfDimensions<1 || numberOfDimensions>3 )
  {
    printF("FourierTransform::backwardTransform:ERROR -- invalid numberOfDimensions=%d",numberOfDimensions);
    OV_ABORT("FourierTransform::backwardTransform:ERROR -- invalid numberOfDimensions");
  
  }

  IndexBox & fftBox = dbase.get<IndexBox>("fftBox");

  int i1a=fftBox.base(dir1), i1b=fftBox.bound(dir1);
  int i2a=fftBox.base(dir2), i2b=fftBox.bound(dir2);
  int i3a=fftBox.base(dir3), i3b=fftBox.bound(dir3);
  
  // if( numberOfDimensions==1 )
  //   i2b=i2a;

  if( debug & 2 )
    printf("backwardTransform: myid=%d: [i1a,i1b][i2a,i2b][i3a,i3b]=[%d,%d][%d,%d][%d,%d]\n",myid,i1a,i1b,i2a,i2b,i3a,i3b);
  
  // ----- Fill in data -----
  #define zl(i1,i2,i3) pzl[ (i1) + n1*( (i2) + n2*(i3) ) ]

  int n1=i1b-i1a+1;  // check me 
  int n2=i2b-i2a+1;
  int n3=i3b-i3a+1;
  int i=0;
  for( int i3=0; i3<n3; i3++ )
  {
    for( int i2=0; i2<n2; i2++ )
    {
      for( int i1=0; i1<n1; i1++ )
      {
	// zl(i1,i2,i3) = work[i] + work[i+1]*I; i+=2;
        work[i] = std::real(zl(i1,i2,i3));  i++;
        work[i] = std::imag(zl(i1,i2,i3));  i++;
	if( adjustForOneDimension==1 )
	{
          // --- add data for a 1D FFT on a 2D array - x-direction
	  work[i] = 0.;  i++;
	  work[i] = 0.;  i++;
	}
	
      }
    }
  }

  
  if( adjustForOneDimension==2 )
  {
    // --- add data for a 1D FFT on a 2D array - y-direction

    // In 1D we have two points in the y-direction since fftmpi does not support 1D FFTs
    // we just duplicate the data -> implies zero out fourier transform variables in 2nd spot
    // printF(" 1D: i=%d =? dataSize/2=%d, dataSize=%d\n",i,dataSize/2,dataSize);
    assert( i== dataSize/2 );
    const int i0=i;
    // for( int j=0; j<dataSize/2; j++ ){  work[i]=work[i-i0]; i++; } // 
    for( int j=0; j<dataSize/2; j++ ){  work[i]=0.; i++; } // 
  }

  // printF(" 1D: i=%d =? dataSize=%d\n",i,dataSize);
  assert( i==dataSize );

// for( int i=0; i<dataSize; i++ ){ work[i]=f[i]; } // copy f to work

  if( numberOfDimensions<=2 )
  {
    FFT2d *& pfft2d = dbase.get<FFT2d*>("pfft2d");
    if( pfft2d==NULL )
    {
      printF("FourierTransform::backward:ERROR: The FFT has not been initialized yet!\n");
      OV_ABORT("ERROR");
    }
    FFT2d & fft = *pfft2d;

    fft.compute(work,work,-1);  // backward FFT
  }
  else if( numberOfDimensions==3 )
  {
    FFT3d *& pfft3d = dbase.get<FFT3d*>("pfft3d");
    if( pfft3d==NULL )
    {
      printF("FourierTransform::backward:ERROR: The FFT has not been initialized yet!\n");
      OV_ABORT("ERROR");
    }
    FFT3d & fft = *pfft3d;

    if( usePhysicalTransforms )
      fft.compute(work,work,-1);  // backward FFT
    else
      fft.compute(work,work,+1);  // forward FFT

  }
  else
  {
    OV_ABORT("FourierTransform::backwardTransform:ERROR: numberOfDimensions!");
  }
  
  // LocalReal *pzloc = (LocalReal*)zl;
  // for( int i=0; i<dataSize; i++ ){ pzloc[i]=work[i]; } //  copy work to fHat

      
  i=0;
  for( int i3=i3a; i3<=i3b; i3++ )
  {
    for( int i2=i2a; i2<=i2b; i2++ )
    {
      for( int i1=i1a; i1<=i1b; i1++ )
      {
	u(i1,i2,i3,mc) = work[i]; i+=2;  // real part 
	if( adjustForOneDimension==1 )
	  i+=2;
      }
    }
  }


  // if( false && np==1 )
  // {
  //   printF("FourierTransform:AFTER backward complex: zl:\n");
  //   for( int i3=0; i3<n3; i3++ )
  //   {
  //     if( numberOfDimensions==3 )
  // 	printF("---------------- i3=%d ---------------\n",i3);
  //     for( int i2=0; i2<n2; i2++ )
  //     {
  // 	for( int i1=0; i1<n1; i1++ )
  // 	{
  // 	  printF("[%9.2e,%9.2e] ",std::real(zl(i1,i2,i3)),std::imag(zl(i1,i2,i3)));
  // 	}
  // 	printF("\n");
  //     }
  //   }
  // }
    


  delete [] work;

  cpuTime+=getCPU()-time;

  return 0;
  

  // *********************************************************************************
  // ---- END FFTMPI ---
  // *********************************************************************************
  #endif



  OV_ABORT("FourierTransform::backwardTransform: finish me");
    
  return 0;
}

// ============================================================================================
// \brief Update periodic boundaries (and parallel ghost).
/// \note This version takes a realMappedGridFunction
///
/// \param C (input) : components to update 
// ============================================================================================
int FourierTransform::periodicUpdate( realMappedGridFunction & u, const Range & C /* = nullRange */ )
{
  const LocalReal time = getCPU();
  
  // The GhostBoundaryUpdate class is used to update parallel ghost and periodic ghost 
  if( !dbase.has_key("ghostBoundaryUpdate") )
  {
    dbase.put<GhostBoundaryUpdate>("ghostBoundaryUpdate");
    GhostBoundaryUpdate & ghostBoundaryUpdate = dbase.get<GhostBoundaryUpdate>("ghostBoundaryUpdate");

    const int & debug = dbase.get<int>("debug");
    ghostBoundaryUpdate.setDebug(debug);

    const aString & debugFileName = dbase.get<aString>("debugFileName");
    aString name = "gbu" + debugFileName;
    ghostBoundaryUpdate.setDebugFileName( name );

    int *face = dbase.get<int[2]>("face");
    int sidefft = face[0];
    int axisfft = face[1];

    MappedGrid & mg = *u.getMappedGrid();
    ghostBoundaryUpdate.initialize( u, sidefft,axisfft );
  }

  GhostBoundaryUpdate & ghostBoundaryUpdate = dbase.get<GhostBoundaryUpdate>("ghostBoundaryUpdate");
  ghostBoundaryUpdate.updateGhostBoundaries( u,C  );

  cpuTime+=getCPU()-time;

  return 0;
}

// ============================================================================================
// \brief Update periodic boundaries (and parallel ghost).
/// \note This version takes a realMappedGridFunction
/// \param C (input) : components to update 
// ============================================================================================
int FourierTransform::periodicUpdate( RealArray & u, const Range & C /* = nullRange */ )
{
   const LocalReal time = getCPU();
   
  // The GhostBoundaryUpdate class is used to update parallel ghost and periodic ghost 
  if( !dbase.has_key("ghostBoundaryUpdate") )
  {
    printF("FourierTransform::periodicUpdate:ERROR: GhostBoundaryUpdate has not been initialized\n"
	   "  You should call the intialize routine that takes a realMappedGridFunction\n");
    OV_ABORT("ERROR");
  }

  GhostBoundaryUpdate & ghostBoundaryUpdate = dbase.get<GhostBoundaryUpdate>("ghostBoundaryUpdate");
  ghostBoundaryUpdate.updateGhostBoundaries( u,C  );

  cpuTime+=getCPU()-time;  

  return 0;
}




// ============================================================================================
// \brief Update periodic (and parallel) boundaries
/// \param gridIndexRange(0:1,0:2) : grid index range 
/// \param dimension(0:1,0:2) : dimension array, indicates which ghost to assign
/// \param x (input) : holds the parallel distribution
/// \param uLocal (input/output) : update periodic boundaries 
// ============================================================================================
int FourierTransform::periodicUpdate( const IntegerArray & gridIndexRange, const IntegerArray & dimension,
				      const IntegerArray & indexRange, const IntegerArray & isPeriodic, 
				      realArray & x, RealArray & u, const int mc  )
{

  if( true )
  {
    // *new* way 

    // The GhostBoundaryUpdate class is used to update parallel ghost and periodic ghost 
    if( !dbase.has_key("ghostBoundaryUpdate") )
    {
      dbase.put<GhostBoundaryUpdate>("ghostBoundaryUpdate");
      GhostBoundaryUpdate & ghostBoundaryUpdate = dbase.get<GhostBoundaryUpdate>("ghostBoundaryUpdate");

      const int & debug = dbase.get<int>("debug");
      ghostBoundaryUpdate.setDebug(debug);

      const aString & debugFileName = dbase.get<aString>("debugFileName");
      aString name = "gbu" + debugFileName;
      ghostBoundaryUpdate.setDebugFileName( name );

      int *face = dbase.get<int[2]>("face");
      int sidefft = face[0];
      int axisfft = face[1];

     const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
     ghostBoundaryUpdate.initialize( x, numberOfDimensions, gridIndexRange,dimension,indexRange,isPeriodic, sidefft,axisfft );

    }

    GhostBoundaryUpdate & ghostBoundaryUpdate = dbase.get<GhostBoundaryUpdate>("ghostBoundaryUpdate");
    ghostBoundaryUpdate.updateGhostBoundaries( u,mc );
  }

  else
  {
    // **OLD WAY ***


    // Note: see periodicUpdate in gf/mappedGridFunction.C 
    // uses ParallelUtility::copy(u,I, u,Ip,nd );
  
    // copyArray( realSerialArray & dest, int destProcessor, realSerialArray & src, int srcProcessor )

    const int myid = max(0,Communication_Manager::My_Process_Number);
    const int np   = max(1,Communication_Manager::Number_Of_Processors);

    const int & debug = dbase.get<int>("debug");
    const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");

    if( debug & 1 )
      printF("FourierTransform::periodicUpdate: numberOfDimensions=%d, mc=%d\n",
	     numberOfDimensions,mc);

    IndexBox & fftBox = dbase.get<IndexBox>("fftBox");

    const int i1a=fftBox.base(0), i1b=fftBox.bound(0);
    const int i2a=fftBox.base(1), i2b=fftBox.bound(1);
    const int i3a=fftBox.base(2), i3b=fftBox.bound(2);
  
    const int n1=i1b-i1a+1;  
    const int n2=i2b-i2a+1;
    const int n3=i3b-i3a+1;

    if( debug & 2 )
      printf("periodicUpdate: myid=%d: [i1a,i1b][i2a,i2b][i3a,i3b]=[%d,%d][%d,%d][%d,%d]\n",myid,i1a,i1b,i2a,i2b,i3a,i3b);

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];

    IntegerArray gid(2,3);
    // gid(0,0) = i1a; gid(1,0) = i1b+1;
    // gid(0,1) = i2a; gid(1,1) = i2b+1;
    // gid(0,2) = i3a; gid(1,2) = i3b+1;
  
    gid = gridIndexRange;
  
    if( np==1 )
    {

      // --------- SERIAL PERIODIC UPDATE ------
      int side=0, axis=0;
      if( n1>1 )
      {
	getBoundaryIndex( gid,side,axis,I1,I2,I3 );
	// right = left 
	u(I1+n1,I2,I3,mc) = u(I1,I2,I3,mc);
      }
  
      if( n2>1 )
      {
	side=0, axis=1;
	getBoundaryIndex( gid,side,axis,I1,I2,I3 );
	// top = bottom
	// Index J1 = Range(I1.getBase(),I1.getBound()+1); // add a point 
	u(I1,I2+n2,I3,mc) = u(I1,I2,I3,mc);
      }
  
      if( numberOfDimensions==3 && n3>1 )
      {
	side=0, axis=2;
	getBoundaryIndex( gid,side,axis,I1,I2,I3 );
	// back = front
	// Index J2 = Range(I2.getBase(),I2.getBound()+1); // add a point 
	u(I1,I2,I3+n3,mc) = u(I1,I2,I3,mc);
      }
  
    }
    else
    {
      // --------- PARALLEL PERIODIC UPDATE ------

      printF("PARALLEL PERIODIC UPDATE -- FINISH ME ...\n");
      // OV_GET_SERIAL_ARRAY(real,u,uLocal);
  
      Index Jv[3], &J1=Jv[0], &J2=Jv[1], &J3=Jv[2];
      Index Kv[3], &K1=Kv[0], &K2=Kv[1], &K3=Kv[2];

      int side=0, axis=0;
      if( n1>1 )
      {
	getBoundaryIndex( gid,side,axis,I1,I2,I3 );
	// right = left 
	// u(I1+n1,I2,I3,mc) = u(I1,I2,I3,mc);

	// J1=I1; J2=I2; J3=I3;
	// int includeGhost=1;
	// bool ok1 = ParallelUtility::getLocalArrayBounds(u,uLocal,J1,J2,J3,includeGhost);    

	// K1=I1; K2=I2; K3=I3;
	// int includeGhost=1;
	// bool ok2 = ParallelUtility::getLocalArrayBounds(u,uLocal,K1,K2,K3,includeGhost);    
      
	// realSerialArray src, dest;
	// if( ok1 )
	// 	src.redim(J1,J2,J3);
	// if( ok2 )
	// 	dest.redim(K1,K2,K3);

	// realSerialArray dest(J1,J2,J3);
	// int destProcessor=0, srcProcessor=0;
	// src = u(I1,I2,I3,mc);
	// copyArray( dest, destProcessor, src, srcProcessor );
	// u(I1+n1,I2,I3,mc)=dest;
      
      }


    }
  
  }
  
  return 0;
}

