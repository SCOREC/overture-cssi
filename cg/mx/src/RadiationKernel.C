#include "RadiationKernel.h"
#include "display.h"

//       subroutine bcperq21(p,f,ploc,c,len,dt,n,m,md,ns,ord,fold,phi,amc,
//      &                    fftsave,bcinit)
// c
// c  this routine uses an adams-moulton formula to compute -
// c  in Fourier variables - 21-pole approximation to the
// c  planar kernel
// c
// c     (d/dt - c beta_j w) phi_j = c alpha_j w^2 phat 
// c
// c     fhat = sum phi_j
// c
// c     w = k*scl , k=1, ... 
// c      
// c  double precision: p(n,m) - m fields to which the operator should be applied
// c
// c  double precision: f(n,m) - the m results
// c
// c  double precision: ploc(n) - workspace 
// c
// c  double precision: c - the wave speed
// c
// c  double precision: len - the period
// c
// c  double precision: dt - the time step
// c
// c  integer: n the number of grid points - most efficient if it has small
// c           prime factors, preferably even
// c
// c  integer: m the number of fields
// c
// c  integer: md the maximum mode used in the bc md < n/2
// c
// c  integer: ns>=2*n+15 
// c
// c  integer: ord - time-stepping order - note that the stability domain for
// c                 Adams-Moulton methods gets small if this is too big
// c 
// c  complex*16: fold(0:ord-2,md,21,m) - stored values for time-stepping
// c
// c  complex*16: phi(md,21,m) - the auxiliary functions computed here
// c
// c  double precision: amc(-1:ord-2) - Adams-Moulton coefficients (computed here)
// c                              use amcof.f
// c
// c  double precision: fftsave(ns) - used by fftpack - link to rffti,rfftf,rfftb
// c
// c  integer bcinit: initialize to zero 
// c

#define bcperq21d EXTERN_C_NAME(bcperq21d)
#define bcperq31d EXTERN_C_NAME(bcperq31d)
#define bcper3dq21 EXTERN_C_NAME(bcper3dq21)
#define bccyld EXTERN_C_NAME(bccyld)
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
}



RadiationKernel::~RadiationKernel()
{
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



// ============================================================================
/// \brief: initialize the RadiationKernel.
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
      // ** 2020/06/29
      amc = new double[numberOfModes1];
      // amc = new double[orderOfTimeStepping+1];
    
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
//  Assign the radiation kernel. 
//
//  This routine assumes that the input values run from u(0),...,u(numberOfGridPoints-1)
//
// 
//   The output hu(i) is defined for all i (with periodic images assigned too)
// ========================================================================================
{
  real time=getCPU();
  
  const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");

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
	// printF("RadiationKernel: period1=%g, period2=%g\n",period1,period2);

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

