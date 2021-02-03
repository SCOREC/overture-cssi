// This file automatically generated from trad.bC with bpp.
// ---------------------------------------------------------------------------
//
//  Test the "exact" radiation boundary conditions
//
// Solve the second order wave equation
//       u_tt + c^2 ( u_xx + u_yy ) =0 
// with radiation boundary conditions
//
//   
// ---------------------------------------------------------------------------

#include "Overture.h"
#include "Ogshow.h"  
#include "MappedGridOperators.h"
#include "SquareMapping.h"
#include "BoxMapping.h"
#include "PlotStuff.h"
#include "ShowFileReader.h"
#include "display.h"
#include "ParallelUtility.h"
#include "gridFunctionNorms.h"

#include "RadiationKernel.h"

#include "RadiationBoundaryCondition.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++)                                       for(i2=I2Base; i2<=I2Bound; i2++)                                     for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3)                                       I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++)                                       for(i2=I2Base; i2<=I2Bound; i2++)                                     for(i1=I1Base; i1<=I1Bound; i1++)


#define exmax exmax_
#define wpulse wpulse_
#define wdpulse wdpulse_
#define radEval radeval_

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
                          const int & gridIndexRange,real& u1,real&  u, real & xy, const int &boundaryCondition, 
                          const int &md1a,const int &md1b,
                          real& hu,real& hux,real& huxx,real& huxxx, const int &ipar,const real& rpar, const int &ierr );

}

static bool useNewBC=true; // false;

static int rside=0, raxis=0;  // apply radiation BC to this side

static real dx[3];
static int orderOfAccuracy=2;


static double wt,wx,wy;
static const int nsources=1;
static double xs[nsources], ys[nsources], tau[nsources], var[nsources], amp[nsources];
      
static double period= 1.;  // period in y
static double gisTime;

// static  RadiationKernel radiationKernelv;

// static  RadiationKernel radiationKernel;
// static  RadiationKernel radiationKernel1,radiationKernel2,radiationKernel3;

static const int numberOfTimeLevels=10;
static RealArray uxb,uxxb,uxxxb, hub, huxb;



//==================================================================================================
// Evaluate Tom Hagstom's exact solution defined as an integral of Guassian sources
// 
//
//==================================================================================================



// ramp1 has 1-derivative zero at 0 and 1
#define ramp1(t)  (t)*(t)*( -(t)/3.+.5 )*6.

// ramp(0)=0  ramp(1)=1 -- three derivatives zero at 0 and 1
#define ramp(t)    ( -84*pow(t,5.)+35*pow(t,4.)-20*pow(t,7.)+70*pow(t,6.) )

// fourth order dissipation 2D:
#define FD4_2D(u,i1,i2,i3) (    -( u(i1-2,i2,i3)+u(i1+2,i2,i3)+u(i1,i2-2,i3)+u(i1,i2+2,i3) )   +4.*( u(i1-1,i2,i3)+u(i1+1,i2,i3)+u(i1,i2-1,i3)+u(i1,i2+1,i3) ) -12.*u(i1,i2,i3) )

// fourth order dissipation 3D:
#define FD4_3D(u,i1,i2,i3) (    -( u(i1-2,i2,i3)+u(i1+2,i2,i3)+u(i1,i2-2,i3)+u(i1,i2+2,i3)+u(i1,i2,i3-2)+u(i1,i2,i3+2) )   +4.*( u(i1-1,i2,i3)+u(i1+1,i2,i3)+u(i1,i2-1,i3)+u(i1,i2+1,i3)+u(i1,i2,i3-1)+u(i1,i2,i3+1) ) -18.*u(i1,i2,i3) )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  


  // Pulse parameters:
static real alpha=30.; // 200.;
static real pulsePow=1.; // 10.; // 20
static real a0=1.;  // 3.;
// static real xPulse=-1.;
static real xPulse=0.;
static real yPulse=0.;
static real zPulse=0.;
static real ampPulse=5.;

enum InitialConditionOptionEnum
{
    smoothPulse,
    pulse,
    gaussianIntegral
};


// Assign the initial conditions
void
getInitialConditions( InitialConditionOptionEnum icOption, realMappedGridFunction *u, real t, real dt,
                                            real c )
{
    printf("get initial conditions: icOption=%d (0=smoothPulse,1=pulse,2=gaussianIntegral) ampPulse=%g, alpha=%g, a0=%g, pulsePow=%g\n",
       	 icOption,ampPulse,alpha,a0,pulsePow);
    
    MappedGrid & mg = *( u[0].getMappedGrid() );

  // Pulse parameters:
    

// define U0(x,y,t) exp( - alpha*( SQR((x)-(xPulse-c*dt)) + SQR((y)-yPulse) ) )
// #define U0(x,y,t) exp( - alpha*( SQR((x)-(xPulse+c*(t))) ) )
// define U0(x,y,t) exp( - alpha*( pow( a0*( (x)-(xPulse+c*(t)) ),20.) ) )
#define U0(x,y,t) ampPulse*exp( - alpha*( pow( a0*( (x)*(x) + (y)*(y) ),pulsePow) ) )
#define U03d(x,y,z,t) ampPulse*exp( - alpha*( pow( a0*( (x-xPulse)*(x-xPulse) + (y-yPulse)*(y-yPulse) + (z-zPulse)*(z-zPulse) ),pulsePow) ) )

#define U0(x,y,t) ampPulse*exp( - alpha*( pow( a0*( (x-xPulse)*(x-xPulse) + (y-yPulse)*(y-yPulse) ),pulsePow) ) )

    mg.update(MappedGrid::THEvertex);  // build the array of vertices
    realArray & vertex = mg.vertex();

    OV_GET_SERIAL_ARRAY(real,vertex,xLocal);
    OV_GET_SERIAL_ARRAY(real,u[0],u0Local);
    OV_GET_SERIAL_ARRAY(real,u[1],u1Local);
    


    Index I1,I2,I3;
    int i1,i2,i3;
    getIndex(mg.dimension(),I1,I2,I3);
    int includeGhost=1;
    bool ok = ParallelUtility::getLocalArrayBounds(u[0],u0Local,I1,I2,I3,includeGhost);     

    if( mg.numberOfDimensions()==2 )
    {
        if( icOption==smoothPulse )
        {
            u1Local(I1,I2,I3)=U0(xLocal(I1,I2,I3,0),xLocal(I1,I2,I3,1),-dt);
            u0Local(I1,I2,I3)=U0(xLocal(I1,I2,I3,0),xLocal(I1,I2,I3,1),0.);
        }
        else
        {
            real t=0.;
            real w; 
            FOR_3(i1,i2,i3,I1,I2,I3) 
            {
      	real x=xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1);
      	
            {
      //   const int nsources=1;
      //   double xs[nsources], ys[nsources], tau[nsources], var[nsources], amp[nsources];
      //   xs[0]=0.;
      //   ys[0]=1.e-8*1./3.;  // should not be on a grid point
      //   tau[0]=-.95;
      //   var[0]=30.;
      //   amp[0]=1.;
      //   double period= 1.;  // period in y
                gisTime=t;
                wpulse(w,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,gisTime);
            }
      	u0Local(i1,i2,i3)=w; 

            {
      //   const int nsources=1;
      //   double xs[nsources], ys[nsources], tau[nsources], var[nsources], amp[nsources];
      //   xs[0]=0.;
      //   ys[0]=1.e-8*1./3.;  // should not be on a grid point
      //   tau[0]=-.95;
      //   var[0]=30.;
      //   amp[0]=1.;
      //   double period= 1.;  // period in y
                gisTime=t-dt;
                wpulse(w,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,gisTime);
            }
      	u1Local(i1,i2,i3)=w;
            }
        }
    }
    else
    {
        u1Local(I1,I2,I3)=U03d(xLocal(I1,I2,I3,0),xLocal(I1,I2,I3,1),xLocal(I1,I2,I3,2),-dt);
        u0Local(I1,I2,I3)=U03d(xLocal(I1,I2,I3,0),xLocal(I1,I2,I3,1),xLocal(I1,I2,I3,2),0.);
    }
    
    printf("done initial conditions\n");

}


static ShowFileReader *referenceShowFileReader=NULL;
static bool compareToReferenceShowFile=false;

static void 
getErrors( InitialConditionOptionEnum icOption, 
                      realMappedGridFunction & u, realMappedGridFunction & err, 
                      real t )
// ==============================================================================================================
// Compute the error compared to a reference solution saved in a show file
// ==============================================================================================================
{
    
    MappedGrid & mg = *( u.getMappedGrid() );
    
    Index I1,I2,I3;
    int i1,i2,i3;

    mg.update(MappedGrid::THEvertex);  // build the array of vertices
    realArray & vertex = mg.vertex();
    OV_GET_SERIAL_ARRAY(real,vertex,xLocal);
    OV_GET_SERIAL_ARRAY(real,u,uLocal);
    OV_GET_SERIAL_ARRAY(real,err,errLocal);

    getIndex(mg.dimension(),I1,I2,I3);
    int includeGhost=1;
    bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3,includeGhost);     

    if( mg.numberOfDimensions()==2 )
    {
        if( icOption==gaussianIntegral )
        {
            real w;
            FOR_3(i1,i2,i3,I1,I2,I3) 
            {
      	real x=xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1);
      	
            {
      //   const int nsources=1;
      //   double xs[nsources], ys[nsources], tau[nsources], var[nsources], amp[nsources];
      //   xs[0]=0.;
      //   ys[0]=1.e-8*1./3.;  // should not be on a grid point
      //   tau[0]=-.95;
      //   var[0]=30.;
      //   amp[0]=1.;
      //   double period= 1.;  // period in y
                gisTime=t;
                wpulse(w,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,gisTime);
            }
      	errLocal(i1,i2,i3)=uLocal(i1,i2,i3)-w;

            }
        }
    }
    else
    {
    }
    

}

static void 
setBoundaryConditions( InitialConditionOptionEnum icOption, 
                                              realMappedGridFunction & u, real t, real dt,
                                              realMappedGridFunction & u2, int currentTimeLevel )
// ==============================================================================================================
//  Assign DIRICHLET boundary conditions
//
//  u2: u(t-dt)
// ==============================================================================================================
{
    
    MappedGrid & mg = *( u.getMappedGrid() );
    
    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    int i1,i2,i3;

    mg.update(MappedGrid::THEvertex);  // build the array of vertices
    realArray & vertex = mg.vertex();

    OV_GET_SERIAL_ARRAY(real,vertex,xLocal);
    OV_GET_SERIAL_ARRAY(real,u,uLocal);

    const int extra=2;
    
    for( int side=0; side<=1; side++ )
    {
        for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
        {
            if( mg.boundaryCondition(side,axis)>0 )
            {
      	getBoundaryIndex(mg.gridIndexRange(),side,axis,I1,I2,I3,extra);
      	
                int ex[2]={0,0}; // 
                ex[side]=2;
      	
                Iv[axis]=Range(Iv[axis].getBase()-ex[0],Iv[axis].getBound()+ex[1]);  // assign ghost lines too
      	
      	int includeGhost=1;
      	bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3,includeGhost);     

      	if( icOption==gaussianIntegral && !(side==rside && axis==raxis) ) 
      	{
        	  real w; 
        	  FOR_3(i1,i2,i3,I1,I2,I3) 
        	  {
          	    real x=xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1);
      	
                    {
          //   const int nsources=1;
          //   double xs[nsources], ys[nsources], tau[nsources], var[nsources], amp[nsources];
          //   xs[0]=0.;
          //   ys[0]=1.e-8*1./3.;  // should not be on a grid point
          //   tau[0]=-.95;
          //   var[0]=30.;
          //   amp[0]=1.;
          //   double period= 1.;  // period in y
                        gisTime=t;
                        wpulse(w,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,gisTime);
                    }
          	    uLocal(i1,i2,i3)=w;

        	  }
      	}

            } // end if bc > 0 
        }
    }
    

}


// fourth order dissipation 2D:
#define FD4_2D(u,i1,i2,i3) (    -( u(i1-2,i2,i3)+u(i1+2,i2,i3)+u(i1,i2-2,i3)+u(i1,i2+2,i3) )   +4.*( u(i1-1,i2,i3)+u(i1+1,i2,i3)+u(i1,i2-1,i3)+u(i1,i2+1,i3) ) -12.*u(i1,i2,i3) )



// ================================================================================
//   Test the radiation kernel -- serial versus parallel versions
// ================================================================================
real testParallelKernel( int numSteps, int numComp, int numberOfDimensions, int numberOfPoles, int orderOfTimeStepping,
                   			 int nx, int ny, int nz, int rside, int raxis, int debug, int check )
{
    const int myid = max(0,Communication_Manager::My_Process_Number);
    const int np   = max(1,Communication_Manager::Number_Of_Processors);

    real maxDiff=0.;

  // -------- create a SQUARE or BOX grid ------------------
    real xa=-2., xb=2., ya=0., yb=1., za=0., zb=1.;

    Mapping *pmap=NULL;

    if( numberOfDimensions==2 )
        pmap = new SquareMapping(xa,xb,ya,yb);
    else
        pmap = new BoxMapping(xa,xb,ya,yb,za,zb);
    pmap->incrementReferenceCount();
    
    Mapping & map = *pmap;
    
    int numberOfGhostPoints=2; // orderOfAccuracy/2+1;

  // printF("Grid points: nx=%d, ny=%d, nz=%d\n",nx,ny,nz);
    
    int nxv[3]={nx,ny,nz};  //
    int side,axis;
    for( int axis=0; axis<map.getDomainDimension(); axis++ )
    {
        map.setGridDimensions(axis,nxv[axis]);
        if( axis !=raxis )
            map.setIsPeriodic(axis,Mapping::derivativePeriodic);
    }
    
    MappedGrid mg(map);
    for( int axis=0; axis<map.getDomainDimension(); axis++ )
    {
        for( int side=Start; side<=End; side++ )
            mg.setNumberOfGhostPoints(side,axis,numberOfGhostPoints);
    }
    mg.update(MappedGrid::THEmask | MappedGrid::THEcenter );
    bool isRectangular = mg.isRectangular();

    mg.getDeltaX(dx);
  // printF(" *** nx=%i, ny=%i, nz=%d \n",nxv[0],nxv[1]);

    
//  int numberOfFields=1;

    int numberOfGridPoints1 = raxis==0 ? ny-1 : raxis==1 ? nx-1 : nx-1;
    int numberOfGridPoints2 = raxis==0 ? nz-1 : raxis==1 ? nz-1 : ny-1;

    const int maxModes1=numberOfGridPoints1/2-1;  // this is the max allowed
    int numberOfModes1=min(maxModes1,max(8,numberOfGridPoints1/4));  // this is the max allowed

    const int maxModes2=numberOfGridPoints2/2-1;  // this is the max allowed
    int numberOfModes2=min(maxModes2,max(8,numberOfGridPoints2/4));  // this is the max allowed

// TEST
  // numberOfModes1-=2;
  // numberOfModes2-=2;
    

    real c=1.;

    int currentTimeLevel=0;

    Range N=numComp;
    Range all;
    realMappedGridFunction up(mg,all,all,all,N);  // for parallel transforms
    OV_GET_SERIAL_ARRAY(real,up,upLocal);
    upLocal=0.;

    RealArray hpLocal;
    hpLocal.redim(upLocal);
        
    
    Index Jv[3], &J1=Jv[0], &J2=Jv[1], &J3=Jv[2];
    getIndex(mg.dimension(),J1,J2,J3);
    Jv[raxis]=mg.gridIndexRange(rside,raxis);

  // real period1=yb-ya;
  // real period2=zb-za;
        
    if( (debug & 2 ) || check==0 )
    {
        printF(" ---------------------------------------------------------------\n"
         	   " ------------------- TEST RADIATION KERNEL ---------------------\n"
         	   " ------- Compare serial to parallel : np=%d     ----------------\n"
         	   " ------- numComp=%d numSteps=%d                   --------------\n"
         	   " ------- numberOfPoles=%d, orderOfTimeStepping=%d --------------\n"
         	   " ------- numberOfDimensions=%d, nx=%d, ny=%d, nz=%d ------------\n"
         	   " ------- numberOfModes1=%d, numberOfModes2=%d     --------------\n"
         	   " ------- rside=%d, raxis=%d                       --------------\n"
         	   " ------- Jv=[%d,%d][%d,%d][%d,%d]  (boundary index) ------------\n"
         	   " ---------------------------------------------------------------\n",
         	   np,numComp,numSteps,numberOfPoles,orderOfTimeStepping,numberOfDimensions,nx,ny,nz,
         	   numberOfModes1,numberOfModes2,rside,raxis,
         	   J1.getBase(),J1.getBound(),J2.getBase(),J2.getBound(),J3.getBase(),J3.getBound()
            );
    }
    

    RadiationKernel radiationKernel;
    radiationKernel.setUseFourierTransformClass(true);

    radiationKernel.setDebug( debug );
        
  // --- new init using a mappedGridFunction ---
    radiationKernel.initialize(up,
                       			     rside,raxis,
                       			     numComp, 
                       			     numberOfModes1, numberOfModes2,
                       			     c, orderOfTimeStepping, numberOfPoles );


    if( numberOfDimensions==2 || numberOfDimensions==3 )
    {
    // --- TEST KERNEL : parallel versus serial ---
    // For 3D : 2D solution on boundary at i1=0:
    // add ghost to arrays
        int numGhost=2;
        Range R2=Range(0,0), R3=Range(0,0);
        if( numberOfDimensions>=2 )
            R2 = Range(-numGhost,nx+numGhost);
        if( numberOfDimensions>=3 )
            R3 = Range(-numGhost,ny+numGhost);

    //     RealArray u3(I2,I3), hu3(I2,I3), hup3(I2,I3);
        RealArray u3, hu3, hup3;
        if( numberOfDimensions==2 )
        {
            if( raxis==0 )
            {
      	u3.redim(J2,N); hu3.redim(J2,N); hup3.redim(J2,N);
            }
            else
            {
      	u3.redim(J1,N); hu3.redim(J1,N); hup3.redim(J1,N);
            }
      	
        }
        else
        {
      // 	u3.redim(R2,R3); hu3.redim(R2,R3); hup3.redim(R2,R3);
            if( raxis==0 )
            {
      	u3.redim(J2,J3,N); hu3.redim(J2,J3,N); hup3.redim(J2,J3,N);
            }
            else if( raxis==1 )
            {
      	u3.redim(J1,J3,N); hu3.redim(J1,J3,N); hup3.redim(J1,J3,N);
            }
            else
            {
      	u3.redim(J1,J2,N); hu3.redim(J1,J2,N); hup3.redim(J1,J2,N);
            }
      	
        }
            
        #define MYFUNC(x,y,z,n) (1.+n) + .3*(1.+n)*(x) + .25*(1.+n)*(y) + .5*(1.+n)*(z) 	                        + .1*(1.+n)*(x)*(y)  + .3*(1.+n)*(y)*(z) + .4*(1.+n)*(x)*(z)

        u3=1.;
        int j1,j2,j3;
        FOR_3D(j1,j2,j3,J1,J2,J3) 
        {
            real x = xa + (xb-xa)*j1/ny;
            real y = ya + (yb-ya)*j2/ny;
            real z = za + (zb-za)*j3/nz;
            for( int n=0; n<numComp; n++ )
            {
      	if( numberOfDimensions==2 )
      	{
        	  if( raxis==0 )
        	  {
          	    u3(j2,n) = MYFUNC(x,y,z,n);
        	  }
        	  else if( raxis==1 )
        	  {
          	    u3(j1,n) = MYFUNC(x,y,z,n);
        	  }
        	  else
        	  {
          	    u3(j1,n) = MYFUNC(x,y,z,n);
        	  }
      	}
      	else
      	{
        	  if( raxis==0 )
        	  {
          	    u3(j2,j3,n) = MYFUNC(x,y,z,n);
        	  }
        	  else if( raxis==1 )
        	  {
          	    u3(j1,j3,n) = MYFUNC(x,y,z,n);
        	  }
        	  else
        	  {
          	    u3(j1,j2,n) = MYFUNC(x,y,z,n);
        	  }
      	}
      	
            }
        }

    // int j1,j2,j3;
        int includeGhost=1;
        bool ok = ParallelUtility::getLocalArrayBounds(up,upLocal,J1,J2,J3,includeGhost);
        if( ok )
        {
            FOR_3(j1,j2,j3,J1,J2,J3)
            {
      	real x = xa + (xb-xa)*j1/ny;
      	real y = ya + (yb-ya)*j2/ny;
      	real z = za + (zb-za)*j3/nz;
      	for( int n=0; n<numComp; n++ )
      	{
        	  upLocal(j1,j2,j3,n) = MYFUNC(x,y,z,n);
      	}
      	
            }
        }
        
    // ::display(u3,"u3 serial (start)");
    // ::display(upLocal(J1,J2,J3),"upLocal(J1,J2,J3) parallel (start)");
            
        real dt=.01;
            

        for( int step=0; step<numSteps; step++ )
        {
            
      // --- old way ---
            radiationKernel.evaluateKernel( dt, u3, hu3 );

            if( debug & 2 && myid==0  )
      	::display(hu3,"hu3 serial","%14.4e ");
        
      // --- new way ---
            if( false )
            {
      	radiationKernel.evaluateKernelParallelV1( dt, u3, hup3, upLocal, hpLocal );
            }
            else
            {
      	radiationKernel.evaluateKernelParallel( dt, upLocal, hpLocal );
            }
            
        }
    //% end for step 

        if( debug & 2  )
        {
      // Index Ib1,Ib2,Ib3;
      // getBoundaryIndex(mg.gridIndexRange(),rside,raxis,Ib1,Ib2,Ib3);
      // ::display(hpLocal(Ib1,Ib2,Ib3),"hpLocal parallel","%14.4e ");
            if( np==1 )
            {
      	RealArray temp(J1,J2,J3,N);
      	temp = hpLocal(J1,J2,J3,N);
      	if( raxis==0 )
        	  temp.reshape(J2,J3,N);
      	else if( raxis==1 )
        	  temp.reshape(J1,J3,N);
      	else
        	  temp.reshape(J1,J2,N);

      	::display(temp,"hpLocal","%14.4e ");
            }
            else
            {
      	if( myid==0 )
      	{
        	  ::display(hpLocal,sPrintF("hpLocal parallel myid=%d",myid),"%14.4e ");
      	}
      	
            }
            
        }
        
    // int includeGhost=1;
    // bool ok = ParallelUtility::getLocalArrayBounds(up,upLocal,J1,J2,J3,includeGhost);
        if( ok )
        {
            FOR_3(j1,j2,j3,J1,J2,J3) 
            {
      	for( int n=0; n<numComp; n++ )
      	{
        	  if( numberOfDimensions==2 )
        	  {
          	    if( raxis==0 )
          	    {
            	      maxDiff = max(maxDiff,fabs(hpLocal(j1,j2,j3,n)-hu3(j2,n)));
          	    }
          	    else if( raxis==1 )
          	    {
            	      maxDiff = max(maxDiff,fabs(hpLocal(j1,j2,j3,n)-hu3(j1,n)));
          	    }
          	    else
          	    {
            	      maxDiff = max(maxDiff,fabs(hpLocal(j1,j2,j3,n)-hu3(j1,n)));
          	    }
        	  }
        	  else
        	  {
          	    if( raxis==0 )
          	    {
            	      maxDiff = max(maxDiff,fabs(hpLocal(j1,j2,j3,n)-hu3(j2,j3,n)));
          	    }
          	    else if( raxis==1 )
          	    {
            	      maxDiff = max(maxDiff,fabs(hpLocal(j1,j2,j3,n)-hu3(j1,j3,n)));
          	    }
          	    else
          	    {
            	      maxDiff = max(maxDiff,fabs(hpLocal(j1,j2,j3,n)-hu3(j1,j2,n)));
          	    }
        	  }
      	}
      	
            }
        }
        
        maxDiff = ParallelUtility::getMaxValue( maxDiff );

    // real maxHu = max(fabs(hu3));
        real maxUp = maxNorm(up);
            
    // real maxDiff = max(fabs(hu3-hup3))/max(1.,maxHu);
    // printF(" 3D: myid=%d, maxHu=%9.2e, max |hu-hup|/max(hu) =%8.2e\n",myid,maxHu,maxDiff);

        printF("\n Serial/Parallel: maxUp=%9.2e, max |hu-hup|/max(up) =%8.2e \n\n",maxUp,maxDiff);
        if( check==0 )
            printF("DONE: test Kernel in %d dimensions\n",numberOfDimensions);

        fflush(0);
            
    }

    if(  pmap->decrementReferenceCount()==0 )
    {
        delete pmap;
    }
    

    return maxDiff;
}

// ==============================================================================
//  OLD test of the Radiationkernel
// ==============================================================================
int testRaditionKernel( int numberOfDimensions, int nx, int ny, int nz,
                  			int numberOfPoles,  int orderOfTimeStepping,
                  			real tFinal, real tPlot, real cfl, real ad4, int debug  )
{
    
    real xa=-2.,xb=2.;
    real ya=0., yb=1.;
    real  za=0., zb=1.;

    Mapping *pmap=NULL;

    if( numberOfDimensions==2 )
    {
        pmap = new SquareMapping(xa,xb,ya,yb);
    }
    else
    {
        pmap = new BoxMapping(xa,xb,ya,yb,za,zb);
    }
    pmap->incrementReferenceCount();
    
    Mapping & map = *pmap;
    
    int numberOfGhostPoints=2; // orderOfAccuracy/2+1;

  // printF("Grid points: nx=%d, ny=%d, nz=%d\n",nx,ny,nz);
    
    int nxv[3]={nx,ny,nz};  //
    int side,axis;
    for( int axis=0; axis<map.getDomainDimension(); axis++ )
    {
        map.setGridDimensions(axis,nxv[axis]);
        if( axis !=raxis )
            map.setIsPeriodic(axis,Mapping::derivativePeriodic);
    }
    
  // const int axisp1 = (raxis+1) % numberOfDimensions;
  // const int axisp2 = (raxis+2) % numberOfDimensions;
  // map.setIsPeriodic(axisp1,Mapping::derivativePeriodic);
  // if( numberOfDimensions==3 )
  //   map.setIsPeriodic(axisp2,Mapping::derivativePeriodic);
    
    MappedGrid mg(map);
    for( int axis=0; axis<map.getDomainDimension(); axis++ )
    {
        for( int side=Start; side<=End; side++ )
            mg.setNumberOfGhostPoints(side,axis,numberOfGhostPoints);
    }
    mg.update(MappedGrid::THEmask | MappedGrid::THEcenter );
    bool isRectangular = mg.isRectangular();

    mg.getDeltaX(dx);
    printF(" *** nx=%i, ny=%i, dx=%e dy=%e cfl=%6.4f ad4=%6.4f \n",nxv[0],nxv[1],dx[0],dx[1],cfl,ad4);

    
    int numberOfGridPoints=ny-1;
    int numberOfFields=1;
    const int maxModes=numberOfGridPoints/2-1;  // this is the max allowed
    int numberOfModes=min(maxModes,max(8,numberOfGridPoints/4));  // this is the max allowed
    real period=1.;
    real c=1.;
    
    int currentTimeLevel=0;


  // solution on the boundary:
  // RealArray ub(ny-1,numberOfTimeLevels); 
  // ub=0.;
  // uxb.redim(ny-1,numberOfTimeLevels);
  // uxb=0.;
  // uxxb.redim(ny-1,numberOfTimeLevels);
  // uxxb=0.;
  // uxxxb.redim(ny-1,numberOfTimeLevels);
  // uxxxb=0.;
    


    if( true ) // option == "testKernel" )
    {
    // ---- OLD TEST KERNEL ----
        
        Range all;
        realMappedGridFunction up(mg,all,all,all);  // for parallel transforms
        OV_GET_SERIAL_ARRAY(real,up,upLocal);
        upLocal=0.;

        RealArray hpLocal;
        hpLocal.redim(upLocal);
        
        real xa=-2.,xb=2.;
        real ya=0., yb=1.;
        real za=0., zb=1.;
    
        real dy=(yb-ya)/(ny-1);

        real dt=dy*cfl;
    
        printF(" ---- ny=%d, cfl=%g, dt=%g ----\n",ny,cfl,dt);

        int numberOfModes2 = floor( (nz-1)/2 );
        
        int numberOfGridPoints1=ny-1;   // -1 to avoid periodic images 
        int numberOfGridPoints2=nz-1;

        int nx=ny;
        Range I1=nx, I2=numberOfGridPoints1, I3=numberOfGridPoints2; // excludes periodic images

        Index Jv[3], &J1=Jv[0], &J2=Jv[1], &J3=Jv[2];
        getIndex(mg.dimension(),J1,J2,J3);
        Jv[raxis]=mg.gridIndexRange(rside,raxis);

        real period1=yb-ya;
        real period2=zb-za;
        
        printF(" ---------------------------------------------------------------\n"
         	   " ------------------- TEST RADIATION KERNEL ---------------------\n"
         	   " -------------- Gaussian Integral Exact Solution----------------\n"
         	   " ------- numberOfPoles=%d, orderOfTimeStepping=%d --------------\n"
         	   " ------- numberOfDimensions=%d, nx=%d, ny=%d, nz=%d ------------\n"
         	   " ------- rside=%d, raxis=%d                       --------------\n"
         	   " ------- Jv=[%d,%d][%d,%d][%d,%d]  (boundary index) ------------\n"
                      " ---------------------------------------------------------------\n",
         	   numberOfPoles,orderOfTimeStepping,numberOfDimensions,nx,ny,nz,rside,raxis,
         	   J1.getBase(),J1.getBound(),J2.getBase(),J2.getBound(),J3.getBase(),J3.getBound()
            );



        RadiationKernel radiationKernel;
        radiationKernel.setUseFourierTransformClass(true);

        radiationKernel.setDebug( debug );
        
    // --- new init using a mappedGridFunction ---
        radiationKernel.initialize(up,
                         			       rside,raxis,
                         			       numberOfFields, 
                         			       numberOfModes, numberOfModes2,
                         			       c, orderOfTimeStepping, numberOfPoles );



        Range R2=I2, R3=I3;                            // ** check me Oct 24, 2020
        RealArray u3(I2,I3), hu3(I2,I3), hup3(I2,I3);  // ** check me Oct 24, 2020

        RealArray u(R2),hu(R2),ut(R2),ux(R2);
        RealArray hup(R2);


    // solution on the boundary:
        const int numberOfTimeLevels=10;
        RealArray ub(R2,numberOfTimeLevels),  utb(R2,numberOfTimeLevels), uxtb(R2,numberOfTimeLevels);
        RealArray uxb(R2,numberOfTimeLevels), uxxb(R2,numberOfTimeLevels);
        RealArray ugb(R2,numberOfTimeLevels);
    
    // RealArray ub(ny-1,numberOfTimeLevels),  utb(ny-1,numberOfTimeLevels), uxtb(ny-1,numberOfTimeLevels);
    // RealArray uxb(ny-1,numberOfTimeLevels), uxxb(ny-1,numberOfTimeLevels);
    // RealArray ugb(ny-1,numberOfTimeLevels);
    
        ub=0.;
        uxb=0.;
        uxxb=0.;
    
        utb=0.;
        uxtb=0.;
        ugb=0.;
    

        Range I=ny-1;

        real eps=1.e-5; // **** what should this be ?
    

        real t=0.;

        int numberOfSteps=int( tFinal/dt+1.5 );
        int nPlot = int( tPlot/dt+1.5 );
    
        int currentTimeLevel=0;
        for( int step=0; step<numberOfSteps; step++ )
        {
            t=step*dt;
            printF("\n ================= step=%d t=%9.3e ny=%d, numberOfPoles=%d, orderOfTimeStepping=%d ===============\n",
           	     step,t,ny,numberOfPoles,orderOfTimeStepping);

            for( int i=0; i<ny-1; i++ )
            {
      	real x=xa; 
        
      	real y=ya+i*dy; 
        
      	real uex,uey,uhz;
        
        // Get the true solution
            {
                gisTime=t;
                exmax(wt,wx,wy,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,gisTime);
                uex = wy;
                uey =-wx;
                uhz= wt;
            }
      	u(i)=uey;
        
        // Approx. the time-derivative of the true solution by differences 
      	real tp=t-eps;
            {
                gisTime=tp;
                exmax(wt,wx,wy,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,gisTime);
                uex = wy;
                uey =-wx;
                uhz= wt;
            }
      	ut(i)=uey; 
      	tp=t+eps;
            {
                gisTime=tp;
                exmax(wt,wx,wy,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,x,y,gisTime);
                uex = wy;
                uey =-wx;
                uhz= wt;
            }
      	ut(i)=(uey-ut(i))/(2.*eps);
        
        // Approx. the x-derivative of the true solution by differences 
      	real xp=x-eps;
            {
                gisTime=t;
                exmax(wt,wx,wy,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,xp,y,gisTime);
                uex = wy;
                uey =-wx;
                uhz= wt;
            }
      	ux(i)=uey;
      	xp=x+eps;
            {
                gisTime=t;
                exmax(wt,wx,wy,nsources,xs[0],ys[0],tau[0],var[0],amp[0],period,xp,y,gisTime);
                uex = wy;
                uey =-wx;
                uhz= wt;
            }
      	ux(i)=(uey-ux(i))/(2.*eps);
        

            }
    
            if( numberOfDimensions==2 )
            {
	// radiationKernel.evaluateKernelParallelOld( dt, u, hup, upLocal,hpLocal  );
      	radiationKernel.evaluateKernelParallel( dt, upLocal,hpLocal  );
	// ::display(hup,"hup parallel","%16.8e ");

      	radiationKernel.evaluateKernel( dt, u, hu );
	// ::display(hu,"hu serial","%16.8e ");
        
      	real maxHu = max(fabs(hu));
            
      	real maxDiff = max(fabs(hu-hup))/max(1.,maxHu);
      	printF(" maxHu=%9.2e, max |hu-hup|/max(hu) =%8.2e\n",maxHu,maxDiff);
            }
            else if( numberOfDimensions==3 )
            {
        // ----- TEST KERNEL IN 3D ----

      	for( int i3=0; i3<numberOfGridPoints2; i3++ )
      	{
        	  u3(I2,i3)=u(I2);   // independent of i3 for now 
      	}
      	
                radiationKernel.evaluateKernel( dt, u3, hu3 );

        // copy back into 2D case 
                const int i3=0;
                hu(I2) = hu3(I2,i3);   // independent of i3 for now 

            }
            else
            {
                OV_ABORT("error: numberOfDimensions");
            }
            
            int i0=ny/2;

            const int m=currentTimeLevel;
            int mm1= (m-1+numberOfTimeLevels) % numberOfTimeLevels;
            int mm2= (m-2+numberOfTimeLevels) % numberOfTimeLevels;
            int mm3= (m-3+numberOfTimeLevels) % numberOfTimeLevels;
            int mm4= (m-4+numberOfTimeLevels) % numberOfTimeLevels;
        
            ub(I,m)=u(I);
      // One sided 4th-order approx to ut on the boundary:
      // utb(I,m) = -( (-25./12.)*ub(I,m)+4.*ub(I,mm1)-3.*ub(I,mm2)+4./3.*ub(I,mm3)-.25*ub(I,mm4) )/dt;
      // uxb(I,m) = utb(I,m)+hu(I);

      // fourth-order Adams-Moulton coefficients:
            const real am41=(9./24.), am42=(19./24.), am43=(-5./24.), am44=(1./24.);

      //  ut + un + H(u) = 0
      //  ut = -( un + H(u) ) = G 
      //  [u(n+1)-u(n)]/dt = am41*G(n+1) + am42*G(n) + am43*G(n-1) + am44*G(n-2) 

      // Solving for un(n+1) gives
      //
      //   am41*un(n+1) = -am41*H(u(n+1)) +  am42*G(n) + am43*G(n-1) + am44*G(n-2)  - [u(n+1)-u(n)]/dt
            real sigma=-1.;  // 1-2*side
        
            uxb(I,m) = sigma*( -hu(I) + ( am42*ugb(I,mm1) + am43*ugb(I,mm2) + am44*ugb(I,mm3) 
                            				    - (ub(I,m)-ub(I,mm1))/dt )/am41 );
        
            real uxpc=uxb(i0,m);
        

            utb(I,m)= - ( sigma*uxb(I,m) + hu(I) ); 

            real utpc=utb(i0,m);
        

      // One sided 4th-order approx to ut on the boundary:
            utb(I,m) = -( (-25./12.)*ub(I,m)+4.*ub(I,mm1)-3.*ub(I,mm2)+4./3.*ub(I,mm3)-.25*ub(I,mm4) )/dt;
            uxb(I,m) = utb(I,m)+hu(I);

            ugb(I,m)= - ( sigma*uxb(I,m) + hu(I) );  // save G -- this is just utb
        
        
            if( ( step % nPlot) == 0 )
            {
      	printF(" step=%i t=%8.2e, dt=%8.2e, u=%8.2e ut=%8.2e utb-ut=%8.2e uxb-ux=%8.2e  "
             	       "hu=%8.2e ut-ux+hu=%8.2e\n", 
             	       step,t,dt,u(i0),ut(i0),utb(i0,m)-ut(i0),uxb(i0,m)-ux(i0),hu(i0),ut(i0)-ux(i0)+hu(i0));
        
      	printF("   ux=%9.3e uxpc=%9.3e ut=%9.3e utpc=%9.3e \n",ux(i0),uxpc,ut(i0),utpc);
            
            }
        
            currentTimeLevel=(currentTimeLevel+1) % numberOfTimeLevels;
    
        } // end for step
        
        printF("DONE: test Kernel\n");
        OV_ABORT("STOP HERE FOR NOW");
        
            
    } // end test kernel 
    

    return 0;


} // test RadiatoionKernel 




// ******************************************************************************
// ************************** Test Radiation BC's *******************************
// ******************************************************************************

int 
main(int argc, char *argv[])
{

    int debug=0;
    bool useOpt=true; // false;

    int numberOfDimensions=2;
    int numComp=1;
    int numSteps=1; // for testRadiationKernel
    
    aString option = "testKernel";
    int useParallelVersion=false;

  // These are for the exact solution:
    xs[0]=0.;
    ys[0]=1.e-8*1./3.;  // should not be on a grid point
    tau[0]=-.95;
    var[0]=30.;
    amp[0]=1.;


    real tFinal=1.; 
    real tPlot=.2;  // plot this often

    real cfl=.25, nu=0.;
    real ad4=0.;  // coeff of the artificial dissipation.

    int powerPow=4;
    enum sigmaProfileEnum
    {
        ramp1Profile,
        ramp3Profile,
        powerProfile
    } sigmaProfile=powerProfile; // ramp3Profile;
    
    InitialConditionOptionEnum icOption=smoothPulse;
  // InitialConditionOptionEnum icOption=gaussianIntegral;

//  real sx=-1.,sy=-1.;

  // For fourth-order we use order=5 Adams-Moulton time-stepping for the Radiation Kernel Kernel
    int orderOfTimeStepping=5; // 6; // 5; // 4;
    int numberOfPoles=21;
    
    
          // "   -sx=[]    : coeff of absorbing layer\n"
          // "   -power3   : use power law x^3 profile for sigma\n"
          // "   -power4,5,6   : use power law x^4 profile for sigma\n";

    Overture::start(argc,argv);  // initialize Overture
    Optimization_Manager::setForceVSG_Update(Off);
    Overture::turnOnMemoryChecking(true);

    const int myid = max(0,Communication_Manager::My_Process_Number);
    const int np   = max(1,Communication_Manager::Number_Of_Processors);


    printF("Usage: `trad [options]' \n"
                    "   -noplot                               \n"
                    "   -option=[checkParallel|testKernel|testRad]  \n"
                    "   -icOption=[sp|p|gi]  : sp=smoothPulse, gi=Gaussian Integral \n"
                    "   -noopt    : do not used optimized version  \n"
                    "   -order=[2/4] : order of accuracy\n"
                    "   -rside=[0/1] : left or right side for BC\n"
                    "   -nx=[]    : number of grid points on unit square\n"
                    "   -ny=[]    : number of grid points on unit square\n"
                    "   -nz=[]    : number of grid points on unit square\n"
                    "   -nd=[2/3] : number of dimensions\n"
                    "   -numComp<i> : number of components\n"
                    "   -numberOfPoles=[21|31] : number of poles in Kernel expansion\n"
                    "   -orderOfTimeStepping=[1,2,...] : order of time-stepping for Kernel.\n"
                    "   -debug[]  : debug parameter\n"
                    "   -tFinal=[]: final time\n"
        	  "   -cfl=[]   : cfl to use\n"
        	  "   -tPlot=[] : times between plots\n"
                    "   -check=1 : run regression checks\n"
                    "   -useParallelVersion=[0|1] : if true, use new parallel version.\n"
        );
    

  // --- Number of grid points in each direction ---
    int nx=9, ny=9, nz=9;

    int check=0;  // set to 1 to perform a sequence of tests

    aString buff;
    aString commandFileName="";
    bool plotOption=true;  // by default we plot interactively
    int len=0;
    if( argc > 1 )
    {
        for( int i=1; i<argc; i++ )
        {
            aString line=argv[i];

            if( len=line.matches("-check=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&check);
      	printF("Setting check=%i\n",check);
            }

            else if( line=="-noplot" )
                plotOption=FALSE; 
            else if( line=="-noopt" )
                useOpt=false;

            else if( len=line.matches("-option=") )
            {
      	option = line(len,line.length()-1);
                printF("Setting option=[%s]\n",(const char*)option);
            }
            
            else if( len=line.matches("-icOption=") )
            {
      	aString opt = line(len,line.length()-1);
      	if( opt == "sp" )
        	  icOption=smoothPulse;
      	else if( opt == "p" )
        	  icOption=pulse;
      	else if( opt == "gi" )
        	  icOption=gaussianIntegral;
      	
                printF(" Setting icOption=%d (0=smoothPulse,1=pulse,2=GaussianIntegral))\n",(int) icOption);
            }
            else if( len=line.matches("-numComp=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&numComp);
                printF(" Setting numComp=%i\n",numComp);
            }

            else if( len=line.matches("-numSteps=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&numSteps);
                printF(" Setting numSteps=%i\n",numSteps);
            }

      // else if( len=line.matches("-sx=") )
      // {
      // 	sScanF(line(len,line.length()-1),"%e",&sx);
      //   sy=sx;
      //   printF(" Setting sx=sy=%e\n",sx);
      // }
            else if( len=line.matches("-tFinal=") )
            {
      	sScanF(line(len,line.length()-1),"%e",&tFinal);
                printF(" Setting tFinal=%e\n",tFinal);
            }
            else if( len=line.matches("-tPlot=") )
            {
      	sScanF(line(len,line.length()-1),"%e",&tPlot);
                printF(" Setting tPlot=%e\n",tPlot);
            }
            else if( len=line.matches("-alpha=") )
            {
      	sScanF(line(len,line.length()-1),"%e",&alpha);
                printF(" Setting alpha=%e\n",alpha);
            }
            else if( len=line.matches("-xPulse=") )
            {
      	sScanF(line(len,line.length()-1),"%e",&xPulse);
                printF(" Setting xPulse=%e\n",xPulse);
            }
            else if( len=line.matches("-cfl=") )
            {
      	sScanF(line(len,line.length()-1),"%e",&cfl);
                printF(" Setting cfl=%e\n",cfl);
            }
            else if( len=line.matches("-nx=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&nx);
                printF(" Setting nx=%i\n",nx);
            }
            else if( len=line.matches("-ny=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&ny);
                printF(" Setting ny=%i\n",ny);
            }
            else if( len=line.matches("-nz=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&nz);
                printF(" Setting nz=%i\n",nz);
            }
            else if( len=line.matches("-numberOfPoles=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&numberOfPoles);
                printF(" Setting numberOfPoles=%i\n",numberOfPoles);
            }
            else if( len=line.matches("-orderOfTimeStepping=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&orderOfTimeStepping);
                printF(" Setting orderOfTimeStepping=%i\n",orderOfTimeStepping);
            }
            else if( len=line.matches("-debug=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&debug);
                printF(" Setting debug=%i\n",debug);
            }
            else if( len=line.matches("-nd=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&numberOfDimensions);
                printF(" Setting numberOfDimensions=%i\n",numberOfDimensions);
            }
            else if( len=line.matches("-order=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&orderOfAccuracy);
                printF(" Setting orderOfAccuracy=%i\n",orderOfAccuracy);
            }
            else if( len=line.matches("-rside=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&rside);
                printF(" Setting rside=%i\n",rside);
            }
            else if( len=line.matches("-raxis=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&raxis);
                printF(" Setting raxis=%i\n",raxis);
            }

            else if( len=line.matches("-useParallelVersion=") )
            {
      	sScanF(line(len,line.length()-1),"%i",&useParallelVersion);
                printF(" Setting useParallelVersion=%i\n",useParallelVersion);
            }

      // else if( line.matches("-power3") )
      // {
      // 	sigmaProfile=powerProfile;
      // 	powerPow=3;
      //   printF(" Setting sigmaProfile=powerProfile (power=%i)\n",powerPow);
      // }
      // else if( line.matches("-power4") )
      // {
      // 	sigmaProfile=powerProfile;
      // 	powerPow=4;
      //   printF(" Setting sigmaProfile=powerProfile (power=%i)\n",powerPow);
      // }
      // else if( line.matches("-power5") )
      // {
      // 	sigmaProfile=powerProfile;
      // 	powerPow=5;
      //   printF(" Setting sigmaProfile=powerProfile (power=%i)\n",powerPow);
      // }
      // else if( line.matches("-power6") )
      // {
      // 	sigmaProfile=powerProfile;
      // 	powerPow=6;
      //   printF(" Setting sigmaProfile=powerProfile (power=%i)\n",powerPow);
      // }
            else
      	commandFileName=line;
        }
    }




  // **************************************************************************
  // **************** CHECK SERIAL VERSUS PARALLEL ****************************
  // **************************************************************************
    if( option == "checkParallel" )
    {
    // **NEW WAY **
        real maxErr=0.;

        int numChecks= check==0 ? 1 : 10;
        for( int icheck=0; icheck<numChecks; icheck++ )
        {
            if( numChecks>1 )
            {
        // --- 3D checks ---
      	if( icheck==0 )
      	{
        	  numberOfDimensions=2; nx=7; ny=7; nz=1; rside=0; raxis=0;
      	}
      	else if( icheck==1 )
      	{
        	  numberOfDimensions=2; nx=9; ny=7; nz=1; rside=1; raxis=0; 
      	}
      	else if( icheck==2 )
      	{
        	  numberOfDimensions=2; nx=9; ny=17; nz=1; rside=0; raxis=1;
      	}
      	else if( icheck==3 )
      	{
        	  numberOfDimensions=2; nx=17; ny=17; nz=1;  rside=1; raxis=1;
      	}
        // --- 3D checks ---
      	else if( icheck==4 )
      	{
        	  numberOfDimensions=3; nx=13; ny=13; nz=13; rside=0;  raxis=0;
      	}
      	else if( icheck==5 )
      	{
        	  numberOfDimensions=3; nx=7; ny=9; nz=13; rside=1;  raxis=0;
      	}
      	else if( icheck==6 )
      	{
        	  numberOfDimensions=3; nx=9; ny=9; nz=9; rside=0;  raxis=1;
      	}
      	else if( icheck==7 )
      	{
        	  numberOfDimensions=3; nx=17; ny=17; nz=17; rside=1;  raxis=1;
      	}
      	else if( icheck==8 )
      	{
        	  numberOfDimensions=3; nx=9; ny=17; nz=13; rside=0;  raxis=2;
      	}
      	else if( icheck==9 )
      	{
        	  numberOfDimensions=3; nx=7; ny=13; nz=7; rside=1;  raxis=2;
      	}
      	else
      	{
        	  OV_ABORT("error");
      	}
      	

            }
            else
            {
	// printF(" tft: >>> nd=%d, ndfft=%d, nz=%d, ny=%d, nz=%d\n",nd,ndfft,nx,ny,nz);
            }

            real err = testParallelKernel( numSteps,numComp, numberOfDimensions, numberOfPoles,orderOfTimeStepping,
                              				      nx, ny, nz, rside, raxis, debug, check );

            printF(" check %d: numberOfDimensions=%d, numComp=%d, nx=%d, ny=%d, nz=%d, [rside,raxis]=[%d,%d] err=%9.2e\n",
           	     icheck,numberOfDimensions,numComp,nx,ny,nz,rside,raxis,err);

            maxErr = max(maxErr,err);
            
        } // end for icheck 
        
        if( check )
        {
            if( maxErr > REAL_EPSILON*1.e5 )
            {
      	printF("\n================== FAILURE in checks: maxErr =%9.2e ===============\n",maxErr);
            }
            else
            {
      	printF("\n================== SUCCESS: np=%d, numSteps=%d, numComp=%d, Done all %d checks: maxErr =%9.2e ===============\n",
             	       np,numSteps,numComp,numChecks,maxErr);
            }
        
        }
        
    // OV_ABORT("STOP HERE FOR NOW");
        Overture::finish();          
        return 0;

    }
    




  // -------- create a SQUARE or BOX grid ------------------
    real xa=-2., xb=0., ya=0., yb=1., za=-.5, zb=.5;

    if( rside==0 && raxis==0 )
    {
        xa=-2., xb=0.;
    }
    else if( rside==1 && raxis==0 )
    {
        xa=0, xb=2.;
    }
    else
    {
        Overture::abort("error");
    }
    if( icOption==smoothPulse )
    {
        xa=-1; xb=1.;
        ya=-2; yb=2.;
    }
    
    Mapping *pmap=NULL;

    if( numberOfDimensions==2 )
    {
        pmap = new SquareMapping(xa,xb,ya,yb);
    }
    else
    {
        pmap = new BoxMapping(xa,xb,ya,yb,za,zb);
    }
    pmap->incrementReferenceCount();
    
    Mapping & map = *pmap;
    
    int numberOfGhostPoints=2; // orderOfAccuracy/2+1;

  // printF("Grid points: nx=%d, ny=%d, nz=%d\n",nx,ny,nz);
    
    int nxv[3]={nx,ny,nz};  //
    int side,axis;
    for( int axis=0; axis<map.getDomainDimension(); axis++ )
    {
        map.setGridDimensions(axis,nxv[axis]);
        if( axis !=raxis )
            map.setIsPeriodic(axis,Mapping::derivativePeriodic);
    }
    

    MappedGrid mg(map);
    for( int axis=0; axis<map.getDomainDimension(); axis++ )
    {
        for( int side=Start; side<=End; side++ )
            mg.setNumberOfGhostPoints(side,axis,numberOfGhostPoints);
    }
    mg.update(MappedGrid::THEmask | MappedGrid::THEcenter );
    bool isRectangular = mg.isRectangular();

    mg.getDeltaX(dx);
    printF(" *** nx=%i, ny=%i, dx=%e dy=%e cfl=%6.4f ad4=%6.4f \n",nxv[0],nxv[1],dx[0],dx[1],cfl,ad4);

    
    int numberOfGridPoints=ny-1;
    int numberOfFields=1;
    const int maxModes=numberOfGridPoints/2-1;  // this is the max allowed
    int numberOfModes=min(maxModes,max(8,numberOfGridPoints/4));  // this is the max allowed
    real period=1.;
    real c=1.;
    
    int currentTimeLevel=0;





    if( option == "testKernel" )
    {
    // ---- OLD TEST KERNEL ----
        testRaditionKernel( numberOfDimensions, nx, ny, nz,
                  			numberOfPoles, orderOfTimeStepping,
                  			tFinal, tPlot, cfl, ad4, debug );
            
    } // end test kernel 
    

    
  // *********************************************************************************
  // *********************************************************************************
  // *********************************************************************************

    bool testWaveSolver=true;
    
    if( testWaveSolver )
    {
    
        PlotStuff ps(plotOption,"wave equation");
        PlotStuffParameters psp;
    
    // By default start saving the command file called "wave.cmd"
        aString logFile="trad.cmd";
        ps.saveCommandFile(logFile);
        printF("User commands are being saved in the file [%s]\n",logFile);
    
        if( commandFileName!="" )
            ps.readCommandFile(commandFileName);
    
    
        aString nameOfShowFile="trad.show";
        Ogshow show;
        bool showFileIsOpen=false;
        
        FILE *debugFile=NULL;
        if( debug>0 )
        {
            debugFile= fopen("trad.debug","w");
        }
    
    
        
    
        MappedGridOperators operators(mg);                        // operators for a CompositeGrid
        if( orderOfAccuracy==4 )
            operators.setOrderOfAccuracy(4);                          // for fourth order
        BoundaryConditionParameters bcParams;
    
        Range all;
        realMappedGridFunction u[2];
    // realMappedGridFunction vx[2]; // , wx[2], vy[2], wy[2], vz[2], wz[2];
        realMappedGridFunction uLap(mg),uLapSq(mg);
        
        realMappedGridFunction err(mg,all,all,all);
        err=0.;
    
        u[0].updateToMatchGrid(mg);
        u[0].setOperators(operators);                                 
        u[0].setName("u");                                              // name the grid function
        u[0]=0.;
        u[1]=u[0];
    
        const int nc = 2; 
        
        realMappedGridFunction uvw(mg,all,all,all,nc);
        uvw.setName("u",0);
        uvw.setName("err",1);
        
    
  //   vx[0].updateToMatchGrid(mg); vx[1].updateToMatchGrid(mg);  
  //   vx[0].setOperators(operators); vx[1].setOperators(operators);   vx[0]=0.; vx[1]=0.;
    
    
        const realArray & center = mg.center();
        const realArray & x = center(all,all,all,0);
        const realArray & y = center(all,all,all,1);
        
    // realMappedGridFunction laplacian(mg);  // holds laplacian
        
        real t=0;
        real cSquared=c*c;
            
        
    // estimate the time step -- for now approximate by a first order wave equation.
        real dt;
        
    
  //   const int numberOfDimensions=mg.numberOfDimensions();
    
        if( numberOfDimensions==2 )
            dt=cfl*1./( c*sqrt( 1./(dx[0]*dx[0])+1./(dx[1]*dx[1]) ) );  // only valid for rectangular grids.
        else
            dt=cfl*1./( c*sqrt( 1./(dx[0]*dx[0])+1./(dx[1]*dx[1])+1./(dx[2]*dx[2]) ) ); 
    
        real dt0=dt;
        int numStepsToPlot = int(tPlot/dt+.999999);
    // int numSteps = int(tFinal/dt+1.);
    // numSteps = (numSteps/numStepsToPlot)*numStepsToPlot;
    
        numSteps= int(tFinal/tPlot+.999999)*numStepsToPlot;
        dt=tFinal/numSteps;
        assert( dt<= dt0 );
        
        printF("--- tPlot=%f, tPlot/dt=%f numStepsToPlot=%i numSteps=%i \n",tPlot,tPlot/dt,numStepsToPlot,numSteps);
        
    
        real dtSquared=dt*dt;
        
    // Assign the initial conditions
    
        getInitialConditions( icOption, u, t, dt, c );
    // setBoundaryConditions( icOption, u, t, dt );  // this is needed to init radiation BC?
    
    
        Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
        Index Jv[3], &J1=Jv[0], &J2=Jv[1], &J3=Jv[2];
    
        getIndex(mg.gridIndexRange(),I1,I2,I3);
        getIndex(mg.gridIndexRange(),J1,J2,J3);
    
        
        mg.update(MappedGrid::THEmask);
        intArray & mask = mg.mask();
        OV_GET_SERIAL_ARRAY(int,mask,maskLocal);
        
        int includeGhost=1;
        bool ok = ParallelUtility::getLocalArrayBounds(mask,maskLocal,I1,I2,I3,includeGhost);     
    
        printf("myid=%d: maskLocal bounds: [%d,%d][%d,%d]\n",myid,maskLocal.getBase(0),maskLocal.getBound(0),
           	 maskLocal.getBase(1),maskLocal.getBound(1));
    
        RadiationBoundaryCondition radiationBoundaryCondition(orderOfAccuracy);
        radiationBoundaryCondition.setDebug(debug);
        radiationBoundaryCondition.useParallelVersion( useParallelVersion );
        
        int nc1=0, nc2=0;      // component range
        printF("Setting NLRBC for (rside,raxis)=(%d,%d)\n",rside,raxis);
        
        radiationBoundaryCondition.setNumberOfPoles(numberOfPoles);
    
        radiationBoundaryCondition.initialize(u[0],rside,raxis,nc1,nc2,c);
    
    // radiationBoundaryCondition.initialize(mg,rside,raxis,nc1,nc2,c);
    
    
    
    // Build a dialog menu for changing parameters
        GUIState gui;
        DialogData & dialog=gui;
    
        dialog.setWindowTitle("Radiation Conditions");
        dialog.setExitCommand("finish", "finish");
    
        dialog.setOptionMenuColumns(1);
    
  //    aString accuracyLabel[] = {"second order", "fourth order", "" };
  //    dialog.addOptionMenu("accuracy:", accuracyLabel, accuracyLabel, (orderOfAccuracy==2 ? 0 : 1) );
    
  //    aString initialConditionLabel[] = {"smooth pulse", "pulse", "" };
  //    dialog.addOptionMenu("Initial Condition:",initialConditionLabel,initialConditionLabel,(int)icOption );
    
        aString pbLabels[] = {"continue",
                                                    "contour",
                                                    "grid",
                                                    "exit",
                      			""};
        int numRows=2;
        dialog.setPushButtons( pbLabels, pbLabels, numRows ); 
    
        const int numberOfComponents = uvw.getComponentBound(0)-uvw.getComponentBase(0)+1;
    // create a new menu with options for choosing a component.
        aString *cmd = new aString[numberOfComponents+1];
        aString *label = new aString[numberOfComponents+1];
        for( int m=0; m<numberOfComponents; m++ )
        {
            label[m]=uvw.getName(m);
            cmd[m]="plot:"+uvw.getName(m);
    
        }
        cmd[numberOfComponents]="";
        label[numberOfComponents]="";
            
        dialog.addOptionMenu("plot component:", cmd,label,0);
        delete [] cmd;
        delete [] label;
    
        bool saveShowFile=false;
    
        aString tbCommands[] = {"save show file",
                                                        "compare to reference show file",
                        			  ""};
        int tbState[10];
        tbState[0] = saveShowFile==true;
        tbState[1] = compareToReferenceShowFile==true;
        int numColumns=1;
        dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 
    
    // ----- Text strings ------
        const int numberOfTextStrings=20;
        aString textCommands[numberOfTextStrings];
        aString textLabels[numberOfTextStrings];
        aString textStrings[numberOfTextStrings];
    
        int nt=0;
        textCommands[nt] = "cfl";  textLabels[nt]=textCommands[nt];
        sPrintF(textStrings[nt], "%g",cfl);  nt++; 
        textCommands[nt] = "tFinal";  textLabels[nt]=textCommands[nt];
        sPrintF(textStrings[nt], "%g",tFinal);  nt++; 
        textCommands[nt] = "tPlot";  textLabels[nt]=textCommands[nt];
        sPrintF(textStrings[nt], "%g",tPlot);  nt++; 
  //    textCommands[nt] = "tShow";  textLabels[nt]=textCommands[nt];
  //    sPrintF(textStrings[nt], "%g",tShow);  nt++; 
  //    textCommands[nt] = "artificial dissipation";  textLabels[nt]=textCommands[nt];
  //    sPrintF(textStrings[nt], "%g",ad4);  nt++; 
  //    textCommands[nt] = "pulse params";  textLabels[nt]=textCommands[nt];
  //    sPrintF(textStrings[nt], "%g %g %g (alpha,a0,pulsePow)",alpha,a0,pulsePow);  nt++; 
    // null strings terminal list
        textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
        dialog.setTextBoxes(textCommands, textLabels, textStrings);
    
        
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
        ps.pushGUI(gui);
        aString answer, line;
        int grid;
        int current=0;
        int step=0;
    
        real time0=getCPU();
        for(;;) 
        {
            real maxU =0.;
      // maxU=max(fabs(u[current]));
            maxU=maxNorm(u[current]);
            getErrors( icOption,u[current],err,t);
    
    
            real errMax = maxNorm(err);
    
            real cpu0= getCPU()-time0;
            cpu0=ParallelUtility::getMaxValue(cpu0);
            
            printF("completed step %6i, t=%8.2e cfl=%3.2f max(u)=%8.2e, err=%8.2e (cpu =%8.2e)\n",
                                step,t,cfl,maxU,errMax,cpu0);

      // real time1=getCPU();
            ps.erase();
            psp.set(GI_TOP_LABEL,sPrintF(buff,"Wave equation, RBC, t=%5.3f, order=%i",t,orderOfAccuracy));
    
    
            OV_GET_SERIAL_ARRAY(real,u[current],ucLocal);
            OV_GET_SERIAL_ARRAY(real,uvw,uvwLocal);
            OV_GET_SERIAL_ARRAY(real,err,errLocal);
            
    
      //uvw(all,all,all,0)=u[current];
      //uvw(all,all,all,1)=err;
            
            uvwLocal(all,all,all,0)=ucLocal;
            uvwLocal(all,all,all,1)=errLocal;
            
    
            
            PlotIt::contour(ps,uvw,psp);
    
            ps.getAnswer(answer,"");      
    
      // time0=getCPU()-time1;  // subtract off time for waiting
            if( answer=="exit" || answer=="finish" )
            {
                break;
            }
            else if( answer.matches("contour") )
            {
                ps.erase();
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
                psp.set(GI_TOP_LABEL,sPrintF(buff,"Wave equation, t=%5.3f",t));
    
    
                PlotIt::contour(ps,uvw,psp);
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
            }
            else if( answer.matches("grid") )
            {
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
                PlotIt::plot(ps,mg,psp);                          // plot the grid
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
            }
            else if( len=answer.matches("cfl") )
            {
                real cflOld=cfl;
    
                sScanF(answer(len,answer.length()-1),"%e",&cfl);
                cout << " cfl=" << cfl << endl;
                dialog.setTextLabel("cfl",sPrintF(line, "%g",cfl));
    
                dt=dt*cfl/cflOld;
            }
            else if( len=answer.matches("tFinal") )
            {
                sScanF(answer(len,answer.length()-1),"%e",&tFinal);
                cout << " tFinal=" << tFinal << endl;
                dialog.setTextLabel("tFinal",sPrintF(line, "%g",tFinal));
            }
            else if( len=answer.matches("tPlot") )
            {
                sScanF(answer(len,answer.length()-1),"%e",&tPlot);
                cout << " tPlot=" << tPlot << endl;
                dialog.setTextLabel("tPlot",sPrintF(line, "%g",tPlot));
            }
            else if( len=answer.matches("plot:") )
            {
        // plot a new component
                aString name = answer(len,answer.length()-1);
                int component=-1;
                for( int n=0; n<numberOfComponents; n++ )
                {
          	if( name==uvw.getName(n) )
          	{
            	  component=n;
            	  break;
          	}
                }
                if( component==-1 )
                {
          	printF("ERROR: unknown component name =[%s]\n",(const char*)name);
          	component=0;
                }
                printF(" Set the component to plot to %i\n",component);
                
                dialog.getOptionMenu(0).setCurrentChoice(component);
                psp.set(GI_COMPONENT_FOR_CONTOURS,component);
            }
            else if( len=answer.matches("save show file") )
            {
                int value;
                sScanF(answer(len,answer.length()-1),"%i",&value); saveShowFile=value;
                dialog.setToggleState("save show file",saveShowFile==true);
                if( saveShowFile )
          	printF("Save show file %s \n",(const char*)nameOfShowFile);
            }
            else if( len=answer.matches("compare to reference show file") )
            {
                int value;
                sScanF(answer(len,answer.length()-1),"%i",&value); compareToReferenceShowFile=value;
                dialog.setToggleState("compare to reference show file",compareToReferenceShowFile==true);
            }
            
            else if( answer=="continue" )
            {
                if( saveShowFile && !showFileIsOpen )
                {
                    showFileIsOpen=true;
        	  show.open( nameOfShowFile );                               // create a show file
        	  show.saveGeneralComment(sPrintF(buff,"Wave equation, RBC, t=%5.3f",t));
  	// show.setFlushFrequency(2);                                 
                }
    
                if( saveShowFile  )  
                {
        	  show.startFrame();                                         // start a new frame
        	  show.saveComment(0,sPrintF(buff,"Wave equation"));
        	  show.saveComment(1,sPrintF(buff,"Wave equation, RBC, t=%5.3f",t));
        	  show.saveSolution( u[current] );                                        // save the current grid function
                }
                
    
  //        t=0.;
  //        getInitialConditions( icOption, u, t, dt, c );
    
                const int numberOfTimeSteps=int( tFinal/dt+.5);
                const int plotSteps = (int)max(1.,tPlot/dt+.5);
                for( int i=0; i<plotSteps; i++ )                    // take some time steps
                {
                    current=step %2;
                    int next = (step+1) %2;
    
                    if( debug & 2 )
        	  {
          	    display(u[current],sPrintF("\n u at step=%i, time t=%9.3e",step,t),debugFile,"%9.2e ");
        	  }
          	
    
          	int ex=0;
    
          	realMappedGridFunction & u1 = u[current];
          	realMappedGridFunction & u2 = u[next];
            
      	OV_GET_SERIAL_ARRAY(real,u1,u1Local);
      	OV_GET_SERIAL_ARRAY(real,u2,u2Local);
      	OV_GET_SERIAL_ARRAY(real,uLap,uLapLocal);
      	OV_GET_SERIAL_ARRAY(real,uLapSq,uLapSqLocal);
          	
            	  
          	real *u1gp= u1Local.Array_Descriptor.Array_View_Pointer3;
          	real *u2gp= u2Local.Array_Descriptor.Array_View_Pointer3;
          	real *lapp= uLapLocal.Array_Descriptor.Array_View_Pointer3;
          	const int uDim0=u1.getRawDataSize(0);
          	const int uDim1=u1.getRawDataSize(1);
          	const int d1=uDim0, d2=d1*uDim1; 
      	#define U1G(i1,i2,i3) u1gp[(i1)+(i2)*d1+(i3)*d2]
      	#define U2G(i1,i2,i3) u2gp[(i1)+(i2)*d1+(i3)*d2]
    
      	#define LAP(i1,i2,i3) lapp[(i1)+(i2)*d1+(i3)*d2]
    
    
      	getIndex(mg.gridIndexRange(),I1,I2,I3);
      	int includeGhost=1;
      	bool ok = ParallelUtility::getLocalArrayBounds(u1,u1Local,I1,I2,I3,includeGhost);     
    
    
      	const real cdtsq=(c*c)*(dt*dt);
    
          	if( ok )
          	{
            	  if( orderOfAccuracy==2 )
            	  {
    
              	    getIndex(mg.gridIndexRange(),J1,J2,J3);
              	    ok = ParallelUtility::getLocalArrayBounds(u1,u1Local,J1,J2,J3,includeGhost);     
    
              	    operators.derivative(MappedGridOperators::laplacianOperator,u1,uLap,J1,J2,J3);  
              	    if( ad4>0. )
              	    {
                	      real ad4dt=ad4*dt;
    
                	      if( mg.numberOfDimensions()==2 )
                	      {
                		int i1,i2,i3;
                		FOR_3(i1,i2,i3,I1,I2,I3) // loop over all points
                		{
  		  // add a 'fourth' order dissipation  ad4 h^4 dt (u.xxxx).t to (c*dt)^2*laplacian(u)
                  		  LAP(i1,i2,i3)=cdtsq*LAP(i1,i2,i3)+ad4dt*( FD4_2D(U1G,i1,i2,i3)-FD4_2D(U2G,i1,i2,i3) );
                		}
                	      }
                	      u2Local(I1,I2,I3,ex)=2.*u1Local(I1,I2,I3,ex)-u2Local(I1,I2,I3,ex) + uLapLocal(I1,I2,I3,ex);
              	    }
              	    else
              	    {
                	      u2Local(I1,I2,I3,ex)=2.*u1Local(I1,I2,I3,ex)-u2Local(I1,I2,I3,ex) + cdtsq*uLapLocal(I1,I2,I3,ex);
              	    }
    
  	    // u2(I1,I2,I3,ex)=2.*u1(I1,I2,I3,ex)-u2(I1,I2,I3,ex) +(dtSquared)*( u1.laplacian(I1,I2,I3,ex)(I1,I2,I3,ex) );
    
            	  }
            	  else if( orderOfAccuracy==4 )
            	  {
              	    real cdtsq12=cdtsq*cdtsq/12.;
            	  
  	    // evaluate the Lapacian squared to 2nd order
              	    operators.setOrderOfAccuracy(2);
              	    getIndex(mg.gridIndexRange(),J1,J2,J3,1);
              	    ok = ParallelUtility::getLocalArrayBounds(u1,u1Local,J1,J2,J3,includeGhost);     
    
              	    operators.derivative(MappedGridOperators::laplacianOperator,u1Local,uLapLocal,J1,J2,J3);  
              	    operators.derivative(MappedGridOperators::laplacianOperator,uLapLocal,uLapSqLocal,I1,I2,I3);
    
              	    operators.setOrderOfAccuracy(4); 
    
  	    // evaluate the Laplacian to fourth order
              	    operators.derivative(MappedGridOperators::laplacianOperator,u1,uLap,J1,J2,J3); 
            	  
              	    u2Local(I1,I2,I3,ex)=2.*u1Local(I1,I2,I3,ex)-u2Local(I1,I2,I3,ex) + cdtsq*uLapLocal(I1,I2,I3,ex) + cdtsq12*uLapSqLocal(I1,I2,I3,ex); 
            	  
            	  }
            	  else
            	  {
              	    Overture::abort("error");
            	  }
    
          	} // end if ok 
    
    
          	t+=dt;
                    currentTimeLevel=(currentTimeLevel+1) % numberOfTimeLevels;
  	// Apply radiation boundary conditions
    
    
          // This will assign Dirichlet BC's:
          	setBoundaryConditions( icOption, u2, t, dt, u1, currentTimeLevel );
    
          	if( false )
          	{
  	  // ::display(uLocal(Ig1,Ig2,Ig3,0),"u-ghost After radEval","%6.3f ");
            	  ::display(u2Local,"u2 setBoundaryConditions","%6.3f ");
          	}
    
          // This will assign the radiation BC
          	if( true )
            	  radiationBoundaryCondition.assignBoundaryConditions( u2,t,dt,u1 );
    
          	
  	  // u2.applyBoundaryCondition(0,BCTypes::dirichlet,BCTypes::allBoundaries,0.);
    
          	u2.finishBoundaryConditions();
          	
          	
          	step++;
                }
                
                current=step %2;
    
            }
            else
            {
                cout << "Unknown command = [" << answer << "]\n";
                ps.stopReadingCommandFile();
                  
            }
        }
    
        ps.popGUI();  // pop dialog
    
        if( saveShowFile && !showFileIsOpen )
        {
            showFileIsOpen=true;
            show.open( nameOfShowFile );                               // create a show file
            show.saveGeneralComment(sPrintF(buff,"Wave equation, RBC, t=%5.3f",t));
      // show.setFlushFrequency(2);                                 
        }
        if( showFileIsOpen )
          show.close();
        
        if( debugFile!=NULL )
            fclose(debugFile);
    }

    
    Overture::finish();          
    return 0;
}
