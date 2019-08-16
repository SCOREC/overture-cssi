// ====================================================================================
///  \file texact.C
///  \brief test program for exact solutions to Maxwell's equations
// ===================================================================================


#include "Overture.h"
#include "PlotStuff.h"
#include "DispersiveMaterialParameters.h"
#include "display.h"


#include "SphereExactSolutions.h"

#include "PlaneInterfaceExactSolution.h"

typedef ::real LocalReal;
typedef ::real OV_real;

int computeS( const LocalReal epsHat[2], const LocalReal muHat[2], const LocalReal br,
	      LocalReal & sr, LocalReal & si, LocalReal & kr, LocalReal & ki  );

int 
main(int argc, char *argv[]) 
{
  Overture::start(argc,argv);  // initialize Overture and A++/P++

  printF("Usage: texact -option=[0|1|2|3] -g=<grid>  -dm=[none|gdm]...\n"
         "  option =  0 : pec sphere\n"
         "            1 : dispersive sphere\n"
         "            2 : plane interface 2d\n"
         "            3 : plane interface 3d\n"
          );

  const int pecSphere=0, dispersiveSphere=1, planeInterface2D=2, planeInterface3D=3;
  int option = pecSphere;  
  option = dispersiveSphere;

  option=planeInterface2D;

  aString nameOfGridFile="";

  int plotOption=true;
  const int none=0, gdm=1;
  int dispersionModel =none;
  
  // char buff[180];
  int len=0;
  if( argc > 1 )
  { // look at arguments for "-noplot" or "-cfl=<value>"
    aString line;
    for( int i=1; i<argc; i++ )
    {
      line=argv[i];
      if( line=="-noplot" || line=="noplot" )
        plotOption=false;
      else if( len=line.matches("-option=") )
      {
        sScanF(line(len,line.length()-1),"%i",&option);
	printF("option = %i\n",option);
      }
      else if( len=line.matches("-g=") )
      {
        nameOfGridFile=line(len,line.length()-1);
	printF("nameOfGridFile=[%s]\n",(const char*)nameOfGridFile);
      }
      else if( len=line.matches("-dm=gdm") )
      {
	dispersionModel=gdm;
        printF("Setting dispersionModel=gdm\n");
      }
     
    }
  }





  if( nameOfGridFile == "" )
  {
    // -- default grid names ---
    if( option==pecSphere )
      nameOfGridFile="sphereInABoxGride2.order2.hdf";
    else if( option==dispersiveSphere )
      nameOfGridFile="solidSphereInABoxe2.order2.hdf";
    else if( option==planeInterface2D )
      nameOfGridFile="twoSquaresInterfacee8.order2.hdf";
    else if( option==planeInterface3D )
      nameOfGridFile= "twoBoxesInterfacee4.order2.hdf";
    else
    {
      OV_ABORT("ERROR: unknown option");
    }
  }
  


  CompositeGrid cg;
  getFromADataBase(cg,nameOfGridFile); // read in the grid with the default loadBalancer

  cg.update(MappedGrid::THEmask | GridCollection::THEdomain); // create CG's for each domain
    

  real a=1.;  // radius

  real k=1;  // wave number 
  real kr = twoPi*k, ki=0.;
  real sr=0., si = -kr;  // s = sr + I*si 
  
  if( option==dispersiveSphere )
  {
    sr=.1;  // TEST complex s 
  }

  real sc[2]={sr,si};   // real and imaginary parts

  real kc[2]={kr,ki};   // real and imaginary parts
  
  // Outside: 
  real eps0[2]={1.,0.1};   // real and imaginary parts
  real mu0[2] ={1.,0.};   // real and imaginary parts

  // Inside: 
  // real eps1[2]={4.,0.};   // real and imaginary parts
  real eps1[2]={.25,0.};   // real and imaginary parts
  real mu1[2] ={1.,0.};   // real and imaginary parts


  if( option==pecSphere || option==dispersiveSphere )
  {
    SphereExactSolutions ses;
    

    // --- Scattering from a PEC sphere currently assumes beta = real
    // 
    //      beta^2 = (m k)^2 = - s^2 epsHat muHat = br^2  = REAL 
    // Then s^2 = kr/sqrt( -epsHat*muHat )
    if( option==pecSphere )
    {
      real betar = kr;

      computeS( eps0, mu0, betar, sr,si, kr,ki );
      sc[0]=sr; sc[1]=si;
      kc[0]=kr; kc[1]=ki;
    
      printF("PEC sphere: s=[%g,%g]\n",sr,si);

    }
  
    
    ses.initialize( cg.numberOfDomains(), a, sc,kc, eps0, mu0, eps1, mu1 );

    if( true )
    {
      real x[3] ={1.,1.,1.}; // 
      real Ev[6], Hv[6];
      ses.eval( x,Ev,Hv );
      printF("\n RESULT: x=(%g,%g,%g) Ev=(%g,%g,%g,%g,%g,%g)\n",x[0],x[1],x[2],Ev[0],Ev[1],Ev[2],Ev[3],Ev[4],Ev[5]);
      printF("                        Hv=(%g,%g,%g,%g,%g,%g)\n",               Hv[0],Hv[1],Hv[2],Hv[3],Hv[4],Hv[5]);
    }
    if( true )
    {
      real x[3] ={.5,.5,.5}; // 
      real Ev[6], Hv[6];
      ses.eval( x,Ev,Hv );
      printF("\n RESULT: x=(%g,%g,%g) Ev=(%g,%g,%g,%g,%g,%g)\n",x[0],x[1],x[2],Ev[0],Ev[1],Ev[2],Ev[3],Ev[4],Ev[5]);
      printF("                        Hv=(%g,%g,%g,%g,%g,%g)\n",               Hv[0],Hv[1],Hv[2],Hv[3],Hv[4],Hv[5]);
  
    }
  
    // --- check the solution ---
    ses.check();

    Range all;
    int nc=6;  // real and imaginary parts of [Ex,Ey,Ez] 
    realCompositeGridFunction u(cg,all,all,all,nc);
    u=0.;
  
    u.setName("E");
    u.setName("Exr",0);
    u.setName("Eyr",1);
    u.setName("Ezr",2);

    u.setName("Exi",3);
    u.setName("Eyi",4);
    u.setName("Ezi",5);
  
    u=0.;

    real time0=getCPU();
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      ses.eval( u[grid], cg.domainNumber(grid),grid  );
    }
    real time=getCPU()-time0;
    printF("Time to evaluate on the grid = %10.2e (s)\n",time);

    // GenericGraphicsInterface & gi = *Overture::getGraphicsInterface("texact",plotOption,argc,argv);

    bool openGraphicsWindow=TRUE;
    PlotStuff gi(openGraphicsWindow,"texact");  // create a PlotStuff object
    PlotStuffParameters psp;                      // This object is used to change plotting parameters
  

    psp.set(GI_TOP_LABEL,sPrintF("Dielectric Sphere"));
    // PlotIt::plot(gi,cg,psp);

    PlotIt::contour(gi,u,psp);


  }
  else
  {
    // --------------- PLANE MATERIAL INTERFACE ---------------
    
     PlaneInterfaceExactSolution pes;
   
     DispersiveMaterialParameters dmp1, dmp2;
     real kvr[3]={0.,0.,0.};  //  set below 
     real kvi[3]={0.,0.,0.};  // 

     // We need av.kv = 0
     real av[3]={0.,0.,0};  // 
     if( cg.numberOfDimensions()==2 )
     {
       kvr[0]=1.;
       kvr[1]=.5;
       real kNorm=sqrt( SQR(kvr[0]) + SQR(kvr[1]) + SQR(kvr[2]) );
       av[0]=-kvr[1]/kNorm;  // do this for now
       av[1]= kvr[0]/kNorm;
     }
     else
     {
       kvr[0]=1.;
       kvr[1]=.5;
       kvr[2]=.5;
       real bv[3]= {1.,-1.,1.};
// TEST rotated 
       // kvr[0]=1.;
       // kvr[1]=1.;
       // kvr[2]=.0;
       // bv[0]=0.;
       // bv[1]=1;
       // bv[2]=1;
       // a = k X b
       real kDotB= kvr[0]*bv[0]+kvr[1]*bv[1]+kvr[2]*bv[2];
       assert( fabs(kDotB)>1.e-10 );
       
       av[0]= (kvr[1]*bv[2]-kvr[2]*bv[1])/kDotB;
       av[1]= (kvr[2]*bv[0]-kvr[0]*bv[2])/kDotB;
       av[2]= (kvr[0]*bv[1]-kvr[1]*bv[0])/kDotB;
     }
     

     real eps1=1., mu1=1.;
     real eps2=1., mu2=1.;
     eps2=2.;
      
  
     dmp1.setEpsInf( eps1 );
     dmp2.setEpsInf( eps2 );

     if( dispersionModel==gdm )
     {
       real a0=1., a1=.1, b0=1.5, b1=.05;
       int eqn=0;
       if( false )
       {
         int numPolarizationVectors1=1;
	 dmp1.setNumberOfPolarizationVectors( numPolarizationVectors1 );
	 dmp1.setParameters( eqn,a0,a1,b0,b1 );       
       }
       
       dmp2.setParameters( eqn,a0,a1,b0,b1 );       
       int numPolarizationVectors2=1;
       dmp2.setNumberOfPolarizationVectors( numPolarizationVectors2 );
       // a0=1., a1=.1, b0=1.5, b1=.05;
       eqn=0;
       dmp2.setParameters( eqn,a0,a1,b0,b1 );       
     }
     


     pes.initialize( cg, dmp1,dmp2,av,kvr,kvi );


     Range all;
     int nc= cg.numberOfDimensions()==2 ? 3 : 6;  //  [Ex,Ey,Hz] or [Ex,Ey,Ez, Hx,Hy,Hz ] 
     realCompositeGridFunction u(cg,all,all,all,nc);
     u=0.;
  
     u.setName("E");
     u.setName("Ex",0);
     u.setName("Ey",1);
     if( cg.numberOfDimensions()==2 )
       u.setName("Hz",2);
     else
     {
       u.setName("Ez",2);
       u.setName("Hx",3);
       u.setName("Hy",4);
       u.setName("Hz",5);
     }
     
     int npvMax=4;
     realCompositeGridFunction pv(cg,all,all,all,npvMax);

     


    int plotOption=1;
    // GenericGraphicsInterface & gi = *Overture::getGraphicsInterface("texact",plotOption,argc,argv);

    bool openGraphicsWindow=TRUE;
    PlotStuff gi(openGraphicsWindow,"texact");  // create a PlotStuff object
    PlotStuffParameters psp;                      // This object is used to change plotting parameters
  

    real dt=.1;
    const int numSteps=11;
    for( int step=0; step<numSteps; step++ )
    {
      real t = dt*step;
     
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
	Index I1,I2,I3;
	getIndex(cg[grid].dimension(),I1,I2,I3);

	bool computeMagneticField=true;
	pes.eval( t, cg, grid, u[grid], pv[grid],I1,I2,I3, 0, computeMagneticField);

      }
     
      gi.erase();
      psp.set(GI_TOP_LABEL,sPrintF("Plane Material Interface t=%g",t));
      // PlotIt::plot(gi,cg,psp);
      if( step==0 || step==(numSteps-1) )
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);      
      else
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);      

      PlotIt::contour(gi,u,psp);

      gi.redraw(TRUE);
      
    }
     

  }
  




  // realCompositeGridFunction v(cg,all,all,all,1);
  // v.setName("E",0);
  // v=0.;
  // PlotIt::contour(gi,v,psp);


  Overture::finish(); 
  return 0;
  

  
}



// Include complex down here to minimize name conflicts
#include <complex>


// ====================================================================================
/// \brief Given a real value for beta (PEC sphere case), compute complex k and s
// 
// Compute (sr,si) = beta/sqrt(eps*mu) (Complex) 
//   beta^2 = - s^2 eps*mu
//   beta = k*m
// 
//   m = sqrt(eps*mu) 
//   s = +- beta/sqrt( -eps*mu )
//   k = beta/m; 
// ====================================================================================
int computeS( const LocalReal epsHat[2], const LocalReal muHat[2], const LocalReal beta,
	      LocalReal & sr, LocalReal & si,  LocalReal & kr, LocalReal & ki  )
{
  std::complex<LocalReal> m,eps,mu,s,k;
   
  eps = epsHat[0] + 1i*epsHat[1]; // complex eps 
  mu  = muHat[0]  + 1i*muHat[1];
  
  s = beta/sqrt(-eps*mu);  
  if( std::imag(s)>0. ){ s=-s; } // choose root with si < 0 
   
  sr = std::real(s);
  si = std::imag(s);

  m = sqrt( eps*mu );  // complex index of refraction 
  k = beta/m;
  
  kr = std::real(k);
  ki = std::imag(k);
}
