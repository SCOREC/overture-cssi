// ***********************************************************
// *** This macro defines the trigonometric TZ functions *****
// ***********************************************************
#beginMacro  defineTrigonometricTZMacro()     

// const int nc = numberOfComponents + int(useChargeDensity);  // include charge density in TZ
const int nc = numberOfComponentsForTZ;

RealArray fx(nc),fy(nc),fz(nc),ft(nc);
RealArray gx(nc),gy(nc),gz(nc),gt(nc);
gx=0.;
gy=0.;
gz=0.;
gt=0.;
RealArray amplitude(nc), cc(nc);
amplitude=1.;
cc=0.;


fx=omega[0];
fy = numberOfDimensions>1 ? omega[1] : 0.;
fz = numberOfDimensions>2 ? omega[2] : 0.;
ft = omega[3];



if( numberOfDimensions==2  )
{   

  if( omega[0] != omega[1] ) 
  {
    printF("Cgmx:Trig TZ: invalid values for omega: omega[0]=%9.3e, omega[1]=%9.3e \n",omega[0],omega[1]);
    printF("Expecting omega[0]==omega[1] (for divergence free field) \n");

    OV_ABORT("Invalid values for omega[0..1]");
  }

  const int uc=ex, vc=ey, wc=hz;


  if( solveForAllFields )
  {
    // u=cos(pi x) cos( pi y )* .5
    // v=sin(pi x) sin( pi y )* .5 
    // w=cos(    ) sin(      )* .5
    assert( omega[0]==omega[1] );

    gx(ey)=.5/omega[0];   // shift by pi/2 to turn cos() into sin()
    gy(ey)=.5/omega[1];

    amplitude(ex)=.5;  cc(ex)=.0;
    amplitude(ey)=.5;  cc(ey)=.0;

    gy(hz)=.5/omega[1]; // turn off for testing symmetry
    cc(hz)=.0;



    gx(hy)=.5/omega[0];   // shift by pi/2 to turn cos() into sin()
    gy(hy)=.5/omega[1];

    amplitude(hx)=.5;  cc(hx)=.0;
    amplitude(hy)=.5;  cc(hy)=.0;

    gy(ez)=.5/omega[1]; // turn off for testing symmetry
    cc(ez)=.0;

  }
  else if( !useChargeDensity )
  {
    // rho=0 : create div(E)=0 TZ function
    
    // u=cos(pi x) cos( pi y )* .5
    // v=sin(pi x) sin( pi y )* .5 
    // w=cos(    ) sin(      )* .5
    assert( omega[0]==omega[1] );

    gx(vc)=.5/omega[0];   // shift by pi/2 to turn cos() into sin()
    gy(vc)=.5/omega[1];

    amplitude(uc)=.5;  cc(uc)=.0;
    amplitude(vc)=.5;  cc(vc)=.0;

    gy(wc)=.5/omega[1]; // turn off for testing symmetry
    cc(wc)=.0;
  }
  else
  {
    // rho= cos(pi x) cos( pi y )* omega[0]*pi 
    // 
    // u=sin(pi x) cos( pi y )* .5
    // v=cos(pi x) sin( pi y )* .5 
    // w=cos(    ) sin(      )* .5
    assert( omega[0]==omega[1] );

    assert( rc==numberOfComponents );

    gx(uc)=.5/omega[0];   // shift by pi/2 to turn cos() into sin()
    gy(vc)=.5/omega[1];

    amplitude(rc)=omega[0]*Pi;  cc(rc)=.0;
    amplitude(uc)=.5;           cc(uc)=.0;
    amplitude(vc)=.5;           cc(vc)=.0;

    gy(wc)=.5/omega[1]; // turn off for testing symmetry
    cc(wc)=.0;

  }
  
}
else if( numberOfDimensions==3 )
{
  if( solveForElectricField )
  {
    const int uc=ex, vc=ey, wc=ez;

    // u=   cos(pi x) cos( pi y ) cos( pi z)  // **** fix ***
    // v=.5 sin(pi x) sin( pi y ) cos( pi z)
    // w=.5 sin(pi x) cos( pi y ) sin( pi z)
	
    if( omega[0]==omega[1] && omega[0]==omega[2] )
    {
      gx(vc)=.5/omega[0];
      gy(vc)=.5/omega[1];
      amplitude(vc)=.5;
	
      gx(wc)=.5/omega[0];
      gz(wc)=.5/omega[2];
      amplitude(wc)=.5;
    }
    else if( omega[0]==omega[2] && omega[1]==0 )
    {
      // pseudo 2D case
      gx(wc)=.5/omega[0];   // shift by pi/2 to turn cos() into sin()
      gz(wc)=.5/omega[2];

      amplitude(uc)=.5;  cc(uc)=.0;
      amplitude(wc)=.5;  cc(wc)=.0;
    }
    else if( omega[0]==omega[1] && omega[2]==0 )
    {
      // pseudo 2D case

      // u=cos(pi x) cos( pi y ) cos( 0*pi z)* .5
      // v=sin(pi x) sin( pi y ) cos( 0*pi z)* .5 
      // w=cos(pi x) cos( pi y ) cos( 0*pi z)* .5


      gx(vc)=.5/omega[0];   // shift by pi/2 to turn cos() into sin()
      gy(vc)=.5/omega[1];

      amplitude(uc)=.5;  cc(uc)=.0;
      amplitude(vc)=.5;  cc(vc)=.0;
      amplitude(wc)=.5;  cc(wc)=.0;
    }
    else
    {
      printF("Cgmx: invalid values for omega: omega[0]=%9.3e, omega[1]=%9.3e, omega[2]=%9.3e\n",
	     omega[0],omega[1],omega[2]);
      printF("Expecting all equal values for omega[0]==omega[2] && omega[1]==0 \n"
             " or omega[0]==omega[1] && omega[2]==0 (for divergence free field)\n");
      OV_ABORT("Invalid values for omega[0..2]");
    }
  }
  
  if( solveForMagneticField )
    {
      const int uc=hx, vc=hy, wc=hz;

    // u=   cos(pi x) cos( pi y ) cos( pi z)  // **** fix ***
    // v=.5 sin(pi x) sin( pi y ) cos( pi z)
    // w=.5 sin(pi x) cos( pi y ) sin( pi z)
	
    if( omega[0]==omega[1] && omega[0]==omega[2] )
    {
      gx(vc)=.5/omega[0];
      gy(vc)=.5/omega[1];
      amplitude(vc)=.5;
	
      gx(wc)=.5/omega[0];
      gz(wc)=.5/omega[2];
      amplitude(wc)=.5;
    }
    else if( omega[0]==omega[2] && omega[1]==0 )
    {
      // pseudo 2D case
      gx(wc)=.5/omega[0];   // shift by pi/2 to turn cos() into sin()
      gz(wc)=.5/omega[2];

      amplitude(uc)=.5;  cc(uc)=.0;
      amplitude(wc)=.5;  cc(wc)=.0;
    }
    else
    {
      Overture::abort("Invalid values for omega[0..2]");
    }
  }
	
}
if( method==sosup )
{
  // Set the TZ function for (ext,eyt,...) equal to the time derivative of (ex,ey,...)

  // time dependence for time-derivatives of E:
  const int numberOfFieldComponents=3;
  for( int n=ex, nt=ext; n<ex+numberOfFieldComponents; n++,nt++ )
  {
    fx(nt)=fx(n); fy(nt)=fy(n); fz(nt)=fz(n); ft(nt)=ft(n);
    gx(nt)=gx(n); gy(nt)=gy(n); gz(nt)=gz(n); 
    gt(nt)=.5/ft(n);  // shift phase by pi/2 to turn cos(Pi*ft*(t-gt)) into sin(Pi*ft*(t-gt))
    amplitude(nt)=-amplitude(n)*ft(n)*Pi;  // amplitude
    cc(nt)=0.; 
  }
  

//   fx(eyt)=fx(ey); gx(eyt)=gx(ey); amplitude(eyt)=-amplitude(ey);  cc(eyt)=0.; ft(eyt)=ft(ey); gt(eyt)=.5/omega[3]; 

//   fx(hzt)=fx(hz); gx(hzt)=gx(hz); amplitude(hzt)=-amplitude(hz);  cc(hzt)=0.; ft(hzt)=ft(hz); gt(hzt)=.5/omega[3]; 

}

int lastTZComponent = hz+1;  // keep track of the last TZ component assigned above so we can fill in MLA components
if( dispersionModel != noDispersion )
{
  // -- dispersion components:  
  int numberOfPolarizationComponents=maxNumberOfPolarizationVectors;  // ** CHECK ME ***
  
  lastTZComponent += numberOfPolarizationComponents*numberOfDimensions;
  if( method==bamx )
    lastTZComponent += numberOfPolarizationComponents*numberOfDimensions;

  
  if( method==nfdtd )
  {
     printF("\n >>>> set TRIG TZ for P maxNumberOfPolarizationVectors=%i <<<<\n\n",maxNumberOfPolarizationVectors);

    for( int iv=0; iv<maxNumberOfPolarizationVectors; iv++ )
    {
      const int pc= iv*numberOfDimensions;
      // Note: Corner extrapolation may assume that div(E)=0 
      for( int dir=0; dir<numberOfDimensions; dir++ )
      {
        int ptz = pxc + pc + dir;
        amplitude(ptz)= 1 + 1./(2.+iv);
        // make P be divergence free like E
        gx(ptz) = gx(ex+dir);
        gy(ptz) = gy(ex+dir);
        ft(ptz) = omega[3]*( 1. + 1./(3. + iv) );
      }
    }
  }
  else if( method==bamx )
  {
     printF("\n >>>> set TRIG TZ for P BAMX  <<<<\n\n");
     // BA Maxwell -- assign polarization components
  
    const int numPolarizationTerms = 2*maxNumberOfPolarizationComponents;  // note "2*" we save p and p.t 
    for( int m=0; m<numPolarizationTerms; m++ )
    {
      int pc = hz+m+1; // polarization component index in TZ functions
      amplitude(pc) = 1./( 1 + 2.*pc/numPolarizationTerms );
      fx(pc) = omega[0]*( (m+1.)/(m+3.) );
      fy(pc) = omega[1]*( (m+2.)/(m+3.) );
      gx(pc) = 0.5/omega[0]*( 2*(m-2.)/(m+4.) );
      gy(pc) = .25/omega[1]*(   (m-2.)/(m+4.) );
      ft(pc) = omega[3]*( 1. + (1.+m)/(3. + m) );
        
    }
  }
 
}

if( dispersionModel != noDispersion &&  nonlinearModel==multilevelAtomic )
{
  // Nonlinear model : multilevelAtomic 
  printF("\n >>> INIT TRIG TZ: Fill in nonlinear TZ variables starting at lastTZComponent=%d <<<\n",lastTZComponent);


  const int & maxNumberOfNonlinearVectors = parameters.dbase.get<int>("maxNumberOfNonlinearVectors");
  for( int m=0; m<maxNumberOfNonlinearVectors; m++ )
  {
    int na = lastTZComponent+m; 
    amplitude(na) = (1.+m)/( 2.+m );
    fx(na) = omega[0]*( 1 + (1.0*m/maxNumberOfNonlinearVectors) );
    fy(na) = omega[1]*( 1 + (0.5*m/maxNumberOfNonlinearVectors) );;
    gx(na) = (1./omega[0])*( (m-1.)/maxNumberOfNonlinearVectors );
    gy(na) = (1./omega[1])*( (m-2.)/maxNumberOfNonlinearVectors );
    cc(na) = 1./(2.+m);
    
  }
   
  
  
}


tz = new OGTrigFunction(fx,fy,fz,ft);
    
((OGTrigFunction*)tz)->setShifts(gx,gy,gz,gt);
((OGTrigFunction*)tz)->setAmplitudes(amplitude);
((OGTrigFunction*)tz)->setConstants(cc);

#endMacro
