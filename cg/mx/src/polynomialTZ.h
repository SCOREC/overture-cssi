
// ===============================================================================================
//  Macro: initialize3DPolyTW is a bpp macro that sets up the E and H polynomial tw function
// ===============================================================================================
#beginMacro initialize3DPolyTW(ux,uy,uz)

// -----------------------------------------------------------------------------
// --------------------- DEFINE POLYNOMIAL TZ SOLUTIONS ------------------------
// -----------------------------------------------------------------------------

// Always include linear terms in TZ if degreSpace>=1 *wdh* Sept 18, 2016 
if( degreeSpace >=1 )
{
  spatialCoefficientsForTZ(0,0,0,ux)=1.;      // u=1 + x + y + z
  spatialCoefficientsForTZ(1,0,0,ux)=1.;
  spatialCoefficientsForTZ(0,1,0,ux)=1.;
  spatialCoefficientsForTZ(0,0,1,ux)=1.;

  spatialCoefficientsForTZ(0,0,0,uy)= 2.;      // v=2+x-2y+z
  spatialCoefficientsForTZ(1,0,0,uy)= 1.;
  spatialCoefficientsForTZ(0,1,0,uy)=-2.;
  spatialCoefficientsForTZ(0,0,1,uy)= 1.;
    
  spatialCoefficientsForTZ(1,0,0,uz)=-1.;      // w=-x+y+z
  spatialCoefficientsForTZ(0,1,0,uz)= 1.;
  spatialCoefficientsForTZ(0,0,1,uz)= 1.;

  // eps and mu should remain positive 
  if( epsc>=0 )
  {
    spatialCoefficientsForTZ(1,0,0,epsc)=eps*.01;  // x
    spatialCoefficientsForTZ(0,1,0,epsc)=eps*.02;  // y
    spatialCoefficientsForTZ(0,0,1,epsc)=eps*.12;  // z
  }
  if( muc>=0 )
  {
    spatialCoefficientsForTZ(1,0,0,muc )=mu*.015;   // x
    spatialCoefficientsForTZ(0,1,0,muc )=mu*.0125;  // y
    spatialCoefficientsForTZ(0,0,1,muc )=mu*.095;   // z
  }
  
}

if( degreeSpace==2 )
{
  spatialCoefficientsForTZ(2,0,0,ux)=1.;      // u=x^2 + 2xy + y^2 + xz  - .25*yz -.5*z^2
  spatialCoefficientsForTZ(1,1,0,ux)=2.;
  spatialCoefficientsForTZ(0,2,0,ux)=1.;
  spatialCoefficientsForTZ(1,0,1,ux)=1.;
  spatialCoefficientsForTZ(0,1,1,ux)=-.25;
  spatialCoefficientsForTZ(0,0,2,ux)=-.5;
      
  spatialCoefficientsForTZ(2,0,0,uy)= 1.;      // v=x^2 -2xy - y^2 + 3yz + .25*xz +.5*z^2
  spatialCoefficientsForTZ(1,1,0,uy)=-2.;
  spatialCoefficientsForTZ(0,2,0,uy)=-1.;
  spatialCoefficientsForTZ(0,1,1,uy)=+3.;
  spatialCoefficientsForTZ(1,0,1,uy)=.25;
  spatialCoefficientsForTZ(0,0,2,uy)=.5;
      
  spatialCoefficientsForTZ(2,0,0,uz)= 1.;      // w=x^2 + y^2 - 2 z^2 + .25*xy 
  spatialCoefficientsForTZ(0,2,0,uz)= 1.;
  spatialCoefficientsForTZ(0,0,2,uz)=-2.;
  spatialCoefficientsForTZ(1,1,0,uz)=.25;

  // eps and mu should remain positive 
  if( epsc>=0 )
  {
    spatialCoefficientsForTZ(1,0,0,epsc)=eps*.01;  // x
    spatialCoefficientsForTZ(0,1,0,epsc)=eps*.02;  // y
    spatialCoefficientsForTZ(0,0,1,epsc)=eps*.12;  // z
    spatialCoefficientsForTZ(2,0,0,epsc)=eps*.1;   // x^2
    spatialCoefficientsForTZ(0,2,0,epsc)=eps*.15;  // y^2        
    spatialCoefficientsForTZ(0,0,2,epsc)=eps*.11;  // z^2        
  }
  if( muc>=0 )
  {
    spatialCoefficientsForTZ(1,0,0,muc )=mu*.015;   // x
    spatialCoefficientsForTZ(0,1,0,muc )=mu*.0125;  // y
    spatialCoefficientsForTZ(0,0,1,muc )=mu*.095;   // z
    spatialCoefficientsForTZ(2,0,0,muc )=mu*.125;   // x^2
    spatialCoefficientsForTZ(0,2,0,muc )=mu*.15;    // y^2
    spatialCoefficientsForTZ(0,0,2,muc )=mu*.13;    // z^2
  }
  
}
else if( degreeSpace==0 )
{
  spatialCoefficientsForTZ(0,0,0,ux)=1.; // -1.; 
  spatialCoefficientsForTZ(0,0,0,uy)=1.; //-.5;
  spatialCoefficientsForTZ(0,0,0,uz)=1.; //.75; 
}
else if( degreeSpace==3 )
{
  spatialCoefficientsForTZ(2,0,0,ux)=1.;      // u=x^2 + 2xy + y^2 + xz 
  spatialCoefficientsForTZ(1,1,0,ux)=2.;    //        + .125( x^3 + y^3 + z^3 ) -.75*x*y^2 + x^2*z +.4yz
  spatialCoefficientsForTZ(0,2,0,ux)=1.;
  spatialCoefficientsForTZ(1,0,1,ux)=1.;
      
  spatialCoefficientsForTZ(3,0,0,ux)=.125; 
  spatialCoefficientsForTZ(0,3,0,ux)=.125; 
  spatialCoefficientsForTZ(0,0,3,ux)=.125; 
  spatialCoefficientsForTZ(1,2,0,ux)=-.75;
  spatialCoefficientsForTZ(2,0,1,ux)=+1.; 
  spatialCoefficientsForTZ(0,1,1,ux)=.4; 


  spatialCoefficientsForTZ(2,0,0,uy)= 1.;      // v=x^2 -2xy - y^2 + 3yz 
  spatialCoefficientsForTZ(1,1,0,uy)=-2.;      //    + .25( x^3 + y^3 + z^3 ) -.375*x^2 y  -.375*y*z^2  
  spatialCoefficientsForTZ(0,2,0,uy)=-1.;
  spatialCoefficientsForTZ(0,1,1,uy)=+3.;
      
  spatialCoefficientsForTZ(3,0,0,uy)=.25; 
  spatialCoefficientsForTZ(0,3,0,uy)=.25; 
  spatialCoefficientsForTZ(0,0,3,uy)=.25; 
  spatialCoefficientsForTZ(2,1,0,uy)=-3.*.125; 
  spatialCoefficientsForTZ(0,1,2,uy)=-3.*.125; 
      
      
  spatialCoefficientsForTZ(2,0,0,uz)= 1.;      // w=x^2 + y^2 - 2 z^2 
  spatialCoefficientsForTZ(0,2,0,uz)= 1.;      //      + .25x^3 -.2y^3 +.125 z^3 - x z^2 -.6*xy^2
  spatialCoefficientsForTZ(0,0,2,uz)=-2.;
      
  spatialCoefficientsForTZ(3,0,0,uz)=.25; 
  spatialCoefficientsForTZ(0,3,0,uz)=-.2; 
  spatialCoefficientsForTZ(0,0,3,uz)=.125; 
  spatialCoefficientsForTZ(1,0,2,uz)=-1.;
  spatialCoefficientsForTZ(1,2,0,uz)=-.6;

}
else if( degreeSpace==4 )
{
  spatialCoefficientsForTZ(2,0,0,ux)=1.;      // u=x^2 + 2xy + y^2 + xz
  spatialCoefficientsForTZ(1,1,0,ux)=2.;
  spatialCoefficientsForTZ(0,2,0,ux)=1.;
  spatialCoefficientsForTZ(1,0,1,ux)=1.;
  spatialCoefficientsForTZ(3,0,0,ux)=.5;      // + .5*x^3

  spatialCoefficientsForTZ(4,0,0,ux)=.125;    // + .125*x^4 + .125*y^4 + .125*z^4  -.5*xz^3
  spatialCoefficientsForTZ(0,4,0,ux)=.125;    
  spatialCoefficientsForTZ(0,0,4,ux)=.125; 
  spatialCoefficientsForTZ(1,0,3,ux)=-.5; 
  spatialCoefficientsForTZ(0,1,3,ux)=.25;    // + .25*y*z^3 -.25*y^2*z^2 +.25*y^3z
  spatialCoefficientsForTZ(0,2,2,ux)=-.25; 
  spatialCoefficientsForTZ(0,3,1,ux)=.25; 
      
      
  spatialCoefficientsForTZ(2,0,0,uy)= 1.;      // v=x^2 -2xy - y^2 + 3yz
  spatialCoefficientsForTZ(1,1,0,uy)=-2.;
  spatialCoefficientsForTZ(0,2,0,uy)=-1.;
  spatialCoefficientsForTZ(0,1,1,uy)=+3.;
      
  spatialCoefficientsForTZ(2,1,0,uy)=-1.5;     // -1.5x^2*y
      
  spatialCoefficientsForTZ(4,0,0,uy)=.25; 
  spatialCoefficientsForTZ(0,4,0,uy)=.25; 
  spatialCoefficientsForTZ(0,0,4,uy)=.25; 
  spatialCoefficientsForTZ(3,1,0,uy)=-.5; 
  spatialCoefficientsForTZ(1,0,3,uy)=.25;    // + .25*x*z^3 -.25*x^2*z^2 +.25*x^3z
  spatialCoefficientsForTZ(2,0,2,uy)=-.25; 
  spatialCoefficientsForTZ(3,0,1,uy)=.25; 
      
      
  spatialCoefficientsForTZ(2,0,0,uz)= 1.;      // w=x^2 + y^2 - 2 z^2
  spatialCoefficientsForTZ(0,2,0,uz)= 1.;
  spatialCoefficientsForTZ(0,0,2,uz)=-2.;
      
  spatialCoefficientsForTZ(4,0,0,uz)=.25; 
  spatialCoefficientsForTZ(0,4,0,uz)=-.2; 
  spatialCoefficientsForTZ(0,0,4,uz)=.125; 
  spatialCoefficientsForTZ(0,3,1,uz)=-1.;
  spatialCoefficientsForTZ(1,3,0,uz)=.25;    // + .25*x*y^3 -.25*x^2*y^2 +.25*x^3y
  spatialCoefficientsForTZ(2,2,0,uz)=-.25; 
  spatialCoefficientsForTZ(3,1,0,uz)=.25; 
}
else if( degreeSpace>=5 )
{
  if( true || degreeSpace!=5 ) printF(" ****WARNING***** using a TZ function with degree=5 in space *****\n");
	  
  spatialCoefficientsForTZ(2,0,0,ux)=1.;      // u=x^2 + 2xy + y^2 + xz
  spatialCoefficientsForTZ(1,1,0,ux)=2.;
  spatialCoefficientsForTZ(0,2,0,ux)=1.;
  spatialCoefficientsForTZ(1,0,1,ux)=1.;
    
  spatialCoefficientsForTZ(4,0,0,ux)=.125;    // + .125*x^4 + .125*y^4 + .125*z^4  -.5*xz^3
  spatialCoefficientsForTZ(0,4,0,ux)=.125;    
  spatialCoefficientsForTZ(0,0,4,ux)=.125; 
  spatialCoefficientsForTZ(1,0,3,ux)=-.5; 
  spatialCoefficientsForTZ(0,1,3,ux)=.25;    // + .25*y*z^3 -.25*y^2*z^2 +.25*y^3z
  spatialCoefficientsForTZ(0,2,2,ux)=-.25; 
  spatialCoefficientsForTZ(0,3,1,ux)=.25; 
    
  spatialCoefficientsForTZ(0,5,0,ux)=.125;   // y^5
    
    
  spatialCoefficientsForTZ(2,0,0,uy)= 1.;      // v=x^2 -2xy - y^2 + 3yz
  spatialCoefficientsForTZ(1,1,0,uy)=-2.;
  spatialCoefficientsForTZ(0,2,0,uy)=-1.;
  spatialCoefficientsForTZ(0,1,1,uy)=+3.;
    
  spatialCoefficientsForTZ(4,0,0,uy)=.25; 
  spatialCoefficientsForTZ(0,4,0,uy)=.25; 
  spatialCoefficientsForTZ(0,0,4,uy)=.25; 
  spatialCoefficientsForTZ(3,1,0,uy)=-.5; 
  spatialCoefficientsForTZ(1,0,3,uy)=.25;    // + .25*x*z^3 -.25*x^2*z^2 +.25*x^3z
  spatialCoefficientsForTZ(2,0,2,uy)=-.25; 
  spatialCoefficientsForTZ(3,0,1,uy)=.25; 
    
  // spatialCoefficientsForTZ(5,0,0,uy)=.125;  // x^5
    
    
  spatialCoefficientsForTZ(2,0,0,uz)= 1.;      // w=x^2 + y^2 - 2 z^2
  spatialCoefficientsForTZ(0,2,0,uz)= 1.;
  spatialCoefficientsForTZ(0,0,2,uz)=-2.;
    
  spatialCoefficientsForTZ(4,0,0,uz)=.25; 
  spatialCoefficientsForTZ(0,4,0,uz)=-.2; 
  spatialCoefficientsForTZ(0,0,4,uz)=.125; 
  spatialCoefficientsForTZ(0,3,1,uz)=-1.;
  spatialCoefficientsForTZ(1,3,0,uz)=.25;    // + .25*x*y^3 -.25*x^2*y^2 +.25*x^3y
  spatialCoefficientsForTZ(2,2,0,uz)=-.25; 
  spatialCoefficientsForTZ(3,1,0,uz)=.25; 
    
  // spatialCoefficientsForTZ(5,0,0,uz)=.125;
}
else
{
  printF("Maxwell:: not implemented for degree in space =%i \n",degreeSpace);
  OV_ABORT("error");
}

// *** end the initialize3DPolyTW bpp macro 
#endMacro
  

// ===============================================================================================
// *** Macro: This macro defines the polynomial TZ functions 
// ===============================================================================================
#beginMacro definePolynomialTZMacro()

tz = new OGPolyFunction(degreeSpace,numberOfDimensions,numberOfComponentsForTZ,degreeTime);

const int ndp=max(max(5,degreeSpace+1),degreeTime+1);
    
printF("\n $$$$$$$ assignInitialConditions: build OGPolyFunction: numCompTz=%i degreeSpace=%i, degreeTime=%i ndp=%i $$$$\n",
       numberOfComponentsForTZ,degreeSpace,degreeTime,ndp);
if( pxc>=0 )
{
  printF(" pxc=%i pyc=%i pzc=%i\n",pxc,pyc,pzc);
}


RealArray spatialCoefficientsForTZ(ndp,ndp,ndp,numberOfComponentsForTZ);  
spatialCoefficientsForTZ=0.;
RealArray timeCoefficientsForTZ(ndp,numberOfComponentsForTZ);      
timeCoefficientsForTZ=0.;

// Default coefficients for eps, mu, sigmaE and sigmaH:
if( epsc>=0 )
{
  assert( epsc>=0 && muc>=0 && sigmaEc>=0 && sigmaHc>=0 );
  printF(" *** numberOfComponentsForTZ=%i, rc=%i, epsc,muc,sigmaEc,sigmaHc=%i,%i,%i,%i, eps,mu=%e,%e\n",
         numberOfComponentsForTZ,rc,epsc,muc,sigmaEc,sigmaHc,eps,mu);

  spatialCoefficientsForTZ(0,0,0,epsc)=eps;
  spatialCoefficientsForTZ(0,0,0,muc )=mu; 
  spatialCoefficientsForTZ(0,0,0,sigmaEc)=0.;  
  spatialCoefficientsForTZ(0,0,0,sigmaHc)=0.;
}


    
if( numberOfDimensions==2 )
{
  // ************************************************************************
  // ********************** TWO SPACE DIMENSIONS ****************************
  // ************************************************************************

  printF("TZ: solveForAllFields=%d, numberOfComponents=%d\n",solveForAllFields,numberOfComponents);

    
    

  if( degreeSpace==0 )
  {
    if( !solveForAllFields )
    {
      spatialCoefficientsForTZ(0,0,0,ex)=1.;      // u=1
      spatialCoefficientsForTZ(0,0,0,ey)= 2.;      // v=2
      spatialCoefficientsForTZ(0,0,0,hz)=-1.;      // w=-1
    }
    else
    {
      for( int c=0; c<numberOfComponents; c++ )
      {
        spatialCoefficientsForTZ(0,0,0,c)= 1. + c*.5;
      }
      
    }
    
    // -- dispersion components: 
    if( dispersionModel != noDispersion )
    {
      if( method==nfdtd )
      {
	for( int iv=0; iv<maxNumberOfPolarizationVectors; iv++ )
	{
	  const int pc= iv*numberOfDimensions;
	  spatialCoefficientsForTZ(0,0,0,pxc+pc)=1.; 
	  spatialCoefficientsForTZ(0,0,0,pyc+pc)=2.; 
	}
      }
      else if( method==bamx )
      {
	 // done below
      }
      
    }
    
  }
  else if( degreeSpace==1 )
  {
    if( !solveForAllFields )
    {
      spatialCoefficientsForTZ(0,0,0,ex)=1.;      // u=1+x+y
      spatialCoefficientsForTZ(1,0,0,ex)=1.;
      spatialCoefficientsForTZ(0,1,0,ex)=1.;

      spatialCoefficientsForTZ(0,0,0,ey)= 2.;      // v=2+x-y
      spatialCoefficientsForTZ(1,0,0,ey)= 1.;
      spatialCoefficientsForTZ(0,1,0,ey)=-1.;

      spatialCoefficientsForTZ(0,0,0,hz)=-1.;      // w=-1+x + y
      spatialCoefficientsForTZ(1,0,0,hz)= 1.;
      spatialCoefficientsForTZ(0,1,0,hz)= 1.;
    }
    else
    {
      for( int c=0; c<numberOfComponents; c++ )
      {
        spatialCoefficientsForTZ(0,0,0,c)= 1. + c;
        spatialCoefficientsForTZ(1,0,0,c)= .5 + c*.5;
        spatialCoefficientsForTZ(0,1,0,c)= .4 + c*.25;
      }
      // We need div(Ex,Ey,Ez)=0 , div(Hx,Hy,Hz)=0 
      spatialCoefficientsForTZ(0,1,0,ey)=- spatialCoefficientsForTZ(1,0,0,ex);
      spatialCoefficientsForTZ(0,1,0,hy)=- spatialCoefficientsForTZ(1,0,0,hx);

    }
    
    // eps and mu should remain positive but do this for now:
    if( epsc>=0 )
    {
      spatialCoefficientsForTZ(1,0,0,epsc)=eps*.01;  // x*eps*.01
      spatialCoefficientsForTZ(0,1,0,epsc)=eps*.02;  // y*eps*.02 

      spatialCoefficientsForTZ(1,0,0,muc )=mu*.015;   // x
      spatialCoefficientsForTZ(0,1,0,muc )=mu*.0125;  // y
    }
    
    // -- dispersion components: 
    if( dispersionModel != noDispersion )
    {
      for( int iv=0; iv<maxNumberOfPolarizationVectors; iv++ )
      {
        const int pc= iv*numberOfDimensions;
        spatialCoefficientsForTZ(0,0,0,pxc+pc)=1.; 
        spatialCoefficientsForTZ(1,0,0,pxc+pc)=.5;
        spatialCoefficientsForTZ(0,1,0,pxc+pc)=.4;

        spatialCoefficientsForTZ(0,0,0,pyc+pc)=2.; 
        spatialCoefficientsForTZ(1,0,0,pyc+pc)=-.5;
        spatialCoefficientsForTZ(0,1,0,pyc+pc)=.5;
      }
    }


  }
  else if( degreeSpace==2 )
  {
    if( !solveForAllFields )
    {

      spatialCoefficientsForTZ(2,0,0,ex)=1.;      // u=x^2 + 2xy + y^2 
      spatialCoefficientsForTZ(1,1,0,ex)=2.;
      spatialCoefficientsForTZ(0,2,0,ex)=1.;

      spatialCoefficientsForTZ(2,0,0,ey)= 1.;      // v=x^2 -2xy - y^2 
      spatialCoefficientsForTZ(1,1,0,ey)=-2.;
      spatialCoefficientsForTZ(0,2,0,ey)=-1.;

      spatialCoefficientsForTZ(2,0,0,hz)= 1.;      // w=x^2 + y^2 -1 +.5 xy
      spatialCoefficientsForTZ(0,2,0,hz)= 1.;
      spatialCoefficientsForTZ(0,0,0,hz)=-1.; 
      spatialCoefficientsForTZ(1,1,0,hz)= .5;
    }
    else
    {
      // We need div(Ex,Ey,Ez)=0 , div(Hx,Hy,Hz)=0 
      spatialCoefficientsForTZ(2,0,0,ex)=1.;      // Ex =x^2 + 2xy + y^2 
      spatialCoefficientsForTZ(1,1,0,ex)=2.;
      spatialCoefficientsForTZ(0,2,0,ex)=1.;

      spatialCoefficientsForTZ(2,0,0,ey)= 1.;      // Ey =x^2 -2xy - y^2 
      spatialCoefficientsForTZ(1,1,0,ey)=-2.;
      spatialCoefficientsForTZ(0,2,0,ey)=-1.;

      spatialCoefficientsForTZ(2,0,0,ez)= .5;    
      spatialCoefficientsForTZ(0,2,0,ez)=.25;
      spatialCoefficientsForTZ(0,0,0,ez)=-.5; 
      spatialCoefficientsForTZ(1,1,0,ez)= .3;

      spatialCoefficientsForTZ(2,0,0,hx)=.5;     
      spatialCoefficientsForTZ(1,1,0,hx)=1.;
      spatialCoefficientsForTZ(0,2,0,hx)=.5;

      spatialCoefficientsForTZ(2,0,0,hy)= .5;    
      spatialCoefficientsForTZ(1,1,0,hy)=-1.;
      spatialCoefficientsForTZ(0,2,0,hy)=-.5;

      spatialCoefficientsForTZ(2,0,0,hz)= 1.;      // Hz=x^2 + y^2 -1 +.5 xy
      spatialCoefficientsForTZ(0,2,0,hz)= 1.;
      spatialCoefficientsForTZ(0,0,0,hz)=-1.; 
      spatialCoefficientsForTZ(1,1,0,hz)= .5;

    }
    
    
    // eps and mu should remain positive 
    if( epsc>=0 )
    {
      spatialCoefficientsForTZ(1,0,0,epsc)=eps*.01;  // x
      spatialCoefficientsForTZ(0,1,0,epsc)=eps*.02;  // y
      spatialCoefficientsForTZ(2,0,0,epsc)=eps*.1;   // x^2
      spatialCoefficientsForTZ(0,2,0,epsc)=eps*.15;  // y^2        

      spatialCoefficientsForTZ(1,0,0,muc )=mu*.015;   // x
      spatialCoefficientsForTZ(0,1,0,muc )=mu*.0125;  // y
      spatialCoefficientsForTZ(2,0,0,muc )=mu*.125;   // x^2
      spatialCoefficientsForTZ(0,2,0,muc )=mu*.15;    // y^2
    }
    
    // -- dispersion components: 
    if( dispersionModel != noDispersion )
    {
      printF("\n >>>> set TZ for P maxNumberOfPolarizationVectors=%i <<<<\n\n",maxNumberOfPolarizationVectors);

      if( method==nfdtd )
      {
	for( int iv=0; iv<maxNumberOfPolarizationVectors; iv++ )
	{
	  const int pc= iv*numberOfDimensions;
	  // Corner extrapolation may assume that div(E)=0 
	  spatialCoefficientsForTZ(2,0,0,pxc+pc)=1.*.5;      // px=(x^2 + 2xy + y^2)*.5
	  spatialCoefficientsForTZ(1,1,0,pxc+pc)=2.*.5;
	  spatialCoefficientsForTZ(0,2,0,pxc+pc)=1.*.5;

	  spatialCoefficientsForTZ(2,0,0,pyc+pc)= 1.*.5;      // py=(x^2 -2xy - y^2)*.5
	  spatialCoefficientsForTZ(1,1,0,pyc+pc)=-2.*.5;
	  spatialCoefficientsForTZ(0,2,0,pyc+pc)=-1.*.5;
	}
      }
      else if( method==bamx )
      {
	// done below
      }
    }

  }
  else if( degreeSpace==3 )
  {

    spatialCoefficientsForTZ(2,0,0,ex)=1.;      // u=x^2 + 2xy + y^2 + .5*y^3 + .25*x^2*y + .2*x^3  - .3*x*y^2
    spatialCoefficientsForTZ(1,1,0,ex)=2.;
    spatialCoefficientsForTZ(0,2,0,ex)=1.;
    spatialCoefficientsForTZ(0,3,0,ex)=.5;
    spatialCoefficientsForTZ(2,1,0,ex)=.25;
    spatialCoefficientsForTZ(3,0,0,0,ex)=.2;
    spatialCoefficientsForTZ(1,2,0,0,ex)=-.3;

    spatialCoefficientsForTZ(2,0,0,ey)= 1.;      // v=x^2 -2xy - y^2 -.5*x^3 -.25*x*y^2  -.6*x^2*y + .1*y^3
    spatialCoefficientsForTZ(1,1,0,ey)=-2.;
    spatialCoefficientsForTZ(0,2,0,ey)=-1.;
    spatialCoefficientsForTZ(3,0,0,ey)=-.5;
    spatialCoefficientsForTZ(1,2,0,ey)=-.25;
    spatialCoefficientsForTZ(2,1,0,ey)=-.6;
    spatialCoefficientsForTZ(0,3,0,ey)= .1;

    spatialCoefficientsForTZ(2,0,0,hz)= 1.;      // w=x^2 + y^2 -1 +.5 xy + .25*x^3 - .25*y^3
    spatialCoefficientsForTZ(0,2,0,hz)= 1.;
    spatialCoefficientsForTZ(0,0,0,hz)=-1.; 
    spatialCoefficientsForTZ(1,1,0,hz)= .5;
    spatialCoefficientsForTZ(3,0,0,hz)= .25;
    spatialCoefficientsForTZ(0,3,0,hz)=-.25;

    if( solveForAllFields )
    {
      spatialCoefficientsForTZ(2,0,0,hx)  = 1.*.5;   // ** fix me -- change more from Ex,Ey   
      spatialCoefficientsForTZ(1,1,0,hx)  = 2.*.5; 
      spatialCoefficientsForTZ(0,2,0,hx)  = 1.*.5;
      spatialCoefficientsForTZ(0,3,0,hx)  = .5*.5;
      spatialCoefficientsForTZ(2,1,0,hx)  =.25*.5;
      spatialCoefficientsForTZ(3,0,0,0,hx)= .2*.5;
      spatialCoefficientsForTZ(1,2,0,0,hx)=-.3*.5;

      spatialCoefficientsForTZ(2,0,0,hy) =  1.*.5;    
      spatialCoefficientsForTZ(1,1,0,hy) = -2.*.5;
      spatialCoefficientsForTZ(0,2,0,hy) = -1.*.5;
      spatialCoefficientsForTZ(3,0,0,hy) = -.5*.5;
      spatialCoefficientsForTZ(1,2,0,hy) =-.25*.5;
      spatialCoefficientsForTZ(2,1,0,hy) = -.6*.5;
      spatialCoefficientsForTZ(0,3,0,hy) =  .1*.5;

      spatialCoefficientsForTZ(2,0,0,ez)= .8;    
      spatialCoefficientsForTZ(0,2,0,ez)= .4;
      spatialCoefficientsForTZ(0,0,0,ez)=-.1; 
      spatialCoefficientsForTZ(1,1,0,ez)= .1;
      spatialCoefficientsForTZ(3,0,0,ez)= .15;
      spatialCoefficientsForTZ(0,3,0,ez)=-.15;

    }
    

    // -- dispersion components: 
    // ** FINISH ME **
    if( dispersionModel != noDispersion )
    {
      printF("\n >>>> set TZ for P maxNumberOfPolarizationVectors=%i ** FINISH ME ** <<<<\n\n",
         maxNumberOfPolarizationVectors);

      if( method==nfdtd )
      {
	for( int iv=0; iv<maxNumberOfPolarizationVectors; iv++ )
	{
	  const int pc= iv*numberOfDimensions;
	  // Corner extrapolation may assume that div(E)=0 
	  spatialCoefficientsForTZ(2,0,0,pxc+pc)=1.*.5;      // px=(x^2 + 2xy + y^2)*.5
	  spatialCoefficientsForTZ(1,1,0,pxc+pc)=2.*.5;
	  spatialCoefficientsForTZ(0,2,0,pxc+pc)=1.*.5;

	  spatialCoefficientsForTZ(2,0,0,pyc+pc)= 1.*.5;      // py=(x^2 -2xy - y^2)*.5
	  spatialCoefficientsForTZ(1,1,0,pyc+pc)=-2.*.5;
	  spatialCoefficientsForTZ(0,2,0,pyc+pc)=-1.*.5;
	}
      }
      else if( method==bamx )
      {
	// done below
      }

    }


  }
  else if( degreeSpace==4 || degreeSpace==5 )
  {
    if( degreeSpace!=4 ) printF(" ****WARNING***** using a TZ function with degree=4 in space *****\n");
	  
    spatialCoefficientsForTZ(2,0,0,hz)= 1.;      // p=x^2 + y^2 -1 +.5 xy + x^4 + y^4 
    spatialCoefficientsForTZ(0,2,0,hz)= 1.;
    spatialCoefficientsForTZ(0,0,0,hz)=-1.; 
    spatialCoefficientsForTZ(1,1,0,hz)= .5;

    spatialCoefficientsForTZ(4,0,0,hz)= 1.;     
    spatialCoefficientsForTZ(0,4,0,hz)= 1.;     
    spatialCoefficientsForTZ(2,2,0,hz)= -.3;


    spatialCoefficientsForTZ(2,0,0,ex)=1.;      // u=x^2 + 2xy + y^2 + .2*x^4 + .5*y^4 + xy^3
    spatialCoefficientsForTZ(1,1,0,ex)=2.;
    spatialCoefficientsForTZ(0,2,0,ex)=1.;

    spatialCoefficientsForTZ(4,0,0,ex)=.2;   
    spatialCoefficientsForTZ(0,4,0,ex)=.5;   
    spatialCoefficientsForTZ(1,3,0,ex)=1.;   


    spatialCoefficientsForTZ(2,0,0,ey)= 1.;      // v=x^2 -2xy - y^2 +.125*x^4 -.25*y^4 -.8*x^3 y
    spatialCoefficientsForTZ(1,1,0,ey)=-2.;
    spatialCoefficientsForTZ(0,2,0,ey)=-1.;

    spatialCoefficientsForTZ(4,0,0,ey)=.125;
    spatialCoefficientsForTZ(0,4,0,ey)=-.25;
    spatialCoefficientsForTZ(3,1,0,ey)=-.8;

    if( solveForAllFields )
    {

      spatialCoefficientsForTZ(2,0,0,hx)=1.*.3;      // u=x^2 + 2xy + y^2 + .2*x^4 + .5*y^4 + xy^3
      spatialCoefficientsForTZ(1,1,0,hx)=2.*.3;
      spatialCoefficientsForTZ(0,2,0,hx)=1.*.3;
      spatialCoefficientsForTZ(4,0,0,hx)=.2*.3;   
      spatialCoefficientsForTZ(0,4,0,hx)=.5*.3;   
      spatialCoefficientsForTZ(1,3,0,hx)=1.*.3;   

      spatialCoefficientsForTZ(2,0,0,hy)= 1.*.3;      // v=x^2 -2xy - y^2 +.125*x^4 -.25*y^4 -.8*x^3 y
      spatialCoefficientsForTZ(1,1,0,hy)=-2.*.3;
      spatialCoefficientsForTZ(0,2,0,hy)=-1.*.3;
      spatialCoefficientsForTZ(4,0,0,hy)=.125*.3;
      spatialCoefficientsForTZ(0,4,0,hy)=-.25*.3;
      spatialCoefficientsForTZ(3,1,0,hy)=-.8*.3;

      spatialCoefficientsForTZ(2,0,0,ez)= .6;      // p=x^2 + y^2 -1 +.5 xy + x^4 + y^4 
      spatialCoefficientsForTZ(0,2,0,ez)= .4;
      spatialCoefficientsForTZ(0,0,0,ez)=-.3; 
      spatialCoefficientsForTZ(1,1,0,ez)= .2;
      spatialCoefficientsForTZ(4,0,0,ez)= .3;     
      spatialCoefficientsForTZ(0,4,0,ez)=.25;     
      spatialCoefficientsForTZ(2,2,0,ez)= -.2;



    }
    

    // -- dispersion components: 
    // ** FINISH ME **
    if( dispersionModel != noDispersion )
    {
      printF("\n >>>> set TZ for P maxNumberOfPolarizationVectors=%i ** FINISH ME ** <<<<\n\n",
	     maxNumberOfPolarizationVectors);

      if( method==nfdtd )
      {
	for( int iv=0; iv<maxNumberOfPolarizationVectors; iv++ )
	{
	  const int pc= iv*numberOfDimensions;
	  // Corner extrapolation may assume that div(E)=0 
	  spatialCoefficientsForTZ(2,0,0,pxc+pc)=1.*.5;      // px=(x^2 + 2xy + y^2)*.5
	  spatialCoefficientsForTZ(1,1,0,pxc+pc)=2.*.5;
	  spatialCoefficientsForTZ(0,2,0,pxc+pc)=1.*.5;

	  spatialCoefficientsForTZ(2,0,0,pyc+pc)= 1.*.5;      // py=(x^2 -2xy - y^2)*.5
	  spatialCoefficientsForTZ(1,1,0,pyc+pc)=-2.*.5;
	  spatialCoefficientsForTZ(0,2,0,pyc+pc)=-1.*.5;
	}
      }
      else if( method==bamx )
      {
	// done below
      }

    }

  }
  else if( degreeSpace>=6 )
  {
    assert( !solveForAllFields ); // finish me 

    if( degreeSpace!=6 ) printF(" ****WARNING***** using a TZ function with degree=4 in space *****\n");
	  
    spatialCoefficientsForTZ(1,0,0,hz)= 1.;
    spatialCoefficientsForTZ(0,0,0,hz)= 1.;

    spatialCoefficientsForTZ(2,0,0,hz)= 1.;      // p=x^2 + y^2 -1 +.5 xy + x^4 + y^4 
    spatialCoefficientsForTZ(0,2,0,hz)= 1.;
    spatialCoefficientsForTZ(0,0,0,hz)=-1.; 
    spatialCoefficientsForTZ(1,1,0,hz)= .5;

    spatialCoefficientsForTZ(4,0,0,hz)= .2;     
    spatialCoefficientsForTZ(0,4,0,hz)= .4;     
    spatialCoefficientsForTZ(2,2,0,hz)= -.3;

    spatialCoefficientsForTZ(3,2,0,hz)= .4;  
    spatialCoefficientsForTZ(2,3,0,hz)= .8;  
    spatialCoefficientsForTZ(3,3,0,hz)= .7;  

    spatialCoefficientsForTZ(5,1,0,hz)= .25;  
    spatialCoefficientsForTZ(1,5,0,hz)=-.25;  

    spatialCoefficientsForTZ(6,0,0,hz)= .2;     
    spatialCoefficientsForTZ(0,6,0,hz)=-.2;     

    //    spatialCoefficientsForTZ=0.; // ************************************************
	  
    //spatialCoefficientsForTZ(2,4,0,hz)= 1.;
    //spatialCoefficientsForTZ(4,2,0,hz)= 1.;

    spatialCoefficientsForTZ(2,0,0,ex)=1.;      // u=x^2 + 2xy + y^2 + .2*x^4 + .5*y^4 + xy^3
    spatialCoefficientsForTZ(1,1,0,ex)=2.;
    spatialCoefficientsForTZ(0,2,0,ex)=1.;

    spatialCoefficientsForTZ(4,0,0,ex)=.2;   
    spatialCoefficientsForTZ(0,4,0,ex)=.5;   
    spatialCoefficientsForTZ(1,3,0,ex)=1.;   

    spatialCoefficientsForTZ(3,2,0,ex)=.1;      // .1*x^3*y^2

    spatialCoefficientsForTZ(4,2,0,ex)=.3;      // .3 x^4 y^2 ** III
    spatialCoefficientsForTZ(3,3,0,ex)=.4;      // .4 x^3 y^3 ** IV 

    spatialCoefficientsForTZ(6,0,0,ex)=.1;      //  + .1*x^6 +.25*y^6 -.6*x*y^5
    spatialCoefficientsForTZ(0,6,0,ex)=.25;
    spatialCoefficientsForTZ(1,5,0,ex)=-.6;


    spatialCoefficientsForTZ(2,0,0,ey)= 1.;      // v=x^2 -2xy - y^2 +.125*x^4 -.25*y^4 -.8*x^3 y
    spatialCoefficientsForTZ(1,1,0,ey)=-2.;
    spatialCoefficientsForTZ(0,2,0,ey)=-1.;

    spatialCoefficientsForTZ(2,3,0,ey)=-.1;      // -.1*x^2*y^3

    spatialCoefficientsForTZ(3,3,0,ey)=-.4;     //-.4 x^3 y^3 ** III 
    spatialCoefficientsForTZ(2,4,0,ey)=-.3;      //-.3 x^2 y^4 ** IV

    spatialCoefficientsForTZ(4,0,0,ey)=.125;
    spatialCoefficientsForTZ(0,4,0,ey)=-.25;
    spatialCoefficientsForTZ(3,1,0,ey)=-.8;

    spatialCoefficientsForTZ(6,0,0,ey)=.3;    //   .3*x^6 +.1*y^6  + .6*x^5*y 
    spatialCoefficientsForTZ(0,6,0,ey)=.1;
    spatialCoefficientsForTZ(5,1,0,ey)=-.6;

    // -- dispersion components: 
    // ** FINISH ME **
    if( dispersionModel != noDispersion )
    {
      if( method==nfdtd )
      {
	printF("\n >>>> set TZ for P maxNumberOfPolarizationVectors=%i ** FINISH ME ** <<<<\n\n",
	       maxNumberOfPolarizationVectors);

	for( int iv=0; iv<maxNumberOfPolarizationVectors; iv++ )
	{
	  const int pc= iv*numberOfDimensions;
	  // Corner extrapolation may assume that div(E)=0 
	  spatialCoefficientsForTZ(2,0,0,pxc+pc)=1.*.5;      // px=(x^2 + 2xy + y^2)*.5
	  spatialCoefficientsForTZ(1,1,0,pxc+pc)=2.*.5;
	  spatialCoefficientsForTZ(0,2,0,pxc+pc)=1.*.5;

	  spatialCoefficientsForTZ(2,0,0,pyc+pc)= 1.*.5;      // py=(x^2 -2xy - y^2)*.5
	  spatialCoefficientsForTZ(1,1,0,pyc+pc)=-2.*.5;
	  spatialCoefficientsForTZ(0,2,0,pyc+pc)=-1.*.5;
	}
      }
      else if( method==bamx )
      {
	// done below
      }
    }
    
  }
  else
  {
    printF("Maxwell:: not implemented for degree in space =%i \n",degreeSpace);
    Overture::abort("error");
  }
}
// *****************************************************************
// ******************* Three Dimensions ****************************
// *****************************************************************
else if( numberOfDimensions==3 )
{
  // ** finish me -- make the E and H poly's be different
  printF("*** initTZ functions: solveForElectricField=%i solveForMagneticField=%i\n",
         solveForElectricField,solveForMagneticField);

  if ( solveForElectricField )
  {
    initialize3DPolyTW(ex,ey,ez);
  }
  
  if ( solveForMagneticField )
  {
    initialize3DPolyTW(hx,hy,hz);
  }

  if( pxc>=0 ) 
  {
    initialize3DPolyTW(pxc,pyc,pzc);
  }
  
  
}
else
{
  OV_ABORT("ERROR:unimplemented number of dimensions");
}


if( dispersionModel != noDispersion &&  method==bamx )
{
  // BA Maxwell -- assign polarization components
  const int degreeSpaceZ = numberOfDimensions==2 ? degreeSpace : 0;
  const int numPolarizationTerms = 2*maxNumberOfPolarizationComponents;  // note "2*" we save p and p.t 
  for( int m=0; m<numPolarizationTerms; m++ )
  {
    int pc = hz+m+1; // polarization component index in TZ functions
    for( int iz=0; iz<=degreeSpaceZ; iz++ )
    {
      for( int iy=0; iy<=degreeSpace; iy++ )
      {
	for( int ix=0; ix<=degreeSpace; ix++ )
	{
	  // printF("*** initTZ functions: in P pc=%d, m=%d, ix=%d, iy=%d\n",pc,m,ix,iy);
	  spatialCoefficientsForTZ(ix,iy,iz,pc)=(ix+.5*iy+.3*iz+ (2.*m)/maxNumberOfPolarizationComponents)/(degreeSpace*5. + 1.);
	}
      }
      
    }
  }
}

int lastTZComponent = hz+1;  // keep track of the last TZ component assigned above so we can fill in MLA components
if( dispersionModel != noDispersion )
{
   
  int numberOfPolarizationComponents=maxNumberOfPolarizationVectors;  // ** CHECK ME ***
  
  lastTZComponent += numberOfPolarizationComponents*numberOfDimensions;
  if( method==bamx )
    lastTZComponent += numberOfPolarizationComponents*numberOfDimensions;
}

if( dispersionModel != noDispersion &&  nonlinearModel==multilevelAtomic )
{
  // Nonlinear model : multilevelAtomic 
  // printF("\n >>> INIT TZ: Fill in nonlinear TZ variables starting at lastTZComponent=%d <<<\n\n",lastTZComponent);


  const int degreeSpaceZ = numberOfDimensions==2 ? degreeSpace : 0;

  const int & maxNumberOfNonlinearVectors = parameters.dbase.get<int>("maxNumberOfNonlinearVectors");
  for( int m=0; m<maxNumberOfNonlinearVectors; m++ )
  {
    int na = lastTZComponent+m; 

    // printF("POLY-TZ: multilevelAtomic: set TZ for N_%d (TZ component=%d)\n",m,na);
    
    for( int iz=0; iz<=degreeSpaceZ; iz++ )
    {
      for( int iy=0; iy<=degreeSpace; iy++ )
      {
	for( int ix=0; ix<=degreeSpace; ix++ )
	{
	  // printF("*** initTZ functions: in P pc=%d, m=%d, ix=%d, iy=%d\n",pc,m,ix,iy);
	  spatialCoefficientsForTZ(ix,iy,iz,na)=(ix+.25*iy+.35*iz+ (1.5*m)/maxNumberOfNonlinearVectors)/(degreeSpace*5. + 1.);
	}
      }
      
    }
  }
}


for( int n=0; n<numberOfComponentsForTZ; n++ )
{
  for( int i=0; i<ndp; i++ )
    timeCoefficientsForTZ(i,n)= i<=degreeTime ? 1./(i+1) : 0. ;

}
  
if( method==sosup )
{
  // Set the TZ function for (ext,eyt,...) equal to the time derivative of (ex,ey,...)
  const int numberOfFieldComponents=3;  // 2D: (ex,ey,hz),  3D: (ex,ey,ez)
  for( int n=ex, nt=ext; n<ex+numberOfFieldComponents; n++, nt++ )
  {
    for( int i1=0; i1<ndp; i1++ )for( int i2=0; i2<ndp; i2++ )for( int i3=0; i3<ndp; i3++ )
    {
      spatialCoefficientsForTZ(i1,i2,i3,nt)=spatialCoefficientsForTZ(i1,i2,i3,n);
    }
    // E =   a0 + a1*t + a2*t^2 + ...  = [a0,a1,a2,...
    // E_t =      a1   +2*a2*t + 3*a3*t^2  = [a1,2*a2,3*a3
    for( int i=0; i<ndp; i++ )
      timeCoefficientsForTZ(i,nt)= i<degreeTime ? real(i+1.)/(i+2.) : 0. ;

  }
}

    

// Make eps, mu, .. constant in time : 
Range all;
if( rc>=0 )
{
  timeCoefficientsForTZ(all,rc)=0.;
  timeCoefficientsForTZ(0,rc)=1.;
}

if( epsc>=0 )
{
  timeCoefficientsForTZ(all,epsc   )=0.; timeCoefficientsForTZ(0,epsc)   =1.;
  timeCoefficientsForTZ(all,muc    )=0.; timeCoefficientsForTZ(0,muc)    =1.;
  timeCoefficientsForTZ(all,sigmaEc)=0.; timeCoefficientsForTZ(0,sigmaEc)=1.;
  timeCoefficientsForTZ(all,sigmaHc)=0.; timeCoefficientsForTZ(0,sigmaHc)=1.;
}

// printF("TZ: rc=%i epsc=%i\n",rc,epsc);

// ::display(timeCoefficientsForTZ,"timeCoefficientsForTZ","%6.2f ");
// ::display(spatialCoefficientsForTZ,"spatialCoefficientsForTZ","%6.2f ");
    
((OGPolyFunction*)tz)->setCoefficients( spatialCoefficientsForTZ,timeCoefficientsForTZ ); 

// real epsEx = ((OGPolyFunction*)tz)->gd(0,0,0,0,.0,0.,0.,epsc,0.);
// printF(" ********** epsEx = %e *********\n",epsEx);

#endMacro // polynomial TZ macro 

