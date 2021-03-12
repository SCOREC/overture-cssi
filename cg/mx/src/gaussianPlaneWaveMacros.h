#beginMacro gpwCommentMacro()
  // This file contains macros that define the Gaussian Plane wave
  //
  // NOTE: This file is included in C++ AND FORTRAN so DO NOT ADD COMMENTS OUTSIDE OF MACRO DEFINITIONS
  //       since one or the other compiler will not recognize them
  //        
#endMacro

#beginMacro setGaussianPlaneWaveInitialConditions()
{
  // =============== Macro: GAUSSIAN PLANE WAVE INITIAL CONDITIONS ==================
  real k0GaussianPlaneWave = parameters.dbase.get<real>("k0GaussianPlaneWave");
  k0GaussianPlaneWave *= twoPi;  // scale by 2*pi

  // (Hz).t = (1/mu) *[ (Ex).y - (Ey).x ]
  if( debug & 2 )
    printF("Setting initial condition to be a Gaussian plane wave, kx,ky,kz=%e %e %e, (x0,y0,z0)=(%e,%e,%e) k0=%g\n",
                kx,ky,kz,x0GaussianPlaneWave,y0GaussianPlaneWave,z0GaussianPlaneWave,k0GaussianPlaneWave);
    
  realSerialArray xei,xhi;

  if( numberOfDimensions==2 && !solveForAllFields )
  {
      
   
    // *new* way June 16 -- use pwc[]
    xei=kx*(xe-x0GaussianPlaneWave) + ky*(ye-y0GaussianPlaneWave) -cc*tE;
    RealArray gpw(I1,I2,I3);
    if( k0GaussianPlaneWave != 0. )
    {
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei))) * cos( k0GaussianPlaneWave*xei );
    }
    else
    {
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei)));
    }

    uLocal(Ie1,Ie2,Ie3,ex)=gpw(Ie1,Ie2,Ie3) * pwc[0];
    uLocal(Ie1,Ie2,Ie3,ey)=gpw(Ie1,Ie2,Ie3) * pwc[1];
    uLocal(Ie1,Ie2,Ie3,hz)=gpw(Ie1,Ie2,Ie3) * pwc[5];

    if( method==nfdtd  || method==bamx )
    {
      // Set values at previous time 
      xei+=cc*dt;
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei)));
      umLocal(Ie1,Ie2,Ie3,ex)=gpw(Ie1,Ie2,Ie3) * pwc[0];
      umLocal(Ie1,Ie2,Ie3,ey)=gpw(Ie1,Ie2,Ie3) * pwc[1];
      umLocal(Ie1,Ie2,Ie3,hz)=gpw(Ie1,Ie2,Ie3) * pwc[5];
    }

    if( saveExtraForcingLevels )
    {
      OV_ABORT("finish me: Gaussian Plane Wave: saveExtraForcingLevels");
    }
   
    
    
  }
  else 
  {
    // Use pwc[] coefficients so we can assign TEz and TMz modes 
    // *new* May 16, 2020
    if( numberOfDimensions==2 )
      xei=kx*(xe-x0GaussianPlaneWave) + ky*(ye-y0GaussianPlaneWave) -cc*tE;
    else 
      xei=kx*(xe-x0GaussianPlaneWave) + ky*(ye-y0GaussianPlaneWave) + kz*(ze-z0GaussianPlaneWave)  -cc*tE;

    RealArray gpw(I1,I2,I3);
    if( k0GaussianPlaneWave != 0. )
    {
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei))) * cos( k0GaussianPlaneWave*xei );
    }
    else
    {
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei)));
    }
     

    uLocal(Ie1,Ie2,Ie3,ex)=gpw(Ie1,Ie2,Ie3) * pwc[0];
    uLocal(Ie1,Ie2,Ie3,ey)=gpw(Ie1,Ie2,Ie3) * pwc[1];
    uLocal(Ie1,Ie2,Ie3,ez)=gpw(Ie1,Ie2,Ie3) * pwc[2];

    if( solveForAllFields )
    {
      uLocal(Ie1,Ie2,Ie3,hx)=gpw(Ie1,Ie2,Ie3) * pwc[3];
      uLocal(Ie1,Ie2,Ie3,hy)=gpw(Ie1,Ie2,Ie3) * pwc[4];
      uLocal(Ie1,Ie2,Ie3,hz)=gpw(Ie1,Ie2,Ie3) * pwc[5];
    }
    
    if( method==nfdtd  || method==bamx )
    {
      // Set values at previous time 
      xei+=cc*dt;
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei)));
      umLocal(Ie1,Ie2,Ie3,ex)=gpw(Ie1,Ie2,Ie3) * pwc[0];
      umLocal(Ie1,Ie2,Ie3,ey)=gpw(Ie1,Ie2,Ie3) * pwc[1];
      umLocal(Ie1,Ie2,Ie3,ez)=gpw(Ie1,Ie2,Ie3) * pwc[2];

      if( solveForAllFields )
      {
        umLocal(Ie1,Ie2,Ie3,hx)=gpw(Ie1,Ie2,Ie3) * pwc[3];
        umLocal(Ie1,Ie2,Ie3,hy)=gpw(Ie1,Ie2,Ie3) * pwc[4];
        umLocal(Ie1,Ie2,Ie3,hz)=gpw(Ie1,Ie2,Ie3) * pwc[5];
      }
      
    }
    
  }
}

#endMacro


#beginMacro getErrorsGaussianPlaneWave()
{
  // =============================================================================================
  // MACRO: Get errors for a Gaussian plane wave
  // =============================================================================================
      
  printF("GET ERRORS for GaussianPlaneWave solveForAllFields=%d, t=%9.3e\n",solveForAllFields,tE);

  real k0GaussianPlaneWave = parameters.dbase.get<real>("k0GaussianPlaneWave");
  k0GaussianPlaneWave *= twoPi;  // scale by 2*pi
  if( numberOfDimensions==2 && !solveForAllFields )
  {
 
    // ** new way June 16, 2020 -- use plane-wave-coefficients to scale
    OV_GET_SERIAL_ARRAY(real,(*cgerrp)[grid],errLocal);
   
    realSerialArray xei(Ie1,Ie2,Ie3);
    xei= ( kx*(xe(Ie1,Ie2,Ie3)-x0GaussianPlaneWave) +
           ky*(ye(Ie1,Ie2,Ie3)-y0GaussianPlaneWave) 
           -cc*tE );
    RealArray gpw(Ie1,Ie2,Ie3);
    if( k0GaussianPlaneWave!=0. )
    {
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei))) * cos( k0GaussianPlaneWave*xei );
    }
    else
    {
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei)));
    }

    errLocal(Ie1,Ie2,Ie3,ex) = uLocal(Ie1,Ie2,Ie3,ex) - gpw(Ie1,Ie2,Ie3) * pwc[0];
    errLocal(Ie1,Ie2,Ie3,ey) = uLocal(Ie1,Ie2,Ie3,ey) - gpw(Ie1,Ie2,Ie3) * pwc[1];
    errLocal(Ie1,Ie2,Ie3,hz) = uLocal(Ie1,Ie2,Ie3,hz) - gpw(Ie1,Ie2,Ie3) * pwc[5];

  
  }
  else
  {
    // *new* May 16, 2020
    OV_GET_SERIAL_ARRAY(real,(*cgerrp)[grid],errLocal);
     
    realSerialArray xei(Ie1,Ie2,Ie3);
    if( numberOfDimensions==2 )
      xei=  kx*(xe(Ie1,Ie2,Ie3)-x0GaussianPlaneWave) +
            ky*(ye(Ie1,Ie2,Ie3)-y0GaussianPlaneWave) 
            -cc*tE;
    else
      xei=  kx*(xe(Ie1,Ie2,Ie3)-x0GaussianPlaneWave) +
            ky*(ye(Ie1,Ie2,Ie3)-y0GaussianPlaneWave) +
            kz*(ze(Ie1,Ie2,Ie3)-z0GaussianPlaneWave)
            -cc*tE;

    RealArray gpw(Ie1,Ie2,Ie3);
    if( k0GaussianPlaneWave!=0. )
    {
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei))) * cos( k0GaussianPlaneWave*xei );
    }
    else
    {
      gpw(Ie1,Ie2,Ie3) = exp(-betaGaussianPlaneWave*((xei)*(xei)));
    }
    

    errLocal(Ie1,Ie2,Ie3,ex) = uLocal(Ie1,Ie2,Ie3,ex) - gpw(Ie1,Ie2,Ie3) * pwc[0];
    errLocal(Ie1,Ie2,Ie3,ey) = uLocal(Ie1,Ie2,Ie3,ey) - gpw(Ie1,Ie2,Ie3) * pwc[1];
    errLocal(Ie1,Ie2,Ie3,ez) = uLocal(Ie1,Ie2,Ie3,ez) - gpw(Ie1,Ie2,Ie3) * pwc[2];

    if( solveForAllFields )
    {
      errLocal(Ie1,Ie2,Ie3,hx) = uLocal(Ie1,Ie2,Ie3,hx) - gpw(Ie1,Ie2,Ie3) * pwc[3];
      errLocal(Ie1,Ie2,Ie3,hy) = uLocal(Ie1,Ie2,Ie3,hy) - gpw(Ie1,Ie2,Ie3) * pwc[4];
      errLocal(Ie1,Ie2,Ie3,hz) = uLocal(Ie1,Ie2,Ie3,hz) - gpw(Ie1,Ie2,Ie3) * pwc[5];

    }
    

  }
  
}
#endMacro

#beginMacro getGaussianPlaneWave(OP, u,t,x,y)
  !===============================================================================================
  ! Macro:
  !   Define a Gaussian Plane Wave incident field in 2D   
  !===============================================================================================

  xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t)
  expxi = amp*exp(-betaGP*xi*xi ) * cos( k0GaussianPlaneWave*xi )
  if( solveForAllFields.eq.0 )then
    ! *new* way June 16, 2020 *wdh*
    u(i1,i2,i3,ex) = u(i1,i2,i3,ex) OP expxi*pwc(0)
    u(i1,i2,i3,ey) = u(i1,i2,i3,ey) OP expxi*pwc(1)
    u(i1,i2,i3,hz) = u(i1,i2,i3,hz) OP expxi*pwc(5)

    ! *old way*
    ! uhz = expxi
    ! uex = uhz*(-ky/(eps*cc))
    ! uey = uhz*( kx/(eps*cc))
    ! u(i1,i2,i3,ex) = u(i1,i2,i3,ex) OP uex
    ! u(i1,i2,i3,ey) = u(i1,i2,i3,ey) OP uey
    ! u(i1,i2,i3,hz) = u(i1,i2,i3,hz) OP uhz
  else

    u(i1,i2,i3,ex) = u(i1,i2,i3,ex) OP expxi*pwc(0)
    u(i1,i2,i3,ey) = u(i1,i2,i3,ey) OP expxi*pwc(1)
    u(i1,i2,i3,ez) = u(i1,i2,i3,ez) OP expxi*pwc(2)

    u(i1,i2,i3,hx) = u(i1,i2,i3,hx) OP expxi*pwc(3)
    u(i1,i2,i3,hy) = u(i1,i2,i3,hy) OP expxi*pwc(4)
    u(i1,i2,i3,hz) = u(i1,i2,i3,hz) OP expxi*pwc(5)

  end if    

#endMacro 

#beginMacro getGaussianPlaneWave3d(OP, u,t,x,y,z)
  !===============================================================================================
  ! Macro:
  !   Define a Gaussian Plane Wave incident field in 3D   
  !   Changed: May 16, 2020 -- use pwc() coefficients 
  !===============================================================================================

  xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - cc*(t)
  expxi = amp*exp(-betaGP*xi*xi ) * cos( k0GaussianPlaneWave*xi )

  u(i1,i2,i3,ex) = u(i1,i2,i3,ex) OP expxi*pwc(0)
  u(i1,i2,i3,ey) = u(i1,i2,i3,ey) OP expxi*pwc(1)
  u(i1,i2,i3,ez) = u(i1,i2,i3,ez) OP expxi*pwc(2)

  if( solveForAllFields.eq.1 )then
    u(i1,i2,i3,hx) = u(i1,i2,i3,hx) OP expxi*pwc(3)
    u(i1,i2,i3,hy) = u(i1,i2,i3,hy) OP expxi*pwc(4)
    u(i1,i2,i3,hz) = u(i1,i2,i3,hz) OP expxi*pwc(5)
  end if

#endMacro 

