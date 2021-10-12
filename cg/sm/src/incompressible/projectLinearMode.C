// This file automatically generated from projectLinearMode.bC with bpp.
#include "Cgsm.h"
#include "SmParameters.h"
#include "display.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)  

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)  

// =============================================================================
/// \brief Project out the linear in time and space mode for the second-order system
///       Incompressible Linear Elastic Equations
///
///  The SOS has a linear in time mode:
///       u1 =   C*y*t 
///       u2 = - C*x*t 
// =============================================================================
int Cgsm::
projectLinearMode( int current, real t, real dt )
{

    int & globalStepNumber = parameters.dbase.get<int >("globalStepNumber");

    FILE *& debugFile  =parameters.dbase.get<FILE* >("debugFile");
    FILE *& logFile    =parameters.dbase.get<FILE* >("logFile");
    FILE *& pDebugFile =parameters.dbase.get<FILE* >("pDebugFile");
    
    const int numberOfDimensions           = cg.numberOfDimensions();
    const int numberOfComponentGrids       = cg.numberOfComponentGrids();
    const int & numberOfComponents         = parameters.dbase.get<int >("numberOfComponents");
    const int & projectLinearModeFrequency = parameters.dbase.get<int>("projectLinearModeFrequency");  
      
  // const int & uc                     = parameters.dbase.get<int >("uc");
  // const int & vc                     = parameters.dbase.get<int >("vc");
  // const int & wc                     = parameters.dbase.get<int >("wc");
  // const int & rc                     = parameters.dbase.get<int >("rc");
  // const int & tc                     = parameters.dbase.get<int >("tc");
  // const int & pc                     = parameters.dbase.get<int >("pc");
  // const int & orderOfAccuracyInSpace = parameters.dbase.get<int>("orderOfAccuracy");
  // const int & orderOfTimeAccuracy    = parameters.dbase.get<int>("orderOfTimeAccuracy");
  // const int & numberOfCorrections    = parameters.dbase.get<int>("numberOfCorrections"); 
  // const int & skipLastPressureSolve  = parameters.dbase.get<int>("skipLastPressureSolve");  // For ILE predictor-corrector schemes  

    const int & u1c                   = parameters.dbase.get<int >("u1c");
    const int & u2c                   = parameters.dbase.get<int >("u2c");
    const int & u3c                   = parameters.dbase.get<int >("u3c");

    if( globalStepNumber<=projectLinearModeFrequency ||  debug() & 4 )
        printF("\n>>>>>> projectLinearMode called at globalStepNumber=%d, t=%9.3e\n",globalStepNumber,t);

    SmParameters::TimeSteppingMethodSm & timeSteppingMethodSm = 
                                                                      parameters.dbase.get<SmParameters::TimeSteppingMethodSm>("timeSteppingMethodSm");
    RealArray & timing = parameters.dbase.get<RealArray >("timing");

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    Range C=numberOfComponents;
    const int prev = (current-1+numberOfTimeLevels) % numberOfTimeLevels;
    const int next = (current+1)                    % numberOfTimeLevels;


  // --------------- START LOOP OVER GRIDS -------------------
    Real xDotX=0., yDotY=0., xDotU2=0., yDotU1=0.;
    Real oneDotU1=0., oneDotU2=0., oneDotOne=0., oneDotX=0., oneDotY=0., u1DotU1=0., u2DotU2=0.;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        real time0=getCPU();

        MappedGrid & mg = cg[grid];
        getIndex(mg.gridIndexRange(),I1,I2,I3);

        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
        OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);

        realArray & u = gf[current].u[grid];
        OV_GET_SERIAL_ARRAY(real,u ,uLocal);

    // Find the least squares solution to
    //     || u1 - (a1*t + b1*y*t) ||
    //     || u2 - (a2*t + b2*x*t) ||

        FOR_3D(i1,i2,i3,I1,I2,I3) 
        {
            if( maskLocal(i1,i2,i3)>0 )
            {
                  oneDotOne += 1.;

                  oneDotX   += xLocal(i1,i2,i3,0);         
                  oneDotY   += xLocal(i1,i2,i3,1);

                  oneDotU1  += uLocal(i1,i2,i3,u1c); 
                  oneDotU2  += uLocal(i1,i2,i3,u2c);          

                  xDotX     += SQR(xLocal(i1,i2,i3,0));
                  yDotY     += SQR(xLocal(i1,i2,i3,1));
      
                  yDotU1    += xLocal(i1,i2,i3,1)*uLocal(i1,i2,i3,u1c); 
                  xDotU2    += xLocal(i1,i2,i3,0)*uLocal(i1,i2,i3,u2c); 

                  u1DotU1   += SQR(uLocal(i1,i2,i3,u1c)); 
                  u2DotU2   += SQR(uLocal(i1,i2,i3,u2c));          
            }   
        }


    } // end for grid 

  // For a1,b1 : solve the least squares problem 
  // Solve [ 1 y_1] [a1] = [u1_1]
  //       [ 1 y_2] [b1]   [u1_2]
  //       [ 1 y_3]        [u1_3]
  //       [ ...  ]        [...]
  // 
  // A^T A = [ sum 1    sum y_i ]
  //         [ sum y_i  sum y_i^2]
  // 
    Real det, a11,a12, a21, a22, r1,r2, a1,b1, a2,b2;
    a11 = oneDotOne; a12 = oneDotY;
    a21 = oneDotY;   a22 = yDotY; 
    r1  = oneDotU1;  r2  = yDotU1;
    det = a11*a22 - a12*a21;
    a1 = ( a22*r1 - a12*r2)/det;
    b1 = (-a21*r1 + a11*r2)/det;

    a11 = oneDotOne; a12 = oneDotX;
    a21 = oneDotX;   a22 = xDotX;
    r1  = oneDotU2;  r2  = xDotU2;
    det = a11*a22 - a12*a21;
    a2 = ( a22*r1 - a12*r2)/det;
    b2 = (-a21*r1 + a11*r2)/det;  


  // Real a1 = oneDotU1/oneDotOne;
  // Real a2 = oneDotU1/oneDotOne;

  // Real b1 = + yDotU1/yDoty; 
  // Real b2 = + xDotU2/xDotx; 

    a1 /= t; b1 /= t;
    a2 /= t; b2 /= t;
    if( globalStepNumber<= 3*projectLinearModeFrequency )
    {
        printF(">>> projectLinearMode: u1 ~= (a1 + b1*y)*t,  u2 ~= (a2 + b2*x)*t \n"
                      "       a1=%9.3e, a2=%9.3e, b1=%9.3e, b2=%9.3e\n",a1,a2,b1,b2);
    }

  // ------ NOW ACTUALLY PROJECT --------
    const Real tm=t-dt, tn=t+dt; 
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        real time0=getCPU();

        MappedGrid & mg = cg[grid];
        getIndex(mg.dimension(),I1,I2,I3);

        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
        OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);

        realArray & um = gf[prev].u[grid];
        realArray & u  = gf[current].u[grid];
        realArray & un = gf[next].u[grid];

        OV_GET_SERIAL_ARRAY(real,um,umLocal);
        OV_GET_SERIAL_ARRAY(real,u ,uLocal);
        OV_GET_SERIAL_ARRAY(real,un,unLocal);

        FOR_3D(i1,i2,i3,I1,I2,I3) 
        {
            const Real x = xLocal(i1,i2,i3,0), y=xLocal(i1,i2,i3,1); 

            umLocal(i1,i2,i3,u1c) -= a1*tm + b1*tm*y;
            umLocal(i1,i2,i3,u2c) -= a2*tm + b2*tm*x;

            uLocal(i1,i2,i3,u1c)  -= a1*t + b1*t*y;
            uLocal(i1,i2,i3,u2c)  -= a2*t + b2*t*x;

            unLocal(i1,i2,i3,u1c) -= a1*tn + b1*tn*y;
            unLocal(i1,i2,i3,u2c) -= a2*tn + b2*tn*x;

        }


    } // end for grid 
  


    return 0;
}
