// This file automatically generated from getChampParameters.bC with bpp.
// =========================================================================================
// Evaluate the parameters for the CHAMP iteration for partitioned schemes.
// =========================================================================================
#include "Cgad.h"
#include "Interface.h"  

// =========================================================================================
//   Here is the start at a routine that estimates the Champ coupling parameters pL and pR
// =========================================================================================


// =========================================================================================
/// \brief Evaluate the parameters for the CHAMP iteration for partitioned schemes.
///
/// Input
///    dx : grid spacing in normal direction for computing lambdaDL = DL*dt/dx^2 
///  
/// Output:
///      champParameters(0,side,axis,grid)=pL;     // optimized Scwartz Parameter for side 1 : Sl = pL/dx
///      champParameters(1,side,axis,grid)=pR;     // optimized Scwartz Parameter for side 2 : Sr = pR/dx 
///      champParameters(2,side,axis,grid)=KLR;  // K1/K2
///      champParameters(3,side,axis,grid)=DLR;   // D1/D2 
/// 
///  pL    = champParameters(0,side,axis,grid);    // optimized Scwartz Parameter for side 1
///  pR    = champParameters(1,side,axis,grid);    // optimized Scwartz Parameter for side 2
///  KLR = champParameters(2,side,axis,grid);    // K1/K2
///  DLR  = champParameters(3,side,axis,grid);    // D1/D2   
///  Sl    = champParameters(4,side,axis,grid);  
///  dxs   = champParameters(5,side,axis,grid);
///  DL    = champParameters(6,side,axis,grid);
///  KL    = champParameters(7,side,axis,grid);
///  DR    = champParameters(8,side,axis,grid);
///  KR    = champParameters(9,side,axis,grid);      
// =========================================================================================
int getChampParameters( int grid, int side, int axis, int grid2, int side2, int axis2,
                                                real dt, real dx, Parameters & parameters, RealArray & champParameters )
{

    const int orderOfAccuracy = parameters.dbase.get<int>("orderOfAccuracy"); 

  // Retrieve the parameters from the opposite side of the interface:
    Parameters & paramRight = getInterfaceParameters( grid, side, axis, parameters);

  // This only works if the opposite side is a Cgad object *** FIX ME ***
    std::vector<real> & kappaLeft  = parameters.dbase.get<std::vector<real> >("kappa");
    std::vector<real> & kappaRight = paramRight.dbase.get<std::vector<real> >("kappa");
    const Real & kThermalLeft  = parameters.dbase.get<Real>("thermalConductivity");
    const Real & kThermalRight = paramRight.dbase.get<Real>("thermalConductivity");


    const Real DL = kappaLeft[0];
    const Real KL = kThermalLeft;
    const Real DR = kappaRight[0];
    const Real KR = kThermalRight; 

    Real KLR = KL/KR; // kThermalLeft/kThermalRight;
    Real DLR = DL/DR; // kappaLeft[0]/kappaRight[0];

    if( dt==0 || dx==0 || isnan(dx) )
    {
        if( dx==0 )
            printF("getChampPars: ERROR dx=0!\n");
        if( isnan(dx) )
            printF("getChampPars: ERROR dx=nan!\n");
        if( dt==0 )
            printF("getChampPars: ERROR dt=0!\n");
        OV_ABORT("ERROR");
    }

    Real lambdaD = DL*dt/(dx*dx); 
    printF(">>> getChampPars: DL=%g, dt=%g, dx=%g, lambdaD=%g. \n"
                  "    On input: pL = champParameters(0,side,axis,grid)=%g, pR=champParameters(1,side,axis,grid)=%g\n",
              DL,dt,dx,lambdaD,champParameters(0,side,axis,grid),champParameters(1,side,axis,grid));

    Real pL=1., pR=1.;  // these are set below

  // do this for now 
  // Real pL =1.;  // default value 
  // Real pR =1.; 
  // printF("***** fillChampBC: (side,axis,grid)=(%d,%d,%d) : KLR=Kleft/Kright=%g, DLR=Dleft/Dright=%g, Sl=%g, Sr=%g\n ****",side,axis,grid,KLR,DLR,Sl,Sr);

    if( champParameters(0,side,axis,grid)<0  )
    {
    // ---- Guess the Champ parameters pL and pR -----
      // -- Here is a table of known parameters computed elsewhere -- 
            const int lamDt=0, DLRt=1, KLRt=2, pLt=3, pRt=4; // entries in the t=table
      // table(0,i) = lambdaD;   // DL*dt/dx^2 
      // table(1,i) = DLR;       // DL/DR
      // table(2,i) = KLR;       // KL/KR
      // table(3,i) = pL;    
      // table(4,i) = pR;
            const int maxNumberInTable=20, numPerRow=5;
            Real pTable[numPerRow*maxNumberInTable];
            #define table(i,j) pTable[i+5*(j)]
      // Values computed with ~/Dropbox/AMP/champ4/champMatlab/timeStepping/optimalParameters/getTimePvalue
            int j=0; 
            if( orderOfAccuracy==2 )
            {
            }
            else if( orderOfAccuracy==4 )
            {
        // -------> Optimal pL, pR: BDF4: lambdaDL=14.4, DLR=1.125, KLR=2, -pL=6.9355e-01 -pR=1.0180e+00 |A|_TS=7.823e-01.
        // -------> Optimal pL, pR: BDF4: lambdaDL=57.6, DLR=1.125, KLR=2, -pL=1.5700e-01 -pR=1.0383e+00 |A|_TS=7.844e-01.
                table(lamDt,j) = 14.4; table(DLRt,j)=1.125; table(KLRt,j)=2.; table(pLt,j)=6.9355e-01; table(pRt,j)=1.0180e+00;  j++; 
                table(lamDt,j) = 57.6; table(DLRt,j)=1.125; table(KLRt,j)=2.; table(pLt,j)=1.5696e-01; table(pRt,j)=1.0395e+00;  j++; 
        // four disks: dt=.01, dx=1/40
                table(lamDt,j)=  1.6000e+01; table(DLRt,j)=  8.3333e-01; table(KLRt,j)=  1.2500e+00; table(pLt,j)=  5.1933e-01; table(pRt,j)=  8.7768e-01;  j++; // BDF4 |A|=7.92e-01
                table(lamDt,j)=  1.6000e+01; table(DLRt,j)=  7.6923e-01; table(KLRt,j)=  1.6667e+00; table(pLt,j)=  4.4610e-01; table(pRt,j)=  1.0624e+00;  j++; // BDF4 |A|=7.83e-01
                table(lamDt,j)=  1.6000e+01; table(DLRt,j)=  7.1429e-01; table(KLRt,j)=  2.5000e+00; table(pLt,j)=  3.9217e-01; table(pRt,j)=  1.3311e+00;  j++; // BDF4 |A|=7.61e-01
                table(lamDt,j)=  1.6000e+01; table(DLRt,j)=  6.6667e-01; table(KLRt,j)=  5.0000e+00; table(pLt,j)=  3.4456e-01; table(pRt,j)=  1.8885e+00;  j++; // BDF4 |A|=7.13e-01
        // G4 : dt=.01, dx=1/160 = .00625
                table(lamDt,j)=  2.5600e+02; table(DLRt,j)=  8.3333e-01; table(KLRt,j)=  1.2500e+00; table(pLt,j)=  1.3679e-01; table(pRt,j)=  1.7935e-01;  j++; // BDF4 |A|=8.71e-01
                table(lamDt,j)=  2.5600e+02; table(DLRt,j)=  7.6923e-01; table(KLRt,j)=  1.6667e+00; table(pLt,j)=  8.1162e-02; table(pRt,j)=  3.1286e-01;  j++; // BDF4 |A|=8.29e-01
                table(lamDt,j)=  2.5600e+02; table(DLRt,j)=  7.1429e-01; table(KLRt,j)=  2.5000e+00; table(pLt,j)=  7.2620e-02; table(pRt,j)=  3.2703e-01;  j++; // BDF4 |A|=7.91e-01
                table(lamDt,j)=  2.5600e+02; table(DLRt,j)=  6.6667e-01; table(KLRt,j)=  5.0000e+00; table(pLt,j)=  5.3259e-02; table(pRt,j)=  3.7317e-01;  j++; // BDF4 |A|=7.11e-01
            }
            else
            {
                OV_ABORT("ERROR: getChampParameters: unexpected order of accuracy");
            }
            const int numberInTable=j; 
            assert( numberInTable < maxNumberInTable );
      // assert( j==numberInTable );
      // Also look to see if the table holds the oppossite side parameters:
      //   [DL,KL] <-> [DR,KR]
            const real lamDR = lambdaD/DLR; // dt*DR/dx^2
            const real DRL = 1./DLR;
            const real KRL = 1./KLR;
      // --- Look for close values in the table ---
            const Real tol=.1; // Relative error in lambdaD, DLR and KLR 
            bool found=false, oppositeSideFound=false;
            for( int j=0; j<numberInTable; j++ )
            {
              if( fabs(lambdaD-table(lamDt,j))< tol*table(lamDt,j) &&
                      fabs(DLR    -table(DLRt,j)) < tol*table(DLRt,j) &&
                      fabs(KLR    -table(KLRt,j)) < tol*table(KLRt,j) )
              {
                  pL = table(pLt,j); pR = table(pRt,j); found=true;
              }
       // look for [DL,KL] <-> [DR,KR]
              if( fabs(lamDR-table(lamDt,j))< tol*table(lamDt,j) &&
                      fabs(DRL  -table(DLRt,j)) < tol*table(DLRt,j) &&
                      fabs(KRL  -table(KLRt,j)) < tol*table(KLRt,j) )
              {
                  pR = table(pLt,j); pL = table(pRt,j); found=true;  oppositeSideFound=true; // NOTE flip pL and pR
              }
            }
            if( found )
            {
                if( !oppositeSideFound )
                    printF("\n>>getChampCoupling: Values FOUND in the table: lambdaD=%g, DLR=%g, KLR=%g. Using pL=%g, pR=%g\n\n",lambdaD,DLR,KLR,pL,pR);
                else
                    printF("\n>>getChampCoupling: OPPOSITE SIDE values FOUND in the table: lambdaD=%g, DLR=%g, KLR=%g. Using pL=%g, pR=%g\n\n",lambdaD,DLR,KLR,pL,pR);
            }
            else
            {
                pL=1.; pR=1.; // default 
                if( orderOfAccuracy==2 )
                {
                    printF("\n>>getChampCoupling: FINISH ME FOR GUESSING pL and PR for orderOfAccuracy=2\n");
                }
                else if( orderOfAccuracy==4 )
                {
                    if( lambdaD>=100 )
                    {
            // Here is a guess for CHAMP4 : *fix me* for CHAMP2
                        if( DLR>=.5 && DLR <=2. )
                        { // DLR is close to 1 -- estimate values from graph of (pL,pR) versus KLR (DLR=1) in champ4 paper
                            if( KLR >1 )
                            {
                                pL=.05; pR=.2; 
                            }
                            else
                            {
                                pL=.2; pR=.05; 
                            }
                        }
                        if( KLR>=.5 && KLR<=2. )
                        {
              // KLR is close to 1 -- estimate  values from graph of (pL,pR) versus DLR (KLR=1) in champ4 paper
                            if( DLR >1 )
                            {
                                pL=.5; pR=.06; 
                            }
                            else
                            {
                                pL=.04; pR=.15; 
                            }        
                        }
                    }
                }
                else
                {
                    OV_ABORT("ERROR: getChampParameters: unexpected order of accuracy");
                }    
                printF("\ngetChampCoupling: Values NOT FOUND in the table: lambdaD=%g, DLR=%g, KLR=%g. Using: pL=%g, pR=%g\n\n",lambdaD,DLR,KLR,pL,pR);
            }
        champParameters(0,side,axis,grid)=pL;
        champParameters(1,side,axis,grid)=pR;

        RealArray & champParametersRight = paramRight.dbase.get<RealArray>("champParameters");
        if( champParametersRight(0,side,axis,grid)<0 )
        { // Set pL and pR on right side too: 
            champParametersRight(0,side2,axis2,grid2)=pR;  // flip pL <-> pR for opposite side 
            champParametersRight(1,side2,axis2,grid2)=pL;
        }

        // OV_ABORT("STOP HERE FOR NOW");

    // use default.
    // This value may have been set in CgAd: setupPdeParameters
    // champParameters(0,side,axis,grid)=pL;     // optimized Scwartz Parameter for side 1 : Sl = pL/dx
    }
  // champParameters(1,side,axis,grid)=pR;     // optimized Scwartz Parameter for side 2 : Sr = pR/dx   ** FIX ME : get from other side if needed

    champParameters(2,side,axis,grid)=KLR;   // K1/K2
    champParameters(3,side,axis,grid)=DLR;   // D1/D2         

    champParameters(6,side,axis,grid) = DL; // kappaLeft[0];
    champParameters(7,side,axis,grid) = KL; // kThermalLeft;
    champParameters(8,side,axis,grid) = DR; // kappaRight[0];
    champParameters(9,side,axis,grid) = KR; // kThermalRight; 

    return 0;
}