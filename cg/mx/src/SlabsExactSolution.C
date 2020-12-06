// This file automatically generated from SlabsExactSolution.bC with bpp.
#include "SlabsExactSolution.h"

#include "DispersiveMaterialParameters.h"

#include "ParallelUtility.h"

// ===============================================================================
// Class to define an exact solution to Maxwell's equations
// 
//     Scattering from one or more SLABS
//
//     -------------------------------------------------------------------
//     |                 |     |       |      |                          |
//     |         0       |  1  |   2   |   3  |           4              |
//     |                 |     |       |      |                          |
//     -------------------------------------------------------------------
///
// ===============================================================================


    

#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++)                                       for(i2=I2Base; i2<=I2Bound; i2++)                                     for(i1=I1Base; i1<=I1Bound; i1++)


// - 
// - // ===================================================================================
// - /// \brief Evaluate the frequency space solution (complex valued) at a single point.
// - /// \param x[3] (input): point to evaluate the solution
// - /// 
// - /// \param E[6] (output):  real and imaginary parts of E : [Exr,Eyr,Ezr,Exi,Eyi,Ezi]
// - /// \param H[6] (output):  optionally compute real and imaginary parts of H : [Hxr,Hyr,Hzr,Hxi,Hyi,Hzi]
// - ///
// - /// \note: This may not be so efficient.
// - // ===================================================================================
// - int SlabsExactSolution::
// - eval( real xv[3], real *Ev, real *Hv /*= NULL */  )
// - {
// -   OV_ABORT("eval: FINISH ME");
// -     
// -   return 0;
// - }
// - 
// - 

// - 


#include <complex>
typedef ::real LocalReal;
typedef ::real OV_real;

#include "ComplexArray.h"

// ===============================================================================
/// \brief  Constructor for th class that defines exact solutions to Maxwell's equations for a sphere
// ===============================================================================
SlabsExactSolution::
SlabsExactSolution( )
{

    dbase.put<int>("numberOfDimensions")=2;
    dbase.put<int>("numberOfDomains")=3;
    dbase.put<int>("numScatteringCoeff")=-1;
    dbase.put<int>("solveForAllFields")=0;

    dbase.put<int>("initialized")=0;

    dbase.put<int>("scatteringCase")=0;   // There are four scattering cases, 2 polarizations, forward/backward 

    dbase.put<std::complex<LocalReal>* >("Vvc")=NULL;   
    dbase.put<std::complex<LocalReal>* >("kvc")=NULL;   
    dbase.put<std::complex<LocalReal>* >("rtc")=NULL;   
    dbase.put<std::complex<LocalReal> >("s");   


}



// ===============================================================================
/// \brief destructor for the class that defines exact solutions to Maxwell's equations for a sphere
// ===============================================================================
SlabsExactSolution::
~SlabsExactSolution()
{
    std::complex<LocalReal> *& Vvc = dbase.get<std::complex<LocalReal>* >("Vvc"); 
    delete [] Vvc;

    std::complex<LocalReal> *& kvc = dbase.get<std::complex<LocalReal>* >("kvc"); 
    delete [] kvc;

    std::complex<LocalReal> *& rtc = dbase.get<std::complex<LocalReal>* >("rtc"); 
    delete [] rtc;

}

// There are four scattering cases, 2 polarizations, forward/backward 
int SlabsExactSolution::
setScatteringCase( int scatCase )
{
    if( scatCase>=0 && scatCase<4 )
    {
        dbase.get<int>("scatteringCase")=scatCase;
    }
    else
    {
        printF("SlabsExactSolution::setScatteringCase:ERROR: invalid scatCase=%d (must be 0,1,2 or 3)\n",scatCase);
    }
    
}


// - // ===============================================================================
// - /// \brief Initialize the slabs solution.
// - // ===============================================================================
// - int SlabsExactSolution::
// - initializeOld( )
// - {
// -   int & initialized = dbase.get<int>("initialized");
// -   if( initialized )
// -     return 0;
// - 
// -   initialized=1;
// - 
// -   printF("\n ------------------- SlabsExactSolution::initialize ------------------------\n\n");
// - 
// - 
// - 
// - 
// -   const int numDomains=3;
// -   const int forward=0, backward=1;
// - 
// -   LocalReal sr, si;
// -   RealArray kvI(3); // , kv(3,2,2,numDomains), Vv(6,2,2,numDomains);
// - 
// -   std::complex<LocalReal> I(0.0,1.0);
// - 
// - 
// -   LocalReal & a = dbase.get<LocalReal>("a");
// -   LocalReal & b = dbase.get<LocalReal>("b");
// - 
// -   std::complex<LocalReal> *& Vvc = dbase.get<std::complex<LocalReal>* >("Vvc"); 
// -   Vvc = new std::complex<LocalReal> [6*2*2*numDomains];  
// - 
// -   std::complex<LocalReal> *& kvc = dbase.get<std::complex<LocalReal>* >("kvc"); 
// -   kvc = new std::complex<LocalReal> [3*2*2*numDomains];  
// - 
// -   std::complex<LocalReal> *& rtc = dbase.get<std::complex<LocalReal>* >("rtc"); 
// -   rtc = new std::complex<LocalReal> [8*4];  
// - 
// -   std::complex<LocalReal> & s = dbase.get<std::complex<LocalReal> >("s"); 
// - 
// - //  std::complex<LocalReal> kvc[3*2*2*numDomains];
// - //  std::complex<LocalReal> Vvc[6*2*2*numDomains];
// - #define kv(i,j,k,l) kvc[ (i)+3*( (j)+2*( (k) + 2*(l) ) ) ]
// - #define Vv(i,j,k,l) Vvc[ (i)+6*( (j)+2*( (k) + 2*(l) ) ) ]
// -    
// -   // store complex scattering coeff here: 
// - //   std::complex<LocalReal> rtc[8*4]
// -   std::complex<LocalReal> qvc[6];
// - #define rtv(i,j) rtc[(i)+8*(j)] 
// - #define qv(i) qvc[(i)]
// - 
// - 
// -   // File written by bagdm/oblique/baSlab.m 
// - #include "slabSolution.h"
// - 
// -   s = sr + I*si;
// -   
// -   // std::complex<LocalReal> s(sr,si);
// - 
// - 
// -   printF(" kvI=[%g,%g,%g]*(2*pi)\n",kvI(0)/twoPi,kvI(1)/twoPi,kvI(2)/twoPi);
// -   printF("\n ===== SLAB: kx and eigenvectors in the %d domains =====\n",numDomains);
// -   
// -   for( int d=0; d<numDomains; d++)
// -   {
// -     for( int i=0; i<2; i++ )  // polarizations
// -     {
// -       for( int j=0; j<6; j++ )
// -         qv(j) = Vv(j,i,backward,d);
// - 
// -       fprintf(stdout,"MED%d: reflected   %d: kx=%5.2f+%5.2fI Vv=[%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f]+ [%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f] I\n",    
// -               d,i,std::real(kv(0,i,backward,d)),imag(kv(0,i,backward,d)),    
// -               std::real(qv(0)),std::real(qv(1)),std::real(qv(2)),std::real(qv(3)),std::real(qv(4)),std::real(qv(5)),     
// -               imag(qv(0)),imag(qv(1)),imag(qv(2)),imag(qv(3)),imag(qv(4)),imag(qv(5)));
// -     }
// -     for( int i=0; i<2; i++ )  // polarizations
// -     {
// -        for( int j=0; j<6; j++ )
// -         qv(j) = Vv(j,i,forward,d);
// -        
// -        fprintf(stdout, "MED%d: transmitted %d: kx=%5.2f+%5.2fI Vv=[%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f]+ [%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f] I\n",    
// - 	       d,i,std::real(kv(0,i,forward,d)),imag(kv(0,i,forward,d)),     
// -               std::real(qv(0)),std::real(qv(1)),std::real(qv(2)),std::real(qv(3)),std::real(qv(4)),std::real(qv(5)),     
// -               imag(qv(0)),imag(qv(1)),imag(qv(2)),imag(qv(3)),imag(qv(4)),imag(qv(5)));
// -     }
// -   } // end for d (domain)
// -   
// -   for( int icase=0; icase<4; icase++ ) //  4 cases 
// -   {
// -     if( icase<2 )
// -       fprintf(stdout,"\n +++++ reflection and transmission coefficients : forward : polarization %d +++++\n",icase);
// -     else
// -       fprintf(stdout,"\n +++++ reflection and transmission coefficients : backward: polarization %d +++++\n",icase-2);
// - 
// -     std::complex<LocalReal> x0,x1,x2,x3,x4,x5,x6,x7;
// -     x0 = rtv(0,icase);
// -     x1 = rtv(1,icase);
// -     x2 = rtv(2,icase);
// -     x3 = rtv(3,icase);
// -     x4 = rtv(4,icase);
// -     x5 = rtv(5,icase);
// -     x6 = rtv(6,icase);
// -     x7 = rtv(7,icase);
// - 
// -     fprintf(stdout," rL1   =[%6.3f + %6.3f I], rL2   =[%6.3f + %6.3f I]\n",std::real(x0),imag(x0),std::real(x1),imag(x1));
// -     fprintf(stdout," rM1   =[%6.3f + %6.3f I], rM2   =[%6.3f + %6.3f I]\n",std::real(x2),imag(x2),std::real(x3),imag(x3));
// -     fprintf(stdout," tauM1 =[%6.3f + %6.3f I], tauM2 =[%6.3f + %6.3f I]\n",std::real(x4),imag(x4),std::real(x5),imag(x5));
// -     fprintf(stdout," tauR1 =[%6.3f + %6.3f I], tauR2 =[%6.3f + %6.3f I]\n",std::real(x6),imag(x6),std::real(x7),imag(x7));
// -   }
// -   
// - 
// - 	
// - 
// -   // ::display(Vv,"eigenvectors");
// - 
// - 
// -   printF("\n ------------------- SlabsExactSolution:: END initialize ------------------------\n\n");
// - 
// -   return 0;
// - }

// lapack routines
#ifdef OV_USE_DOUBLE
    #define GESV EXTERN_C_NAME(zgesv)
    #define GEEV EXTERN_C_NAME(zgeev)
#else
    #define GESV EXTERN_C_NAME(cgesv)
    #define GEEV EXTERN_C_NAME(cgeev)
#endif

extern "C"
{
    /* Solve  A*X = B (complex) */
    void GESV( int & N, int & NRHS, std::complex<LocalReal>  & a, const int & lda, int & ipvt, std::complex<LocalReal> & b, int & LDB, int & info );

    void GEEV( char *jobvl,
           	     char* jobvr,
           	     int & n,
           	     std::complex<LocalReal> & a,
           	     const int & lda,
           	     std::complex<LocalReal> & w,
           	     std::complex<LocalReal> &vl,
           	     int & ldvl,
           	     std::complex<LocalReal> &vr,
           	     int & ldvr,
           	     std::complex<LocalReal> & work,
           	     int & lwork,
           	     LocalReal & rwork,
           	     int & info );
}


//==============================================================================
//===================  COMPUTE THE SLAB SOLUION ================================
///
/// \brief Determine the exact solution from the scattering of a plane wave from one or more slabs
///
/// \notes See the matlab file getReflectionAndTransmission.m
//==============================================================================
int SlabsExactSolution::
initialize( CompositeGrid & cg, int numberOfDomains,
          	    std::vector<DispersiveMaterialParameters> & dispersiveMaterialParameters,
          	    const LocalReal & omega, const RealArray & kvI, const int solveForAllFields )
{

    int & initialized = dbase.get<int>("initialized");
    if( initialized )
    {
        printF("SlabsExactSolution::initialize - WARNING: already initialized! Nothing to be done\n");
        return 0;
    }
    
    initialized=1;

    dbase.get<int>("numberOfDomains")=numberOfDomains;
    dbase.get<int>("solveForAllFields")=solveForAllFields;

    printf("SlabsExactSolution::initialize *new* numberOfDomains=%d, omega=%g, kvI=[%g,%g,%g]\n",
           	     numberOfDomains,omega,kvI(0),kvI(1),kvI(2)); 

  // Keep a pointer to the vector of DispersiveMaterialParameters
    dbase.put<std::vector<DispersiveMaterialParameters>* >("pDispersiveMaterialParameters") = &dispersiveMaterialParameters;
    

    int idebug=3;
    const int forward=0, backward=1;

    const int ex=0, ey=1, ez=2, hx=3, hy=4, hz=5;


  // --- find the positions of the interfaces -----
    RealArray interfacePosition(numberOfDomains-1);
    assert( numberOfDomains==cg.numberOfComponentGrids() );
    for( int grid=0; grid<cg.numberOfComponentGrids()-1; grid++ )
    {
        MappedGrid & mg = cg[grid];
        const IntegerArray & gid = mg.gridIndexRange();
        
        const bool isRectangular=mg.isRectangular();
        if( !isRectangular )
            mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

        OV_GET_SERIAL_ARRAY_CONDITIONAL(real,mg.center(),xLocal,!isRectangular); // *wdh* added conditional, Nov 21, 2020

        LocalReal dvx[3]={1.,1.,1.}, xab[2][3]={{0.,0.,0.},{0.,0.,0.}};
        if( isRectangular )
        {
            mg.getRectangularGridParameters( dvx, xab );

            interfacePosition(grid)= xab[1][0];  // right-most x-position
            printF("SES: grid=%d: xInterface=%9.3e\n",grid, interfacePosition(grid));
            if( grid>0 && interfacePosition(grid) < interfacePosition(grid-1) )
            {
      	printF("SlabsExactSolution::ERROR: x-coordinates if slab interfaces are not increasing\n");
      	OV_ABORT("SlabsExactSolution::fix me");
            }
            

        }    
        else
        {
            interfacePosition(grid)=xLocal(gid(1,0),0,0,0);  // fix me for parallel
            #ifdef USE_PPP
                OV_ABORT("finish me");
            #endif 
        }
    }

    display(interfacePosition,"interfacePosition","%12.9f ");
  // OV_ABORT("stop here for now");
    



    std::complex<LocalReal> I(0.0,1.0);

  // kv{m}(1:3,numTransmitted,forward) = [kxt; kvI(2:3)]; 
  //   Vv{m}(1:6,numTransmitted,forward)
             		 
    std::complex<LocalReal> *& Vvc = dbase.get<std::complex<LocalReal>* >("Vvc"); 
    delete [] Vvc;
    Vvc = new  std::complex<LocalReal> [6*2*2*numberOfDomains];  
    #define Vv(i,polar,fb,domain) Vvc[(i)+6*( (polar) + 2*( (fb)+ 2*( domain ) ) ) ]

    std::complex<LocalReal> *& kvc = dbase.get<std::complex<LocalReal>* >("kvc"); 
    delete [] kvc;
    kvc = new  std::complex<LocalReal> [3*2*2*numberOfDomains]; 
    #define kv(i,polar,fb,domain) kvc[(i)+3*( (polar) + 2*( (fb)+ 2*( domain ) ) ) ]

        
  // ---- store reflection/transmission coefficients here : there are 4 cases  ---
  // rtc(i,icase)   i = 4*(numberOfDomains-1), icase=0,1,2,3
  // numScatteringCoeff =  2(left) + 4*(middle domains) + 2(right) = 4*(numDomains-1) 
    std::complex<LocalReal> *& rtc = dbase.get<std::complex<LocalReal>* >("rtc"); 
    delete [] rtc;
    const int numCases=4;  // forward/backward, 2 polar
    int & numScatteringCoeff = dbase.get<int>("numScatteringCoeff");
    numScatteringCoeff = 4*(numberOfDomains-1);
    rtc = new std::complex<LocalReal> [numScatteringCoeff*numCases];  
    #define rtv(i,j) rtc[(i)+numScatteringCoeff*(j)]
    
    std::complex<LocalReal> & s = dbase.get<std::complex<LocalReal> >("s"); 

    s = -I*omega;



    RealArray L0(6,6); L0=0.; //  curl matrix without kx 
//  RealArray Lx(6,6); Lx=0.; //  curl matrix 1's in the positions for kx
        
    std::complex<LocalReal> Lxc[6*6];
    #define Lx(i,j) Lxc[(i)+6*(j)]

    LocalReal kx = 0, ky=kvI(1), kz=kvI(2); 
    L0(ex,hy)=-kz; L0(ex,hz)= ky;
    L0(ey,hz)=-kx; L0(ey,hx)= kz;
    L0(ez,hx)=-ky; L0(ez,hy)= kx;
    L0(hx,ey)= kz; L0(hx,ez)=-ky;
    L0(hy,ez)= kx; L0(hy,ex)=-kz;
    L0(hz,ex)= ky; L0(hz,ey)=-kx;
    
    kx = 1; ky=0; kz=0;
    Lx(ex,hy)=-kz; Lx(ex,hz)= ky;
    Lx(ey,hz)=-kx; Lx(ey,hx)= kz;
    Lx(ez,hx)=-ky; Lx(ez,hy)= kx;
    Lx(hx,ey)= kz; Lx(hx,ez)=-ky;
    Lx(hy,ez)= kx; Lx(hy,ex)=-kz;
    Lx(hz,ex)= ky; Lx(hz,ey)=-kx;
    
    kx = kvI(0); ky=kvI(1); kz=kvI(2); // reset

    

    RealArray K0(6,6);

    std::complex<LocalReal> Kc[6*6];
    #define K(i,j) Kc[(i)+6*(j)]

    std::complex<LocalReal> Bc[6*6];
    #define B(i,j) Bc[(i)+6*(j)]

    std::complex<LocalReal> Ac[6*6];
    #define A(i,j) Ac[(i)+6*(j)]



    int n=6, nrhs=6, lda=6, ldb=6, info;
    IntegerArray ipvt(6);

    std::complex<LocalReal> lambdac[6];
    #define lambda(i) lambdac[i]

    int lwork=8*n;  // >= 2*n
    std::complex<LocalReal> work[lwork];
    LocalReal rwork[2*6];
          
    std::complex<LocalReal> vlc[6*6];
    #define vl(i,j)  vrc[(i)+6*(j)]
    std::complex<LocalReal> vrc[6*6];
    #define V(i,j)  vrc[(i)+6*(j)]

    ComplexArray Ks(6,6);
    
    for( int m=0; m<numberOfDomains; m++ )
    {
        DispersiveMaterialParameters & dmp = dispersiveMaterialParameters[m];
        printF("\n ############# SlabsExactSolution::initialize: domain=%d: name=%s\n",m,(const char*)dmp.getMaterialName());
        dmp.getBianisotropicMaterialMatrix( K0 );
        ::display(K0,"K0","%5.2f ");
        
    // Here is the complex version  *new* *check me*
        Complex ss;
        ss=s;
        dmp.getK( Ks, ss );
        Ks.display();
    // OV_ABORT("stop here for now");
        
    // RealArray & gdmPar = dmp.getBianisotropicGDMParameters();
    // display(gdmPar,"gdmPar","%5.2f "); 

    // Ks = getMaterialMatrix( medium{m},s );

    // %  is*K(s) + kx*Lx + L0   : generalized eigenvalue problem for kx
    // %   (- is*K(s) - L0 ) v = kx Lx v
    // %   =>  B v = kx Lx v 
    // %   =>  B^{-1} Lx v = (1/kx) v   if B is invertible
    // %       A v = (1/kx) v 
    // B = zeros(6,6);
    // B = - ( 1i*s*Ks + L0 );
    // % B 
    // % inv(B)
    // A = inv(B)*Lx;
        
    // [V,D] = eig(A); % V = right eigenvectors 
    // lam = diag(D);    

        for( int j=0; j<6; j++ )
        {
            for( int i=0; i<6; i++ )
            {
	// K(i,j) = K0(i,j);  // finish for dispersive case 
      	K(i,j) = Ks(i,j);  // finish for dispersive case 
            }
        }
        

        for( int j=0; j<6; j++ )
        {
            for( int i=0; i<6; i++ )
            {
      	B(i,j) = -( I*s*K(i,j) + L0(i,j) );
                A(i,j) = Lx(i,j);  // RHS to solve below -- answer returned in A 
            }
        }
        
    // Solve B*A = Lx   for A given B 

        GESV( n, nrhs, B(0,0), lda, ipvt(0), A(0,0), ldb,info );
        if( info==0 )
        {
            printF("info=0: succesful return from  Linear-solver GESV\n");
        }
        else
        {
            printF("ERROR RETURN FROM Linear-solver GESV: info=%d: < argument -info is invalid, >0: matrix is singular\n",info);
            OV_ABORT("error");
        }
        
    // "N"  = do not compute left eigenvectors 
    // "V" = do compute right eigenvectors
        GEEV("N", "V", n, A(0,0), lda, lambda(0), vl(0,0), lda, V(0,0), lda, work[0], lwork, rwork[0], info);


        if( info !=0 )
        {
            printF("ERROR return info=%d from eigenvalue routine DGEEV\n",info);
            OV_ABORT("error");
        }
        for( int i=0; i<6; i++ )
        {
            printF(" lambda(%d)= %12.5e + %12.5e I\n",i,std::real(lambda(i)),imag(lambda(i)));
        }
        

        const LocalReal tol = 100.*sqrt(REAL_EPSILON);  //  tolerance for computing eigenvalues -- double roots may have large errors

    
        int numTransmitted=0; 
        int numReflected  =0; 
        for( int i=0; i<6; i++ )
        {
            
            std::complex<LocalReal> kcInverse = lambda(i); 
            if( abs(kcInverse) > tol )  //  skip zero eigenvalues 
            {
                std::complex<LocalReal> kxt = 1./kcInverse; 
      	if( idebug>1 )
      	{
                    printF(" kxt(%d)=(%12.4e,%12.4e), qr=[%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f]+I [%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f]\n",
              		  i,std::real(kxt),std::imag(kxt),
              		  std::real(V(0,i)),std::real(V(1,i)),std::real(V(2,i)),std::real(V(3,i)),std::real(V(4,i)),std::real(V(5,i)),
              		  std::imag(V(0,i)),std::imag(V(2,i)),std::imag(V(2,i)),std::imag(V(3,i)),std::imag(V(4,i)),std::imag(V(5,i)));
      	}
      	
                if( std::real(kxt) > 0. )
      	{
          // --- transmitted wave ---
                    kv(0,numTransmitted,forward,m) = kxt;      // kx component of wave vector 
            	  kv(1,numTransmitted,forward,m) = kvI(1);   // tangential components match incident 
                    kv(2,numTransmitted,forward,m) = kvI(2);
        	  for( int j=0; j<6; j++ ){ Vv(j,numTransmitted,forward,m) = V(j,i); } // 

                    numTransmitted=numTransmitted+1;
      	}
                else
      	{
          // --- reflected wave ---
                    kv(0,numReflected,backward,m) = kxt;
                    kv(1,numReflected,backward,m) = kvI(1);
                    kv(2,numReflected,backward,m) = kvI(2);
                    for( int j=0; j<6; j++ ){ Vv(j,numReflected,backward,m) = V(j,i); } // 

                    numReflected=numReflected+1;
      	}
            }
            else
            {
                if( idebug>1 )
        	  printF(" lambda(%d)=1/kxt = (%12.4e,%12.4e) (skipped).\n",i,std::real(kcInverse),imag(kcInverse)); 
            }
        }
    // end for i 
            
        if( idebug>1 )
            printF("\n +++++ MEDIUM %d: There were %d reflected and %d transmitted waves.\n",m,numReflected,numTransmitted); 

        if( numTransmitted != 2 )
        {
            printF("ERROR: MEDIUM %d : expected 2 transmitted waves (?!)\n",m); 
            OV_ABORT("This should not normally happen!");
        }
    
    // -- re-order eigenvectors in the first and last media 
    //    the TEz mode (Ex,Ey,Hz) comes before TMz mode (Ez,Hx,Hy) (if these exists)
    //    Here we assume that Matlab has generated the TEz and TMz modes as eigenvectors in isotopic media.

        if( m==0 || m==(numberOfDomains-1) )
        {
            std::complex<LocalReal> temp;
            for( int dir=0; dir<2; dir++ ) //  forward,backward
            {
      	if( abs(Vv(hz,1,dir,m)) > abs(Vv(hz,0,dir,m)) ) //  choose hz component since ex or ey components may be zero depending on the polarization
      	{
          // flip polarizations TEz <--> TMz 
                    if( idebug>1 ) printF(" FLIP eigenvectors (polarization) m=%d, dir=%d in first or last domain (assumed isotropic)\n",m,dir); 
        	  for( int j=0; j<3; j++ )
        	  {
          	    temp = kv(j,0,dir,m);
          	    kv(j,0,dir,m) = kv(j,1,dir,m);
          	    kv(j,1,dir,m) = temp;
        	  }  
        	  for( int j=0; j<6; j++ )
        	  {
          	    temp = Vv(j,0,dir,m);
          	    Vv(j,0,dir,m) = Vv(j,1,dir,m);
          	    Vv(j,1,dir,m) = temp;
        	  }  
        	  
      	}
            }
        }





    }
  // end for m (domain)

    if( idebug>1 )
    {
        std::complex<LocalReal> pqv[6];
        #define qv(i) pqv[(i)]

        printF("\n ----------------- getReflectionAndTransmission: Summary -------------------\n");
        for( int m=0; m<numberOfDomains; m++ )
        {
            printF(" --- domain %d : %s\n",m,(const char*)dispersiveMaterialParameters[m].getMaterialName());
        }
        
        for( int m=0; m<numberOfDomains; m++ )
        {
            for( int i=0; i<2; i++ ) 
            {
      	for( int j=0; j<6; j++ ){ qv(j) = Vv(j,i,backward,m); } // copy for easier printing 
                printF("MED%d: reflected   %d: kv=%5.2f+%5.2fI qv=[%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f]+ [%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f] I\n",
             	       m,i,std::real(kv(0,i,backward,m)),std::imag(kv(0,i,backward,m)), 
                            std::real(qv(0)),std::real(qv(1)),std::real(qv(2)),std::real(qv(3)),std::real(qv(4)),std::real(qv(5)), 
                            std::imag(qv(0)),std::imag(qv(1)),std::imag(qv(2)),std::imag(qv(3)),std::imag(qv(4)),std::imag(qv(5)));
            }
            for( int i=0; i<2; i++ ) 
            {
      	for( int j=0; j<6; j++ ){ qv(j) = Vv(j,i,forward,m); } // 
                printF("MED%d: transmitted %d: kv=%5.2f+%5.2fI qv=[%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f]+ [%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f] I\n",
             	       m,i,std::real(kv(0,i,forward,m)),std::imag(kv(0,i,forward,m)), 
                            std::real(qv(0)),std::real(qv(1)),std::real(qv(2)),std::real(qv(3)),std::real(qv(4)),std::real(qv(5)), 
                            std::imag(qv(0)),std::imag(qv(1)),std::imag(qv(2)),std::imag(qv(3)),std::imag(qv(4)),std::imag(qv(5)));
            }
        }
    }
    
  //
  // ----------- Find the relection and transmission coefficients --------
  //
  // Unknowns : x= [rL1, rL2, rM1, rM2, tauM1, tauM2, tauR1, tauR2]
  // 
  // [ tv.Ev ] =0 and [tv.Hv ] =0 : 
  // x=0: 
  //      SUM rL_j tv.qvr1(:,j) - SUM_j rM_j tv.qvr2(:,j) - SUM_j tauM_j tv.qvt2(:,j) - SUM_j tauR_j tv.qvt3(:,j) = - tv . qvI 
  // 
  // x=L: 
  //     


    RealArray tv(3,2);  //  tangent vectors
    tv(0,0)=0; tv(1,0)=1; tv(2,0)=0; 
    tv(0,1)=0; tv(1,1)=0; tv(2,1)=1; 
    
    
  // number of unknowns is (numDomains-1)*4
    const int nrt=(numberOfDomains-1)*4;
    std::complex<LocalReal> AAc[nrt*nrt];
    #define A(i,j) AAc[(i)+nrt*(j)]

    std::complex<LocalReal> rhsc[nrt*4]; // holds 4 right-hand-sides
    #define rhs(i,j) rhsc[(i)+nrt*(j)]
    
    for( int i=0; i<nrt; i++ )
    {
        for( int j=0; j<nrt; j++ )
        {
            A(i,j)=0.;
        }
    }
    

// ==========================================================================
// Macro:   dotProd
//
//  isign (input) : 1 or -1 to form plus or minus of the dot-product 
//  xi (input) : interface location
//  vc (input) : ex or hx to take dot-product with E or H
//  polar (input) : 0 or 1 -- polarization  
//  fb    (input) : 0 or 1, forward or backward
// 
//  val (output)  : dotProduct 
// ==========================================================================
    

  // LocalReal a=-.5, b=.5;

  // ------------------------------------------------------------------------------------------
  // ----------------------- FIRST INTERFACE:  domain=0 | domain=1    ------------------------
  // ------------------------------------------------------------------------------------------
  //  Match Tangential components of: 
  //    qvI + SUM_j r_{j,d0} exp() V_{j,d0}^-  =  SUM_j r_j exp() V_{j,d1}^- + SUM_j t_j exp() V_{j,d1}^+
  // 
    int domain=0;
    int dLeft=domain, dRight=dLeft+1;
    LocalReal xi = interfacePosition(domain);   // interface position
    for( int eh=0; eh<=1; eh++ ) // E or H 
    {
        int ec = eh==0 ? ex : hx;  // apply tv.E or tv.H 
        
        for( int itan=0; itan<2; itan++ ) // two tangents
        {
            int i = itan + 2*eh + 4*domain;                     // equation i: row in matrix 
                A(i,0) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,0) += tv(ii,itan)*Vv(ec+ii,0,backward,dLeft);
                }
                A(i,0) *= +1.*exp( I*kv(0,0,backward,dLeft)*xi );
                A(i,1) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,1) += tv(ii,itan)*Vv(ec+ii,1,backward,dLeft);
                }
                A(i,1) *= +1.*exp( I*kv(0,1,backward,dLeft)*xi );

                A(i,2) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,2) += tv(ii,itan)*Vv(ec+ii,0,backward,dRight);
                }
                A(i,2) *= -1.*exp( I*kv(0,0,backward,dRight)*xi );
                A(i,3) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,3) += tv(ii,itan)*Vv(ec+ii,1,backward,dRight);
                }
                A(i,3) *= -1.*exp( I*kv(0,1,backward,dRight)*xi );
                        
                A(i,4) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,4) += tv(ii,itan)*Vv(ec+ii,0,forward,dRight);
                }
                A(i,4) *= -1.*exp( I*kv(0,0,forward,dRight)*xi );
                A(i,5) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,5) += tv(ii,itan)*Vv(ec+ii,1,forward,dRight);
                }
                A(i,5) *= -1.*exp( I*kv(0,1,forward,dRight)*xi );
      
      // A(i,6)=0.;                                          // coeff of tR1
      // A(i,7)=0.;                                          // coeff of tR2
              

      // we solve for 4 rhs at once           
                rhs(i,0) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    rhs(i,0) += tv(ii,itan)*Vv(ec+ii,0,forward,dLeft);
                }
                rhs(i,0) *= -1.*exp( I*kv(0,0,forward,dLeft)*xi );
                rhs(i,1) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    rhs(i,1) += tv(ii,itan)*Vv(ec+ii,1,forward,dLeft);
                }
                rhs(i,1) *= -1.*exp( I*kv(0,1,forward,dLeft)*xi );
            rhs(i,2)=0.;                                            // backward, polarization 0
            rhs(i,3)=0.;                                            // backward, polarization 1
          
        }

    } // end for eh 
    
  // ------------------------------------------------------------------------------------------
  // ------------------ INTERFACES BETWEEN MIDDLE SLABS :  domain | domain+1  -----------------
  // ------------------------------------------------------------------------------------------
    for( int domain=1; domain<numberOfDomains-2; domain++ )
    {
        int dLeft=domain, dRight=domain+1;
        xi = interfacePosition(domain);   // interface position

        for( int eh=0; eh<=1; eh++ ) // E or H 
        {
            int ec = eh==0 ? ex : hx;  // apply tv.E or tv.H 
        
            for( int itan=0; itan<2; itan++ ) // two tangents
            {
      	const int i = itan + 2*eh + 4*(domain);               // equation i: ow in matrix  
                const int j =      - 2    + 4*(domain);                // first unknown in this equation 
                A(i,j) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j) += tv(ii,itan)*Vv(ec+ii,0,backward,dLeft);
                }
                A(i,j) *= +1.*exp( I*kv(0,0,backward,dLeft)*xi );
                A(i,j+1) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+1) += tv(ii,itan)*Vv(ec+ii,1,backward,dLeft);
                }
                A(i,j+1) *= +1.*exp( I*kv(0,1,backward,dLeft)*xi );

                A(i,j+2) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+2) += tv(ii,itan)*Vv(ec+ii,0,forward,dLeft);
                }
                A(i,j+2) *= +1.*exp( I*kv(0,0,forward,dLeft)*xi );
                A(i,j+3) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+3) += tv(ii,itan)*Vv(ec+ii,1,forward,dLeft);
                }
                A(i,j+3) *= +1.*exp( I*kv(0,1,forward,dLeft)*xi );

                A(i,j+4) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+4) += tv(ii,itan)*Vv(ec+ii,0,backward,dRight);
                }
                A(i,j+4) *= -1.*exp( I*kv(0,0,backward,dRight)*xi );
                A(i,j+5) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+5) += tv(ii,itan)*Vv(ec+ii,1,backward,dRight);
                }
                A(i,j+5) *= -1.*exp( I*kv(0,1,backward,dRight)*xi );
                        
                A(i,j+6) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+6) += tv(ii,itan)*Vv(ec+ii,0,forward,dRight);
                }
                A(i,j+6) *= -1.*exp( I*kv(0,0,forward,dRight)*xi );
                A(i,j+7) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+7) += tv(ii,itan)*Vv(ec+ii,1,forward,dRight);
                }
                A(i,j+7) *= -1.*exp( I*kv(0,1,forward,dRight)*xi );
  

	// we solve for 4 rhs at once -- rhs are zero for middle domains           
      	rhs(i,0)=0.;
      	rhs(i,1)=0.;
      	rhs(i,2)=0.;                                    
      	rhs(i,3)=0.;                                    
          
            }

        } // end for eh 

    } // end for domain 
    
    
  // ------------------------------------------------------------------------------------------
  // ----------------------- LAST INTERFACE    domain=numDomains-2 | domain=numDomains-1  ------------------------
  // ------------------------------------------------------------------------------------------
    domain = numberOfDomains-2;
    dLeft=domain, dRight=dLeft+1;
    xi = interfacePosition(domain); // interface position
    for( int eh=0; eh<=1; eh++ )  // E or H 
    {
        int ec = eh==0 ? ex : hx;  // apply tv.E or tv.H 
        
        for( int itan=0; itan<2; itan++ ) // two tangents
        {
            const int i = itan + 2*eh + 4*(domain);                // row in matrix 
            const int j =      - 2    + 4*(domain);                // first unknown in this equation 
      // A(i,0) = 0.;                                        // coeff of rL1
      // A(i,1) = 0.;                                        // coeff of rL2 
            
                A(i,j) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j) += tv(ii,itan)*Vv(ec+ii,0,backward,dLeft);
                }
                A(i,j) *= 1.*exp( I*kv(0,0,backward,dLeft)*xi );
                A(i,j+1) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+1) += tv(ii,itan)*Vv(ec+ii,1,backward,dLeft);
                }
                A(i,j+1) *= 1.*exp( I*kv(0,1,backward,dLeft)*xi );

                A(i,j+2) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+2) += tv(ii,itan)*Vv(ec+ii,0,forward,dLeft);
                }
                A(i,j+2) *= 1.*exp( I*kv(0,0,forward,dLeft)*xi );
                A(i,j+3) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+3) += tv(ii,itan)*Vv(ec+ii,1,forward,dLeft);
                }
                A(i,j+3) *= 1.*exp( I*kv(0,1,forward,dLeft)*xi );

                A(i,j+4) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+4) += tv(ii,itan)*Vv(ec+ii,0,forward,dRight);
                }
                A(i,j+4) *= -1.*exp( I*kv(0,0,forward,dRight)*xi );
                A(i,j+5) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    A(i,j+5) += tv(ii,itan)*Vv(ec+ii,1,forward,dRight);
                }
                A(i,j+5) *= -1.*exp( I*kv(0,1,forward,dRight)*xi );

      // we solve for 4 rhs at once           
            rhs(i,0)=0.;                                            // forward, polarization 0
            rhs(i,1)=0.;                                            // forward, polarization 1
                rhs(i,2) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    rhs(i,2) += tv(ii,itan)*Vv(ec+ii,0,backward,dRight);
                }
                rhs(i,2) *= 1.*exp( I*kv(0,0,backward,dRight)*xi );
                rhs(i,3) = 0.;
                for( int ii=0; ii<3; ii++ )
                {
                    rhs(i,3) += tv(ii,itan)*Vv(ec+ii,1,backward,dRight);
                }
                rhs(i,3) *= 1.*exp( I*kv(0,1,backward,dRight)*xi );
          
        }

    } // end for eh 
    
      
    if( true )
    {
        printF(" ---- A(i,j)----\n");
        for( int i=0; i<nrt; i++ )
        {
            printF(" i=%2d:",i);
            for( int j=0; j<nrt; j++ )
      	printF(" %6.2f, ",std::real(A(i,j)));
            printF("\n");
        }
    }
    


  // Solve A*rt = rhs   for A given B 

    if( true )
    {
        int n=nrt, nrhs=4, lda=nrt, ldb=nrt, info;
        ipvt.redim(n);

        GESV( n, nrhs, A(0,0), lda, ipvt(0), rhs(0,0), ldb,info );
    
        if( info==0 )
        {
            printF("info=0: succesful return from  Linear-solver GESV (reflection and transmission voeff)\n");
        }
        else
        {
            printF("ERROR RETURN FROM Linear-solver GESV: info=%d: < argument -info is invalid, >0: matrix is singular\n",info);
            OV_ABORT("error");
        }
    }
    
  // -- store scattering coefficients --- 
    for( int i=0; i<numScatteringCoeff; i++ )
    {
        for( int icase=0; icase<numCases; icase++ )
        {
          rtv(i,icase) = rhs(i,icase);
        }
    }


    if( idebug>1 )
    {
        for( int j=0; j<numCases; j++ ) // 4 cases 
        {
            if( j<2 )
                printF("\n +++++ Case %d: reflection and transmission coefficients : forward : polarization %d +++++\n",j,j);
            else
                printF("\n +++++ Case %d: reflection and transmission coefficients : backward: polarization %d +++++\n",j,j-2);

      // Left domain: 
            int domain=0;
            printF(" d=%d: rL1   =[%6.3f + %6.3f I], rL2   =[%6.3f + %6.3f I]\n",domain,std::real(rhs(0,j)),imag(rhs(0,j)),std::real(rhs(1,j)),imag(rhs(1,j)));

            for( int domain=1; domain<numberOfDomains-1; domain++ )
            {
      	int ic = 2 +4*(domain-1);
                printF(" d=%d: rM1   =[%6.3f + %6.3f I], rM2   =[%6.3f + %6.3f I]\n",domain,std::real(rhs(ic  ,j)),imag(rhs(ic  ,j)),std::real(rhs(ic+1,j)),imag(rhs(ic+1,j)));
      	printF(" d=%d: tauM1 =[%6.3f + %6.3f I], tauM2 =[%6.3f + %6.3f I]\n",domain,std::real(rhs(ic+2,j)),imag(rhs(ic+2,j)),std::real(rhs(ic+3,j)),imag(rhs(ic+3,j)));
            }
            
      // Right domain: 
            domain=numberOfDomains-1;
            printF(" d=%d: tauR1 =[%6.3f + %6.3f I], tauR2 =[%6.3f + %6.3f I]\n",domain,std::real(rhs(6,j)),imag(rhs(6,j)),std::real(rhs(7,j)),imag(rhs(7,j)));
        }
    }
    

  

  // OV_ABORT("SlabsExactSolution::initialize *new*  stop here for now");
      

    return 0;
}



// - //==============================================================================
// - int SlabsExactSolution::
// - SlabsExactSolution::evalTest()
// - {
// -   printF("SlabsExactSolution::evalTest: ENTERING EVAL TEST ****\n");
// - 
// -   
// -   int & numDomains         = dbase.get<int>("numberOfDomains");
// -   int & numScatteringCoeff = dbase.get<int>("numScatteringCoeff");
// - 
// -   numDomains=3;
// -   numScatteringCoeff=8;
// -   
// -   const int forward=0, backward=1;
// - 
// - 
// -   LocalReal sr, si;
// -   RealArray kvI(3); // , kv(3,2,2,numDomains), Vv(6,2,2,numDomains);
// - 
// -   std::complex<LocalReal> I(0.0,1.0);
// - 
// -   LocalReal & a = dbase.get<LocalReal>("a");
// -   LocalReal & b = dbase.get<LocalReal>("b");
// - 
// -   std::complex<LocalReal> *& Vvc = dbase.get<std::complex<LocalReal>* >("Vvc"); 
// -   std::complex<LocalReal> *& kvc = dbase.get<std::complex<LocalReal>* >("kvc"); 
// -   std::complex<LocalReal> *& rtc = dbase.get<std::complex<LocalReal>* >("rtc"); 
// -   std::complex<LocalReal> & s    = dbase.get<std::complex<LocalReal> >("s"); 
// - 
// - #define kv(i,j,k,l) kvc[ (i)+3*( (j)+2*( (k) + 2*(l) ) ) ]
// - #define Vv(i,j,k,l) Vvc[ (i)+6*( (j)+2*( (k) + 2*(l) ) ) ]
// - 
// - /* ---  
// - 
// -   dbase.put<std::complex<LocalReal>* >("Vvc");   
// -   dbase.put<std::complex<LocalReal>* >("kvc");   
// -   dbase.put<std::complex<LocalReal>* >("rtc");   
// -   dbase.put<std::complex<LocalReal> >("s");   
// - 
// -   std::complex<LocalReal> *& Vvc = dbase.get<std::complex<LocalReal>* >("Vvc"); 
// -   Vvc = new std::complex<LocalReal> [6*2*2*numDomains];  // ** delete me ***
// - 
// -   std::complex<LocalReal> *& kvc = dbase.get<std::complex<LocalReal>* >("kvc"); 
// -   kvc = new std::complex<LocalReal> [3*2*2*numDomains];  // ** delete me ***
// - 
// -   std::complex<LocalReal> *& rtc = dbase.get<std::complex<LocalReal>* >("rtc"); 
// -   rtc = new std::complex<LocalReal> [8*4];  // ** delete me ***
// - 
// -   std::complex<LocalReal> & s = dbase.get<std::complex<LocalReal> >("s"); 
// - 
// - //  std::complex<LocalReal> kvc[3*2*2*numDomains];
// - //  std::complex<LocalReal> Vvc[6*2*2*numDomains];
// - #define kv(i,j,k,l) kvc[ (i)+3*( (j)+2*( (k) + 2*(l) ) ) ]
// - #define Vv(i,j,k,l) Vvc[ (i)+6*( (j)+2*( (k) + 2*(l) ) ) ]
// -    
// -   // store complex scattering coeff here: 
// - //   std::complex<LocalReal> rtc[8*4]
// -   std::complex<LocalReal> qvc[6];
// - #define rtv(i,j) rtc[(i)+8*(j)] 
// - #define qv(i) qvc[(i)]
// - 
// - 
// -   // File written by bagdm/oblique/baSlab.m 
// - #include "slabSolution.h"
// - 
// -   s = sr + I*si;
// -   
// -   // std::complex<LocalReal> s(sr,si);
// - 
// - 
// -   printF(" kvI=[%g,%g,%g]*(2*pi)\n",kvI(0)/twoPi,kvI(1)/twoPi,kvI(2)/twoPi);
// -   printF("\n ===== SLAB: kx and eigenvectors in the %d domains =====\n",numDomains);
// -   
// -   for( int d=0; d<numDomains; d++)
// -   {
// -     for( int i=0; i<2; i++ )  // polarizations
// -     {
// -       for( int j=0; j<6; j++ )
// -         qv(j) = Vv(j,i,backward,d);
// - 
// -       fprintf(stdout,"MED%d: reflected   %d: kx=%5.2f+%5.2fI Vv=[%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f]+ [%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f] I\n",    
// -               d,i,std::real(kv(0,i,backward,d)),imag(kv(0,i,backward,d)),    
// -               std::real(qv(0)),std::real(qv(1)),std::real(qv(2)),std::real(qv(3)),std::real(qv(4)),std::real(qv(5)),     
// -               imag(qv(0)),imag(qv(1)),imag(qv(2)),imag(qv(3)),imag(qv(4)),imag(qv(5)));
// -     }
// -     for( int i=0; i<2; i++ )  // polarizations
// -     {
// -        for( int j=0; j<6; j++ )
// -         qv(j) = Vv(j,i,forward,d);
// -        
// -        fprintf(stdout, "MED%d: transmitted %d: kx=%5.2f+%5.2fI Vv=[%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f]+ [%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f] I\n",    
// - 	       d,i,std::real(kv(0,i,forward,d)),imag(kv(0,i,forward,d)),     
// -               std::real(qv(0)),std::real(qv(1)),std::real(qv(2)),std::real(qv(3)),std::real(qv(4)),std::real(qv(5)),     
// -               imag(qv(0)),imag(qv(1)),imag(qv(2)),imag(qv(3)),imag(qv(4)),imag(qv(5)));
// -     }
// -   } // end for d (domain)
// -   
// -   for( int icase=0; icase<4; icase++ ) //  4 cases 
// -   {
// -     if( icase<2 )
// -       fprintf(stdout,"\n +++++ reflection and transmission coefficients : forward : polarization %d +++++\n",icase);
// -     else
// -       fprintf(stdout,"\n +++++ reflection and transmission coefficients : backward: polarization %d +++++\n",icase-2);
// - 
// -     std::complex<LocalReal> x0,x1,x2,x3,x4,x5,x6,x7;
// -     x0 = rtv(0,icase);
// -     x1 = rtv(1,icase);
// -     x2 = rtv(2,icase);
// -     x3 = rtv(3,icase);
// -     x4 = rtv(4,icase);
// -     x5 = rtv(5,icase);
// -     x6 = rtv(6,icase);
// -     x7 = rtv(7,icase);
// - 
// -     fprintf(stdout," rL1   =[%6.3f + %6.3f I], rL2   =[%6.3f + %6.3f I]\n",std::real(x0),imag(x0),std::real(x1),imag(x1));
// -     fprintf(stdout," rM1   =[%6.3f + %6.3f I], rM2   =[%6.3f + %6.3f I]\n",std::real(x2),imag(x2),std::real(x3),imag(x3));
// -     fprintf(stdout," tauM1 =[%6.3f + %6.3f I], tauM2 =[%6.3f + %6.3f I]\n",std::real(x4),imag(x4),std::real(x5),imag(x5));
// -     fprintf(stdout," tauR1 =[%6.3f + %6.3f I], tauR2 =[%6.3f + %6.3f I]\n",std::real(x6),imag(x6),std::real(x7),imag(x7));
// -   }
// -   
// -   --- */
// - 	
// - 
// -   // ::display(Vv,"eigenvectors");
// - 
// -   LocalReal x=-1., y=0, z=0;
// -   LocalReal t=.1;
// - 
// -   Range M=6;
// -   RealArray q(M);
// - 
// - 
// -   int icase=0;  // forward-scat, TEz-polarization
// -   const int pf=0, pb=0;    // polarization for forward and backward INCIDENT states (depends on icase)
// -   LocalReal a0=1., a2=0.;  // INCIDENT field amplitudes (depends on icase)
// - 
// - 
// -   std::complex<LocalReal> rL1,rL2,tL1,tL2, rM1,rM2,tM1,tM2, rR1, rR2, tR1,tR2;
// - 
// -   rL1=rtv(0,icase);
// -   rL2=rtv(1,icase);
// - 
// -   rM1=rtv(2,icase);
// -   rM2=rtv(3,icase);
// -   tM1=rtv(4,icase);
// -   tM2=rtv(5,icase);
// - 
// -   tR1=rtv(6,icase);
// -   tR2=rtv(7,icase);
// - 
// -   const int d0=0, d1=1, d2=2;  // 3 domains 
// -    
// -   if( x<=a )
// -   {
// -     for( int m=0; m<6; m++ )
// -     {
// -       q(m) = std::real(
// - 	a0 *exp(I*(kv(0,pf,forward ,d0)*x + kv(1,pf,forward ,d0)*y) + s*t )*Vv(m,pf,forward ,d0)  +    
// - 	rL1*exp(I*(kv(0, 0,backward,d0)*x + kv(1, 0,backward,d0)*y) + s*t )*Vv(m, 0,backward,d0) +    
// - 	rL2*exp(I*(kv(0, 1,backward,d0)*x + kv(1, 1,backward,d0)*y) + s*t )*Vv(m, 1,backward,d0)
// - 	); 
// -     }
// -   }
// -   else if( x <= b )
// -   {
// -     for( int m=0; m<6; m++ )
// -     {
// -       q(m) = std::real(
// - 	rM1*exp(I*(kv(0, 0,backward,d1)*x + kv(1, 0,backward,d1)*y) + s*t )*Vv(m,0,backward,d1) +    
// - 	rM2*exp(I*(kv(0, 1,backward,d1)*x + kv(1, 1,backward,d1)*y) + s*t )*Vv(m,1,backward,d1) +    
// - 	tM1*exp(I*(kv(0, 0,forward ,d1)*x + kv(1, 0,forward ,d1)*y) + s*t )*Vv(m,0,forward ,d1) +    
// - 	tM2*exp(I*(kv(0, 1,forward ,d1)*x + kv(1, 1,forward ,d1)*y) + s*t )*Vv(m,1,forward ,d1) ); 
// -     }
// -     
// -   }
// -   else
// -   {
// -     for( int m=0; m<6; m++ )
// -     {
// -       q(m) = std::real(
// - 	tR1*exp(I*(kv(0, 0,forward ,d2)*x + kv(1, 0,forward ,d2)*y) + s*t )*Vv(m, 0,forward ,d2) +     
// - 	tR2*exp(I*(kv(0, 1,forward ,d2)*x + kv(1, 1,forward ,d2)*y) + s*t )*Vv(m, 1,forward ,d2) +    
// - 	a2 *exp(I*(kv(0,pb,backward,d2)*x + kv(1,pb,backward,d2)*y) + s*t )*Vv(m,pb,backward,d2)   );
// -     }
// -     
// -   }
// -   
// - 
// -   ::display(q,"q");
// -   
// -       
// -   return 0;
// - 
// - }



// ==========================================================================================
/// \brief  Evaluate the solution and save in an array.
///
/// \param numberOfTimeDerivatives (input) : evaluate this many time-derivatives of the solution.
/// \param computeMagneticField (input): if true return the magnetic field in 3D (in 2D the magnetic field is always computed). 
// ==========================================================================================
int SlabsExactSolution::
eval(LocalReal t, CompositeGrid & cg, int grid, 
          realArray & ua, realArray & pv,
          const Index & I1a, const Index &I2a, const Index &I3a, 
          int numberOfTimeDerivatives /* = 0 */,
          bool computeMagneticField /* = false */ )
{

    const int & numberOfDomains    = dbase.get<int>("numberOfDomains");
    const int & numScatteringCoeff = dbase.get<int>("numScatteringCoeff");
    const int & solveForAllFields  = dbase.get<int>("solveForAllFields");
    const int & scatteringCase     = dbase.get<int>("scatteringCase");


  // domain number for this grid: 
//   const int myDomain = cg.domainNumber(grid);
    const int myDomain = grid; // **** TEMP ****


    if( t <= 0. )
        printF("--SlabsExactSolution--  eval on grid=%i, domain=%i, at t=%9.3e (numberOfDomains=%d) \n",grid,myDomain,t,numberOfDomains);

    MappedGrid & mg = cg[grid];
    const int numberOfDimensions = cg.numberOfDimensions();
    
    
    std::vector<DispersiveMaterialParameters> & dmpVector = *dbase.get<std::vector<DispersiveMaterialParameters>* >("pDispersiveMaterialParameters");

//    parameters.dbase.put<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");

    DispersiveMaterialParameters & dmp    = dmpVector[myDomain];
    const int numberOfPolarizationVectors = dmp.numberOfPolarizationVectors;
    RealArray & gdmPar                    = dmp.modelParameters;
    if( false && dmp.isDispersiveMaterial() )
    {
    //   OV_ABORT("SlabsExactSolution::eval -- finish me for dispersive materials");
        printF("SES:eval: grid=%d: numberOfPolarizationVectors=%d\n",grid,numberOfPolarizationVectors);
        ::display(gdmPar,"gdmPar");
    }
    
  // ** finish me for BA ***
  // const IntegerArray & NpBA1                  = dmp1.dbase.get<IntegerArray>("NpBA");
  // const RealArray & bianisotropicParameters1  = dmp1.dbase.get<RealArray>("bianisotropicParameters");

    
    OV_GET_SERIAL_ARRAY(real,ua,uLocal);

    Index I1=I1a, I2=I2a, I3=I3a;
    bool ok = ParallelUtility::getLocalArrayBounds(ua,uLocal,I1,I2,I3,1);   
    if( !ok ) return 0;  // no points on this processor (NOTE: no communication should be done after this point)

  // -- we optimize for Cartesian grids (we can avoid creating the vertex array)
    const bool isRectangular=mg.isRectangular();
    if( !isRectangular )
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

    OV_GET_SERIAL_ARRAY_CONDITIONAL(real,mg.center(),xLocal,!isRectangular); // *wdh* added conditional, Nov 21, 2020


    LocalReal dvx[3]={1.,1.,1.}, xab[2][3]={{0.,0.,0.},{0.,0.,0.}};
    int iv0[3]={0,0,0}; //
    int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];  // NOTE: iv[0]==i1, iv[1]==i2, iv[2]==i3
    LocalReal xv[3]={0.,0.,0.};
    if( isRectangular )
    {
        mg.getRectangularGridParameters( dvx, xab );
        for( int dir=0; dir<mg.numberOfDimensions(); dir++ )
        {
            iv0[dir]=mg.gridIndexRange(0,dir);
            if( mg.isAllCellCentered() )
      	xab[0][dir]+=.5*dvx[dir];  // offset for cell centered
        }
    }
  // This macro defines the grid points for rectangular grids:
#undef XC
#define XC(iv,axis) (xab[0][axis]+dvx[axis]*(iv[axis]-iv0[axis]))

    
  // -- Store components here: 
    const int ex=0, ey=1, ez=2;
    const int hx=3, hy=4, hz=numberOfDimensions==2 ? 2 :  5;

    if( computeMagneticField && numberOfDimensions==3 && ua.getLength(3)<6 )
    {
        printF(" SlabsExactSolution::ERROR: Not enough space in ua to hold the H field\n");
        OV_ABORT("error");
    }
    

  // --- Get Arrays for the dispersive model ----

    RealArray pLocal;
    if( numberOfPolarizationVectors>0 )
    {
        OV_GET_SERIAL_ARRAY(real, pv,pLoc);
        pLocal.reference(pLoc);
    }

    std::complex<LocalReal> I(0.0,1.0);

    std::complex<LocalReal> *& Vvc = dbase.get<std::complex<LocalReal>* >("Vvc"); 
    std::complex<LocalReal> *& kvc = dbase.get<std::complex<LocalReal>* >("kvc"); 
    std::complex<LocalReal> *& rtc = dbase.get<std::complex<LocalReal>* >("rtc"); 
    std::complex<LocalReal> & s    = dbase.get<std::complex<LocalReal> >("s"); 

    Range M=6;

  // const int numDomains=3;
    const int forward=0, backward=1;

    const int icase=scatteringCase;  
    int pf=0, pb=0;          // polarization for forward and backward INCIDENT states (depends on icase)
    LocalReal a0=0., a2=0.;  // INCIDENT field amplitudes (depends on icase)
    if( scatteringCase==0 )
    { // forward-scat, TEz-polarization
        a0=1.; 
    }
    else if( scatteringCase==1 )
    { // forward-scat, TMz-polarization
        a0=1.; pf=1;
    }
    else if( scatteringCase==2 )
    { // back-scat, TEz-polarization
        a2=1.;
    }
    else if( scatteringCase==3 )
    { // back-scat, TMz-polarization
        a2=1.; pf=1;
    }
    else
    {
        OV_ABORT("SES::eval:ERROR: scatteringCase");
    }
    


    std::complex<LocalReal> rL1,rL2,tL1,tL2, rM1,rM2,tM1,tM2, rR1, rR2, tR1,tR2;

    std::complex<LocalReal> ratio;
    std::complex<LocalReal> qc[6];
    #define q(i) qc[(i)]
    
  // left domain
    if( myDomain==0 )
    {
        rL1=rtv(0,icase);
        rL2=rtv(1,icase);
    }
    
  // middle domains
    int ic = 2 + (myDomain-1)*4;
    if( myDomain>0 )
    {
        rM1=rtv(ic  ,icase);
        rM2=rtv(ic+1,icase);
        tM1=rtv(ic+2,icase);
        tM2=rtv(ic+3,icase);
    }
    
  // right domain 
    if( myDomain == (numberOfDomains-1) )
    {
    // tR1=rtv(6,icase);
    // tR2=rtv(7,icase);
        tR1=rtv(ic  ,icase);
        tR2=rtv(ic+1,icase);
    }
    
  // const int d0=0, d1=1, d2=2;  // 3 domains 

    LocalReal x0[3]={0.,0.,0.};   //     
    LocalReal x,y,z=0.;
    if( numberOfTimeDerivatives==0 )
    {
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
            if( !isRectangular )
            {
      	x= xLocal(i1,i2,i3,0)-x0[0];   // shift point to reference coordinates 
      	y= xLocal(i1,i2,i3,1)-x0[1];
      	if( numberOfDimensions==3 ) z= xLocal(i1,i2,i3,2)-x0[2];
            }
            else
            {
      	x=XC(iv,0)-x0[0];
      	y=XC(iv,1)-x0[1];
      	if( numberOfDimensions==3 ) z=XC(iv,2)-x0[2];
            }

            const int d=myDomain;
      	
            if( myDomain==0 )
            {
      	for( int m=0; m<6; m++ )
      	{
        	  q(m) = (
          	    a0 *exp(I*(kv(0,pf,forward ,d)*x + kv(1,pf,forward ,d)*y + kv(2,pf,forward ,d)*z) + s*t )*Vv(m,pf,forward ,d)  +    
          	    rL1*exp(I*(kv(0, 0,backward,d)*x + kv(1, 0,backward,d)*y + kv(2, 0,backward,d)*z) + s*t )*Vv(m, 0,backward,d) +    
          	    rL2*exp(I*(kv(0, 1,backward,d)*x + kv(1, 1,backward,d)*y + kv(2, 1,backward,d)*z) + s*t )*Vv(m, 1,backward,d)
          	    ); 
      	}

            }
            else if( myDomain < numberOfDomains-1 )
            {
      	for( int m=0; m<6; m++ )
      	{
        	  q(m) = (
          	    rM1*exp(I*(kv(0, 0,backward,d)*x + kv(1, 0,backward,d)*y + kv(2, 0,backward,d)*z) + s*t )*Vv(m,0,backward,d) +    
          	    rM2*exp(I*(kv(0, 1,backward,d)*x + kv(1, 1,backward,d)*y + kv(2, 1,backward,d)*z) + s*t )*Vv(m,1,backward,d) +    
          	    tM1*exp(I*(kv(0, 0,forward ,d)*x + kv(1, 0,forward ,d)*y + kv(2, 0,forward ,d)*z) + s*t )*Vv(m,0,forward ,d) +    
          	    tM2*exp(I*(kv(0, 1,forward ,d)*x + kv(1, 1,forward ,d)*y + kv(2, 1,forward ,d)*z) + s*t )*Vv(m,1,forward ,d) ); 
      	}
            }
            else if( myDomain==(numberOfDomains-1) )
            {
      	for( int m=0; m<6; m++ )
      	{
        	  q(m) = (
          	    tR1*exp(I*(kv(0, 0,forward ,d)*x + kv(1, 0,forward ,d)*y + kv(2, 0,forward ,d)*z) + s*t )*Vv(m, 0,forward ,d) +     
          	    tR2*exp(I*(kv(0, 1,forward ,d)*x + kv(1, 1,forward ,d)*y + kv(2, 1,forward ,d)*z) + s*t )*Vv(m, 1,forward ,d) +    
          	    a2 *exp(I*(kv(0,pb,backward,d)*x + kv(1,pb,backward,d)*y + kv(2,pb,backward,d)*z) + s*t )*Vv(m,pb,backward,d)   );
      	}

            }
            else
            {
      	OV_ABORT("SES: error: myDomain");
            }
            if( solveForAllFields )
            {
      	uLocal(i1,i2,i3,ex) = std::real(q(0));
      	uLocal(i1,i2,i3,ey) = std::real(q(1));
      	uLocal(i1,i2,i3,ez) = std::real(q(2));
      	uLocal(i1,i2,i3,hx) = std::real(q(3));
      	uLocal(i1,i2,i3,hy) = std::real(q(4));
      	uLocal(i1,i2,i3,hz) = std::real(q(5));
            }
            else if( numberOfDimensions==2 )
            {
      	uLocal(i1,i2,i3,ex) = std::real(q(0));
      	uLocal(i1,i2,i3,ey) = std::real(q(1));
      	uLocal(i1,i2,i3,hz) = std::real(q(5));
            }
            else if( numberOfDimensions==3 )
            {
      	uLocal(i1,i2,i3,ex) = std::real(q(0));
      	uLocal(i1,i2,i3,ey) = std::real(q(1));
      	uLocal(i1,i2,i3,ez) = std::real(q(2));
            }

            if( numberOfPolarizationVectors>0 )
            {
        // --- assign polarization vectors ---
      	for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
      	{
        	  const int pc= iv*numberOfDimensions;
                    LocalReal a0 = gdmPar(0,iv), a1=gdmPar(1,iv), b0=gdmPar(2,iv), b1=gdmPar(3,iv);
                    ratio =  (a0 + a1*s)/(b0 + b1*s + s*s);
        	  pLocal(i1,i2,i3,pc  ) = std::real( q(0)*ratio );  // Px 
        	  pLocal(i1,i2,i3,pc+1) = std::real( q(1)*ratio );  // Py 
                    
                    if( solveForAllFields )
        	  {
          	    OV_ABORT("SES: finish me -- dispersive, solveForAllFields");
        	  }
        	  else if( numberOfDimensions==2 )
        	  {
	    // ** fix me for Hz 
        	  }
        	  else if( numberOfDimensions==3 )
        	  {
                        pLocal(i1,i2,i3,pc+2) = std::real( q(2)*ratio );  // Pz 
        	  }
        	  
      	}
            }
      // -------- BA Material -----
      // ** finish me for BA ***  see also BAPlaneInterface.bC
      // if( method==bamx )
      // {
      // 	const IntegerArray & NpBA = mt==0 ? NpBA1 : NpBA2;
      // 	const RealArray & bianisotropicParameters = mt==0 ? bianisotropicParameters1 : bianisotropicParameters2;
      // 	int m=0;
      // 	for( int k1=0; k1<6; k1++ )
      // 	{
      // 	  for( int k2=0; k2<6; k2++ )
      // 	  {
      // 	    int ec=k2;  // *check me**
      // 	    for( int n=0; n<NpBA(k1,k2); n++ )
      // 	    {
      // 	      LocalReal a0 = bianisotropicParameters(0,n,k1,k2);
      // 	      LocalReal a1 = bianisotropicParameters(1,n,k1,k2);
      // 	      LocalReal b0 = bianisotropicParameters(2,n,k1,k2);
      // 	      LocalReal b1 = bianisotropicParameters(3,n,k1,k2);

      // 	      pijm = (a0+a1*s)/(b0+b1*s+s*s) * q(ec);
      // 	      pLocal(i1,i2,i3,m) = std::real(pijm);   // Pijm 
      // 	      m++;   
      // 	      pLocal(i1,i2,i3,m) = std::real(s*pijm); // Q = dP/dt
      // 	      m++;

      // 	    }
      // 	  }
      // 	}
      // }	
            

        }  // end for i1,i2,i3


    }
    else
    {
        OV_ABORT("ERROR: numberOfTimeDerivatives != 0 ");
        
    } // end if number of time derivatives 
        
    

    return 0;

}





// ===============================================================================
/// \brief Check the solution.
// ===============================================================================
int SlabsExactSolution::
SlabsExactSolution::check()
{

  // ------------- CHECK THAT THE EQUATIONS ARE SATISFIED AT POINTS INSIDE AND OUTSIDE ------

    printF("------------ SlabsExactSolution::check: CHECK THE EQUATIONS ------------\n\n");

    LocalReal maxErr=0.;

    printF("\n------------ FINSHED CHECK EQUATIONS ------------\n\n");


    
    return 0;
}


