// This file automatically generated from PlaneInterfaceExactSolution.bC with bpp.
#include "PlaneInterfaceExactSolution.h"

#include "DispersiveMaterialParameters.h"

#include "ParallelUtility.h"

// ===============================================================================
// Class to define exact solutions to Maxwell's equations
// 
//     Scattering from a plane material interface in 2D and 3D , dispersive or not
// ===============================================================================



// #define scatSphere EXTERN_C_NAME(scatsphere)

// extern "C"
// {

// void scatSphere(const int&nd ,
// 	     const int&n1a,const int&n1b,const int&n2a,const int&n2b,const int&n3a,const int&n3b,const int&nd1a,
// 	     const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,const int&nd4a,const int&nd4b,
// 	     const real& xy, real&u, const int&ipar, const real&rpar );

// }


#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++)                                       for(i2=I2Base; i2<=I2Bound; i2++)                                     for(i1=I1Base; i1<=I1Bound; i1++)



#define evalMaterialInterfaceSolution3d EXTERN_C_NAME(evalmaterialinterfacesolution3d)
extern "C"
{

void  evalMaterialInterfaceSolution3d( const real&sr,const real&si, 
                                                                              const real&chiSum1r, const real&chiSum1i, const real&chiSum2r, const real&chiSum2i,
                                                                              const real&av,const real&nv,const real&kv,
                                                                              const real&eps1,const real&eps2,
                                                                              const real&mu1,const real&mu2,
                                                                              const real&krr,const real&kri,
                                                                              const real&ktr,const real&kti,
                                                                              real & arr, real & ari, real & atr, real & ati );

}




// ===============================================================================
/// \brief  Constructor for th class that defines exact solutions to Maxwell's equations for a sphere
// ===============================================================================
PlaneInterfaceExactSolution::
PlaneInterfaceExactSolution()
{

    dbase.put<int>("numberOfDimensions");

    dbase.put<real>("amp");

    dbase.put<real[2]>("s"); 

    dbase.put<real>("eps1"); 
    dbase.put<real>("mu1"); 

    dbase.put<real>("eps2"); 
    dbase.put<real>("mu2"); 

    
    dbase.put<real[3]>("av");

    dbase.put<real[3]>("arr");
    dbase.put<real[3]>("ari");
    dbase.put<real[3]>("atr");
    dbase.put<real[3]>("ati");

    dbase.put<real[3]>("krr");
    dbase.put<real[3]>("kri");
    dbase.put<real[3]>("ktr");
    dbase.put<real[3]>("kti");

    dbase.put<real[3]>("kv");
    dbase.put<real[3]>("kvr");
    dbase.put<real[3]>("kvi");

    dbase.put<real[3]>("kpvr");
    dbase.put<real[3]>("kpvi");


    dbase.put<int>("numberOfPolarizationVectors1");

    dbase.put<int>("numberOfPolarizationVectors2");

    dbase.put<real>("eps1Hatr");
    dbase.put<real>("eps1Hati");

    dbase.put<real>("eps2Hatr");
    dbase.put<real>("eps2Hati");

    dbase.put<real>("chiSum1r");
    dbase.put<real>("chiSum1i");

    dbase.put<real>("chiSum2r");
    dbase.put<real>("chiSum2i");

    dbase.put<real[10]>("chi1r");
    dbase.put<real[10]>("chi1i");

    dbase.put<real[10]>("chi2r");
    dbase.put<real[10]>("chi2i");

    dbase.put<real[3]>("normalPlaneMaterialInterface");
    dbase.put<real[3]>("x0PlaneMaterialInterface");
    
}


// ===============================================================================
/// \brief destructor for the class that defines exact solutions to Maxwell's equations for a sphere
// ===============================================================================
PlaneInterfaceExactSolution::
~PlaneInterfaceExactSolution()
{
}


// ===============================================================================
/// \brief Initialize the plane material interface solution.
// ===============================================================================
int PlaneInterfaceExactSolution::
initialize( CompositeGrid & cg, DispersiveMaterialParameters & dmp1, DispersiveMaterialParameters & dmp2,
          	    real *av_, real *kvr_, real *kvi_ )
{
    int & numberOfDimensions = dbase.get<int>("numberOfDimensions");
    numberOfDimensions=cg.numberOfDimensions();


    assert( cg.numberOfDomains()==2 );

    real & a = dbase.get<real>("amp");
    a =1.;  // amplitude 

  // -- Find the normal to the interface and a point on it --
    int grid=0, side=1, axis=0;
    
    MappedGrid & mg = cg[grid];
    mg.update(MappedGrid::THEcenter | MappedGrid::THEvertex | MappedGrid::THEvertexBoundaryNormal);  
    OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);
    OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
    
  // ** DO THIS FOR NOW **
    real kx=kvr_[0], ky=kvr_[1], kz=kvr_[2];

    real (&av)[3] = dbase.get<real[3]>("av");
    for( int i=0; i<3; i++ ){ av[i]=av_[i]; }  // 
    
    real kNorm = sqrt( SQR(kvr_[0]) + SQR(kvr_[1]) + SQR(kvr_[2]) + SQR(kvi_[0]) + SQR(kvi_[1]) + SQR(kvi_[2]) );
    real aNorm = sqrt( SQR(av[0]) + SQR(av[1]) + SQR(av[2]) );

    real aDotK = av[0]*kx + av[1]*ky + av[2]*kz;
    if( fabs(aDotK) > REAL_EPSILON*100. * max(aNorm,kNorm) )
    {

        printF("PlaneInterfaceExactSolution::initialize:ERROR: a.k != 0\n");
        printF(" av=[%g,%g,%g] kv=[%g,%g,%g]\n",av[0],av[1],av[2],kx,ky,kz);
        OV_ABORT("error");
    }
    

    real (&normalPlaneMaterialInterface)[3] = dbase.get<real[3]>("normalPlaneMaterialInterface");
    real (&x0PlaneMaterialInterface)[3] = dbase.get<real[3]>("x0PlaneMaterialInterface");

    real *nv = normalPlaneMaterialInterface;
    real *x0 = x0PlaneMaterialInterface;
    for( int i=0; i<3; i++ ){ nv[i]=0.; x0[i]=0.; }  // 
    
    const IntegerArray & gid = mg.gridIndexRange();
  // mid point of a face: 
    int i1 = gid(1,0), i2= int( (gid(0,1)+gid(1,1))/2 ), i3= int( (gid(0,2)+gid(1,2))/2 );
    for( int axis=0; axis<cg.numberOfDimensions(); axis++ )
    {
        nv[axis]=normal(i1,i2,i3,axis);
        x0[axis]=xLocal(i1,i2,i3,axis);
    }


    printF("\n ------------------- PlaneInterfaceExactSolution::initialize ------------------------\n\n");
    printF(" av=[%g,%g,%g] x0=[%g,%g,%g] nv=[%g,%g,%g]\n",av[0],av[1],av[2],x0[0],x0[1],x0[2],nv[0],nv[1],nv[2]);

  // Should n.k be always positive ??
    real nDotK = nv[0]*kx + nv[1]*ky + nv[2]*kz;
    printF(" n.k = %g\n",nDotK);
    if( nDotK < 0 )
    {
        printF("\n\n *****************************************************************************\n");
        printF("PlaneInterfaceExactSolution::initialize: WARNING n.k < 0 -- this may not work\n");
        printF(" *****************************************************************************\n\n");
    }
    
    

    
    real & eps1 = dbase.get<real>("eps1");
    real & eps2 = dbase.get<real>("eps2");
    real & mu1 = dbase.get<real>("mu1");
    real & mu2 = dbase.get<real>("mu2");
    

    eps1 = dmp1.getEpsInf();
    eps2 = dmp2.getEpsInf();
    mu1=1.;  // do this for now 
    mu2=1.;
    

    printF(" eps1=%g, eps2=%g, [kx,ky,kz]=[%g,%g,%g] \n",eps1,eps2,kx,ky,kz);



  // ----------------------------------------------------
  // ---- DISPERSIVE PLANE WAVE MATERIAL INTERFACE ------
  // ----------------------------------------------------
  // NOTES:
  //    (1) incident wave number is given --> compute s=sr + I*si 
  //    (2) Given s, compute wave number in right state


    LocalReal (&s)[2] = dbase.get<LocalReal[2]>("s");

    real & sr = s[0], &si = s[1];

    real (&kv)[3] = dbase.get<real[3]>("kv");
    real (&kvr)[3] = dbase.get<real[3]>("kvr");
    real (&kvi)[3] = dbase.get<real[3]>("kvi");

    real (&kpvr)[3] = dbase.get<real[3]>("kpvr");
    real (&kpvi)[3] = dbase.get<real[3]>("kpvi");

    real &kxr =kvr[0], &kxi=kvi[0], &kyr=kvr[1], &kyi=kvi[1];          // Incident wave number (complex)
    real &kxpr =kpvr[0], &kxpi=kpvi[0], &kypr=kpvr[1], &kypi=kvi[1];   // 
    
    int & numberOfPolarizationVectors1 = dbase.get<int>("numberOfPolarizationVectors1");
    int & numberOfPolarizationVectors2 = dbase.get<int>("numberOfPolarizationVectors2");

    real (&chi1r)[10] = dbase.get<real[10]>("chi1r");
    real (&chi1i)[10] = dbase.get<real[10]>("chi1i");

    real (&chi2r)[10] = dbase.get<real[10]>("chi2r");
    real (&chi2i)[10] = dbase.get<real[10]>("chi2i");



  // real chi1r[10],chi1i[10];
  // real chi2r[10],chi2i[10];

    real & chiSum1r = dbase.get<real>("chiSum1r");
    real & chiSum1i = dbase.get<real>("chiSum1i");
    
    real & chiSum2r = dbase.get<real>("chiSum2r");
    real & chiSum2i = dbase.get<real>("chiSum2i");
    

    chiSum1r=0.; chiSum1i=0;
    chiSum2r=0.; chiSum2i=0;
  
      
    const int gridLeft = 0;
    const int gridRight=cg.numberOfComponentGrids()-1;

    real c1=1./sqrt(eps1*mu1);  // incident 
    real c2=1./sqrt(eps2*mu2);  // transmitted

    int domain=0;
    numberOfPolarizationVectors1=dmp1.numberOfPolarizationVectors;
    assert( numberOfPolarizationVectors1<10 );

    kxr=twoPi*kx; kxi=0.; kyr=twoPi*ky; kyi=0.;  // Incident wave number (complex)

    kv[0]=twoPi*kx; kv[1]=twoPi*ky; kv[2]=twoPi*kz;  // is this needed?
    

    const real kk = twoPi*sqrt( kx*kx+ky*ky+kz*kz );   
    dmp1.evaluateDispersionRelation( c1,kk, sr, si, chi1r,chi1i,chiSum1r,chiSum1i ); 
  // si = -si; // reverse the direction NO -- changes chi1 !

                  
  // -- -right domain --
    domain=1;
    numberOfPolarizationVectors2=dmp2.numberOfPolarizationVectors;
    assert( numberOfPolarizationVectors2<10 );
      
    real kr,ki;
    dmp2.evaluateComplexWaveNumber( c2,sr,si, kr,ki, chi2r,chi2i,chiSum2r,chiSum2i );
  //  kxp^2 + kyp^2 = (kr+I*ki)^2 = (kr^2-ki^2) + 2*I*kr*ki 
  // kxp = kxpr + I*kpri = sqrt( (kr+I*ki)^2 - kyp^2 )
    getTransmisionWaveNumber( kr,ki, kxr,kxi, kyr,kyi, kxpr,kxpi, kypr,kypi );
  // // do this for now -- assume normal incidence
  // assert( ky==0. );
  // kxpr=kr; kxpi=ki;
  // kypr=0.; kypi=0.;
      

    real & eps1Hatr = dbase.get<real>("eps1Hatr");
    real & eps1Hati = dbase.get<real>("eps1Hati");
    
    real & eps2Hatr = dbase.get<real>("eps2Hatr");
    real & eps2Hati = dbase.get<real>("eps2Hati");

    eps1Hatr = eps1*(1.+chiSum1r); eps1Hati=eps1*(chiSum1i);
    eps2Hatr = eps2*(1.+chiSum2r); eps2Hati=eps2*(chiSum2i);

    if( true )
    {
        printF(" s=(%16.10e,%16.10e) kx=(%16.10e,%16.10e) ky=(%16.10e,%16.10e) \n"
         	   "   -> k2=(kr,ki)=(%16.10e,%16.10e) \n"
         	   " Right: kxp=(%16.10e,%16.10e) kyp=(%16.10e,%16.10e)\n"
         	   ,sr,si,kxr,kxi,kyr,kyi,kr,ki,kxpr,kxpi,kypr,kypi);
        for( int i=0; i<numberOfPolarizationVectors1; i++ )
            printF("    chi1=(%16.10e,%16.10e) \n",chi1r[i],chi1i[i]);
        for( int i=0; i<numberOfPolarizationVectors2; i++ )
            printF("    chi2=(%16.10e,%16.10e) \n",chi2r[0],chi2i[0]);
            
    }
        
    if( cg.numberOfDimensions()==2 )
    {
    // 2D -- 
    //    sr,si : s= sr + I*si 
    //    kxr,kxi,  kyr,kyi,     : complex wave number on left
    //    kxpr,kxpi,  kypr,kypi, : complex wave number on right (plus)


        
        real x=.5, y=.5;
        real t=0.;

// File generated by overtureFramework/cg/mx/codes/dispersivePlaneWaveInterface.maple
// File generated by DropBox/DMX/codes/dispersivePlaneWaveInterface.maple
real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t530,t531,t532,t533,t534,t535,t536,t537,t538,t539,t540,t541,t542,t543,t544,t545,t546,t547,t548,t549,t550,t551,t552,t553,t554,t555,t556,t557,t558,t559,t560,t561,t562,t563,t564,t565,t566,t567,t568,t569,t570,t571,t572,t573,t574,t575,t576,t577,t578,t579,t580,t581,t582,t583,t584,t585,t586,t587,t588,t589,t590,t591,t592,t593,t594,t595,t596,t597,t598,t599,t600;
// -------------------------------------------------------------------------
// Need: 
//    a = amplitude of the wave, e.g. a=1
//    x,y,t
//    sr,si : s= sr + I*si 
//    alphaP
//    kxr,kxi,  kyr,kyi,     : complex wave number on left
//    kxpr,kxpi,  kypr,kypi, : complex wave number on right (plus)
// -------------------------------------------------------------------------
// Evaluated:                                                              
//    Exr,Eyr   : left state                                                
//    Expr,Eypr : right state                                               
// -------------------------------------------------------------------------
real kNorm = sqrt( kxr*kxr + kxi*kxi + kyr*kyr + kyi*kyi);
real khxr = kxr/kNorm, khxi = kxi/kNorm, khyr=kyr/kNorm, khyi=kyi/kNorm; 
real kpNorm = sqrt( kxpr*kxpr + kxpi*kxpi + kypr*kypr + kypi*kypi);
real khxpr = kxpr/kpNorm, khxpi=kxpi/kpNorm, khypr=kypr/kpNorm, khypi=kypi/kpNorm; 
real kappar,kappai, betar, betai, rr,ri, taur,taui;
real Exr,Eyr,Hzr, Expr,Eypr,Hzpr;
real Exi,Eyi,Hzi, Expi,Eypi,Hzpi;
t1 = khxi * khxpi;
t2 = khxr * khxpr;
t4 = pow(khxpi, 0.2e1);
t5 = pow(khxpr, 0.2e1);
t7 = 0.1e1 / (t4 + t5);
kappar = (t1 + t2) * t7;
t8 = khxi * khxpr;
t9 = khxr * khxpi;
kappai = (t8 - t9) * t7;
t11 = 0.1e1 / mu1;
t12 = mu2 * t11;
t13 = kxi * kxpi;
t15 = kxpr * kxr;
t17 = kxi * kxpr;
t19 = kxpi * kxr;
t21 = khxi * khypi;
t22 = kxi * kypi;
t24 = kxr * kypr;
t26 = khxi * khypr;
t27 = kxi * kypr;
t29 = kxr * kypi;
t33 = khxpi * khyi;
t34 = kxpi * kyi;
t36 = kxpr * kyr;
t38 = khxpi * khyr;
t39 = kxpi * kyr;
t41 = kxpr * kyi;
t45 = t1 * t13 + t1 * t15 + t2 * t13 + t2 * t15 - t8 * t17 + t9 * t17 + t8 * t19 - t9 * t19 + t21 * t22 + t21 * t24 - t26 * t27 + t26 * t29 + t33 * t34 + t33 * t36 - t38 * t39 + t38 * t41;
t46 = khxpr * khyi;
t49 = khxpr * khyr;
t52 = khxr * khypi;
t55 = khxr * khypr;
t58 = khyi * khypi;
t59 = kyi * kypi;
t61 = kypr * kyr;
t63 = khyi * khypr;
t64 = kyi * kypr;
t66 = kypi * kyr;
t68 = khypi * khyr;
t71 = khypr * khyr;
t74 = t55 * t22 + t55 * t24 + t52 * t27 - t52 * t29 + t49 * t34 + t49 * t36 + t46 * t39 - t46 * t41 + t58 * t59 + t58 * t61 + t71 * t59 + t71 * t61 - t63 * t64 + t63 * t66 + t68 * t64 - t68 * t66;
t76 = pow(kxpi, 0.2e1);
t78 = pow(kxpr, 0.2e1);
t80 = khxpi * khypi;
t81 = kxpi * kypi;
t84 = kxpr * kypr;
t87 = khxpi * khypr;
t88 = kxpi * kypr;
t91 = kxpr * kypi;
t96 = khxpr * khypi;
t101 = khxpr * khypr;
t106 = pow(khypi, 0.2e1);
t107 = pow(kypi, 0.2e1);
t109 = pow(kypr, 0.2e1);
t111 = pow(khypr, 0.2e1);
t114 = 0.2e1 * t101 * t81 + 0.2e1 * t101 * t84 + t106 * t107 + t106 * t109 + t111 * t107 + t111 * t109 + t4 * t76 + t4 * t78 + t5 * t76 + t5 * t78 + 0.2e1 * t80 * t81 + 0.2e1 * t80 * t84 - 0.2e1 * t87 * t88 + 0.2e1 * t87 * t91 + 0.2e1 * t96 * t88 - 0.2e1 * t96 * t91;
t115 = 0.1e1 / t114;
betar = t12 * (t45 + t74) * t115;
t133 = t1 * t17 - t1 * t19 + t8 * t13 - t9 * t13 + t8 * t15 - t9 * t15 + t2 * t17 - t2 * t19 + t21 * t27 - t21 * t29 + t26 * t22 + t26 * t24 - t33 * t39 + t33 * t41 - t38 * t34 - t38 * t36;
t150 = -t52 * t22 - t52 * t24 + t55 * t27 - t55 * t29 + t46 * t34 + t46 * t36 - t49 * t39 + t49 * t41 + t58 * t64 - t58 * t66 + t63 * t59 - t68 * t59 + t63 * t61 - t68 * t61 + t71 * t64 - t71 * t66;
betai = t12 * (t133 + t150) * t115;
t153 = pow(betai, 0.2e1);
t154 = pow(betar, 0.2e1);
t155 = pow(kappai, 0.2e1);
t156 = pow(kappar, 0.2e1);
t163 = 0.1e1 / (0.2e1 * kappai * betai + 0.2e1 * kappar * betar + t153 + t154 + t155 + t156);
rr = (t153 + t154 - t155 - t156) * t163;
ri = 0.2e1 * (betai * kappar - betar * kappai) * t163;
taur = 0.2e1 * (betar * t155 + betar * t156 + t153 * kappar + t154 * kappar) * t163;
taui = 0.2e1 * (betai * t155 + betai * t156 + kappai * t153 + t154 * kappai) * t163;
t180 = x * kxi;
t181 = y * kyi;
t182 = t * sr;
t184 = exp(t180 - t181 + t182);
t185 = x * kxr;
t186 = y * kyr;
t187 = t * si;
t188 = t185 - t186 - t187;
t189 = sin(t188);
t190 = t184 * t189;
t191 = khyi * rr;
t193 = khyr * ri;
t195 = cos(t188);
t196 = t184 * t195;
t197 = khyi * ri;
t199 = khyr * rr;
t202 = exp(-t180 - t181 + t182);
t203 = t185 + t186 + t187;
t204 = cos(t203);
t205 = t202 * t204;
t207 = sin(t203);
t208 = t202 * t207;
Exr = -a * (-t208 * khyi + t205 * khyr - t190 * t191 - t190 * t193 + t196 * t197 - t196 * t199);
t212 = khxi * rr;
t214 = khxr * ri;
t216 = khxi * ri;
t218 = khxr * rr;
Eyr = a * (-t208 * khxi + t205 * khxr + t190 * t212 + t190 * t214 - t196 * t216 + t196 * t218);
t226 = exp(-x * kxpi - y * kypi + t182);
t229 = x * kxpr + y * kypr + t187;
t230 = cos(t229);
t231 = t226 * t230;
t232 = khypi * taui;
t234 = khypr * taur;
t236 = sin(t229);
t237 = t226 * t236;
t238 = khypi * taur;
t240 = khypr * taui;
Expr = -a * (-t231 * t232 + t231 * t234 - t237 * t238 - t237 * t240);
t244 = khxpi * taui;
t246 = khxpr * taur;
t248 = khxpi * taur;
t250 = khxpr * taui;
Eypr = a * (-t231 * t244 + t231 * t246 - t237 * t248 - t237 * t250);
Exi = -a * (t205 * khyi + t208 * khyr - t190 * t197 + t190 * t199 - t196 * t191 - t196 * t193);
Eyi = a * (t205 * khxi + t208 * khxr + t190 * t216 - t190 * t218 + t196 * t212 + t196 * t214);
Expi = -a * (t231 * t238 + t231 * t240 - t237 * t232 + t237 * t234);
Eypi = a * (t231 * t248 + t231 * t250 - t237 * t244 + t237 * t246);
t279 = t11 * a;
t280 = t190 * khxi;
t281 = kxi * ri;
t282 = t281 * si;
t284 = kxi * rr;
t285 = t284 * sr;
t287 = kxr * ri;
t288 = t287 * sr;
t290 = kxr * rr;
t291 = t290 * si;
t293 = t190 * khxr;
t294 = t281 * sr;
t296 = t284 * si;
t298 = t287 * si;
t300 = t290 * sr;
t302 = t190 * khyi;
t303 = kyi * ri;
t304 = t303 * si;
t306 = kyi * rr;
t307 = t306 * sr;
t309 = kyr * ri;
t310 = t309 * sr;
t312 = kyr * rr;
t313 = t312 * si;
t315 = -t280 * t282 - t280 * t285 - t280 * t288 + t280 * t291 - t293 * t294 + t293 * t296 + t293 * t298 + t293 * t300 - t302 * t304 - t302 * t307 - t302 * t310 + t302 * t313;
t316 = t190 * khyr;
t317 = t303 * sr;
t319 = t306 * si;
t321 = t309 * si;
t323 = t312 * sr;
t325 = t196 * khxi;
t330 = t196 * khxr;
t335 = -t330 * t282 - t330 * t285 - t330 * t288 + t330 * t291 + t325 * t294 - t325 * t296 - t325 * t298 - t325 * t300 - t316 * t317 + t316 * t319 + t316 * t321 + t316 * t323;
t337 = t196 * khyi;
t342 = t196 * khyr;
t347 = khxi * kxi;
t348 = t347 * si;
t350 = khxi * kxr;
t351 = t350 * sr;
t353 = khxr * kxi;
t354 = t353 * sr;
t356 = khxr * kxr;
t357 = t356 * si;
t359 = t205 * t348 + t205 * t351 + t205 * t354 - t205 * t357 - t342 * t304 - t342 * t307 - t342 * t310 + t342 * t313 + t337 * t317 - t337 * t319 - t337 * t321 - t337 * t323;
t360 = khyi * kyi;
t361 = t360 * si;
t363 = khyi * kyr;
t364 = t363 * sr;
t366 = khyr * kyi;
t367 = t366 * sr;
t369 = khyr * kyr;
t370 = t369 * si;
t372 = t347 * sr;
t374 = t350 * si;
t376 = t353 * si;
t378 = t356 * sr;
t380 = t360 * sr;
t382 = t363 * si;
t384 = t366 * si;
t386 = t369 * sr;
t388 = t205 * t361 + t205 * t364 + t205 * t367 - t205 * t370 - t208 * t372 + t208 * t374 + t208 * t376 + t208 * t378 - t208 * t380 + t208 * t382 + t208 * t384 + t208 * t386;
t391 = pow(si, 0.2e1);
t392 = pow(sr, 0.2e1);
t394 = 0.1e1 / (t391 + t392);
Hzr = t279 * (t315 + t335 + t359 + t388) * t394;
t397 = 0.1e1 / mu2 * a;
t398 = t236 * khypr;
t399 = kypi * si;
t400 = t399 * taur;
t402 = kypi * sr;
t403 = t402 * taui;
t405 = kypr * si;
t406 = t405 * taui;
t408 = kypr * sr;
t409 = t408 * taur;
t411 = t230 * khxpi;
t412 = kxpi * si;
t413 = t412 * taur;
t415 = kxpi * sr;
t416 = t415 * taui;
t418 = kxpr * si;
t419 = t418 * taui;
t421 = kxpr * sr;
t422 = t421 * taur;
t424 = t230 * khxpr;
t425 = t412 * taui;
t427 = t415 * taur;
t429 = t418 * taur;
t431 = t421 * taui;
t433 = t230 * khypi;
t438 = t398 * t400 - t398 * t403 + t398 * t406 + t398 * t409 + t433 * t400 - t433 * t403 + t433 * t406 + t433 * t409 + t411 * t413 - t411 * t416 + t411 * t419 + t411 * t422 + t424 * t425 + t424 * t427 - t424 * t429 + t424 * t431;
t439 = t230 * khypr;
t440 = t399 * taui;
t442 = t402 * taur;
t444 = t405 * taur;
t446 = t408 * taui;
t448 = t236 * khxpi;
t453 = t236 * khxpr;
t458 = t236 * khypi;
t463 = t453 * t413 - t453 * t416 + t453 * t419 + t453 * t422 - t448 * t425 - t448 * t427 + t448 * t429 - t448 * t431 + t439 * t440 + t439 * t442 - t439 * t444 + t439 * t446 - t458 * t440 - t458 * t442 + t458 * t444 - t458 * t446;
Hzpr = t397 * t226 * (t438 + t463) * t394;
t479 = -t280 * t294 + t280 * t296 + t280 * t298 + t280 * t300 + t293 * t282 + t293 * t285 + t293 * t288 - t293 * t291 - t302 * t317 + t302 * t319 + t302 * t321 + t302 * t323;
t492 = -t325 * t282 - t325 * t285 - t325 * t288 + t325 * t291 - t330 * t294 + t330 * t296 + t330 * t298 + t330 * t300 + t316 * t304 + t316 * t307 + t316 * t310 - t316 * t313;
t506 = t205 * t372 - t205 * t374 - t205 * t376 - t205 * t378 - t337 * t304 - t337 * t307 - t337 * t310 + t337 * t313 - t342 * t317 + t342 * t319 + t342 * t321 + t342 * t323;
t519 = t205 * t380 - t205 * t382 - t205 * t384 - t205 * t386 + t208 * t348 + t208 * t351 + t208 * t354 - t208 * t357 + t208 * t361 + t208 * t364 + t208 * t367 - t208 * t370;
Hzi = t279 * (t479 + t492 + t506 + t519) * t394;
t539 = t398 * t440 + t398 * t442 - t398 * t444 + t398 * t446 + t458 * t400 - t458 * t403 + t458 * t406 + t458 * t409 + t411 * t425 + t411 * t427 + t448 * t419 + t448 * t422 + t453 * t425 + t453 * t427 - t453 * t429 + t453 * t431;
t556 = -t439 * t400 + t439 * t403 - t439 * t406 - t439 * t409 - t411 * t429 + t411 * t431 - t424 * t413 + t448 * t413 + t424 * t416 - t448 * t416 - t424 * t419 - t424 * t422 + t433 * t440 + t433 * t442 - t433 * t444 + t433 * t446;
Hzpi = t397 * t226 * (t539 + t556) * t394;


        printF("chiSum1=(%9.3e,%9.3e) chiSum2=(%9.3e,%9.3e) r=(%8.2e,%8.2e) "
         	   "tau=(%8.2e,%8.2e) khy=(%8.2e,%8.2e) khpy=(%8.2e,%8.2e): ",
         	   chiSum1r,chiSum1i,chiSum2r,chiSum2i,rr,ri,taur,taui,khyr,khyi,khypr,khypi);


        checkPlaneMaterialInterfaceJumps( 
            c1,c2,eps1,eps2,mu1,mu2, sr,si, rr,ri, taur,taui, 
            eps1Hatr,eps1Hati, eps2Hatr,eps2Hati,
            chiSum1r,chiSum1i,chiSum2r,chiSum2i,
            kxr,kxi, kyr,kyi, kxpr,kxpi, kypr,kypi );
                        
    // OV_ABORT("stop here for now");
                        
    }

    
    if( cg.numberOfDimensions()==3 )
    {
// -- 3D
//    s=sr + I*si : complex frequency
//    [ax,ay,az] = amplitude vector of incident wave, ax=axr+I*axi etc.
//    [arx,ary,arz] = amplitude vector of reflected wave, arx=arxr+I*arxi etc.
//    [atx,aty,atz] = amplitude vector of transmitted wave, atx=atxr+I*atxi etc.
//    [kx,ky,kz] = incident wave vector, kx=kxr+I*kxi etc. 
//    [krx,kry,krz] = reflected wave vector, krx=krxr+I*krxi etc. 
//    [ktx,kty,ktz] = reflected wave vector, ktx=ktxr+I*ktxi etc. 


    // real kv[3]={twoPi*kx,twoPi*ky,twoPi*kz}; // 
                  
        real (&arr)[3] = dbase.get<real[3]>("arr");
        real (&ari)[3] = dbase.get<real[3]>("ari");
        real (&atr)[3] = dbase.get<real[3]>("atr");
        real (&ati)[3] = dbase.get<real[3]>("ati");

        real (&krr)[3] = dbase.get<real[3]>("krr");
        real (&kri)[3] = dbase.get<real[3]>("kri");
        real (&ktr)[3] = dbase.get<real[3]>("ktr");
        real (&kti)[3] = dbase.get<real[3]>("kti");

        printF("chiSum1=(%9.3e,%9.3e) chiSum2=(%9.3e,%9.3e)\n",
         	   chiSum1r,chiSum1i,chiSum2r,chiSum2i);

        evalMaterialInterfaceSolution3d( sr,si, chiSum1r, chiSum1i, chiSum2r, chiSum2i,
                             				     av[0],nv[0],kv[0], eps1,eps2,mu1,mu2, 
                             				     krr[0],kri[0], ktr[0],kti[0], arr[0],ari[0], atr[0],ati[0] );
        
          if( true )
          {
              printF(" INIT:\n");
              printF(" kr =(%7.3f + %7.3f I,%7.3f + %7.3f I,%7.3f + %7.3f I) \n",krr[0],kri[0],krr[1],kri[1],krr[2],kri[2]);
              printF(" kt =(%7.3f + %7.3f I,%7.3f + %7.3f I,%7.3f + %7.3f I) \n",ktr[0],kti[0],ktr[1],kti[1],ktr[2],kti[2]);

              printF(" ar =(%7.3f + %7.3f I,%7.3f + %7.3f I,%7.3f + %7.3f I) \n",arr[0],ari[0],arr[1],ari[1],arr[2],ari[2]);
              printF(" at =(%7.3f + %7.3f I,%7.3f + %7.3f I,%7.3f + %7.3f I) \n",atr[0],ati[0],atr[1],ati[1],atr[2],ati[2]);

          }

    }
    
    printF("\n ------------------- PlaneInterfaceExactSolution:: END initialize ------------------------\n\n");

    check();
    

    return 0;
}






// ===================================================================================
/// \brief Evaluate the frequency space solution (complex valued) at a single point.
/// \param x[3] (input): point to evaluate the solution
/// 
/// \param E[6] (output):  real and imaginary parts of E : [Exr,Eyr,Ezr,Exi,Eyi,Ezi]
/// \param H[6] (output):  optionally compute real and imaginary parts of H : [Hxr,Hyr,Hzr,Hxi,Hyi,Hzi]
///
/// \note: This may not be so efficient.
// ===================================================================================
int PlaneInterfaceExactSolution::
eval( real xv[3], real *Ev, real *Hv /*= NULL */  )
{
    const int & numberOfDimensions = dbase.get<int>("numberOfDimensions");

    bool computeMagneticField = Hv != NULL;
  // printF("PIES:eval: computeMagneticField=%i, numberOfDimensions=%i\n",(int)computeMagneticField,numberOfDimensions);
    
  // const int numVars = computeMagneticField ? 12 : 6;

    const real & a = dbase.get<real>("amp");

    real (&normalPlaneMaterialInterface)[3] = dbase.get<real[3]>("normalPlaneMaterialInterface");
    real (&x0PlaneMaterialInterface)[3] = dbase.get<real[3]>("x0PlaneMaterialInterface");

    real *nv = normalPlaneMaterialInterface;
    real *x0 = x0PlaneMaterialInterface;

  // Here we assume that the normal points from domain 1 into domain 2
    real nDotX = nv[0]*(xv[0]-x0[0]) + nv[1]*(xv[1]-x0[1]) + nv[2]*(xv[2]-x0[2]);
    
    const int myDomain = nDotX > 0. ? 1 : 0;
    
    LocalReal (&s)[2] = dbase.get<LocalReal[2]>("s");

    real & sr = s[0], &si = s[1];
    real (&av)[3] = dbase.get<real[3]>("av");

    real (&kv)[3] = dbase.get<real[3]>("kv");

    real & eps1 = dbase.get<real>("eps1");
    real & eps2 = dbase.get<real>("eps2");
    real & mu1 = dbase.get<real>("mu1");
    real & mu2 = dbase.get<real>("mu2");

  // -- these are needed to compute P : 
  // int & numberOfPolarizationVectors1 = dbase.get<int>("numberOfPolarizationVectors1");
  // int & numberOfPolarizationVectors2 = dbase.get<int>("numberOfPolarizationVectors2");

  // real (&chi1r)[10] = dbase.get<real[10]>("chi1r");
  // real (&chi1i)[10] = dbase.get<real[10]>("chi1i");

  // real (&chi2r)[10] = dbase.get<real[10]>("chi2r");
  // real (&chi2i)[10] = dbase.get<real[10]>("chi2i");



    real t=0.;
  // const real & a = dbase.get<real>("amp");
    
    real x=xv[0]-x0[0], y=xv[1]-x0[1], z=xv[2]-x0[2]; // *** NOTE: shift point 

    if( numberOfDimensions==2 )
    {

    // -- for 2d eval: 
        real (&kvr)[3] = dbase.get<real[3]>("kvr");
        real (&kvi)[3] = dbase.get<real[3]>("kvi");

        real (&kpvr)[3] = dbase.get<real[3]>("kpvr");
        real (&kpvi)[3] = dbase.get<real[3]>("kpvi");

        real &kxr =kvr[0], &kxi=kvi[0], &kyr=kvr[1], &kyi=kvi[1];   // Incident wave number (complex)
        real &kxpr =kpvr[0], &kxpi=kpvi[0], &kypr=kpvr[1], &kypi=kvi[1];   // 


    // -- un-rotate the point (x,y) to the reference space and evaluate the solution 
        real ct= nv[0], st=nv[1];  // cos(theta), sin(theta)
    // real xa=x-x0[0], ya=y-x0[1];
        real xa=x, ya=y;
        x =  ct*xa + st*ya;
        y = -st*xa + ct*ya;
                    

    // Here are the statements to eval the solution: 
// File generated by overtureFramework/cg/mx/codes/dispersivePlaneWaveInterface.maple
// File generated by DropBox/DMX/codes/dispersivePlaneWaveInterface.maple
real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t530,t531,t532,t533,t534,t535,t536,t537,t538,t539,t540,t541,t542,t543,t544,t545,t546,t547,t548,t549,t550,t551,t552,t553,t554,t555,t556,t557,t558,t559,t560,t561,t562,t563,t564,t565,t566,t567,t568,t569,t570,t571,t572,t573,t574,t575,t576,t577,t578,t579,t580,t581,t582,t583,t584,t585,t586,t587,t588,t589,t590,t591,t592,t593,t594,t595,t596,t597,t598,t599,t600;
// -------------------------------------------------------------------------
// Need: 
//    a = amplitude of the wave, e.g. a=1
//    x,y,t
//    sr,si : s= sr + I*si 
//    alphaP
//    kxr,kxi,  kyr,kyi,     : complex wave number on left
//    kxpr,kxpi,  kypr,kypi, : complex wave number on right (plus)
// -------------------------------------------------------------------------
// Evaluated:                                                              
//    Exr,Eyr   : left state                                                
//    Expr,Eypr : right state                                               
// -------------------------------------------------------------------------
real kNorm = sqrt( kxr*kxr + kxi*kxi + kyr*kyr + kyi*kyi);
real khxr = kxr/kNorm, khxi = kxi/kNorm, khyr=kyr/kNorm, khyi=kyi/kNorm; 
real kpNorm = sqrt( kxpr*kxpr + kxpi*kxpi + kypr*kypr + kypi*kypi);
real khxpr = kxpr/kpNorm, khxpi=kxpi/kpNorm, khypr=kypr/kpNorm, khypi=kypi/kpNorm; 
real kappar,kappai, betar, betai, rr,ri, taur,taui;
real Exr,Eyr,Hzr, Expr,Eypr,Hzpr;
real Exi,Eyi,Hzi, Expi,Eypi,Hzpi;
t1 = khxi * khxpi;
t2 = khxr * khxpr;
t4 = pow(khxpi, 0.2e1);
t5 = pow(khxpr, 0.2e1);
t7 = 0.1e1 / (t4 + t5);
kappar = (t1 + t2) * t7;
t8 = khxi * khxpr;
t9 = khxr * khxpi;
kappai = (t8 - t9) * t7;
t11 = 0.1e1 / mu1;
t12 = mu2 * t11;
t13 = kxi * kxpi;
t15 = kxpr * kxr;
t17 = kxi * kxpr;
t19 = kxpi * kxr;
t21 = khxi * khypi;
t22 = kxi * kypi;
t24 = kxr * kypr;
t26 = khxi * khypr;
t27 = kxi * kypr;
t29 = kxr * kypi;
t33 = khxpi * khyi;
t34 = kxpi * kyi;
t36 = kxpr * kyr;
t38 = khxpi * khyr;
t39 = kxpi * kyr;
t41 = kxpr * kyi;
t45 = t1 * t13 + t1 * t15 + t2 * t13 + t2 * t15 - t8 * t17 + t9 * t17 + t8 * t19 - t9 * t19 + t21 * t22 + t21 * t24 - t26 * t27 + t26 * t29 + t33 * t34 + t33 * t36 - t38 * t39 + t38 * t41;
t46 = khxpr * khyi;
t49 = khxpr * khyr;
t52 = khxr * khypi;
t55 = khxr * khypr;
t58 = khyi * khypi;
t59 = kyi * kypi;
t61 = kypr * kyr;
t63 = khyi * khypr;
t64 = kyi * kypr;
t66 = kypi * kyr;
t68 = khypi * khyr;
t71 = khypr * khyr;
t74 = t55 * t22 + t55 * t24 + t52 * t27 - t52 * t29 + t49 * t34 + t49 * t36 + t46 * t39 - t46 * t41 + t58 * t59 + t58 * t61 + t71 * t59 + t71 * t61 - t63 * t64 + t63 * t66 + t68 * t64 - t68 * t66;
t76 = pow(kxpi, 0.2e1);
t78 = pow(kxpr, 0.2e1);
t80 = khxpi * khypi;
t81 = kxpi * kypi;
t84 = kxpr * kypr;
t87 = khxpi * khypr;
t88 = kxpi * kypr;
t91 = kxpr * kypi;
t96 = khxpr * khypi;
t101 = khxpr * khypr;
t106 = pow(khypi, 0.2e1);
t107 = pow(kypi, 0.2e1);
t109 = pow(kypr, 0.2e1);
t111 = pow(khypr, 0.2e1);
t114 = 0.2e1 * t101 * t81 + 0.2e1 * t101 * t84 + t106 * t107 + t106 * t109 + t111 * t107 + t111 * t109 + t4 * t76 + t4 * t78 + t5 * t76 + t5 * t78 + 0.2e1 * t80 * t81 + 0.2e1 * t80 * t84 - 0.2e1 * t87 * t88 + 0.2e1 * t87 * t91 + 0.2e1 * t96 * t88 - 0.2e1 * t96 * t91;
t115 = 0.1e1 / t114;
betar = t12 * (t45 + t74) * t115;
t133 = t1 * t17 - t1 * t19 + t8 * t13 - t9 * t13 + t8 * t15 - t9 * t15 + t2 * t17 - t2 * t19 + t21 * t27 - t21 * t29 + t26 * t22 + t26 * t24 - t33 * t39 + t33 * t41 - t38 * t34 - t38 * t36;
t150 = -t52 * t22 - t52 * t24 + t55 * t27 - t55 * t29 + t46 * t34 + t46 * t36 - t49 * t39 + t49 * t41 + t58 * t64 - t58 * t66 + t63 * t59 - t68 * t59 + t63 * t61 - t68 * t61 + t71 * t64 - t71 * t66;
betai = t12 * (t133 + t150) * t115;
t153 = pow(betai, 0.2e1);
t154 = pow(betar, 0.2e1);
t155 = pow(kappai, 0.2e1);
t156 = pow(kappar, 0.2e1);
t163 = 0.1e1 / (0.2e1 * kappai * betai + 0.2e1 * kappar * betar + t153 + t154 + t155 + t156);
rr = (t153 + t154 - t155 - t156) * t163;
ri = 0.2e1 * (betai * kappar - betar * kappai) * t163;
taur = 0.2e1 * (betar * t155 + betar * t156 + t153 * kappar + t154 * kappar) * t163;
taui = 0.2e1 * (betai * t155 + betai * t156 + kappai * t153 + t154 * kappai) * t163;
t180 = x * kxi;
t181 = y * kyi;
t182 = t * sr;
t184 = exp(t180 - t181 + t182);
t185 = x * kxr;
t186 = y * kyr;
t187 = t * si;
t188 = t185 - t186 - t187;
t189 = sin(t188);
t190 = t184 * t189;
t191 = khyi * rr;
t193 = khyr * ri;
t195 = cos(t188);
t196 = t184 * t195;
t197 = khyi * ri;
t199 = khyr * rr;
t202 = exp(-t180 - t181 + t182);
t203 = t185 + t186 + t187;
t204 = cos(t203);
t205 = t202 * t204;
t207 = sin(t203);
t208 = t202 * t207;
Exr = -a * (-t208 * khyi + t205 * khyr - t190 * t191 - t190 * t193 + t196 * t197 - t196 * t199);
t212 = khxi * rr;
t214 = khxr * ri;
t216 = khxi * ri;
t218 = khxr * rr;
Eyr = a * (-t208 * khxi + t205 * khxr + t190 * t212 + t190 * t214 - t196 * t216 + t196 * t218);
t226 = exp(-x * kxpi - y * kypi + t182);
t229 = x * kxpr + y * kypr + t187;
t230 = cos(t229);
t231 = t226 * t230;
t232 = khypi * taui;
t234 = khypr * taur;
t236 = sin(t229);
t237 = t226 * t236;
t238 = khypi * taur;
t240 = khypr * taui;
Expr = -a * (-t231 * t232 + t231 * t234 - t237 * t238 - t237 * t240);
t244 = khxpi * taui;
t246 = khxpr * taur;
t248 = khxpi * taur;
t250 = khxpr * taui;
Eypr = a * (-t231 * t244 + t231 * t246 - t237 * t248 - t237 * t250);
Exi = -a * (t205 * khyi + t208 * khyr - t190 * t197 + t190 * t199 - t196 * t191 - t196 * t193);
Eyi = a * (t205 * khxi + t208 * khxr + t190 * t216 - t190 * t218 + t196 * t212 + t196 * t214);
Expi = -a * (t231 * t238 + t231 * t240 - t237 * t232 + t237 * t234);
Eypi = a * (t231 * t248 + t231 * t250 - t237 * t244 + t237 * t246);
t279 = t11 * a;
t280 = t190 * khxi;
t281 = kxi * ri;
t282 = t281 * si;
t284 = kxi * rr;
t285 = t284 * sr;
t287 = kxr * ri;
t288 = t287 * sr;
t290 = kxr * rr;
t291 = t290 * si;
t293 = t190 * khxr;
t294 = t281 * sr;
t296 = t284 * si;
t298 = t287 * si;
t300 = t290 * sr;
t302 = t190 * khyi;
t303 = kyi * ri;
t304 = t303 * si;
t306 = kyi * rr;
t307 = t306 * sr;
t309 = kyr * ri;
t310 = t309 * sr;
t312 = kyr * rr;
t313 = t312 * si;
t315 = -t280 * t282 - t280 * t285 - t280 * t288 + t280 * t291 - t293 * t294 + t293 * t296 + t293 * t298 + t293 * t300 - t302 * t304 - t302 * t307 - t302 * t310 + t302 * t313;
t316 = t190 * khyr;
t317 = t303 * sr;
t319 = t306 * si;
t321 = t309 * si;
t323 = t312 * sr;
t325 = t196 * khxi;
t330 = t196 * khxr;
t335 = -t330 * t282 - t330 * t285 - t330 * t288 + t330 * t291 + t325 * t294 - t325 * t296 - t325 * t298 - t325 * t300 - t316 * t317 + t316 * t319 + t316 * t321 + t316 * t323;
t337 = t196 * khyi;
t342 = t196 * khyr;
t347 = khxi * kxi;
t348 = t347 * si;
t350 = khxi * kxr;
t351 = t350 * sr;
t353 = khxr * kxi;
t354 = t353 * sr;
t356 = khxr * kxr;
t357 = t356 * si;
t359 = t205 * t348 + t205 * t351 + t205 * t354 - t205 * t357 - t342 * t304 - t342 * t307 - t342 * t310 + t342 * t313 + t337 * t317 - t337 * t319 - t337 * t321 - t337 * t323;
t360 = khyi * kyi;
t361 = t360 * si;
t363 = khyi * kyr;
t364 = t363 * sr;
t366 = khyr * kyi;
t367 = t366 * sr;
t369 = khyr * kyr;
t370 = t369 * si;
t372 = t347 * sr;
t374 = t350 * si;
t376 = t353 * si;
t378 = t356 * sr;
t380 = t360 * sr;
t382 = t363 * si;
t384 = t366 * si;
t386 = t369 * sr;
t388 = t205 * t361 + t205 * t364 + t205 * t367 - t205 * t370 - t208 * t372 + t208 * t374 + t208 * t376 + t208 * t378 - t208 * t380 + t208 * t382 + t208 * t384 + t208 * t386;
t391 = pow(si, 0.2e1);
t392 = pow(sr, 0.2e1);
t394 = 0.1e1 / (t391 + t392);
Hzr = t279 * (t315 + t335 + t359 + t388) * t394;
t397 = 0.1e1 / mu2 * a;
t398 = t236 * khypr;
t399 = kypi * si;
t400 = t399 * taur;
t402 = kypi * sr;
t403 = t402 * taui;
t405 = kypr * si;
t406 = t405 * taui;
t408 = kypr * sr;
t409 = t408 * taur;
t411 = t230 * khxpi;
t412 = kxpi * si;
t413 = t412 * taur;
t415 = kxpi * sr;
t416 = t415 * taui;
t418 = kxpr * si;
t419 = t418 * taui;
t421 = kxpr * sr;
t422 = t421 * taur;
t424 = t230 * khxpr;
t425 = t412 * taui;
t427 = t415 * taur;
t429 = t418 * taur;
t431 = t421 * taui;
t433 = t230 * khypi;
t438 = t398 * t400 - t398 * t403 + t398 * t406 + t398 * t409 + t433 * t400 - t433 * t403 + t433 * t406 + t433 * t409 + t411 * t413 - t411 * t416 + t411 * t419 + t411 * t422 + t424 * t425 + t424 * t427 - t424 * t429 + t424 * t431;
t439 = t230 * khypr;
t440 = t399 * taui;
t442 = t402 * taur;
t444 = t405 * taur;
t446 = t408 * taui;
t448 = t236 * khxpi;
t453 = t236 * khxpr;
t458 = t236 * khypi;
t463 = t453 * t413 - t453 * t416 + t453 * t419 + t453 * t422 - t448 * t425 - t448 * t427 + t448 * t429 - t448 * t431 + t439 * t440 + t439 * t442 - t439 * t444 + t439 * t446 - t458 * t440 - t458 * t442 + t458 * t444 - t458 * t446;
Hzpr = t397 * t226 * (t438 + t463) * t394;
t479 = -t280 * t294 + t280 * t296 + t280 * t298 + t280 * t300 + t293 * t282 + t293 * t285 + t293 * t288 - t293 * t291 - t302 * t317 + t302 * t319 + t302 * t321 + t302 * t323;
t492 = -t325 * t282 - t325 * t285 - t325 * t288 + t325 * t291 - t330 * t294 + t330 * t296 + t330 * t298 + t330 * t300 + t316 * t304 + t316 * t307 + t316 * t310 - t316 * t313;
t506 = t205 * t372 - t205 * t374 - t205 * t376 - t205 * t378 - t337 * t304 - t337 * t307 - t337 * t310 + t337 * t313 - t342 * t317 + t342 * t319 + t342 * t321 + t342 * t323;
t519 = t205 * t380 - t205 * t382 - t205 * t384 - t205 * t386 + t208 * t348 + t208 * t351 + t208 * t354 - t208 * t357 + t208 * t361 + t208 * t364 + t208 * t367 - t208 * t370;
Hzi = t279 * (t479 + t492 + t506 + t519) * t394;
t539 = t398 * t440 + t398 * t442 - t398 * t444 + t398 * t446 + t458 * t400 - t458 * t403 + t458 * t406 + t458 * t409 + t411 * t425 + t411 * t427 + t448 * t419 + t448 * t422 + t453 * t425 + t453 * t427 - t453 * t429 + t453 * t431;
t556 = -t439 * t400 + t439 * t403 - t439 * t406 - t439 * t409 - t411 * t429 + t411 * t431 - t424 * t413 + t448 * t413 + t424 * t416 - t448 * t416 - t424 * t419 - t424 * t422 + t433 * t440 + t433 * t442 - t433 * t444 + t433 * t446;
Hzpi = t397 * t226 * (t539 + t556) * t394;


    // ---- rotate the field from the reference space to the rotated space ---
        real Exra, Eyra;
        real Exia, Eyia;
        if( myDomain==0 )
        {
            Exra=Exr, Eyra=Eyr;
            Exia=Exi, Eyia=Eyi;
        }
        else
        {
            Exra=Expr, Eyra=Eypr;
            Exia=Expi, Eyia=Eypi;
        }
        Exr =  ct*Exra - st*Eyra;
        Eyr =  st*Exra + ct*Eyra;
                    
        Exi =  ct*Exia - st*Eyia;
        Eyi =  st*Exia + ct*Eyia;
    
        Ev[0]=Exr; Ev[1]=Eyr; Ev[2]=0.;
        Ev[3]=Exi; Ev[4]=Eyi; Ev[5]=0.;
        
        if( computeMagneticField )
        {
            for( int i=0; i<6; i++ ){ Hv[i]=0.; }  
            if( myDomain==0 )
            {
      	Hv[2]=Hzr;
      	Hv[5]=Hzi;
            }
            else
            {
      	Hv[2]=Hzpr;
      	Hv[5]=Hzpi;
            }
            
        }
        
    }
    else
    {
    // -- for 3d eval: 

        real (&arr)[3] = dbase.get<real[3]>("arr");
        real (&ari)[3] = dbase.get<real[3]>("ari");
        real (&atr)[3] = dbase.get<real[3]>("atr");
        real (&ati)[3] = dbase.get<real[3]>("ati");

        real (&krr)[3] = dbase.get<real[3]>("krr");
        real (&kri)[3] = dbase.get<real[3]>("kri");
        real (&ktr)[3] = dbase.get<real[3]>("ktr");
        real (&kti)[3] = dbase.get<real[3]>("kti");
        const real &ax=av[0], &ay=av[1], &az=av[2]; 

    // These need to be set for the solution evaluation below:
        real axr=ax, axi=0., ayr=ay, ayi=0., azr=az, azi=0.;
        real arxr=arr[0], arxi=ari[0], aryr=arr[1], aryi=ari[1], arzr=arr[2], arzi=ari[2];
        real atxr=atr[0], atxi=ati[0], atyr=atr[1], atyi=ati[1], atzr=atr[2], atzi=ati[2];
                    
    // real kxr=twoPi*kx, kxi=0., kyr=twoPi*ky, kyi=0., kzr=twoPi*kz, kzi=0.;
        real kxr=kv[0], kxi=0., kyr=kv[1], kyi=0., kzr=kv[2], kzi=0.;
        real krxr=krr[0], krxi=kri[0], kryr=krr[1], kryi=kri[1], krzr=krr[2], krzi=kri[2];
        real ktxr=ktr[0], ktxi=kti[0], ktyr=ktr[1], ktyi=kti[1], ktzr=ktr[2], ktzi=kti[2];


        

    // Here are the statements to eval the solution: 
        if( true )
        {
      // *** NOTE: Frequency domain: 
// File generated by Dropbox/DARPA/RPI/adePapers/adegdmi/dispersivePlaneWaveInterface3d.maple
real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150;
// -------------------------------------------------------------------------
// Need: 
//    s=sr + I*si : complex frequency
//    [ax,ay,az] = amplitude vector of incident wave, ax=axr+I*axi etc.
//    [arx,ary,arz] = amplitude vector of reflected wave, arx=arxr+I*arxi etc.
//    [atx,aty,atz] = amplitude vector of transmitted wave, atx=atxr+I*atxi etc.
//    [kx,ky,kz] = incident wave vector, kx=kxr+I*kxi etc. 
//    [krx,kry,krz] = reflected wave vector, krx=krxr+I*krxi etc. 
//    [ktx,kty,ktz] = reflected wave vector, ktx=ktxr+I*ktxi etc. 
// -------------------------------------------------------------------------
// Evaluated:                                                               
//    [Exr,Exi] ,[Eyr,Eyi], [Ezr,Ezi]        : left state                   
//    [Etxr,Etxi] ,[Etyr,Etyi], [Etzr,Etzi]  : right state                  
// -------------------------------------------------------------------------
real Exr,Eyr,Ezr, Exi,Eyi,Ezi;
real Etxr,Etyr,Etzr, Etxi,Etyi,Etzi;
t5 = exp(-kxi * x - y * kyi - z * kzi);
t6 = axr * t5;
t10 = x * kxr + y * kyr + z * kzr;
t11 = cos(t10);
t13 = axi * t5;
t14 = sin(t10);
t20 = exp(-x * krxi - y * kryi - z * krzi);
t21 = arxr * t20;
t25 = x * krxr + y * kryr + z * krzr;
t26 = cos(t25);
t28 = arxi * t20;
t29 = sin(t25);
Exr = t6 * t11 - t13 * t14 + t21 * t26 - t28 * t29;
t31 = ayr * t5;
t33 = ayi * t5;
t35 = aryr * t20;
t37 = aryi * t20;
Eyr = t31 * t11 - t33 * t14 + t35 * t26 - t37 * t29;
t39 = azr * t5;
t41 = azi * t5;
t43 = arzr * t20;
t45 = arzi * t20;
Ezr = t39 * t11 - t41 * t14 + t43 * t26 - t45 * t29;
Exi = t13 * t11 + t6 * t14 + t21 * t29 + t28 * t26;
Eyi = t33 * t11 + t31 * t14 + t37 * t26 + t35 * t29;
Ezi = t41 * t11 + t39 * t14 + t45 * t26 + t43 * t29;
t63 = exp(-x * ktxi - y * ktyi - z * ktzi);
t64 = atxr * t63;
t68 = x * ktxr + y * ktyr + z * ktzr;
t69 = cos(t68);
t71 = atxi * t63;
t72 = sin(t68);
Etxr = t64 * t69 - t71 * t72;
t74 = atyr * t63;
t76 = atyi * t63;
Etyr = t74 * t69 - t76 * t72;
t78 = atzr * t63;
t80 = atzi * t63;
Etzr = t78 * t69 - t80 * t72;
Etxi = t64 * t72 + t71 * t69;
Etyi = t76 * t69 + t74 * t72;
Etzi = t80 * t69 + t78 * t72;

      // #Include "dispersivePlaneWaveInterface3d.h"

            if( myDomain==0 )
            {
      	Ev[0]=Exr; Ev[1]=Eyr; Ev[2]=Ezr;
      	Ev[3]=Exi; Ev[4]=Eyi; Ev[5]=Ezi;
            }
            else
            {
      	Ev[0]=Etxr; Ev[1]=Etyr; Ev[2]=Etzr;
      	Ev[3]=Etxi; Ev[4]=Etyi; Ev[5]=Etzi;
            }
        }
        
        if( computeMagneticField )
        {
      // Here is the computation of H 
      // *** NOTE: Frequency domain: 
// File generated by Dropbox/DARPA/RPI/adePapers/adegdmi/dispersivePlaneWaveInterface3d.maple
real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t530,t531,t532,t533,t534,t535,t536,t537,t538,t539,t540,t541,t542,t543,t544,t545,t546,t547,t548,t549,t550,t551,t552,t553,t554,t555,t556,t557,t558,t559,t560,t561,t562,t563,t564,t565,t566,t567,t568,t569,t570,t571,t572,t573,t574,t575,t576,t577,t578,t579,t580,t581,t582,t583,t584,t585,t586,t587,t588,t589,t590,t591,t592,t593,t594,t595,t596,t597,t598,t599,t600,t601,t602,t603,t604,t605,t606,t607,t608,t609,t610,t611,t612,t613,t614,t615,t616,t617,t618,t619,t620,t621,t622,t623,t624,t625,t626,t627,t628,t629,t630,t631,t632,t633,t634,t635,t636,t637,t638,t639,t640,t641,t642,t643,t644,t645,t646,t647,t648,t649,t650,t651,t652,t653,t654,t655,t656,t657,t658,t659,t660,t661,t662,t663,t664,t665,t666,t667,t668,t669,t670,t671,t672,t673,t674,t675,t676,t677,t678,t679,t680,t681,t682,t683,t684,t685,t686,t687,t688,t689,t690,t691,t692,t693,t694,t695,t696,t697,t698,t699,t700,t701,t702,t703,t704,t705,t706,t707,t708,t709,t710,t711,t712,t713,t714,t715,t716,t717,t718,t719,t720,t721,t722,t723,t724,t725,t726,t727,t728,t729,t730,t731,t732,t733,t734,t735,t736,t737,t738,t739,t740,t741,t742,t743,t744,t745,t746,t747,t748,t749,t750,t751,t752,t753,t754,t755,t756,t757,t758,t759,t760,t761,t762,t763,t764,t765,t766,t767,t768,t769,t770,t771,t772,t773,t774,t775,t776,t777,t778,t779,t780,t781,t782,t783,t784,t785,t786,t787,t788,t789,t790,t791,t792,t793,t794,t795,t796,t797,t798,t799,t800,t801,t802,t803,t804,t805,t806,t807,t808,t809,t810,t811,t812,t813,t814,t815,t816,t817,t818,t819,t820,t821,t822,t823,t824,t825,t826,t827,t828,t829,t830,t831,t832,t833,t834,t835,t836,t837,t838,t839,t840,t841,t842,t843,t844,t845,t846,t847,t848,t849,t850,t851,t852,t853,t854,t855,t856,t857,t858,t859,t860,t861,t862,t863,t864,t865,t866,t867,t868,t869,t870,t871,t872,t873,t874,t875,t876,t877,t878,t879,t880,t881,t882,t883,t884,t885,t886,t887,t888,t889,t890,t891,t892,t893,t894,t895,t896,t897,t898,t899,t900,t901,t902,t903,t904,t905,t906,t907,t908,t909,t910,t911,t912,t913,t914,t915,t916,t917,t918,t919,t920,t921,t922,t923,t924,t925,t926,t927,t928,t929,t930,t931,t932,t933,t934,t935,t936,t937,t938,t939,t940,t941,t942,t943,t944,t945,t946,t947,t948,t949,t950,t951,t952,t953,t954,t955,t956,t957,t958,t959,t960,t961,t962,t963,t964,t965,t966,t967,t968,t969,t970,t971,t972,t973,t974,t975,t976,t977,t978,t979,t980,t981,t982,t983,t984,t985,t986,t987,t988,t989,t990,t991,t992,t993,t994,t995,t996,t997,t998,t999;
// -------------------------------------------------------------------------
// Need: 
//    s=sr + I*si : complex frequency
//    mu1,mu2 : mu on left and right                                         
//    [ax,ay,az] = amplitude vector of incident wave, ax=axr+I*axi etc.
//    [arx,ary,arz] = amplitude vector of reflected wave, arx=arxr+I*arxi etc.
//    [atx,aty,atz] = amplitude vector of transmitted wave, atx=atxr+I*atxi etc.
//    [kx,ky,kz] = incident wave vector, kx=kxr+I*kxi etc. 
//    [krx,kry,krz] = reflected wave vector, krx=krxr+I*krxi etc. 
//    [ktx,kty,ktz] = reflected wave vector, ktx=ktxr+I*ktxi etc. 
// -------------------------------------------------------------------------
// Evaluated:                                                               
//    [Hxr,Hxi] ,[Hyr,Hyi], [Hzr,Hzi]        : left state                   
//    [Htxr,Htxi] ,[Htyr,Htyi], [Htzr,Htzi]  : right state                  
// -------------------------------------------------------------------------
real Hxr,Hyr,Hzr, Hxi,Hyi,Hzi;
real Htxr,Htyr,Htzr, Htxi,Htyi,Htzi;
t1 = 0.1e1 / mu1;
t6 = exp(-kxi * x - y * kyi - z * kzi);
t10 = x * kxr + y * kyr + z * kzr;
t11 = sin(t10);
t12 = t6 * t11;
t13 = ayi * kzi;
t14 = t13 * sr;
t16 = ayi * kzr;
t17 = t16 * si;
t19 = ayr * kzi;
t20 = t19 * si;
t22 = ayr * kzr;
t23 = t22 * sr;
t25 = azi * kyi;
t26 = t25 * sr;
t28 = azi * kyr;
t29 = t28 * si;
t31 = azr * kyi;
t32 = t31 * si;
t34 = azr * kyr;
t35 = t34 * sr;
t37 = cos(t10);
t38 = t37 * t6;
t39 = t13 * si;
t41 = t16 * sr;
t43 = t19 * sr;
t45 = t22 * si;
t47 = t25 * si;
t49 = t28 * sr;
t51 = t31 * sr;
t53 = t34 * si;
t55 = -t12 * t14 + t12 * t17 + t12 * t20 + t12 * t23 + t12 * t26 - t12 * t29 - t12 * t32 - t12 * t35 + t38 * t39 + t38 * t41 + t38 * t43 - t38 * t45 - t38 * t47 - t38 * t49 - t38 * t51 + t38 * t53;
t60 = exp(-x * krxi - y * kryi - z * krzi);
t64 = x * krxr + y * kryr + z * krzr;
t65 = sin(t64);
t66 = t60 * t65;
t67 = aryi * krzi;
t68 = t67 * sr;
t70 = aryi * krzr;
t71 = t70 * si;
t73 = aryr * krzi;
t74 = t73 * si;
t76 = aryr * krzr;
t77 = t76 * sr;
t79 = arzi * kryi;
t80 = t79 * sr;
t82 = arzi * kryr;
t83 = t82 * si;
t85 = arzr * kryi;
t86 = t85 * si;
t88 = arzr * kryr;
t89 = t88 * sr;
t91 = cos(t64);
t92 = t60 * t91;
t93 = t67 * si;
t95 = t70 * sr;
t97 = t73 * sr;
t99 = t76 * si;
t101 = t79 * si;
t103 = t82 * sr;
t105 = t85 * sr;
t107 = t88 * si;
t109 = -t92 * t101 - t92 * t103 - t92 * t105 + t92 * t107 - t66 * t68 + t66 * t71 + t66 * t74 + t66 * t77 + t66 * t80 - t66 * t83 - t66 * t86 - t66 * t89 + t92 * t93 + t92 * t95 + t92 * t97 - t92 * t99;
t112 = pow(si, 0.2e1);
t113 = pow(sr, 0.2e1);
t115 = 0.1e1 / (t112 + t113);
Hxr = -t1 * (t55 + t109) * t115;
t117 = arxr * krzi;
t118 = t117 * si;
t120 = arxr * krzr;
t121 = t120 * sr;
t123 = arzi * krxi;
t124 = t123 * sr;
t126 = arzi * krxr;
t127 = t126 * si;
t129 = arzr * krxi;
t130 = t129 * si;
t132 = arzr * krxr;
t133 = t132 * sr;
t135 = arxi * krzi;
t136 = t135 * si;
t138 = arxi * krzr;
t139 = t138 * sr;
t141 = t117 * sr;
t143 = t120 * si;
t145 = t123 * si;
t147 = t126 * sr;
t149 = t129 * sr;
t151 = t132 * si;
t153 = axi * kzi;
t154 = t153 * sr;
t156 = axi * kzr;
t157 = t156 * si;
t159 = t66 * t118 - t12 * t154 + t12 * t157 + t66 * t121 + t66 * t124 - t66 * t127 - t66 * t130 - t66 * t133 + t92 * t136 + t92 * t139 + t92 * t141 - t92 * t143 - t92 * t145 - t92 * t147 - t92 * t149 + t92 * t151;
t160 = axr * kzi;
t161 = t160 * si;
t163 = axr * kzr;
t164 = t163 * sr;
t166 = azi * kxi;
t167 = t166 * sr;
t169 = azi * kxr;
t170 = t169 * si;
t172 = azr * kxi;
t173 = t172 * si;
t175 = azr * kxr;
t176 = t175 * sr;
t178 = t153 * si;
t180 = t156 * sr;
t182 = t160 * sr;
t184 = t163 * si;
t186 = t166 * si;
t188 = t169 * sr;
t190 = t172 * sr;
t192 = t175 * si;
t194 = t135 * sr;
t196 = t138 * si;
t198 = t12 * t161 + t12 * t164 + t12 * t167 - t12 * t170 - t12 * t173 - t12 * t176 + t38 * t178 + t38 * t180 + t38 * t182 - t38 * t184 - t38 * t186 - t38 * t188 - t38 * t190 + t38 * t192 - t66 * t194 + t66 * t196;
Hyr = t1 * (t159 + t198) * t115;
t201 = arxr * kryr;
t202 = t201 * sr;
t204 = aryi * krxi;
t205 = t204 * sr;
t207 = aryi * krxr;
t208 = t207 * si;
t210 = aryr * krxi;
t211 = t210 * si;
t213 = aryr * krxr;
t214 = t213 * sr;
t216 = arxi * kryi;
t217 = t216 * si;
t219 = arxi * kryr;
t220 = t219 * sr;
t222 = arxr * kryi;
t223 = t222 * sr;
t225 = t201 * si;
t227 = t204 * si;
t229 = t207 * sr;
t231 = t210 * sr;
t233 = t213 * si;
t235 = axi * kyi;
t236 = t235 * sr;
t238 = axi * kyr;
t239 = t238 * si;
t241 = axr * kyi;
t242 = t241 * si;
t244 = -t12 * t236 + t12 * t239 + t12 * t242 + t66 * t202 + t66 * t205 - t66 * t208 - t66 * t211 - t66 * t214 + t92 * t217 + t92 * t220 + t92 * t223 - t92 * t225 - t92 * t227 - t92 * t229 - t92 * t231 + t92 * t233;
t245 = axr * kyr;
t246 = t245 * sr;
t248 = ayi * kxi;
t249 = t248 * sr;
t251 = ayi * kxr;
t252 = t251 * si;
t254 = ayr * kxi;
t255 = t254 * si;
t257 = ayr * kxr;
t258 = t257 * sr;
t260 = t235 * si;
t262 = t238 * sr;
t264 = t241 * sr;
t266 = t245 * si;
t268 = t248 * si;
t270 = t251 * sr;
t272 = t254 * sr;
t274 = t257 * si;
t276 = t216 * sr;
t278 = t219 * si;
t280 = t222 * si;
t282 = t12 * t246 + t12 * t249 - t12 * t252 - t12 * t255 - t12 * t258 + t38 * t260 + t38 * t262 + t38 * t264 - t38 * t266 - t38 * t268 - t38 * t270 - t38 * t272 + t38 * t274 - t66 * t276 + t66 * t278 + t66 * t280;
Hzr = -t1 * (t244 + t282) * t115;
t302 = -t66 * t101 - t66 * t103 - t66 * t105 + t66 * t107 + t38 * t14 - t38 * t17 - t38 * t20 - t38 * t23 - t38 * t26 + t38 * t29 + t38 * t32 + t38 * t35 + t66 * t93 + t66 * t95 + t66 * t97 - t66 * t99;
t319 = t12 * t39 + t12 * t41 + t12 * t43 - t12 * t45 - t12 * t47 - t12 * t49 - t12 * t51 + t12 * t53 + t92 * t68 - t92 * t71 - t92 * t74 - t92 * t77 - t92 * t80 + t92 * t83 + t92 * t86 + t92 * t89;
Hxi = -t1 * (t302 + t319) * t115;
t339 = -t92 * t118 - t92 * t121 + t66 * t136 + t66 * t139 + t66 * t141 - t66 * t143 - t66 * t145 - t66 * t147 - t66 * t149 + t66 * t151 - t38 * t167 + t38 * t170 + t38 * t173 + t38 * t176 + t92 * t194 - t92 * t196;
t356 = t12 * t178 + t12 * t180 + t12 * t182 - t12 * t184 - t12 * t186 - t12 * t188 - t12 * t190 + t12 * t192 - t92 * t124 + t92 * t127 + t92 * t130 + t92 * t133 + t38 * t154 - t38 * t157 - t38 * t161 - t38 * t164;
Hyi = t1 * (t339 + t356) * t115;
t375 = t66 * t217 + t66 * t220 + t66 * t223 - t66 * t225 - t66 * t227 - t66 * t229 - t66 * t231 + t66 * t233 - t38 * t239 - t38 * t242 - t38 * t246 - t38 * t249 + t38 * t252 + t38 * t255 + t38 * t258 + t92 * t276;
t392 = t12 * t260 + t12 * t262 + t12 * t264 - t12 * t266 - t12 * t268 - t12 * t270 - t12 * t272 + t12 * t274 - t92 * t202 - t92 * t205 + t92 * t208 + t92 * t211 + t92 * t214 + t38 * t236 - t92 * t278 - t92 * t280;
Hzi = -t1 * (t375 + t392) * t115;
t401 = exp(-x * ktxi - y * ktyi - z * ktzi);
t402 = 0.1e1 / mu2 * t401;
t406 = x * ktxr + y * ktyr + z * ktzr;
t407 = cos(t406);
t408 = t407 * atyi;
t409 = ktzi * si;
t411 = ktzr * sr;
t413 = t407 * atyr;
t414 = ktzi * sr;
t416 = ktzr * si;
t418 = t407 * atzi;
t419 = ktyi * si;
t421 = ktyr * sr;
t423 = t407 * atzr;
t424 = ktyi * sr;
t426 = ktyr * si;
t428 = sin(t406);
t429 = t428 * atyi;
t432 = t428 * atyr;
t435 = t428 * atzi;
t438 = t428 * atzr;
t441 = t408 * t409 + t408 * t411 + t432 * t409 + t432 * t411 + t413 * t414 - t413 * t416 - t429 * t414 + t429 * t416 - t418 * t419 - t418 * t421 - t438 * t419 - t438 * t421 - t423 * t424 + t423 * t426 + t435 * t424 - t435 * t426;
Htxr = -t402 * t441 * t115;
t444 = t407 * atxi;
t447 = t407 * atxr;
t450 = ktxi * si;
t452 = ktxr * sr;
t454 = ktxi * sr;
t456 = ktxr * si;
t458 = t428 * atxi;
t461 = t428 * atxr;
t468 = t444 * t409 + t461 * t409 + t444 * t411 + t461 * t411 + t447 * t414 - t458 * t414 - t447 * t416 + t458 * t416 - t418 * t450 - t418 * t452 - t423 * t454 + t423 * t456 + t435 * t454 - t435 * t456 - t438 * t450 - t438 * t452;
Htyr = t402 * t468 * t115;
t486 = -t408 * t450 - t408 * t452 - t413 * t454 + t413 * t456 + t444 * t419 + t461 * t419 + t444 * t421 + t461 * t421 + t447 * t424 - t458 * t424 - t447 * t426 + t458 * t426 + t429 * t454 - t429 * t456 - t432 * t450 - t432 * t452;
Htzr = -t402 * t486 * t115;
t505 = t408 * t414 - t408 * t416 - t413 * t409 + t429 * t409 - t413 * t411 + t429 * t411 + t432 * t414 - t432 * t416 - t418 * t424 + t418 * t426 + t423 * t419 - t435 * t419 + t423 * t421 - t435 * t421 - t438 * t424 + t438 * t426;
Htxi = -t402 * t505 * t115;
t524 = -t447 * t409 + t458 * t409 - t447 * t411 + t458 * t411 + t444 * t414 + t461 * t414 - t444 * t416 - t461 * t416 - t418 * t454 + t418 * t456 + t423 * t450 + t423 * t452 - t435 * t450 - t435 * t452 - t438 * t454 + t438 * t456;
Htyi = t402 * t524 * t115;
t542 = -t408 * t454 + t408 * t456 + t413 * t450 + t413 * t452 - t447 * t419 + t458 * t419 - t447 * t421 + t458 * t421 + t444 * t424 + t461 * t424 - t444 * t426 - t461 * t426 - t429 * t450 - t429 * t452 - t432 * t454 + t432 * t456;
Htzi = -t402 * t542 * t115;

      // #Include "dispersivePlaneWaveInterfaceH3d.h"

      // printF(" myDomain=%i: Hz=(%g,%g) Htz=(%g,%g)\n",myDomain,Hzr,Hzi,Htzr,Htzi);
            
            if( myDomain==0 )
            {
      	Hv[0]=Hxr; Hv[1]=Hyr; Hv[2]=Hzr;
      	Hv[3]=Hxi; Hv[4]=Hyi; Hv[5]=Hzi;
            }
            else
            {
      	Hv[0]=Htxr; Hv[1]=Htyr; Hv[2]=Htzr;
      	Hv[3]=Htxi; Hv[4]=Htyi; Hv[5]=Htzi;
            }

        }

        
    }
        
    return 0;
}







// ==========================================================================================
/// \brief  Evaluate the solution and save in an array.
///
/// \param numberOfTimeDerivatives (input) : evaluate this many time-derivatives of the solution.
/// \param computeMagneticField (input): if true return the magnetic field in 3D (in 2D the magnetic field is always computed). 
// ==========================================================================================
int PlaneInterfaceExactSolution::
eval(real t, CompositeGrid & cg, int grid, 
          realArray & ua, realArray & pv,
          const Index & I1a, const Index &I2a, const Index &I3a, 
          int numberOfTimeDerivatives /* = 0 */,
          bool computeMagneticField /* = false */ )
{

  // domain number for this grid: 
    const int myDomain = cg.domainNumber(grid);


    if( t <= 0. )
        printF("--PlaneInterfaceExactSolution--  eval on grid=%i, domain=%i, at t=%9.3e\n",grid,myDomain,t);

    MappedGrid & mg = cg[grid];
    const int numberOfDimensions = cg.numberOfDimensions();
    
    
    OV_GET_SERIAL_ARRAY(real,ua,uLocal);

    Index I1=I1a, I2=I2a, I3=I3a;
    bool ok = ParallelUtility::getLocalArrayBounds(ua,uLocal,I1,I2,I3,1);   
    if( !ok ) return 0;  // no points on this processor (NOTE: no communication should be done after this point)

  // -- we optimize for Cartesian grids (we can avoid creating the vertex array)
    const bool isRectangular=mg.isRectangular();
    if( !isRectangular )
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
    OV_GET_SERIAL_ARRAY(real,mg.center(),xLocal);

    real dvx[3]={1.,1.,1.}, xab[2][3]={{0.,0.,0.},{0.,0.,0.}};
    int iv0[3]={0,0,0}; //
    int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];  // NOTE: iv[0]==i1, iv[1]==i2, iv[2]==i3
    real xv[3]={0.,0.,0.};
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

    
    const real & a = dbase.get<real>("amp");

    LocalReal (&s)[2] = dbase.get<LocalReal[2]>("s");

    real & sr = s[0], &si = s[1];
    real (&av)[3] = dbase.get<real[3]>("av");

    real (&kv)[3] = dbase.get<real[3]>("kv");

    real & eps1 = dbase.get<real>("eps1");
    real & eps2 = dbase.get<real>("eps2");
    real & mu1 = dbase.get<real>("mu1");
    real & mu2 = dbase.get<real>("mu2");

  // -- for 2d eval: 

    real (&kvr)[3] = dbase.get<real[3]>("kvr");
    real (&kvi)[3] = dbase.get<real[3]>("kvi");

    real (&kpvr)[3] = dbase.get<real[3]>("kpvr");
    real (&kpvi)[3] = dbase.get<real[3]>("kpvi");

    real &kxr =kvr[0], &kxi=kvi[0], &kyr=kvr[1], &kyi=kvi[1];   // Incident wave number (complex)
    real &kxpr =kpvr[0], &kxpi=kpvi[0], &kypr=kpvr[1], &kypi=kvi[1];   // 


    int & numberOfPolarizationVectors1 = dbase.get<int>("numberOfPolarizationVectors1");
    int & numberOfPolarizationVectors2 = dbase.get<int>("numberOfPolarizationVectors2");

    real (&chi1r)[10] = dbase.get<real[10]>("chi1r");
    real (&chi1i)[10] = dbase.get<real[10]>("chi1i");

    real (&chi2r)[10] = dbase.get<real[10]>("chi2r");
    real (&chi2i)[10] = dbase.get<real[10]>("chi2i");

  // -- for 3d eval: 

    real (&arr)[3] = dbase.get<real[3]>("arr");
    real (&ari)[3] = dbase.get<real[3]>("ari");
    real (&atr)[3] = dbase.get<real[3]>("atr");
    real (&ati)[3] = dbase.get<real[3]>("ati");

    real (&krr)[3] = dbase.get<real[3]>("krr");
    real (&kri)[3] = dbase.get<real[3]>("kri");
    real (&ktr)[3] = dbase.get<real[3]>("ktr");
    real (&kti)[3] = dbase.get<real[3]>("kti");



    real (&normalPlaneMaterialInterface)[3] = dbase.get<real[3]>("normalPlaneMaterialInterface");
    real (&x0PlaneMaterialInterface)[3] = dbase.get<real[3]>("x0PlaneMaterialInterface");

    real *nv = normalPlaneMaterialInterface;
    real *x0 = x0PlaneMaterialInterface;


  // -- Store components here: 
    const int ex=0, ey=1, ez=2;
    const int hx=3, hy=4, hz=numberOfDimensions==2 ? 2 :  5;

    if( computeMagneticField && numberOfDimensions==3 && ua.getLength(3)<6 )
    {
        printF(" PlaneInterfaceExactSolution::ERROR: Not enough spacein ua to hold the H field\n");
        OV_ABORT("error");
    }
    

  // --- Get Arrays for the dispersive model ----

  // realMappedGridFunction & pCur = getDispersionModelMappedGridFunction( grid,current );
    RealArray pLocal;
    if( (myDomain==0 && numberOfPolarizationVectors1>0) ||
            (myDomain==1 && numberOfPolarizationVectors2>0)  )
    {
        OV_GET_SERIAL_ARRAY(real, pv,pLoc);
        pLocal.reference(pLoc);
    }

        
    real x,y,z=0.;
    if( numberOfTimeDerivatives==0 )
    {
        if( numberOfDimensions==2 )
        {
      // ----------- 2D --------------
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
      	if( !isRectangular )
      	{
        	  x= xLocal(i1,i2,i3,0)-x0[0];   // shift point to reference coordinates 
        	  y= xLocal(i1,i2,i3,1)-x0[1];
      	}
      	else
      	{
        	  x=XC(iv,0)-x0[0];
        	  y=XC(iv,1)-x0[1];
      	}

	// -- un-rotate the point (x,y) to the reference space and evaluate the solution 
      	real ct= nv[0], st=nv[1];  // cos(theta), sin(theta)
	// real xa=x-x0[0], ya=y-x0[1];
      	real xa=x, ya=y;
      	x =  ct*xa + st*ya;
      	y = -st*xa + ct*ya;
                    

	// Here are the statements to eval the solution: 
// File generated by overtureFramework/cg/mx/codes/dispersivePlaneWaveInterface.maple
// File generated by DropBox/DMX/codes/dispersivePlaneWaveInterface.maple
real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t530,t531,t532,t533,t534,t535,t536,t537,t538,t539,t540,t541,t542,t543,t544,t545,t546,t547,t548,t549,t550,t551,t552,t553,t554,t555,t556,t557,t558,t559,t560,t561,t562,t563,t564,t565,t566,t567,t568,t569,t570,t571,t572,t573,t574,t575,t576,t577,t578,t579,t580,t581,t582,t583,t584,t585,t586,t587,t588,t589,t590,t591,t592,t593,t594,t595,t596,t597,t598,t599,t600;
// -------------------------------------------------------------------------
// Need: 
//    a = amplitude of the wave, e.g. a=1
//    x,y,t
//    sr,si : s= sr + I*si 
//    alphaP
//    kxr,kxi,  kyr,kyi,     : complex wave number on left
//    kxpr,kxpi,  kypr,kypi, : complex wave number on right (plus)
// -------------------------------------------------------------------------
// Evaluated:                                                              
//    Exr,Eyr   : left state                                                
//    Expr,Eypr : right state                                               
// -------------------------------------------------------------------------
real kNorm = sqrt( kxr*kxr + kxi*kxi + kyr*kyr + kyi*kyi);
real khxr = kxr/kNorm, khxi = kxi/kNorm, khyr=kyr/kNorm, khyi=kyi/kNorm; 
real kpNorm = sqrt( kxpr*kxpr + kxpi*kxpi + kypr*kypr + kypi*kypi);
real khxpr = kxpr/kpNorm, khxpi=kxpi/kpNorm, khypr=kypr/kpNorm, khypi=kypi/kpNorm; 
real kappar,kappai, betar, betai, rr,ri, taur,taui;
real Exr,Eyr,Hzr, Expr,Eypr,Hzpr;
real Exi,Eyi,Hzi, Expi,Eypi,Hzpi;
t1 = khxi * khxpi;
t2 = khxr * khxpr;
t4 = pow(khxpi, 0.2e1);
t5 = pow(khxpr, 0.2e1);
t7 = 0.1e1 / (t4 + t5);
kappar = (t1 + t2) * t7;
t8 = khxi * khxpr;
t9 = khxr * khxpi;
kappai = (t8 - t9) * t7;
t11 = 0.1e1 / mu1;
t12 = mu2 * t11;
t13 = kxi * kxpi;
t15 = kxpr * kxr;
t17 = kxi * kxpr;
t19 = kxpi * kxr;
t21 = khxi * khypi;
t22 = kxi * kypi;
t24 = kxr * kypr;
t26 = khxi * khypr;
t27 = kxi * kypr;
t29 = kxr * kypi;
t33 = khxpi * khyi;
t34 = kxpi * kyi;
t36 = kxpr * kyr;
t38 = khxpi * khyr;
t39 = kxpi * kyr;
t41 = kxpr * kyi;
t45 = t1 * t13 + t1 * t15 + t2 * t13 + t2 * t15 - t8 * t17 + t9 * t17 + t8 * t19 - t9 * t19 + t21 * t22 + t21 * t24 - t26 * t27 + t26 * t29 + t33 * t34 + t33 * t36 - t38 * t39 + t38 * t41;
t46 = khxpr * khyi;
t49 = khxpr * khyr;
t52 = khxr * khypi;
t55 = khxr * khypr;
t58 = khyi * khypi;
t59 = kyi * kypi;
t61 = kypr * kyr;
t63 = khyi * khypr;
t64 = kyi * kypr;
t66 = kypi * kyr;
t68 = khypi * khyr;
t71 = khypr * khyr;
t74 = t55 * t22 + t55 * t24 + t52 * t27 - t52 * t29 + t49 * t34 + t49 * t36 + t46 * t39 - t46 * t41 + t58 * t59 + t58 * t61 + t71 * t59 + t71 * t61 - t63 * t64 + t63 * t66 + t68 * t64 - t68 * t66;
t76 = pow(kxpi, 0.2e1);
t78 = pow(kxpr, 0.2e1);
t80 = khxpi * khypi;
t81 = kxpi * kypi;
t84 = kxpr * kypr;
t87 = khxpi * khypr;
t88 = kxpi * kypr;
t91 = kxpr * kypi;
t96 = khxpr * khypi;
t101 = khxpr * khypr;
t106 = pow(khypi, 0.2e1);
t107 = pow(kypi, 0.2e1);
t109 = pow(kypr, 0.2e1);
t111 = pow(khypr, 0.2e1);
t114 = 0.2e1 * t101 * t81 + 0.2e1 * t101 * t84 + t106 * t107 + t106 * t109 + t111 * t107 + t111 * t109 + t4 * t76 + t4 * t78 + t5 * t76 + t5 * t78 + 0.2e1 * t80 * t81 + 0.2e1 * t80 * t84 - 0.2e1 * t87 * t88 + 0.2e1 * t87 * t91 + 0.2e1 * t96 * t88 - 0.2e1 * t96 * t91;
t115 = 0.1e1 / t114;
betar = t12 * (t45 + t74) * t115;
t133 = t1 * t17 - t1 * t19 + t8 * t13 - t9 * t13 + t8 * t15 - t9 * t15 + t2 * t17 - t2 * t19 + t21 * t27 - t21 * t29 + t26 * t22 + t26 * t24 - t33 * t39 + t33 * t41 - t38 * t34 - t38 * t36;
t150 = -t52 * t22 - t52 * t24 + t55 * t27 - t55 * t29 + t46 * t34 + t46 * t36 - t49 * t39 + t49 * t41 + t58 * t64 - t58 * t66 + t63 * t59 - t68 * t59 + t63 * t61 - t68 * t61 + t71 * t64 - t71 * t66;
betai = t12 * (t133 + t150) * t115;
t153 = pow(betai, 0.2e1);
t154 = pow(betar, 0.2e1);
t155 = pow(kappai, 0.2e1);
t156 = pow(kappar, 0.2e1);
t163 = 0.1e1 / (0.2e1 * kappai * betai + 0.2e1 * kappar * betar + t153 + t154 + t155 + t156);
rr = (t153 + t154 - t155 - t156) * t163;
ri = 0.2e1 * (betai * kappar - betar * kappai) * t163;
taur = 0.2e1 * (betar * t155 + betar * t156 + t153 * kappar + t154 * kappar) * t163;
taui = 0.2e1 * (betai * t155 + betai * t156 + kappai * t153 + t154 * kappai) * t163;
t180 = x * kxi;
t181 = y * kyi;
t182 = t * sr;
t184 = exp(t180 - t181 + t182);
t185 = x * kxr;
t186 = y * kyr;
t187 = t * si;
t188 = t185 - t186 - t187;
t189 = sin(t188);
t190 = t184 * t189;
t191 = khyi * rr;
t193 = khyr * ri;
t195 = cos(t188);
t196 = t184 * t195;
t197 = khyi * ri;
t199 = khyr * rr;
t202 = exp(-t180 - t181 + t182);
t203 = t185 + t186 + t187;
t204 = cos(t203);
t205 = t202 * t204;
t207 = sin(t203);
t208 = t202 * t207;
Exr = -a * (-t208 * khyi + t205 * khyr - t190 * t191 - t190 * t193 + t196 * t197 - t196 * t199);
t212 = khxi * rr;
t214 = khxr * ri;
t216 = khxi * ri;
t218 = khxr * rr;
Eyr = a * (-t208 * khxi + t205 * khxr + t190 * t212 + t190 * t214 - t196 * t216 + t196 * t218);
t226 = exp(-x * kxpi - y * kypi + t182);
t229 = x * kxpr + y * kypr + t187;
t230 = cos(t229);
t231 = t226 * t230;
t232 = khypi * taui;
t234 = khypr * taur;
t236 = sin(t229);
t237 = t226 * t236;
t238 = khypi * taur;
t240 = khypr * taui;
Expr = -a * (-t231 * t232 + t231 * t234 - t237 * t238 - t237 * t240);
t244 = khxpi * taui;
t246 = khxpr * taur;
t248 = khxpi * taur;
t250 = khxpr * taui;
Eypr = a * (-t231 * t244 + t231 * t246 - t237 * t248 - t237 * t250);
Exi = -a * (t205 * khyi + t208 * khyr - t190 * t197 + t190 * t199 - t196 * t191 - t196 * t193);
Eyi = a * (t205 * khxi + t208 * khxr + t190 * t216 - t190 * t218 + t196 * t212 + t196 * t214);
Expi = -a * (t231 * t238 + t231 * t240 - t237 * t232 + t237 * t234);
Eypi = a * (t231 * t248 + t231 * t250 - t237 * t244 + t237 * t246);
t279 = t11 * a;
t280 = t190 * khxi;
t281 = kxi * ri;
t282 = t281 * si;
t284 = kxi * rr;
t285 = t284 * sr;
t287 = kxr * ri;
t288 = t287 * sr;
t290 = kxr * rr;
t291 = t290 * si;
t293 = t190 * khxr;
t294 = t281 * sr;
t296 = t284 * si;
t298 = t287 * si;
t300 = t290 * sr;
t302 = t190 * khyi;
t303 = kyi * ri;
t304 = t303 * si;
t306 = kyi * rr;
t307 = t306 * sr;
t309 = kyr * ri;
t310 = t309 * sr;
t312 = kyr * rr;
t313 = t312 * si;
t315 = -t280 * t282 - t280 * t285 - t280 * t288 + t280 * t291 - t293 * t294 + t293 * t296 + t293 * t298 + t293 * t300 - t302 * t304 - t302 * t307 - t302 * t310 + t302 * t313;
t316 = t190 * khyr;
t317 = t303 * sr;
t319 = t306 * si;
t321 = t309 * si;
t323 = t312 * sr;
t325 = t196 * khxi;
t330 = t196 * khxr;
t335 = -t330 * t282 - t330 * t285 - t330 * t288 + t330 * t291 + t325 * t294 - t325 * t296 - t325 * t298 - t325 * t300 - t316 * t317 + t316 * t319 + t316 * t321 + t316 * t323;
t337 = t196 * khyi;
t342 = t196 * khyr;
t347 = khxi * kxi;
t348 = t347 * si;
t350 = khxi * kxr;
t351 = t350 * sr;
t353 = khxr * kxi;
t354 = t353 * sr;
t356 = khxr * kxr;
t357 = t356 * si;
t359 = t205 * t348 + t205 * t351 + t205 * t354 - t205 * t357 - t342 * t304 - t342 * t307 - t342 * t310 + t342 * t313 + t337 * t317 - t337 * t319 - t337 * t321 - t337 * t323;
t360 = khyi * kyi;
t361 = t360 * si;
t363 = khyi * kyr;
t364 = t363 * sr;
t366 = khyr * kyi;
t367 = t366 * sr;
t369 = khyr * kyr;
t370 = t369 * si;
t372 = t347 * sr;
t374 = t350 * si;
t376 = t353 * si;
t378 = t356 * sr;
t380 = t360 * sr;
t382 = t363 * si;
t384 = t366 * si;
t386 = t369 * sr;
t388 = t205 * t361 + t205 * t364 + t205 * t367 - t205 * t370 - t208 * t372 + t208 * t374 + t208 * t376 + t208 * t378 - t208 * t380 + t208 * t382 + t208 * t384 + t208 * t386;
t391 = pow(si, 0.2e1);
t392 = pow(sr, 0.2e1);
t394 = 0.1e1 / (t391 + t392);
Hzr = t279 * (t315 + t335 + t359 + t388) * t394;
t397 = 0.1e1 / mu2 * a;
t398 = t236 * khypr;
t399 = kypi * si;
t400 = t399 * taur;
t402 = kypi * sr;
t403 = t402 * taui;
t405 = kypr * si;
t406 = t405 * taui;
t408 = kypr * sr;
t409 = t408 * taur;
t411 = t230 * khxpi;
t412 = kxpi * si;
t413 = t412 * taur;
t415 = kxpi * sr;
t416 = t415 * taui;
t418 = kxpr * si;
t419 = t418 * taui;
t421 = kxpr * sr;
t422 = t421 * taur;
t424 = t230 * khxpr;
t425 = t412 * taui;
t427 = t415 * taur;
t429 = t418 * taur;
t431 = t421 * taui;
t433 = t230 * khypi;
t438 = t398 * t400 - t398 * t403 + t398 * t406 + t398 * t409 + t433 * t400 - t433 * t403 + t433 * t406 + t433 * t409 + t411 * t413 - t411 * t416 + t411 * t419 + t411 * t422 + t424 * t425 + t424 * t427 - t424 * t429 + t424 * t431;
t439 = t230 * khypr;
t440 = t399 * taui;
t442 = t402 * taur;
t444 = t405 * taur;
t446 = t408 * taui;
t448 = t236 * khxpi;
t453 = t236 * khxpr;
t458 = t236 * khypi;
t463 = t453 * t413 - t453 * t416 + t453 * t419 + t453 * t422 - t448 * t425 - t448 * t427 + t448 * t429 - t448 * t431 + t439 * t440 + t439 * t442 - t439 * t444 + t439 * t446 - t458 * t440 - t458 * t442 + t458 * t444 - t458 * t446;
Hzpr = t397 * t226 * (t438 + t463) * t394;
t479 = -t280 * t294 + t280 * t296 + t280 * t298 + t280 * t300 + t293 * t282 + t293 * t285 + t293 * t288 - t293 * t291 - t302 * t317 + t302 * t319 + t302 * t321 + t302 * t323;
t492 = -t325 * t282 - t325 * t285 - t325 * t288 + t325 * t291 - t330 * t294 + t330 * t296 + t330 * t298 + t330 * t300 + t316 * t304 + t316 * t307 + t316 * t310 - t316 * t313;
t506 = t205 * t372 - t205 * t374 - t205 * t376 - t205 * t378 - t337 * t304 - t337 * t307 - t337 * t310 + t337 * t313 - t342 * t317 + t342 * t319 + t342 * t321 + t342 * t323;
t519 = t205 * t380 - t205 * t382 - t205 * t384 - t205 * t386 + t208 * t348 + t208 * t351 + t208 * t354 - t208 * t357 + t208 * t361 + t208 * t364 + t208 * t367 - t208 * t370;
Hzi = t279 * (t479 + t492 + t506 + t519) * t394;
t539 = t398 * t440 + t398 * t442 - t398 * t444 + t398 * t446 + t458 * t400 - t458 * t403 + t458 * t406 + t458 * t409 + t411 * t425 + t411 * t427 + t448 * t419 + t448 * t422 + t453 * t425 + t453 * t427 - t453 * t429 + t453 * t431;
t556 = -t439 * t400 + t439 * t403 - t439 * t406 - t439 * t409 - t411 * t429 + t411 * t431 - t424 * t413 + t448 * t413 + t424 * t416 - t448 * t416 - t424 * t419 - t424 * t422 + t433 * t440 + t433 * t442 - t433 * t444 + t433 * t446;
Hzpi = t397 * t226 * (t539 + t556) * t394;


	// ---- rotate the field from the reference space to the rotated space ---
      	real Exra=Exr, Eyra=Eyr;
      	Exr =  ct*Exra - st*Eyra;
      	Eyr =  st*Exra + ct*Eyra;
                    
      	real Exia=Exi, Eyia=Eyi;
      	Exi =  ct*Exia - st*Eyia;
      	Eyi =  st*Exia + ct*Eyia;
                    
      	real Expra=Expr, Eypra=Eypr;
      	Expr =  ct*Expra - st*Eypra;
      	Eypr =  st*Expra + ct*Eypra;
                    
      	real Expia=Expi, Eypia=Eypi;
      	Expi =  ct*Expia - st*Eypia;
      	Eypi =  st*Expia + ct*Eypia;
                    

                    
      	if( false )
      	{
        	  printF("(i1,i2)=(%3i,%3i): kNorm=%g, kpNorm=%g, kappa=(%g,%g) beta=(%g,%g)\n",
             		 i1,i2,kNorm,kpNorm,kappar,kappai,betar,betai);
        	  printF("    : eps1=%g, eps2=%g, r=(%g,%g) tau=(%g,%g) \n",eps1,eps2,rr,ri,taur,taui);
	  // printF("    : chiSum1=(%g,%g) chiSum2=(%g,%g) \n",chiSum1r,chiSum1i,chiSum2r,chiSum2i);
        	  printF("    : Exr=%g, Eyr=%g, Exi=%g, Eyi=%g Hzr=%g Hzi=%g\n",Exr,Eyr,Exi,Eyi,Hzr,Hzi);
        	  printF("    : Expr=%g, Eypr=%g, Expi=%g, Eypi=%g Hzpr=%g Hzpi=%g\n",Expr,Eypr,Expi,Eypi,Hzpr,Hzpi);

        	  OV_ABORT("finish me");
      	}
                    

      	if( myDomain==0 )
      	{
        	  uLocal(i1,i2,i3,ex) = Exr;
        	  uLocal(i1,i2,i3,ey) = Eyr;
        	  uLocal(i1,i2,i3,hz) = Hzr;

        	  for( int iv=0; iv<numberOfPolarizationVectors1; iv++ )
        	  {
          	    const int pc= iv*numberOfDimensions;
          	    pLocal(i1,i2,i3,pc  ) = eps1*( chi1r[iv]*Exr - chi1i[iv]*Exi );
          	    pLocal(i1,i2,i3,pc+1) = eps1*( chi1r[iv]*Eyr - chi1i[iv]*Eyi );
        	  }

      	}
      	else
      	{
        	  uLocal(i1,i2,i3,ex) = Expr;
        	  uLocal(i1,i2,i3,ey) = Eypr;
        	  uLocal(i1,i2,i3,hz) = Hzpr;
        	  for( int iv=0; iv<numberOfPolarizationVectors2; iv++ )
        	  {
          	    const int pc= iv*numberOfDimensions;
          	    pLocal(i1,i2,i3,pc  ) = eps2*( chi2r[iv]*Expr - chi2i[iv]*Expi );
          	    pLocal(i1,i2,i3,pc+1) = eps2*( chi2r[iv]*Eypr - chi2i[iv]*Eypi );
        	  }
      	}

            }
        }
        else
        {
      // -----------------------------
      // ----------- 3D --------------
      // -----------------------------

            const real &ax=av[0], &ay=av[1], &az=av[2]; 

      // These need to be set for the solution evaluation below:
            real axr=ax, axi=0., ayr=ay, ayi=0., azr=az, azi=0.;
            real arxr=arr[0], arxi=ari[0], aryr=arr[1], aryi=ari[1], arzr=arr[2], arzi=ari[2];
            real atxr=atr[0], atxi=ati[0], atyr=atr[1], atyi=ati[1], atzr=atr[2], atzi=ati[2];
                    
      // real kxr=twoPi*kx, kxi=0., kyr=twoPi*ky, kyi=0., kzr=twoPi*kz, kzi=0.;
            real kxr=kv[0], kxi=0., kyr=kv[1], kyi=0., kzr=kv[2], kzi=0.;
            real krxr=krr[0], krxi=kri[0], kryr=krr[1], kryi=kri[1], krzr=krr[2], krzi=kri[2];
            real ktxr=ktr[0], ktxi=kti[0], ktyr=ktr[1], ktyi=kti[1], ktzr=ktr[2], ktzi=kti[2];


            FOR_3D(i1,i2,i3,I1,I2,I3)
            { 
      	if( !isRectangular )
      	{
        	  x= xLocal(i1,i2,i3,0)-x0[0];
        	  y= xLocal(i1,i2,i3,1)-x0[1];
        	  z= xLocal(i1,i2,i3,2)-x0[2];
      	}
      	else
      	{
        	  x=XC(iv,0)-x0[0];
        	  y=XC(iv,1)-x0[1];
        	  z=XC(iv,2)-x0[2];
      	}

	// Here are the statements to eval the time-dependent solution: 
// File generated by Dropbox/DARPA/RPI/adePapers/adegdmi/dispersivePlaneWaveInterface3d.maple
real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150;
// -------------------------------------------------------------------------
// Need: 
//    s=sr + I*si : complex frequency
//    [ax,ay,az] = amplitude vector of incident wave, ax=axr+I*axi etc.
//    [arx,ary,arz] = amplitude vector of reflected wave, arx=arxr+I*arxi etc.
//    [atx,aty,atz] = amplitude vector of transmitted wave, atx=atxr+I*atxi etc.
//    [kx,ky,kz] = incident wave vector, kx=kxr+I*kxi etc. 
//    [krx,kry,krz] = reflected wave vector, krx=krxr+I*krxi etc. 
//    [ktx,kty,ktz] = reflected wave vector, ktx=ktxr+I*ktxi etc. 
// -------------------------------------------------------------------------
// Evaluated:                                                               
//    [Exr,Exi] ,[Eyr,Eyi], [Ezr,Ezi]        : left state                   
//    [Etxr,Etxi] ,[Etyr,Etyi], [Etzr,Etzi]  : right state                  
// -------------------------------------------------------------------------
real Exr,Eyr,Ezr, Exi,Eyi,Ezi;
real Etxr,Etyr,Etzr, Etxi,Etyi,Etzi;
t2 = exp(t * sr);
t3 = t * si;
t4 = cos(t3);
t5 = t2 * t4;
t10 = exp(-x * krxi - y * kryi - z * krzi);
t14 = x * krxr + y * kryr + z * krzr;
t15 = sin(t14);
t16 = t10 * t15;
t17 = t16 * arxi;
t19 = cos(t14);
t20 = t10 * t19;
t21 = t20 * arxr;
t27 = exp(-x * kxi - y * kyi - z * kzi);
t31 = x * kxr + y * kyr + z * kzr;
t32 = sin(t31);
t33 = t27 * t32;
t34 = t33 * axi;
t36 = cos(t31);
t37 = t27 * t36;
t38 = t37 * axr;
t40 = sin(t3);
t41 = t2 * t40;
t42 = t16 * arxr;
t44 = t20 * arxi;
t46 = t33 * axr;
t48 = t37 * axi;
Exr = -t5 * t17 + t5 * t21 - t5 * t34 + t5 * t38 - t41 * t42 - t41 * t44 - t41 * t46 - t41 * t48;
t50 = t16 * aryi;
t52 = t20 * aryr;
t54 = t33 * ayi;
t56 = t37 * ayr;
t58 = t16 * aryr;
t60 = t20 * aryi;
t62 = t33 * ayr;
t64 = t37 * ayi;
Eyr = -t41 * t58 - t41 * t60 - t41 * t62 - t41 * t64 - t5 * t50 + t5 * t52 - t5 * t54 + t5 * t56;
t66 = t16 * arzi;
t68 = t20 * arzr;
t70 = t33 * azi;
t72 = t37 * azr;
t74 = t16 * arzr;
t76 = t20 * arzi;
t78 = t33 * azr;
t80 = t37 * azi;
Ezr = -t41 * t74 - t41 * t76 - t41 * t78 - t41 * t80 - t5 * t66 + t5 * t68 - t5 * t70 + t5 * t72;
Exi = -t41 * t17 + t41 * t21 - t41 * t34 + t41 * t38 + t5 * t42 + t5 * t44 + t5 * t46 + t5 * t48;
Eyi = -t41 * t50 + t41 * t52 - t41 * t54 + t41 * t56 + t5 * t58 + t5 * t60 + t5 * t62 + t5 * t64;
Ezi = -t41 * t66 + t41 * t68 - t41 * t70 + t41 * t72 + t5 * t74 + t5 * t76 + t5 * t78 + t5 * t80;
t110 = exp(-x * ktxi - y * ktyi - z * ktzi);
t114 = x * ktxr + y * ktyr + z * ktzr;
t115 = cos(t114);
t116 = t110 * t115;
t117 = t5 * atxr;
t119 = t41 * atxi;
t121 = sin(t114);
t122 = t110 * t121;
t123 = t5 * atxi;
t125 = t41 * atxr;
Etxr = t116 * t117 - t116 * t119 - t122 * t123 - t122 * t125;
t127 = t5 * atyr;
t129 = t41 * atyi;
t131 = t5 * atyi;
t133 = t41 * atyr;
Etyr = t116 * t127 - t116 * t129 - t122 * t131 - t122 * t133;
t135 = t5 * atzr;
t137 = t41 * atzi;
t139 = t5 * atzi;
t141 = t41 * atzr;
Etzr = t116 * t135 - t116 * t137 - t122 * t139 - t122 * t141;
Etxi = t116 * t123 + t116 * t125 + t122 * t117 - t122 * t119;
Etyi = t116 * t131 + t116 * t133 + t122 * t127 - t122 * t129;
Etzi = t116 * t139 + t116 * t141 + t122 * t135 - t122 * t137;


      	if( myDomain==0 )
      	{
        	  uLocal(i1,i2,i3,ex) = Exr;
        	  uLocal(i1,i2,i3,ey) = Eyr;
        	  uLocal(i1,i2,i3,ez) = Ezr;

        	  for( int iv=0; iv<numberOfPolarizationVectors1; iv++ )
        	  {
          	    const int pc= iv*numberOfDimensions;
          	    pLocal(i1,i2,i3,pc  ) = eps1*( chi1r[iv]*Exr - chi1i[iv]*Exi );
          	    pLocal(i1,i2,i3,pc+1) = eps1*( chi1r[iv]*Eyr - chi1i[iv]*Eyi );
          	    pLocal(i1,i2,i3,pc+2) = eps1*( chi1r[iv]*Ezr - chi1i[iv]*Ezi );
        	  }

      	}
      	else
      	{
        	  uLocal(i1,i2,i3,ex) = Etxr;
        	  uLocal(i1,i2,i3,ey) = Etyr;
        	  uLocal(i1,i2,i3,ez) = Etzr;
        	  for( int iv=0; iv<numberOfPolarizationVectors2; iv++ )
        	  {
          	    const int pc= iv*numberOfDimensions;
          	    pLocal(i1,i2,i3,pc  ) = eps2*( chi2r[iv]*Etxr - chi2i[iv]*Etxi );
          	    pLocal(i1,i2,i3,pc+1) = eps2*( chi2r[iv]*Etyr - chi2i[iv]*Etyi );
          	    pLocal(i1,i2,i3,pc+2) = eps2*( chi2r[iv]*Etzr - chi2i[iv]*Etzi );
        	  }
      	}
                    

      	if( computeMagneticField )
      	{
  	  // Here are the statements to eval the time-dependent solution: 
// File generated by Dropbox/DARPA/RPI/adePapers/adegdmi/dispersivePlaneWaveInterface3d.maple
real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t530,t531,t532,t533,t534,t535,t536,t537,t538,t539,t540,t541,t542,t543,t544,t545,t546,t547,t548,t549,t550,t551,t552,t553,t554,t555,t556,t557,t558,t559,t560,t561,t562,t563,t564,t565,t566,t567,t568,t569,t570,t571,t572,t573,t574,t575,t576,t577,t578,t579,t580,t581,t582,t583,t584,t585,t586,t587,t588,t589,t590,t591,t592,t593,t594,t595,t596,t597,t598,t599,t600,t601,t602,t603,t604,t605,t606,t607,t608,t609,t610,t611,t612,t613,t614,t615,t616,t617,t618,t619,t620,t621,t622,t623,t624,t625,t626,t627,t628,t629,t630,t631,t632,t633,t634,t635,t636,t637,t638,t639,t640,t641,t642,t643,t644,t645,t646,t647,t648,t649,t650,t651,t652,t653,t654,t655,t656,t657,t658,t659,t660,t661,t662,t663,t664,t665,t666,t667,t668,t669,t670,t671,t672,t673,t674,t675,t676,t677,t678,t679,t680,t681,t682,t683,t684,t685,t686,t687,t688,t689,t690,t691,t692,t693,t694,t695,t696,t697,t698,t699,t700,t701,t702,t703,t704,t705,t706,t707,t708,t709,t710,t711,t712,t713,t714,t715,t716,t717,t718,t719,t720,t721,t722,t723,t724,t725,t726,t727,t728,t729,t730,t731,t732,t733,t734,t735,t736,t737,t738,t739,t740,t741,t742,t743,t744,t745,t746,t747,t748,t749,t750,t751,t752,t753,t754,t755,t756,t757,t758,t759,t760,t761,t762,t763,t764,t765,t766,t767,t768,t769,t770,t771,t772,t773,t774,t775,t776,t777,t778,t779,t780,t781,t782,t783,t784,t785,t786,t787,t788,t789,t790,t791,t792,t793,t794,t795,t796,t797,t798,t799,t800,t801,t802,t803,t804,t805,t806,t807,t808,t809,t810,t811,t812,t813,t814,t815,t816,t817,t818,t819,t820,t821,t822,t823,t824,t825,t826,t827,t828,t829,t830,t831,t832,t833,t834,t835,t836,t837,t838,t839,t840,t841,t842,t843,t844,t845,t846,t847,t848,t849,t850,t851,t852,t853,t854,t855,t856,t857,t858,t859,t860,t861,t862,t863,t864,t865,t866,t867,t868,t869,t870,t871,t872,t873,t874,t875,t876,t877,t878,t879,t880,t881,t882,t883,t884,t885,t886,t887,t888,t889,t890,t891,t892,t893,t894,t895,t896,t897,t898,t899,t900,t901,t902,t903,t904,t905,t906,t907,t908,t909,t910,t911,t912,t913,t914,t915,t916,t917,t918,t919,t920,t921,t922,t923,t924,t925,t926,t927,t928,t929,t930,t931,t932,t933,t934,t935,t936,t937,t938,t939,t940,t941,t942,t943,t944,t945,t946,t947,t948,t949,t950,t951,t952,t953,t954,t955,t956,t957,t958,t959,t960,t961,t962,t963,t964,t965,t966,t967,t968,t969,t970,t971,t972,t973,t974,t975,t976,t977,t978,t979,t980,t981,t982,t983,t984,t985,t986,t987,t988,t989,t990,t991,t992,t993,t994,t995,t996,t997,t998,t999;
// -------------------------------------------------------------------------
// Need: 
//    s=sr + I*si : complex frequency
//    mu1,mu2 : mu on left and right                                         
//    [ax,ay,az] = amplitude vector of incident wave, ax=axr+I*axi etc.
//    [arx,ary,arz] = amplitude vector of reflected wave, arx=arxr+I*arxi etc.
//    [atx,aty,atz] = amplitude vector of transmitted wave, atx=atxr+I*atxi etc.
//    [kx,ky,kz] = incident wave vector, kx=kxr+I*kxi etc. 
//    [krx,kry,krz] = reflected wave vector, krx=krxr+I*krxi etc. 
//    [ktx,kty,ktz] = reflected wave vector, ktx=ktxr+I*ktxi etc. 
// -------------------------------------------------------------------------
// Evaluated:                                                               
//    [Hxr,Hxi] ,[Hyr,Hyi], [Hzr,Hzi]        : left state                   
//    [Htxr,Htxi] ,[Htyr,Htyi], [Htzr,Htzi]  : right state                  
// -------------------------------------------------------------------------
real Hxr,Hyr,Hzr, Hxi,Hyi,Hzi;
real Htxr,Htyr,Htzr, Htxi,Htyi,Htzi;
t3 = exp(t * sr);
t4 = 0.1e1 / mu1 * t3;
t5 = t * si;
t6 = sin(t5);
t11 = exp(-x * krxi - y * kryi - z * krzi);
t12 = t6 * t11;
t16 = x * krxr + y * kryr + z * krzr;
t17 = sin(t16);
t18 = t12 * t17;
t19 = aryi * krzi;
t20 = t19 * si;
t22 = aryi * krzr;
t23 = t22 * sr;
t25 = aryr * krzi;
t26 = t25 * sr;
t28 = aryr * krzr;
t29 = t28 * si;
t31 = arzi * kryi;
t32 = t31 * si;
t34 = arzi * kryr;
t35 = t34 * sr;
t37 = arzr * kryi;
t38 = t37 * sr;
t40 = arzr * kryr;
t41 = t40 * si;
t43 = cos(t16);
t44 = t12 * t43;
t45 = t19 * sr;
t47 = t22 * si;
t49 = t25 * si;
t51 = t28 * sr;
t53 = t31 * sr;
t55 = t34 * si;
t57 = t37 * si;
t59 = t40 * sr;
t61 = t18 * t20 + t18 * t23 + t18 * t26 - t18 * t29 - t18 * t32 - t18 * t35 - t18 * t38 + t18 * t41 + t44 * t45 - t44 * t47 - t44 * t49 - t44 * t51 - t44 * t53 + t44 * t55 + t44 * t57 + t44 * t59;
t66 = exp(-x * kxi - y * kyi - z * kzi);
t67 = t6 * t66;
t71 = x * kxr + y * kyr + z * kzr;
t72 = sin(t71);
t73 = t67 * t72;
t74 = ayi * kzi;
t75 = t74 * si;
t77 = ayi * kzr;
t78 = t77 * sr;
t80 = ayr * kzi;
t81 = t80 * sr;
t83 = ayr * kzr;
t84 = t83 * si;
t86 = azi * kyi;
t87 = t86 * si;
t89 = azi * kyr;
t90 = t89 * sr;
t92 = azr * kyi;
t93 = t92 * sr;
t95 = azr * kyr;
t96 = t95 * si;
t98 = cos(t71);
t99 = t67 * t98;
t100 = t74 * sr;
t102 = t77 * si;
t104 = t80 * si;
t106 = t83 * sr;
t108 = t86 * sr;
t110 = t89 * si;
t112 = t92 * si;
t114 = t95 * sr;
t116 = t99 * t100 - t99 * t102 - t99 * t104 - t99 * t106 - t99 * t108 + t99 * t110 + t99 * t112 + t99 * t114 + t73 * t75 + t73 * t78 + t73 * t81 - t73 * t84 - t73 * t87 - t73 * t90 - t73 * t93 + t73 * t96;
t118 = cos(t5);
t119 = t118 * t11;
t120 = t119 * t17;
t129 = t119 * t43;
t138 = t120 * t45 - t120 * t47 - t120 * t49 - t120 * t51 - t120 * t53 + t120 * t55 + t120 * t57 + t120 * t59 - t129 * t20 - t129 * t23 - t129 * t26 + t129 * t29 + t129 * t32 + t129 * t35 + t129 * t38 - t129 * t41;
t139 = t118 * t66;
t140 = t139 * t72;
t149 = t139 * t98;
t158 = t140 * t100 - t140 * t102 - t140 * t104 - t140 * t106 - t140 * t108 + t140 * t110 + t140 * t112 + t140 * t114 - t149 * t75 - t149 * t78 - t149 * t81 + t149 * t84 + t149 * t87 + t149 * t90 + t149 * t93 - t149 * t96;
t161 = pow(si, 0.2e1);
t162 = pow(sr, 0.2e1);
t164 = 0.1e1 / (t161 + t162);
Hxr = t4 * (t61 + t116 + t138 + t158) * t164;
t166 = azr * kxi;
t167 = t166 * si;
t169 = azr * kxr;
t170 = t169 * sr;
t172 = axi * kzi;
t173 = t172 * si;
t175 = axi * kzr;
t176 = t175 * sr;
t178 = axr * kzi;
t179 = t178 * sr;
t181 = axr * kzr;
t182 = t181 * si;
t184 = azi * kxi;
t185 = t184 * si;
t187 = azi * kxr;
t188 = t187 * sr;
t190 = t166 * sr;
t192 = t169 * si;
t194 = arxi * krzi;
t195 = t194 * si;
t197 = arxi * krzr;
t198 = t197 * sr;
t200 = arxr * krzi;
t201 = t200 * sr;
t203 = arxr * krzr;
t204 = t203 * si;
t206 = arzi * krxi;
t207 = t206 * si;
t209 = arzi * krxr;
t210 = t209 * sr;
t212 = t140 * t167 + t140 * t170 - t149 * t173 - t149 * t176 - t149 * t179 + t149 * t182 + t149 * t185 + t149 * t188 + t149 * t190 - t149 * t192 + t18 * t195 + t18 * t198 + t18 * t201 - t18 * t204 - t18 * t207 - t18 * t210;
t213 = arzr * krxi;
t214 = t213 * sr;
t216 = arzr * krxr;
t217 = t216 * si;
t219 = t213 * si;
t221 = t216 * sr;
t231 = t172 * sr;
t233 = t175 * si;
t235 = t178 * si;
t237 = t181 * sr;
t239 = t73 * t173 + t73 * t176 + t73 * t179 - t18 * t214 + t18 * t217 - t73 * t182 - t73 * t185 - t73 * t188 - t73 * t190 + t73 * t192 + t44 * t219 + t44 * t221 + t99 * t231 - t99 * t233 - t99 * t235 - t99 * t237;
t241 = t184 * sr;
t243 = t187 * si;
t247 = t197 * si;
t249 = t200 * si;
t251 = t203 * sr;
t253 = t206 * sr;
t255 = t209 * si;
t264 = t120 * t219 + t120 * t221 - t120 * t247 - t120 * t249 - t120 * t251 - t120 * t253 + t120 * t255 - t129 * t195 - t129 * t198 - t129 * t201 + t129 * t204 + t129 * t207 + t99 * t167 + t99 * t170 - t99 * t241 + t99 * t243;
t274 = t194 * sr;
t282 = t120 * t274 + t129 * t210 + t129 * t214 - t129 * t217 + t140 * t231 - t140 * t233 - t140 * t235 - t140 * t237 - t140 * t241 + t140 * t243 - t44 * t247 - t44 * t249 - t44 * t251 - t44 * t253 + t44 * t255 + t44 * t274;
Hyr = -t4 * (t212 + t239 + t264 + t282) * t164;
t287 = aryi * krxr;
t288 = t287 * sr;
t290 = aryr * krxi;
t291 = t290 * sr;
t293 = aryr * krxr;
t294 = t293 * si;
t296 = axi * kyi;
t297 = t296 * sr;
t299 = ayr * kxi;
t300 = t299 * sr;
t302 = ayr * kxr;
t303 = t302 * si;
t305 = arxi * kryi;
t306 = t305 * si;
t308 = arxi * kryr;
t309 = t308 * sr;
t311 = arxr * kryi;
t312 = t311 * sr;
t314 = arxr * kryr;
t315 = t314 * si;
t317 = aryi * krxi;
t318 = t317 * si;
t323 = t305 * sr;
t325 = t308 * si;
t327 = t129 * t288 + t129 * t291 - t129 * t294 + t140 * t297 + t149 * t300 - t149 * t303 - t18 * t288 - t18 * t291 + t18 * t294 + t18 * t306 + t18 * t309 + t18 * t312 - t18 * t315 - t18 * t318 + t44 * t323 - t44 * t325;
t328 = t311 * si;
t330 = t314 * sr;
t332 = t317 * sr;
t334 = t287 * si;
t336 = t290 * si;
t338 = t293 * sr;
t340 = t296 * si;
t342 = axi * kyr;
t343 = t342 * sr;
t345 = axr * kyi;
t346 = t345 * sr;
t348 = axr * kyr;
t349 = t348 * si;
t351 = ayi * kxi;
t352 = t351 * si;
t354 = ayi * kxr;
t355 = t354 * sr;
t361 = t120 * t323 - t120 * t325 - t73 * t300 + t73 * t303 - t44 * t328 - t44 * t330 - t44 * t332 + t44 * t334 + t44 * t336 + t44 * t338 + t73 * t340 + t73 * t343 + t73 * t346 - t73 * t349 - t73 * t352 - t73 * t355;
t375 = t342 * si;
t377 = t345 * si;
t379 = t348 * sr;
t381 = t351 * sr;
t383 = -t120 * t328 - t120 * t330 - t120 * t332 + t120 * t334 + t120 * t336 + t120 * t338 - t129 * t306 - t129 * t309 - t129 * t312 + t129 * t315 + t129 * t318 + t99 * t297 - t99 * t375 - t99 * t377 - t99 * t379 - t99 * t381;
t384 = t354 * si;
t386 = t299 * si;
t388 = t302 * sr;
t403 = -t140 * t375 - t140 * t377 - t140 * t379 - t140 * t381 + t140 * t384 + t140 * t386 + t140 * t388 - t149 * t340 - t149 * t343 - t149 * t346 + t149 * t349 + t149 * t352 + t149 * t355 + t99 * t384 + t99 * t386 + t99 * t388;
Hzr = t4 * (t327 + t361 + t383 + t403) * t164;
t423 = -t120 * t32 - t120 * t35 - t120 * t38 + t120 * t41 + t129 * t45 - t129 * t47 - t129 * t49 - t129 * t51 - t129 * t53 + t129 * t55 + t129 * t57 + t129 * t59 + t140 * t75 + t140 * t78 + t140 * t81 - t140 * t84;
t440 = t149 * t100 - t149 * t102 - t149 * t104 - t149 * t106 - t149 * t108 + t149 * t110 + t149 * t112 + t149 * t114 - t140 * t87 - t140 * t90 - t140 * t93 + t140 * t96 - t18 * t45 + t18 * t47 + t18 * t49 + t18 * t51;
t458 = -t73 * t100 + t73 * t102 + t73 * t104 + t73 * t106 + t18 * t53 - t18 * t55 - t18 * t57 - t18 * t59 + t44 * t20 + t44 * t23 + t44 * t26 - t44 * t29 - t44 * t32 - t44 * t35 - t44 * t38 + t44 * t41;
t475 = t73 * t108 - t73 * t110 - t73 * t112 - t73 * t114 + t120 * t20 + t120 * t23 + t120 * t26 - t120 * t29 + t99 * t75 + t99 * t78 + t99 * t81 - t99 * t84 - t99 * t87 - t99 * t90 - t99 * t93 + t99 * t96;
Hxi = -t4 * (t423 + t440 + t458 + t475) * t164;
t496 = t120 * t195 + t120 * t198 + t120 * t201 - t120 * t204 - t120 * t207 - t120 * t210 - t120 * t214 + t120 * t217 - t129 * t247 + t129 * t274 - t99 * t182 - t44 * t207 - t44 * t210 - t44 * t214 + t44 * t217 - t73 * t231;
t513 = t129 * t219 + t129 * t221 - t129 * t249 - t129 * t251 - t129 * t253 + t129 * t255 + t140 * t173 + t140 * t176 + t140 * t179 - t140 * t182 - t140 * t185 - t140 * t188 - t99 * t185 - t99 * t188 - t99 * t190 + t99 * t192;
t531 = -t140 * t190 + t140 * t192 + t149 * t167 + t149 * t170 + t149 * t231 - t149 * t233 - t149 * t235 - t149 * t237 - t149 * t241 + t149 * t243 + t18 * t247 + t18 * t249 + t18 * t251 + t18 * t253 - t18 * t255 - t18 * t274;
t548 = -t73 * t167 - t73 * t170 + t99 * t173 + t99 * t176 + t99 * t179 - t18 * t219 - t18 * t221 + t44 * t195 + t44 * t198 + t44 * t201 - t44 * t204 + t73 * t233 + t73 * t235 + t73 * t237 + t73 * t241 - t73 * t243;
Hyi = t4 * (t496 + t513 + t531 + t548) * t164;
t568 = t120 * t306 - t149 * t375 - t99 * t300 + t99 * t303 + t99 * t340 + t99 * t343 + t99 * t346 - t99 * t349 - t99 * t352 - t99 * t355 + t73 * t377 + t73 * t379 + t73 * t381 - t73 * t384 - t73 * t386 - t73 * t388;
t585 = -t120 * t288 - t120 * t291 + t120 * t294 + t120 * t309 + t120 * t312 - t120 * t315 - t120 * t318 + t129 * t323 - t129 * t325 - t129 * t328 - t129 * t330 - t129 * t332 + t129 * t334 + t129 * t336 + t129 * t338 + t140 * t340;
t603 = -t140 * t300 + t140 * t303 + t140 * t343 + t140 * t346 - t140 * t349 - t140 * t352 - t140 * t355 + t149 * t297 - t149 * t377 - t149 * t379 - t149 * t381 + t149 * t384 + t149 * t386 + t149 * t388 - t18 * t323 + t18 * t325;
t620 = t18 * t328 + t18 * t330 + t18 * t332 - t18 * t334 - t18 * t336 - t18 * t338 - t44 * t288 - t44 * t291 + t44 * t294 - t73 * t297 + t44 * t306 + t44 * t309 + t44 * t312 - t44 * t315 - t44 * t318 + t73 * t375;
Hzi = -t4 * (t568 + t585 + t603 + t620) * t164;
t630 = exp(-x * ktxi - y * ktyi - z * ktzi);
t631 = 0.1e1 / mu2 * t630;
t635 = x * ktxr + y * ktyr + z * ktzr;
t636 = cos(t635);
t637 = t636 * t6;
t638 = atzi * ktyi;
t639 = t638 * sr;
t641 = atzi * ktyr;
t642 = t641 * si;
t644 = atzr * ktyi;
t645 = t644 * si;
t647 = atzr * ktyr;
t648 = t647 * sr;
t650 = sin(t635);
t651 = t650 * t118;
t652 = atyi * ktzi;
t653 = t652 * sr;
t655 = atyi * ktzr;
t656 = t655 * si;
t658 = atyr * ktzi;
t659 = t658 * si;
t661 = atyr * ktzr;
t662 = t661 * sr;
t668 = t650 * t6;
t669 = t652 * si;
t671 = t655 * sr;
t673 = t658 * sr;
t675 = t661 * si;
t677 = t637 * t639 - t637 * t642 - t637 * t645 - t637 * t648 + t651 * t639 - t651 * t642 - t651 * t645 - t651 * t648 - t651 * t653 + t651 * t656 + t651 * t659 + t651 * t662 - t668 * t669 - t668 * t671 - t668 * t673 + t668 * t675;
t678 = t638 * si;
t680 = t641 * sr;
t682 = t644 * sr;
t684 = t647 * si;
t686 = t636 * t118;
t699 = -t637 * t653 + t637 * t656 + t637 * t659 + t637 * t662 + t668 * t678 + t668 * t680 + t668 * t682 - t668 * t684 + t686 * t669 + t686 * t671 + t686 * t673 - t686 * t675 - t686 * t678 - t686 * t680 - t686 * t682 + t686 * t684;
Htxr = -t631 * t3 * (t677 + t699) * t164;
t704 = atxr * ktzi;
t705 = t704 * si;
t707 = atxr * ktzr;
t708 = t707 * sr;
t710 = atzi * ktxi;
t711 = t710 * sr;
t713 = atzi * ktxr;
t714 = t713 * si;
t716 = atzr * ktxi;
t717 = t716 * si;
t719 = atzr * ktxr;
t720 = t719 * sr;
t722 = atxi * ktzi;
t723 = t722 * si;
t725 = atxi * ktzr;
t726 = t725 * sr;
t728 = t704 * sr;
t730 = t707 * si;
t732 = t710 * si;
t734 = t713 * sr;
t736 = t716 * sr;
t738 = t719 * si;
t742 = t651 * t705 + t651 * t708 + t651 * t711 - t651 * t714 - t651 * t717 - t651 * t720 - t668 * t723 - t668 * t726 - t668 * t728 + t668 * t730 + t668 * t732 + t668 * t734 + t668 * t736 - t668 * t738 + t686 * t723 + t686 * t726;
t749 = t722 * sr;
t751 = t725 * si;
t761 = t637 * t705 + t637 * t708 + t637 * t711 - t637 * t714 - t637 * t717 - t637 * t720 - t637 * t749 + t637 * t751 - t651 * t749 + t651 * t751 + t686 * t728 - t686 * t730 - t686 * t732 - t686 * t734 - t686 * t736 + t686 * t738;
Htyr = t631 * t3 * (t742 + t761) * t164;
t765 = atxr * ktyr;
t766 = t765 * si;
t768 = atyi * ktxi;
t769 = t768 * si;
t771 = atyi * ktxr;
t772 = t771 * sr;
t774 = atyr * ktxi;
t775 = t774 * sr;
t777 = atyr * ktxr;
t778 = t777 * si;
t780 = atxi * ktyi;
t781 = t780 * sr;
t783 = atxi * ktyr;
t784 = t783 * si;
t786 = atxr * ktyi;
t787 = t786 * si;
t789 = t765 * sr;
t791 = t768 * sr;
t793 = t771 * si;
t795 = t774 * si;
t797 = t777 * sr;
t802 = -t637 * t781 + t637 * t784 + t637 * t787 + t637 * t789 + t637 * t791 - t637 * t793 - t637 * t795 - t637 * t797 - t651 * t781 + t651 * t784 + t651 * t787 - t686 * t766 - t686 * t769 - t686 * t772 - t686 * t775 + t686 * t778;
t808 = t780 * si;
t810 = t783 * sr;
t812 = t786 * sr;
t822 = t651 * t789 + t651 * t791 - t651 * t793 - t651 * t795 - t651 * t797 + t668 * t766 + t668 * t769 + t668 * t772 + t668 * t775 - t668 * t778 - t668 * t808 - t668 * t810 - t668 * t812 + t686 * t808 + t686 * t810 + t686 * t812;
Htzr = -t631 * t3 * (t802 + t822) * t164;
t843 = t668 * t639 - t668 * t642 - t668 * t645 - t668 * t648 - t651 * t678 - t651 * t680 - t651 * t682 + t651 * t684 - t668 * t653 + t686 * t653 + t668 * t656 - t686 * t656 + t668 * t659 - t686 * t659 + t668 * t662 - t686 * t662;
t860 = t637 * t669 + t637 * t671 + t637 * t673 - t637 * t675 - t637 * t678 - t637 * t680 - t637 * t682 + t637 * t684 - t686 * t639 + t686 * t642 + t686 * t645 + t686 * t648 + t651 * t669 + t651 * t671 + t651 * t673 - t651 * t675;
Htxi = -t631 * t3 * (t843 + t860) * t164;
t881 = t637 * t723 + t637 * t726 + t637 * t728 - t637 * t730 - t637 * t732 - t637 * t734 - t637 * t736 + t637 * t738 - t686 * t705 - t686 * t708 - t686 * t711 + t686 * t714 + t686 * t717 + t686 * t720 + t686 * t749 - t686 * t751;
t898 = t651 * t723 + t651 * t726 + t651 * t728 - t651 * t730 - t651 * t732 - t651 * t734 - t651 * t736 + t651 * t738 + t668 * t705 + t668 * t708 + t668 * t711 - t668 * t714 - t668 * t717 - t668 * t720 - t668 * t749 + t668 * t751;
Htyi = t631 * t3 * (t881 + t898) * t164;
t918 = -t637 * t766 - t637 * t769 - t637 * t772 - t637 * t775 + t637 * t778 + t637 * t808 + t637 * t810 + t637 * t812 + t651 * t808 + t651 * t810 + t651 * t812 - t686 * t789 - t686 * t791 + t686 * t793 + t686 * t795 + t686 * t797;
t935 = -t651 * t766 - t651 * t769 - t651 * t772 - t651 * t775 + t651 * t778 - t668 * t781 + t668 * t784 + t668 * t787 + t668 * t789 + t668 * t791 - t668 * t793 - t668 * t795 - t668 * t797 + t686 * t781 - t686 * t784 - t686 * t787;
Htzi = -t631 * t3 * (t918 + t935) * t164;

        	  if( myDomain==0 )
        	  {
          	    uLocal(i1,i2,i3,hx) = Hxr;
          	    uLocal(i1,i2,i3,hy) = Hyr;
          	    uLocal(i1,i2,i3,hz) = Hzr;
        	  }
        	  else
        	  {
          	    uLocal(i1,i2,i3,hx) = Htxr;
          	    uLocal(i1,i2,i3,hy) = Htyr;
          	    uLocal(i1,i2,i3,hz) = Htzr;
        	  }
      	}
      	

            } // end FOR
            

                
        }
    }
    else
    {
        OV_ABORT("ERROR: numberOfTimeDerivatives != 0 ");
        
    } // end if number of time derivatives 
        
    

    return 0;

}













// Include complex down here to minimize name conflicts
#include <complex>

typedef ::real LocalReal;
typedef ::real OV_real;



// ===============================================================================
/// \brief Check the solution.
// ===============================================================================
int PlaneInterfaceExactSolution::
PlaneInterfaceExactSolution::check()
{

  // ------------- CHECK THAT THE EQUATIONS ARE SATISFIED AT POINTS INSIDE AND OUTSIDE ------

    printF("------------ PlaneInterfaceExactSolution::check: CHECK THE EQUATIONS ------------\n\n");

    LocalReal maxErr=0.;

    const LocalReal & a = dbase.get<LocalReal>("amp");
    LocalReal (&normalPlaneMaterialInterface)[3] = dbase.get<LocalReal[3]>("normalPlaneMaterialInterface");
    LocalReal (&x0PlaneMaterialInterface)[3] = dbase.get<LocalReal[3]>("x0PlaneMaterialInterface");

    LocalReal *nv = normalPlaneMaterialInterface;
    LocalReal *x0 = x0PlaneMaterialInterface;

    LocalReal (&s)[2] = dbase.get<LocalReal[2]>("s");

    LocalReal & sr = s[0], &si = s[1];
    LocalReal (&av)[3] = dbase.get<LocalReal[3]>("av");

    LocalReal (&kv)[3] = dbase.get<LocalReal[3]>("kv");

    LocalReal & eps1 = dbase.get<LocalReal>("eps1");
    LocalReal & eps2 = dbase.get<LocalReal>("eps2");
    LocalReal & mu1 = dbase.get<LocalReal>("mu1");
    LocalReal & mu2 = dbase.get<LocalReal>("mu2");

    LocalReal & eps1Hatr = dbase.get<LocalReal>("eps1Hatr");
    LocalReal & eps1Hati = dbase.get<LocalReal>("eps1Hati");
    
    LocalReal & eps2Hatr = dbase.get<LocalReal>("eps2Hatr");
    LocalReal & eps2Hati = dbase.get<LocalReal>("eps2Hati");


    std::complex<LocalReal> ss(sr,si); // s = -i omega
    std::complex<LocalReal>  eps1Hat, mu1Hat, eps2Hat, mu2Hat, beta1, beta2;

    eps1Hat = eps1Hatr + 1i*eps1Hati;  // complex epsHat
    mu1Hat  = mu1;                     // complex muHat

    eps2Hat = eps2Hatr + 1i*eps2Hati;  // complex epsHat
    mu2Hat  = mu2;                     // complex muHat

    beta1 = sqrt( -ss*ss*eps1Hat*mu1Hat );
    beta2 = sqrt( -ss*ss*eps2Hat*mu2Hat );


    LocalReal deps1  = 10.*pow(REAL_EPSILON,1./2.);   // delta to compute 1st derivatives by differences 
    LocalReal deps2  = 10.*pow(REAL_EPSILON,1./3.);   // delta compute 2nd derivatives by differences 

    LocalReal xx0,yy0,zz0, x[3], xm[3], xp[3], Ev[6], Hv[6];
    std::complex<LocalReal> E[3],Em[3],Ep[3], Ex[3],Ey[3], Ez[3], Exx[3],Eyy[3], Ezz[3], DeltaE[3];
    std::complex<LocalReal> H[3],Hm[3],Hp[3], Hx[3],Hy[3], Hz[3], Hxx[3],Hyy[3], Hzz[3], DeltaH[3];

    std::complex<LocalReal>  epsHat, muHat, beta;

    for( int ip=0; ip<2; ip++ )
    {
        const LocalReal dist=.2;
        if( ip==0 )
        {  // Domain 1 : left 
            for( int i=0; i<3; i++ ){ x[i]=x0[i] -nv[i]*dist; }; 
              beta=beta1; epsHat=eps1Hat; muHat=mu1Hat; 
        }
        else
        { // Domain 2 : right 
            for( int i=0; i<3; i++ ){ x[i]=x0[i] + nv[i]*dist; }; 
            beta=beta2; epsHat=eps2Hat; muHat=mu2Hat; 
        }
        
    // To do: Make a macro to eval and return complex values 

// ----------------------------------------------------------------
// Macro: Compute 2nd derivatives of E and H by differences 
// ----------------------------------------------------------------

// ----------------------------------------------------------------
// Macro: Compute 1st derivatives of E and H by differences 
// ----------------------------------------------------------------

        
    // macro calls 
                eval( x, Ev,Hv );
                for( int i=0; i<3; i++ ) { E[i]= Ev[i]+1i*Ev[i+3]; H[i]= Hv[i]+1i*Hv[i+3]; } //
                xm[0]=x[0]-deps2; xm[1]=x[1]-0; xm[2]=x[2]-0;
                eval( xm, Ev,Hv );
                for( int i=0; i<3; i++ ) { Em[i]=  Ev[i]+1i*Ev[i+3]; Hm[i]= Hv[i]+1i*Hv[i+3]; } //
                xp[0]=x[0]+deps2; xp[1]=x[1]+0; xp[2]=x[2]+0;
                eval( xp, Ev,Hv );
                for( int i=0; i<3; i++ ) { Ep[i]=  Ev[i]+1i*Ev[i+3]; Hp[i]=  Hv[i]+1i*Hv[i+3]; } //
                for( int i=0; i<3; i++ ) { Exx[i]=  (Ep[i]-2.*E[i]+Em[i])/(deps2*deps2); Hxx[i]=  (Hp[i]-2.*H[i]+Hm[i])/(deps2*deps2);} // 
                eval( x, Ev,Hv );
                for( int i=0; i<3; i++ ) { E[i]= Ev[i]+1i*Ev[i+3]; H[i]= Hv[i]+1i*Hv[i+3]; } //
                xm[0]=x[0]-0; xm[1]=x[1]-deps2; xm[2]=x[2]-0;
                eval( xm, Ev,Hv );
                for( int i=0; i<3; i++ ) { Em[i]=  Ev[i]+1i*Ev[i+3]; Hm[i]= Hv[i]+1i*Hv[i+3]; } //
                xp[0]=x[0]+0; xp[1]=x[1]+deps2; xp[2]=x[2]+0;
                eval( xp, Ev,Hv );
                for( int i=0; i<3; i++ ) { Ep[i]=  Ev[i]+1i*Ev[i+3]; Hp[i]=  Hv[i]+1i*Hv[i+3]; } //
                for( int i=0; i<3; i++ ) { Eyy[i]=  (Ep[i]-2.*E[i]+Em[i])/(deps2*deps2); Hyy[i]=  (Hp[i]-2.*H[i]+Hm[i])/(deps2*deps2);} // 
                eval( x, Ev,Hv );
                for( int i=0; i<3; i++ ) { E[i]= Ev[i]+1i*Ev[i+3]; H[i]= Hv[i]+1i*Hv[i+3]; } //
                xm[0]=x[0]-0; xm[1]=x[1]-0; xm[2]=x[2]-deps2;
                eval( xm, Ev,Hv );
                for( int i=0; i<3; i++ ) { Em[i]=  Ev[i]+1i*Ev[i+3]; Hm[i]= Hv[i]+1i*Hv[i+3]; } //
                xp[0]=x[0]+0; xp[1]=x[1]+0; xp[2]=x[2]+deps2;
                eval( xp, Ev,Hv );
                for( int i=0; i<3; i++ ) { Ep[i]=  Ev[i]+1i*Ev[i+3]; Hp[i]=  Hv[i]+1i*Hv[i+3]; } //
                for( int i=0; i<3; i++ ) { Ezz[i]=  (Ep[i]-2.*E[i]+Em[i])/(deps2*deps2); Hzz[i]=  (Hp[i]-2.*H[i]+Hm[i])/(deps2*deps2);} // 
          
    // macro calls
                xm[0]=x[0]-deps1; xm[1]=x[1]-0; xm[2]=x[2]-0;
                eval( xm, Ev,Hv );
                for( int i=0; i<3; i++ ) { Em[i]=  Ev[i]+1i*Ev[i+3]; Hm[i]= Hv[i]+1i*Hv[i+3]; } //
                xp[0]=x[0]+deps1; xp[1]=x[1]+0; xp[2]=x[2]+0;
                eval( xp, Ev,Hv );
                for( int i=0; i<3; i++ ) { Ep[i]=  Ev[i]+1i*Ev[i+3]; Hp[i]=  Hv[i]+1i*Hv[i+3]; } //
                for( int i=0; i<3; i++ ) { Ex[i]=  (Ep[i]-Em[i])/(2.*deps1); Hx[i]=  (Hp[i]-Hm[i])/(2.*deps1);} // 
                xm[0]=x[0]-0; xm[1]=x[1]-deps1; xm[2]=x[2]-0;
                eval( xm, Ev,Hv );
                for( int i=0; i<3; i++ ) { Em[i]=  Ev[i]+1i*Ev[i+3]; Hm[i]= Hv[i]+1i*Hv[i+3]; } //
                xp[0]=x[0]+0; xp[1]=x[1]+deps1; xp[2]=x[2]+0;
                eval( xp, Ev,Hv );
                for( int i=0; i<3; i++ ) { Ep[i]=  Ev[i]+1i*Ev[i+3]; Hp[i]=  Hv[i]+1i*Hv[i+3]; } //
                for( int i=0; i<3; i++ ) { Ey[i]=  (Ep[i]-Em[i])/(2.*deps1); Hy[i]=  (Hp[i]-Hm[i])/(2.*deps1);} // 
                xm[0]=x[0]-0; xm[1]=x[1]-0; xm[2]=x[2]-deps1;
                eval( xm, Ev,Hv );
                for( int i=0; i<3; i++ ) { Em[i]=  Ev[i]+1i*Ev[i+3]; Hm[i]= Hv[i]+1i*Hv[i+3]; } //
                xp[0]=x[0]+0; xp[1]=x[1]+0; xp[2]=x[2]+deps1;
                eval( xp, Ev,Hv );
                for( int i=0; i<3; i++ ) { Ep[i]=  Ev[i]+1i*Ev[i+3]; Hp[i]=  Hv[i]+1i*Hv[i+3]; } //
                for( int i=0; i<3; i++ ) { Ez[i]=  (Ep[i]-Em[i])/(2.*deps1); Hz[i]=  (Hp[i]-Hm[i])/(2.*deps1);} // 
          
        
    // Check E wave equation:   
    //     Delta Ev = s^2 epsHat muHat Ev = - m^2 k^2 Ev = - beta^2 Ev 
        for( int i=0; i<3; i++ ) { DeltaE[i]= Exx[i]+Eyy[i]+Ezz[i]; DeltaH[i]= Hxx[i]+Hyy[i]+Hzz[i]; } //

        LocalReal normDeltaE =  max(abs(DeltaE[0]),abs(DeltaE[1]),abs(DeltaE[2]));
        LocalReal errE = max( abs( DeltaE[0] + beta*beta*E[0] ),
                                                    abs( DeltaE[1] + beta*beta*E[1] ),
                                                    abs( DeltaE[2] + beta*beta*E[2] ) )/normDeltaE;
        
        maxErr = max(maxErr,errE);

        LocalReal normDeltaH =  max(abs(DeltaH[0]),abs(DeltaH[1]),abs(DeltaH[2]));
        LocalReal errH = max( abs( DeltaH[0] + beta*beta*H[0] ),
                                                    abs( DeltaH[1] + beta*beta*H[1] ),
                                                    abs( DeltaH[2] + beta*beta*H[2] ) )/normDeltaH;
        
        maxErr = max(maxErr,errH);


    //  Check E FOS equation:    eps*E_t = curl(H) 
    //  s epsHat E = curl( H )
        std::complex<LocalReal> curlE[3], curlH[3], EtHat[3], HtHat[3];
        LocalReal normCurlH, normCurlE, errEt, errHt;
        curlH[0] = Hy[2] - Hz[1];
        curlH[1] = Hz[0] - Hx[2];
        curlH[2] = Hx[1] - Hy[0];

        normCurlH = max(abs(curlH[0]),abs(curlH[1]),abs(curlH[2])); 
        for( int i=0; i<3; i++ ){ EtHat[i] = ss*epsHat*E[i] - curlH[i]; } // 
        errEt = max( abs(EtHat[0]),abs(EtHat[1]),abs(EtHat[2]) )/normCurlH;

        maxErr = max(maxErr,errEt);


    // Check H FOS equation:
    //     mu*H_t = -curl(E) 
    //  s muHat H = -curl( E )
        curlE[0] = Ey[2] - Ez[1];
        curlE[1] = Ez[0] - Ex[2];
        curlE[2] = Ex[1] - Ey[0];

        normCurlE = max(abs(curlE[0]),abs(curlE[1]),abs(curlE[2])); 
        for( int i=0; i<3; i++ ){ HtHat[i] = ss*muHat*H[i] + curlE[i]; } 
        errHt = max( abs(HtHat[0]),abs(HtHat[1]),abs(HtHat[2]) )/normCurlE;

        maxErr = max(maxErr,errHt);

    // printF("(1) maxErr=%g\n",maxErr);

        printF(" Point: (x,y,z)=(%g,%g,%g) beta%d=(%g,%g) s=(%g,%g) epsHat%d=(%g,%g)\n",
                      x[0],x[1],x[2],ip+1,std::real(beta),std::imag(beta),std::real(ss),std::imag(ss),
                        	                  ip+1,std::real(epsHat),std::imag(epsHat));
        printF(" SOS(E):  |Delta(E)|=%9.2e, | Delta(E) + beta^2*E|/|Delta(E)| =%9.2e\n",normDeltaE,errE);
        printF(" SOS(H):  |Delta(H)|=%9.2e, | Delta(H) + beta^2*H|/|Delta(H)| =%9.2e\n",normDeltaH,errH);
        printF(" FOS(E):  |curl(H) |=%9.2e, | s*eps*E - curl(H)  |/|curl(H) | =%9.2e\n",normCurlH,errEt);
        printF(" FOS(H):  |curl(E) |=%9.2e, | s*mu*H + curl(E)   |/|curl(E) | =%9.2e\n",normCurlE,errHt);

        if( false )
        {
            for( int i=0; i<3; i++ )
            {
                printF(" i=%d: s*eps*E = (%g,%g) curl(H)=(%g,%g) \n",i,
                              std::real(ss*epsHat*E[i]),std::imag(ss*epsHat*E[i]),
                              std::real(curlH[i]),std::imag(curlH[i]));

                printF(" i=%d: DeltaE = (%g,%g) beta^2*E=(%g,%g) \n",i,
                              std::real(DeltaE[i]),std::imag(DeltaE[i]),
                              std::real(beta*beta*E[i]),std::imag(beta*beta*E[i]));
                printF(" i=%d: DeltaH = (%g,%g) beta^2*H=(%g,%g) \n",i,
                              std::real(DeltaH[i]),std::imag(DeltaH[i]),
                              std::real(beta*beta*H[i]),std::imag(beta*beta*H[i]));
          
            }
        }
    // eval( xp, Evp,Hvp );

    } // end for ip


    if( true )
    {
    // --- CHECK INTERFACE JUMP CONDITIONS -------------
    // LocalReal theta,phi;
    // std::complex<LocalReal> Er,Etheta,Ephi,Hr,Htheta,Hphi;
        printF("JUMPS: nv=(%g,%g,%g) \n",nv[0],nv[1],nv[2]);
        
        LocalReal dist=1.e-9;

    // Evaluate left side of the interface 
        for( int i=0; i<3; i++ ){ x[i]=x0[i] -nv[i]*dist; }; 
        eval( x, Ev,Hv );
        for( int i=0; i<3; i++ ) { Em[i]=  Ev[i]+1i*Ev[i+3]; Hm[i]=  Hv[i]+1i*Hv[i+3]; } //

    // Evaluate right side
        for( int i=0; i<3; i++ ){ x[i]=x0[i] + nv[i]*dist; }; 
        eval( x, Ev,Hv );
        for( int i=0; i<3; i++ ) { Ep[i]=  Ev[i]+1i*Ev[i+3]; Hp[i]=  Hv[i]+1i*Hv[i+3]; } //

    // E = jump in E : [E] 
        for( int i=0; i<3; i++ ) { E[i]= Ep[i]-Em[i];  H[i]= Hp[i]-Hm[i]; } //


    // Compute n X E and n X H 
        std::complex<LocalReal> nCrossE[3], nCrossH[3], nDotE;
        nCrossE[0] = nv[1]*E[2] - nv[2]*E[1];
        nCrossE[1] = nv[2]*E[0] - nv[0]*E[2];
        nCrossE[2] = nv[0]*E[1] - nv[1]*E[0];

        nCrossH[0] = nv[1]*H[2] - nv[2]*H[1];
        nCrossH[1] = nv[2]*H[0] - nv[0]*H[2];
        nCrossH[2] = nv[0]*H[1] - nv[1]*H[0];

        nDotE = nv[0]*E[0] + nv[1]*E[1] + nv[2]*E[2];

        printF("JUMPS: abs([n.E]) = %10.2e\n",abs(nDotE));
        printF("JUMPS: [E] =(%10.2e,%10.2e,%10.2e)+ I*(%10.2e,%10.2e,%10.2e) \n",
         	   std::real(E[0]),std::real(E[1]),std::real(E[2]), imag(E[0]), imag(E[1]), imag(E[2]) );
        printF("JUMPS: abs([E])     =(%10.2e,%10.2e,%10.2e)  [H]    =(%10.2e,%10.2e,%10.2e) \n",
             	       abs(E[0]),abs(E[1]),abs(E[2]), abs(H[0]),abs(H[1]),abs(H[2]) );
        printF("JUMPS: abs([n X E]) =(%10.2e,%10.2e,%10.2e) abs([n X H]) =(%10.2e,%10.2e,%10.2e) \n",
         	   abs(nCrossE[0]),abs(nCrossE[1]),abs(nCrossE[2]), abs(nCrossH[0]),abs(nCrossH[1]),abs(nCrossH[2]) );
        maxErr = max(maxErr,abs(nCrossE[0]),abs(nCrossE[1]),abs(nCrossE[2]));
        maxErr = max(maxErr,abs(nCrossH[0]),abs(nCrossH[1]),abs(nCrossH[2]));



    }
    
    
    
    printF("\n Max-error in tests = %9.2e. TESTS %s.\n",maxErr,(maxErr<1.e-5 ? "PASSED" : "***FAILED***"));
    
    printF("\n------------ FINSHED CHECK EQUATIONS ------------\n\n");


    
    return 0;
}



/* -----
// =======================================================================================
///  \brief Compute some derived quanities 
/// \param beta0v[2], beta1v[2] (output) : phase parameters beta= m*k
// =======================================================================================
int PlaneInterfaceExactSolution::
getDispersiveParameters( LocalReal beta0v[2], LocalReal beta1v[2] )
{


  // const int & numDomains = dbase.get<int>("numDomains");

  // const LocalReal &a =dbase.get<LocalReal>("radius");

    LocalReal (&s)[2] = dbase.get<LocalReal[2]>("s");
    LocalReal (&k)[2] = dbase.get<LocalReal[2]>("k");

    LocalReal (&eps0)[2] = dbase.get<LocalReal[2]>("eps0");
    LocalReal (&eps1)[2] = dbase.get<LocalReal[2]>("eps1");
    LocalReal (&mu0)[2] = dbase.get<LocalReal[2]>("mu0");
    LocalReal (&mu1)[2] = dbase.get<LocalReal[2]>("mu1");
    
  // const LocalReal c1 = 1./sqrt(eps0[0]*mu0[0]);
  // const LocalReal c2 = 1./sqrt(eps1[0]*mu1[0]);
    
    LocalReal sr = s[0], si = s[1]; // -twoPi*kx*c1;
    
    std::complex<LocalReal> ss(sr,si); // s = -i omega
    std::complex<LocalReal> epsHat, epsHat0, epsHat1, muHat, muHat0, muHat1, cHat0, cHat1, eta0, eta1, m0,m1, kx;
    std::complex<LocalReal>  beta, beta0, beta1;
    
  // *** WE COULD RETURN OTHER QUANTITIES TOO 

    kx = k[0] + 1i*k[1];  // complex k, wave-number

    epsHat0 = eps0[0] + 1i*eps0[1];  // complex epsHat
    muHat0  = mu0[0]  + 1i*mu0[1];  // complex muHat
    
    epsHat1 = eps1[0] + 1i*eps1[1];  // complex epsHat
    muHat1  = mu1[0]  + 1i*mu1[1];  // complex muHat
    

    cHat0= 1./sqrt(epsHat0*muHat0); 
    m0   = 1./cHat0;              // index of refraction in the outer domain 
    eta0 = sqrt(muHat0/epsHat0);     // impedance
    
    cHat1= 1./sqrt(epsHat1*muHat1); 
    m1   = 1./cHat1;              // index of refraction in the inner domain 
    eta1 = sqrt(muHat1/epsHat1);    // impedance

                           				   
    beta0 = m0*kx;   // complex wave number (outer domain) 

    beta0v[0]=std::real(beta0);
    beta0v[1]=std::imag(beta0);

    beta1v[0]=std::real(beta0);
    beta1v[1]=std::imag(beta0);


    return 0;
}
---- */


// =====================================================================================
/// \brief Utility routine to do some complex arithemetic for the dispersive plane
///    wave material interface.
///
///  Compute kxp=(kxpr,kxpi)  given (kr,ki), and kyp=(kypr,kypi)
///  kxp^2 + kyp^2 = (kr+I*ki)^2 = (kr^2-ki^2) + 2*I*kr*ki 
///  kxp = kxpr + I*kpri = sqrt( (kr+I*ki)^2 - kyp^2 )
// =====================================================================================
void PlaneInterfaceExactSolution::
getTransmisionWaveNumber( const LocalReal & kr,  const LocalReal & ki, 
                                                    const LocalReal & kxr, const LocalReal & kxi, 
                                                    const LocalReal & kyr, const LocalReal & kyi, 
                                                    LocalReal & kxpr, LocalReal & kxpi, 
                                                    LocalReal & kypr, LocalReal & kypi )
{
  // No jump in tangential field: kyp=ky : 
    kypr=kyr;
    kypi=kyi;

  // std::complex<LocalReal> I(0.0,1.0); 
    std::complex<LocalReal> ky(kyr,kyi);
    std::complex<LocalReal> k(kr,ki);
    std::complex<LocalReal> kxp,kyp(kypr,kypi);

  // cout << "kyp=" << kyp << endl;

    kxp = std::sqrt( k*k - kyp*kyp );

    kxpr= std::real(kxp);
    kxpi= std::imag(kxp);
    
  // printF("--getTransmisionWaveNumber--- kx=(%g,%g) ky=(%g,%g) (kr,ki)=(%g,%g) kxp=(%g,%g) kyp=(%g,%g)\n",
  //            kxr,kxi,kyr,kyi,kr,ki,kxpr,kxpi,kypr,kypi);

}


// ---------------------------------------------------------------------------------------
// Check routine : 
//   Check the jump:
//       eps1Hat*khyat*(1-r) = eps2Hat*kyHatp*tau
// ---------------------------------------------------------------------------------------
void PlaneInterfaceExactSolution::
checkPlaneMaterialInterfaceJumps( 
                                                    const LocalReal & c1, const LocalReal & c2,
                                                    const LocalReal & eps1, const LocalReal & eps2,
                                                    const LocalReal & mu1, const LocalReal & mu2,

                                                    const LocalReal & sr, const LocalReal & si,
                                                    const LocalReal & rr, const LocalReal & ri, 
                                                    const LocalReal & taur, const LocalReal & taui, 

                                                    const LocalReal & eps1Hatr, const LocalReal & eps1Hati,
                                                    const LocalReal & eps2Hatr, const LocalReal & eps2Hati,

                                                    const LocalReal & psiSum1r, const LocalReal & psiSum1i,
                                                    const LocalReal & psiSum2r, const LocalReal & psiSum2i,
                                                    const LocalReal & kxr, const LocalReal & kxi,
                                                    const LocalReal & kyr, const LocalReal & kyi,
                                                    const LocalReal & kxpr, const LocalReal & kxpi,
                                                    const LocalReal & kypr, const LocalReal & kypi

                                                        )
{

  // std::complex<LocalReal> I(0.0,1.0); 
    std::complex<LocalReal> psiSum1(psiSum1r,psiSum1i);
    std::complex<LocalReal> psiSum2(psiSum2r,psiSum2i);

    std::complex<LocalReal> eps1Hat(eps1Hatr,eps1Hati);
    std::complex<LocalReal> eps2Hat(eps2Hatr,eps2Hati);
    std::complex<LocalReal> s(sr,si);
    std::complex<LocalReal> kx(kxr,kxi);
    std::complex<LocalReal> ky(kyr,kyi);
    std::complex<LocalReal> kxp(kxpr,kxpi);
    std::complex<LocalReal> kyp(kypr,kypi);
    std::complex<LocalReal> khx,khy, khpx,khpy;
    std::complex<LocalReal> r(rr,ri), tau(taur,taui);
    std::complex<LocalReal> dr1,dr2,jump;
    
    LocalReal kNorm = sqrt(kxr*kxr + kxi*kxi + kyr*kyr + kyi*kyi);
    khx= (kxr + 1i*kxi)/kNorm;
    khy= (kyr + 1i*kyi)/kNorm;
    
    LocalReal kpNorm = sqrt(kxpr*kxpr + kxpi*kxpi + kypr*kypr + kypi*kypi);
    khpx=(kxpr + 1i*kxpi)/kpNorm;
    khpy=(kypr + 1i*kypi)/kpNorm;
    
    printF("\n\n ** s=(%g,%g) kx=(%g,%g) ky=(%g,%g) c1=%g eps1=%g mu1=%g \n",sr,si,kxr,kxi,kyr,kyi,c1,eps1,mu1);
    
    LocalReal maxErr=0.;

    jump = khx*(1.+r) - tau*khpx;
    printF("khx*(1.+r) - tau*khpx                 =(%12.4e,%12.4e)\n",std::real(jump), std::imag(jump));
    maxErr=max(maxErr,abs(jump));
    
    jump =eps1Hat*khy*(1.-r) - eps2Hat*tau*khpy;
    printF(" [epsHat khy(1-r)- epsHat*tau*khy'    =(%12.4e,%12.4e)\n",std::real(jump), std::imag(jump));
    maxErr=max(maxErr,abs(jump));
    
  // dispersion relations: 
    dr1 = s*s + c1*c1*(kx*kx  +ky*ky  ) + s*s*psiSum1;
    dr2 = s*s + c2*c2*(kxp*kxp+kyp*kyp) + s*s*psiSum2;

    printF(" dispersion-relation1                 =(%12.4e,%12.4e)\n",std::real(dr1), std::imag(dr1));
    printF(" dispersion-relation2                 =(%12.4e,%12.4e)\n",std::real(dr2), std::imag(dr2));
    maxErr=max(maxErr,abs(dr1));
    maxErr=max(maxErr,abs(dr2));

    jump = (1.-r)*( ky*khy + kx*khx )/mu1 - tau*( kyp*khpy+kxp*khpx )/mu2;
    printF(" (1-r)*( kSq )/mu1 - tau*( kpSq )/mu2'=(%12.4e,%12.4e)\n",std::real(jump), std::imag(jump));
    maxErr=max(maxErr,abs(jump));

    jump = (1.-r)*( ky*ky + kx*kx )/(kNorm*mu1) - tau*( kyp*kyp+kxp*kxp )/(kpNorm*mu2);
    printF(" (1-r)*( kSq )/mu1-tau*( kpSq )/mu2'  =(%12.4e,%12.4e)\n",std::real(jump), std::imag(jump));
    maxErr=max(maxErr,abs(jump));

    jump = (1.-r)*( s*s*(1.+psiSum1)/(c1*c1) )/(kNorm*mu1) - tau*( s*s*(1.+psiSum2)/(c2*c2) )/(kpNorm*mu2);
    printF(" (1-r)*( kSq )/mu1-tau*( kpSq )/mu2'  =(%12.4e,%12.4e)\n",std::real(jump), std::imag(jump));
    maxErr=max(maxErr,abs(jump));

    jump = (1.-r)*( s*s*(1.+psiSum1)*(mu1*eps1) )/(kNorm*mu1) - tau*( s*s*(1.+psiSum2)*(mu2*eps2) )/(kpNorm*mu2);
    printF(" (1-r)*( kSq )/mu1-tau*( kpSq )/mu2'  =(%12.4e,%12.4e)\n",std::real(jump), std::imag(jump));
    maxErr=max(maxErr,abs(jump));

    jump = (1.-r)*( (1.+psiSum1)*(mu1*eps1) )/(kNorm*mu1) - tau*( (1.+psiSum2)*(mu2*eps2) )/(kpNorm*mu2);
    printF(" (1-r)*( kSq )/mu1-tau*( kpSq )/mu2'  =(%12.4e,%12.4e)\n",std::real(jump), std::imag(jump));
    maxErr=max(maxErr,abs(jump));

    jump = (1.-r)*( (1.+psiSum1)*(mu1*eps1) )*khy/(mu1) - tau*( (1.+psiSum2)*(mu2*eps2) )*khpy/(mu2);
    printF(" (1-r)*( kSq )/mu1-tau*( kpSq )/mu2'  =(%12.4e,%12.4e)\n",std::real(jump), std::imag(jump));
    maxErr=max(maxErr,abs(jump));

    jump = (1.-r)*( (1.+psiSum1)*(eps1) )*khy - tau*( (1.+psiSum2)*(eps2) )*khpy;
    printF(" (1-r)*( kSq )/mu1-tau*( kpSq )/mu2'  =(%12.4e,%12.4e)\n",std::real(jump), std::imag(jump));
    maxErr=max(maxErr,abs(jump));

    if( maxErr< REAL_EPSILON*1000. )
    {
        printF("ALL CHECKS PASSED!  maxErr=%9.3e\n",maxErr);
    }
    else
    {
        printF("**ERROR** SOME CHECKS FAILED!  maxErr=%9.3e\n",maxErr);
    }
    

}
