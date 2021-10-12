*
*  cgsm: compute Rayleigh surface waves
*
* Usage: (not all options implemented yet)
*   
*  cgsm [-noplot] rayleigh -g=<name> -tf=<tFinal> -tp=<tPlot> -diss=<> -dissOrder=<> -order=<2/4> -debug=<num> ...
*                       -bc=[d|sf|slip|dirichlet] -bg=<backGround> -pv=[nc|c|g|h] -godunovOrder=[1|2] -go=[run/halt/og]
* 
*  -bc : boundary conditions: -bc=dirichlet, -bc=sf (traction)
*        -bc=ellipseDeform : specified motion of the boundary
*  -ic : initial conditions, ic=gaussianPulse, ic=zero, ic=special
*  -diss : coeff of artificial diffusion 
*  -go : run, halt, og=open graphics
*  -pv : "pde-variation" : nc=non-conservative, c=conservative, g=godunov, h=hemp
*  -en : error norm, "max", "l1" or "l2"
* 
* Examples:
*   cgsm rayleigh -g=rayleighGridf16 -pv=nc -tp=.1 -ic=special -bc=dirichlet [Ok errors e-6 
*   cgsm rayleigh -g=rayleighGridf16 -pv=nc -tp=.1 -ic=special -bc=sf
*   cgsm rayleigh -g=rayleighGridf16 -pv=c -tp=.1 -ic=special -bc=dirichlet    [Ok errors e-6 
*     bottom BC=given: 
*      -->t=5.5000e-01 dt=2.3e-03 maxNorm errors:[3.9599e-05,2.8359e-05,], maxNorm(u):[4.22e-02,6.51e-02,]
* -- godunov:
*   cgsm rayleigh -g=rayleighGridf16 -pv=g -tp=.1 -tf=10. -ic=special -bc=dirichlet   [Ok
*   cgsm rayleigh -g=rayleighGridf16 -pv=g -tp=.1 -tf=10. -ic=special -bc=sf -k0=4
*     bottom BC=given: 
*    -->t=5.0000e-01 dt=1.6e-03 maxNorm errors:[7.70e-05,5.10e-05,3.01e-04,4.91e-05,4.91e-05,5.34e-05,8.98e-06,6.01e-06,]
* --- set default values for parameters ---
* 
$k0=1.; $a0=.5; $b0=-.5; 
$tFinal=10.; $tPlot=.05; $backGround="square"; $cfl=.9; $bc="sf"; $pv="nc"; $ts="me";  $show=" ";
$filter=0; $filterFreq=1; 
$ic="special";  $exponent=10.; $x0=.5; $y0=.5; $z0=.5; $specialOption="RayleighWave"; 
$noplot=""; $grid="rectangle80.ar10"; $mu=1.; $lambda=1.; $godunovOrder=2;
$debug = 0;  $tPlot=.1; $diss=0.; $dissOrder=2; $bc="sf"; $cons=1; $dsf=0.1;
$tz = "poly"; $degreex=2; $degreet=2; $fx=.5; $fy=$fx; $fz=$fx; $ft=$fx;
$order = 2; $go="run"; 
$en="max";
* 
* ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"degreex=i"=>\$degreex, "degreet=i"=>\$degreet,"diss=f"=>\$diss,\
 "dissOrder=i"=>\$dissOrder,"tp=f"=>\$tPlot, "tz=s"=>\$tz, "show=s"=>\$show,"order=i"=>\$order,"debug=i"=>\$debug, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bc=s"=>\$bc,"ic=s"=>\$ic,"go=s"=>\$go,"noplot=s"=>\$noplot,"ts=s"=>\$ts,\
  "mu=f"=>\$mu,"lambda=f"=>\$lambda,"dtMax=f"=>\$dtMax, "cons=i"=>\$cons,"x0=f"=>\$x0,"y0=f"=>\$y0,"dsf=f"=>\$dsf,\
  "pv=s"=>\$pv,"exponent=f"=>\$exponent,"godunovOrder=f"=>\$godunovOrder,"specialOption=s"=>\$specialOption,\
  "en=s"=>\$en,"filter=i"=>\$filter,"filterFreq=i"=>\$filterFreq,"k0=f"=>\$k0,"a0=f"=>\$a0 );
* -------------------------------------------------------------------------------------------------
if( $solver eq "best" ){ $solver="choose best iterative solver"; }
if( $tz eq "poly" ){ $tz="polynomial"; }else{ $tz="trigonometric"; }
if( $order eq "2" ){ $order = "second order accurate"; }else{ $order = "fourth order accurate"; }
if( $pv eq "nc" ){ $pv = "non-conservative"; $cons=0; }
if( $pv eq "c" ){ $pv = "conservative"; $cons=1; }
if( $pv eq "g" ){ $pv = "godunov"; }
if( $pv eq "h" ){ $pv = "hemp"; }
if( $en eq "max" ){ $errorNorm="maximum norm"; }
if( $en eq "l1" ){ $errorNorm="l1 norm"; }
if( $en eq "l2" ){ $errorNorm="l2 norm"; }
*
if( $ts eq "me" ){ $ts = "modifiedEquationTimeStepping"; }
if( $ts eq "fe" ){ $ts = "forwardEuler"; }
if( $ts eq "ie" ){ $ts = "improvedEuler"; }
if( $ts eq "ab" ){ $ts = "adamsBashforth2"; }
* 
if( $bc eq "d" ){ $bc = "all=displacementBC"; }
if( $bc eq "sf" ){ $bc = "all=tractionBC"; }
if( $bc eq "slip" ){ $bc = "all=slipWall"; }
if( $bc eq "dirichlet" ){ $bc = "all=dirichletBoundaryCondition"; }
if( $bc eq "ellipseDeform" ){ $bc = "all=displacementBC , userDefinedBoundaryData\n ellipse deform\n .25 1.\n done"; }
if( $ic eq "gaussianPulse" ){ $ic="gaussianPulseInitialCondition\n Gaussian pulse: 10 2 $exponent $x0 $y0 $z0 \n"; }
if( $ic eq "zero" ){ $ic = "zeroInitialCondition"; }; 
if( $ic eq "special" ){ $ic = "specialInitialCondition"; }
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
* 
* $tFinal=10.; $tPlot=.05; $backGround="rectangle"; 
* $diss=0.; $cfl=.9;
* 
* $grid = "rectangle80.ar10"; $diss=10.; $tPlot=.2; $cfl=.5; 
*
* 
* Note: artificial dissipation is scaled by c^2
*
$grid
* 
* -new: set-up stage: 
linear elasticity
$pv
 continue
* 
* -- set the time-stepping method:
$ts
* 
apply filter $filter
if( $filter eq 1 ){ $cmds = "filter order 4\n filter frequency $filterFreq\n filter iterations 1\n filter coefficient 1. \n  filter stages 1\n explicit filter\n  exit"; }else{ $cmds = "#"; }
$cmds
*
* ----- trig IC's ----
* twilightZoneInitialCondition
* trigonometric
* TZ omega: 2 2 2 2 (fx,fy,fz,ft)
$errorNorm
* -----------------------------
close forcing options
* 
final time $tFinal
times to plot $tPlot
* 
SMPDE:lambda $lambda
SMPDE:mu $mu 
SMPDE:Godunov order of accuracy $godunovOrder
*
boundary conditions
  $bc 
  square(1,1)=tractionBC
  * square(0,1)=dirichletBoundaryCondition
  * square(0,1)=tractionBC
done  
*
debug $debug
*
displacement scale factor $dsf
dissipation $diss
order of dissipation $dissOrder
cfl $cfl
use conservative difference $cons
* 
plot divergence 1
plot vorticity 1
initial conditions options...
specialInitialCondition
Special initial condition option: $specialOption
  $ySurf=0.; $nk=1; $period=2.; $xShift=-.5; 
  $ySurf $period $xShift
#-   $nk
#-    $k0 $a0 $b0
$amp=1.; 
# -- See cgDoc/sm/rayleigh.maple:
# Fourier series for a Gaussian Pulse: f=exp(-(x/w)^2), width=2.000000e-01. Interval [-1,1] 
16   # number of terms
$cmd="";  # holds commands
$pulseFourierCoeff=$amp*8.8622692545e-02; $cmd .= "  0  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*1.7292554358e-01; $cmd .= "  1  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*1.6058751920e-01; $cmd .= "  2  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*1.4194916916e-01; $cmd .= "  3  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*1.1943245159e-01; $cmd .= "  4  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*9.5648962964e-02; $cmd .= "  5  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*7.2913275847e-02; $cmd .= "  6  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*5.2905568066e-02; $cmd .= "  7  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*3.6539666530e-02; $cmd .= "  8  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*2.4021283208e-02; $cmd .= "  9  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*1.5031290003e-02; $cmd .= " 10  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*8.9529205508e-03; $cmd .= " 11  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*5.0757664823e-03; $cmd .= " 12  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*2.7390941662e-03; $cmd .= " 13  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*1.4069566570e-03; $cmd .= " 14  $pulseFourierCoeff  0. \n";
$pulseFourierCoeff=$amp*6.8789618475e-04; $cmd .= " 15  $pulseFourierCoeff  0. \n";
$cmd .="#";
# -- execute the pulse commands: 
$cmd
# 
# -- two wave-numbers: 
#    2
#    2 .1
#    4 .05
close initial conditions options
*
$checkErrors=0;
if( $ic eq "specialInitialCondition" ){ $checkErrors=1; }
check errors $checkErrors
plot errors $checkErrors
*
* For displacement solvers plot velocity and stress: 
if( $pv eq "non-conservative" || $pv eq "conservative" ){ $plotCommands = "plot velocity 1\n plot stress 1"; }else{ $plotCommands="*"; }
$plotCommands
**********************************
show file options...
  OBPSF:compressed
  OBPSF:open
    $show
  * OBPSF:frequency to save 
  OBPSF:frequency to flush 50
exit
***********************************
continue
* 
  contour
    adjust grid for displacement 1
    plot:v
    exit


erase
displacement
exit
