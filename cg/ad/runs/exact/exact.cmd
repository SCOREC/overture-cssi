#
# cgad -- EXACT SOLUTIONS
#
# Usage:
#   
#  cgad [-noplot] exact -g=<name> -known=[squareEigen|diskEigen|annulusEigen|planeInterfaceMode] -tf=<tFinal> -tp=<tPlot> ...
#         -go=<go/halt/og> -kappa=<value> -solver=<yale/best> -order=<2/4> -ts=[adams2|euler|implicit|midPoint|bdf] ...
#         -a=<val> -b=<val> -bc=<d|m|n> -amr=[0|1] -useNewStep=[0|1] -varCoeff=[poly] ...
#         -bc1=[d|m|n|champ] -bc2=[d|m|n|champ]
#
#  bc : d=dirichlet, m=mixed boundary condition
#  useNewStep : 1=use new advance steps time stepping routines
# 
# Examples:
# 
#  cgad exact -g=square16.order2 -known=squareEigen
# 
# --- set default values for parameters ---
# 
$known="squareEigen"; $kx=1.; $ky=1.; $kz=1.;
$n=0; $m=0; $radius=1.; $amp=1.; $bcOpt=0; # for disk eigenfunctions
$tFinal=1.; $cfl=.9; $kappa=.1;  $kThermal=.1; $useChamp=0; $bc1=""; $bc2=""; $bc3=""; $bc4=""; 
$ts="adams PC"; $noplot=""; $go="halt"; $a=0.; $b=0.; $c=0.;
$useNewTimeSteppingStartup=1; # This option should properly assign past tine values at startup
$debug = 0;  $tPlot=.1; $maxIterations=100; $tol=1.e-16; $atol=1.e-16; 
$tz = "poly"; $degreex=2; $degreet=2;  $ft=1.; $bc="d"; $orderInTime=2; 
$order = 2; $useNewStep=0; $varCoeff=""; $bdfOrder=2; $implicitAdvection=0; $dtMax=-1; 
$amr=0; $ratio=2; $levels=2; $bufferZones=2; $amrTol=1.e-2; 
$xPulse=.5; $yPulse=.5; $zPulse=0.;
$flushFrequency=5;
$append=0; # set to "1" to append to an existing show file  
# 
$solver="yale"; $ogesDebug=0; $ksp="bcgs"; $pc="bjacobi"; $subksp="preonly"; $subpc="ilu"; $iluLevels=3;
# $ksp="gmres"; 
#
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"degreex=i"=>\$degreex, "degreet=i"=>\$degreet, "kappa=f"=>\$kappa,\
 "tp=f"=>\$tPlot, "solver=s"=>\$solver, "tz=s"=>\$tz, "show=s"=>\$show,"order=i"=>\$order,"amr=i"=>\$amr, \
 "ts=s"=>\$ts, "noplot=s"=>\$noplot, "go=s"=>\$go,"debug=i"=>\$debug,"a=f"=>\$a,"b=f"=>\$b,"c=f"=>\$c,\
 "bc=s"=>\$bc, "amrTol=f"=>\$amrTol, "useNewStep=i"=>\$useNewStep,"levels=i"=>\$levels,"varCoeff=s"=>\$varCoeff,\
 "bdfOrder=i"=>\$bdfOrder,"implicitAdvection=i"=>\$implicitAdvection,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,\
 "ft=f"=>\$ft,"useChamp=i"=>\$useChamp, "bc1=s"=>\$bc1, "bc2=s"=>\$bc, "bc3=s"=>\bc3, "bc4=s">\$bc4,\
 "orderInTime=i"=>\$orderInTime,"known=s"=>\$known,"cfl=f"=>\$cfl,"append=i"=>\$append, "show=s"=>\$show,\
"flushFrequency=i"=>\$flushFrequency,"m=i"=>\$m,"n=i"=>\$n,"a=f"=>\$a,"amp=f"=>\$amp,"bcOpt=i"=>\$bcOpt,"dtMax=f"=>\$dtMax,\
 "useNewTimeSteppingStartup=i"=>\$useNewTimeSteppingStartup,"radius=f"=>\$radius );
# -------------------------------------------------------------------------------------------------
if( $solver eq "best" ){ $solver="choose best iterative solver"; }
if( $tz eq "poly" ){ $tz="turn on polynomial"; }elsif( $tz eq "trig" ){ $tz="turn on trigonometric"; }\
               else{ $tz="OBTZ:pulse"; }
if( $order eq "2" ){ $order = "second order accurate"; }else{ $order = "fourth order accurate"; }
if( $ts eq "adams2" || $ts eq "pc2" ){ $ts = "adams PC"; }
if( $ts eq "pc4" || ( $ts eq "pc" && $orderInTime==4) ){ $ts="adams PC order 4"; $orderInTime=4; } 
if( $ts eq "euler" ){ $ts = "forward Euler"; }
if( $ts eq "bdf" ){ $ts = "BDF"; $useNewStep=1; }
if( $ts eq "midPoint" ){ $ts = "midpoint"; }
if( $amr eq "0" ){ $amr="turn off adaptive grids"; }else{ $amr="turn on adaptive grids"; }
if( $go eq "go" ){ $go = "movie mode\n finish"; }
if( $go eq "halt" ){ $go = "*"; }
if( $go eq "og" ){ $go = "open graphics"; }
#
# 
$grid
# 
#
$grid
# 
  advection diffusion
# 
  continue
#
  # open graphics
  # turn off TZ: 
  OBTZ:twilight zone flow 0
  if( $bc eq "d" ){ $bcOpt=0; }elsif( $bc eq "n" ){ $bcOpt=1; }
  OBTZ:user defined known solution
    if( $known eq "squareEigen" ){ $cmd="square eigenfunction\n $kx $ky $kz $bcOpt\n done"; }else{ $cmd="#"; }
    # n=angular number (Jn*cos(n*theta), m=radial number, zero number)
    if( $known eq "diskEigen" ){ $cmd="disk eigenfunction\n $n $m $radius $amp $bcOpt\n done"; }
    if( $known eq "annulusEigen" ){ $cmd="annulus eigenfunction\n $n $m $amp $bcOpt\n done"; }
    $sr=-1.; $si=1.; $c1=.01; $c2=1.; # c1*exp(alpha*x) + c2*exp(-alpha*x)
    if( $known eq "planeInterfaceMode" ){ $cmd="plane interface mode\n $amp $sr $si $ky $c1 $c2 $option\n done"; }    
    $cmd
# 
  if( $known eq "planeInterfaceMode" ){ $assignKnown=1; }else{ $assignKnown=0; }
#
  initial conditions options...
  OBIC:known solution
  # open graphics
#
# -- time-stepping method --
  $ts
#
# This option should properly assign past tine values at startup
use new time-stepping startup $useNewTimeSteppingStartup  
#
if( $orderInTime eq 4 ){ $cmd="fourth order accurate in time"; }else{ $cmd="#"; }
$cmd
#
  implicit factor .5 (1=BE,0=FE)
#  implicit factor 0. (1=BE,0=FE)
  if( $useNewStep eq 1 ){ $cmd ="use new advanceSteps versions"; }else{ $cmd="#"; }
  $cmd
  BDF order $bdfOrder
# 
  choose grids for implicit
    all=implicit
  done
# 
  if( $dtMax<0 ){ $dtMax=$tp; }
  dtMax $dtMax
  cfl $cfl
# 
  pde parameters
    kappa $kappa
    # a $a
    # b $b
    # c $c
#
    $cmd="#"; $dct=0; $dcx=2; 
    if( $varCoeff eq "poly" ){ $cmd="OBPDE:user defined coefficients\n polynomial coefficients\n $dct $dcx\n done"; }
    $cmd 
    treat advection implicitly $implicitAdvection
    thermal conductivity $kThermal 
    assign known solution at boundaries $assignKnown
  done
# 
  boundary conditions
    if( $bc eq "d" ){ $bcCommand="all=dirichletBoundaryCondition"; }\
    elsif( $bc eq "n"){ $bcCommand="all=neumannBoundaryCondition"; }\
    elsif( $bc eq "m"){ $bcCommand="all=mixedBoundaryCondition, mixedDerivative(1.*t+1.*t.n=0.)"; }\
    elsif( $bc eq "dm"){ $bcCommand="all=dirichletBoundaryCondition\n square(1,0)=mixedBoundaryCondition"; }else{ $bcCommand="all=neumannBoundaryCondition"; }
 # ****
    # $bcCommand="all=dirichletBoundaryCondition\n square(1,0)=mixedBoundaryCondition";
    $bcCommand
    if( $bc1 eq "champ" ){ $cmd="bcNumber1=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber1=heatFluxInterface"; }else{ $cmd="#"; }
    $cmd 
    if( $bc2 eq "champ" ){ $cmd="bcNumber2=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber2=heatFluxInterface"; }else{ $cmd="#"; }
    $cmd 
    if( $bc3 eq "champ" ){ $cmd="bcNumber3=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber3=heatFluxInterface"; }else{ $cmd="#"; }
    $cmd 
    if( $bc4 eq "champ" ){ $cmd="bcNumber4=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber4=heatFluxInterface"; }else{ $cmd="#"; }
    $cmd             
  done
  #
  boundary conditions...
    # turn on CHAMP CHT interface codnitions: 
    apply champ interface conditions $useChamp
  done  
#
   implicit time step solver options
    $solver
 # parallel bi-conjugate gradient stabilized
#*     lu preconditioner
#
     maximum number of iterations
      $maxIterations
     relative tolerance
       $tol
     absolute tolerance
       1.e-12
     maximum number of iterations
       $maxIterations
     debug 
       $ogesDebug
    exit
#
  show file options
    if( $append eq 0 ){ $cmd="OBPSF:create new show file"; }else{ $cmd="OBPSF:append to old show file"; }
       $cmd
    compressed
     OBPSF:maximum number of parallel sub-files 8
      open
      $show
    frequency to flush
      $flushFrequency
    exit
#
#
 debug $debug
 final time $tFinal
 times to plot $tPlot
 done
 continue
 contour
 exit
 $go

