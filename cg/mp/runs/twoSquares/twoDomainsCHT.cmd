# 
# cgmp example: solve advection-diffusion in two domains 
#
#
echo to terminal 0
# Usage:
#    cgmp [-noplot] aa -g=<name> -nu=<num> -kappa1=<num> -kappa2=<num> -ktc1=<> -ktc2=<> -tf=<tFinal> -tp=<tPlot> ...
#          -solver=<yale/best> -nc=<num> -degreex[1/2]=<num> -degreet[1/2]=<num> -tz=[poly/trig/none] ...
#          -known=[planeInterfaceMode|concentricCylinders]
#          -bg=<backGroundGrid> -ts=<fe/be/im/pc> -coupled=[0|1]
# where
#  -kappa1, -kappa2 : thermal diffusivity ( T_t + .. = kappa*Delta(T) )
#  -ktc1, -ktc2 : thermal conductivity ( heat flux = - ktc1 grad(T)  kappa=k/(rho*C) )
#  -ts = time-stepping method, fe="forward Euler", be="backward Euler", mid="mid-point" im="implicit"
#  -nc : number of correction steps for implicit predictor-corrector
#  -coupled : 1=solve coupled interface equations, 0=solve decoupled 
# 
# Examples:
# 
#  cgmp aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=1. -ktc1=1. -kappa2=.5 -ktc2=.5 -tf=.05 -tp=.01 
# 
#  cgmp aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=1. -ktc1=1. -kappa2=.5 -ktc2=.5 -tf=.05 -tp=.01 -ts=pc
# 
# -- non-matching interfaces: (in progress...)
#  cgmp -noplot aa -g=twoSquaresInterface1to2.order2 -kappa1=1. -ktc1=1. -kappa2=.1 -ktc2=.5 -tf=.01 -tp=.005 -iOmega=1. -coupled=0 -nc=3 -useNewInterfaceTransfer=1 -debug=3 >! junk
#  cgmp -noplot aa -g=twoSquaresInterfacee1.order2   -kappa1=1. -ktc1=1. -kappa2=.1 -ktc2=.5 -tf=.01 -tp=.005 -iOmega=1. -coupled=0 -nc=3  -useNewInterfaceTransfer=1 [matching case for comparison
#   ..refinement patch: 
#  cgmp -noplot aa -g=twoSquaresInterface1RefineLeft.order2 -kappa1=1. -ktc1=1. -kappa2=.1 -ktc2=.5 -tf=.01 -tp=.005 -iOmega=1. -coupled=0 -nc=3 -useNewInterfaceTransfer=1 -debug=3 
#  cgmp -noplot aa -g=twoSquaresInterface1Refine.order2 -kappa1=1. -ktc1=1. -kappa2=.1 -ktc2=.5 -tf=.01 -tp=.005 -iOmega=1. -coupled=0 -nc=3 -useNewInterfaceTransfer=1 -debug=3  [refined both sides
# 
# -- pc: 
# cgmp aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=1. -ktc1=1. -kappa2=.1 -ktc2=.5 -tf=.05 -tp=.005 -ts=pc -nc=1
#
# -- explicit but solve interface equations decoupled: 
# cgmp aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=1. -ktc1=1. -kappa2=.1 -ktc2=.5 -tf=.05 -tp=.005 -iOmega=1. -coupled=0 -nc=3
# cgmp aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=1. -ktc1=1. -kappa2=.1 -ktc2=.5 -tf=.05 -tp=.005 -iOmega=1. -coupled=0 -nc=3 -ts=pc
#
# implicit:
#  cgmp aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=.9 -ktc1=.9 -kappa2=.1 -ktc2=.1 -tf=.5 -tp=.1 -ts=im -nc=6 -coupled=0 -debug=3
#  cgmp aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=.1 -ktc1=.1 -kappa2=.9 -ktc2=.9 -tf=.5 -tp=.1 -ts=im -nc=6 -coupled=0 -debug=3
# 
#  cgmp noplot aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=.9 -kappa2=.2 -tf=.02 -tp=.02 -ts=im -nc=5 >! junka
# 
# -- test implicit and under-relaxed: 
# cgmp noplot aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=1.01 -ktc1=1.01 -kappa2=1. -ktc2=1. -tf=.2 -tp=.1 -ts=im -nc=10 -go=go -iOmega=.64 -coupled=0
# 
# parallel:
#  mpirun -np 1 $cgmpp aa.cmd
#  mpirun -np 1 -dbg=valgrindebug $cgmpp noplot aa.cmd
#
# --- set default values for parameters ---
# 
$grid="twoSquaresInterfacee1.order2.hdf"; $domain1="leftDomain"; $domain2="rightDomain"; 
$tFinal=10.; $tPlot=.1; $show = " "; $debug=0; $cfl=.9; $ghost=0; $show=""; $dtMax=.1; 
$solver="yale";  $go="halt"; $dtMax=.1; $coupled=1; $tz="poly"; $uMin=1; $uMax=-1; 
$ncc=0; $mcc=1;  # for concentric cylinders
$known="";   # [planeInterface|annulusInterface]
# $ts="fe";            # mp solver
# $tsd="forward Euler"; $numberOfCorrections=0;  # domain solver
$ts="pc"; $numberOfCorrections=1;
$degreex1=2; $degreet1=2; $a1=0.; $b1=0.; 
$degreex2=2; $degreet2=2; $a2=0.; $b2=0.; 
$kappa1=.1; $kappa2=.1; 
$ktc1=-1.; $ktc2=-1.;   # by default set ktc equal to kappa
$itol=1.e-4; $iOmega=1.; $useNewInterfaceTransfer=0; $useChamp=0; 
$useMultiStage=0; # if =1 use new multistage algorithm and request data when needed
$p1=1.; $p2=1.;   # Optimized champ parameters for domain1 and domain2
$useNewStep=0; $bdfOrder=2;
#
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal, "bg=s"=>\$backGround,"known=s"=>\$known,\
 "tp=f"=>\$tPlot, "solver=s"=>\$solver, "tz=s"=>\$tz,"degreex1=i"=>\$degreex1, "degreet1=i"=>\$degreet1,\
 "degreex2=i"=>\$degreex2, "degreet2=i"=>\$degreet2,"show=s"=>\$show,"ts=s"=>\$ts,"go=s"=>\$go,\
 "debug=i"=>\$debug,"nc=i"=> \$numberOfCorrections,"iOmega=f"=>\$iOmega,"noplot=s"=>\$noplot,"itol=f"=>\$itol,"iTol=f"=>\$itol,\
 "kappa1=f"=>\$kappa1,"kappa2=f"=>\$kappa2,"ktc1=f"=>\$ktc1,"ktc2=f"=>\$ktc2,"coupled=i"=>\$coupled,\
 "useNewInterfaceTransfer=i"=>\$useNewInterfaceTransfer,"useMultiStage=i"=>\$useMultiStage,"useChamp=i"=>\$useChamp,\
 "uMin=f"=>\$uMin,"uMax=f"=>\$uMax,"dtMax=f"=>\$dtMax,"domain1=s"=>\$domain1,"domain2=s"=>\$domain2,\
 "ncc=i"=>\$ncc, "mcc=i"=>\$mcc,"p1=f"=>\$p1,"p2=f"=>\$p2,"bdfOrder=i"=>\$bdfOrder );
# -------------------------------------------------------------------------------------------------
if( $solver eq "best" ){ $solver="choose best iterative solver"; }
if( $ts eq "fe" ){ $ts="forward Euler";  $tsd="forward Euler"; }
if( $ts eq "be" ){ $ts="backward Euler"; $tsd="backward Euler"; }
if( $ts eq "im" ){ $ts="implicit";       $tsd="implicit";  }
if( $ts eq "pc" ){ $ts="adams PC";       $tsd="adams PC";  }
if( $ts eq "mid"){ $ts="midpoint";       $tsd="forward Euler"; }  
if( $ts eq "bdf" ){ $ts = "BDF"; $useNewStep=1; }
if( $ktc1 < 0. ){ $ktc1=$kappa1; }if( $ktc2 < 0. ){ $ktc2=$kappa2; }
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
# 
if( $tz eq "none" ){ $tz="turn off twilight zone"; }
if( $tz eq "poly" ){ $tz="turn on twilight zone\n turn on polynomial"; }
if( $tz eq "trig" ){ $tz="turn on twilight zone\n turn on trigonometric"; $uMin=-1; $uMax=1.; }
if( $known ne "" ){ $tz="turn off twilight zone"; }
# 
#
# ---- specfiy the composite grid ----
$grid
#
# 
# ------- start new domain ---------- 
printf("domain1=$domain1\n");
$domainName=$domain1; $solverName="solidA"; 
# Initial conditions for cgad are $ic.  Extra commands for cgad are called "$commands"
$commands="#"; 
$amp=.1; $sr=-1.; $si=1.; $ky=1.; $c1=.01; $c2=1.; $option=0; # c1*exp(alpha*x) + c2*exp(-alpha*x), option=0 : left side
if( $known eq "planeInterfaceMode" ){ $ic="OBIC:known solution"; $commands=" OBTZ:user defined known solution\n plane interface mode\n $amp $sr $si $ky $c1 $c2 $kappa1 $ktc1 $kappa2 $ktc2 $option\n done\n OBIC:known solution"; }
# gi.inputString(answer,"Enter n, m, amp, a, b, c, D1,K1, D2, K2, option, for the exact solution");
$ampcc=1.; $acc=.5; $bcc=1.; $ccc=1.5; $option=0; 
if( $known eq "concentricCylinders" ){ $ic="OBIC:known solution"; $commands=" OBTZ:user defined known solution\n concentric cylinders\n $ncc $mcc $ampcc $acc $bcc $ccc $kappa1 $ktc1 $kappa2 $ktc2 $option\n done\n OBIC:known solution"; }   
if( $known eq "planeInterfaceMode" || $known eq "concentricCylinders" ){ $assignKnown=1; }else{ $assignKnown=0; }
# 
$bc = "all=dirichletBoundaryCondition\n bcNumber100=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
$kappa=$kappa1; $ktc=$ktc1; $degreeSpace=$degreex1; $degreeTime=$degreet1; $a=$a1; $b=$b1;
# Change the optimized Champ parameter
$commands .= "\n OBPDE:champ parameters\n -1 -1 -1 $p1\n done";
include $ENV{CG}/mp/cmd/adDomain.h
# pause
# ------- start new domain ---------- 
$domainName=$domain2; $solverName="solidB"; 
$option=1; # right side 
if( $known eq "planeInterfaceMode" ){ $ic="OBIC:known solution"; $commands=" OBTZ:user defined known solution\n plane interface mode\n $amp $sr $si $ky $c1 $c2 $kappa1 $ktc1 $kappa2 $ktc2 $option\n done\n OBIC:known solution"; }
if( $known eq "concentricCylinders" ){ $ic="OBIC:known solution"; $commands=" OBTZ:user defined known solution\n concentric cylinders\n $ncc $mcc $ampcc $acc $bcc $ccc  $kappa1 $ktc1 $kappa2 $ktc2 $option\n done\n OBIC:known solution"; }   
$bc = "all=dirichletBoundaryCondition\n bcNumber100=mixedBoundaryCondition, mixedDerivative(0.*t+1.*t.n=0.)\n bcNumber100=heatFluxInterface";
$kappa=$kappa2; $ktc=$ktc2; $degreeSpace=$degreex2; $degreeTime=$degreet2; $a=$a2; $b=$b2; 
# Change the optimized Champ parameter
$commands .= "\n OBPDE:champ parameters\n -1 -1 -1 $p2\n done";
include $ENV{CG}/mp/cmd/adDomain.h
# 
 continue
# -- set parameters for cgmp ---
#  midpoint
  $ts 
  number of PC corrections $numberOfCorrections
  OBPDE:interface tolerance $itol
  OBPDE:interface omega $iOmega
  OBPDE:solve coupled interface equations $coupled
  OBPDE:use new interface transfer $useNewInterfaceTransfer
  boundary conditions...
    # turn on CHAMP CHT interface codnitions: 
    apply champ interface conditions $useChamp
  done
  # DEFINE THE MULTI_STAGE ALGORITHM --
  if( $useMultiStage eq 1 ){ $cmd="OBPDE:multi-stage\n actions=takeStep,applyBC classNames=Cgad\n done"; }else{ $cmd="#"; }
  $cmd 
#
   $tz
# --
# 
  final time $tFinal
  times to plot $tPlot
  debug $debug
#  Here is the show file that saves the solutions for both domains
  show file options
    compressed
      open
       $show
    frequency to flush
      1
    exit
  continue
#
  contour
    # set min max if $uMin < $uMax 
    printf("uMin=$uMin, uMax=$uMax\n"); 
    if( $uMin < $uMax ){ $cmd="min max $uMin $uMax"; }else{ $cmd="#"; }
     $cmd
    exit
     $cmd
   exit
continue
echo to terminal 1
$go

 movie mode
 finish
