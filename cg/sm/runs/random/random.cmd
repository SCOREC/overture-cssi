#
# cgsm -- RANDOM INITIAL CONDIIONS
#
# Usage:
#   
#  cgsm [-noplot] random -g=<name> -bc1=[d|t] -bc2=[d|t] -tf=<tFinal> -tp=<tPlot> ...
#                    -bcn=[d|sf|mixed] -diss=<> -order=<2/4> -debug=<num> -bg=<backGround> -cons=[0/1] ...
#                    -pv=[nc|c|g|h] -godunovOrder=[1|2] -mu=<> -rho=<> -go=[run/halt/og]
# 
#  Annulus:
#    icase =1 : displacement + displacement
#           2 : traction + disp
#           3 : traction + traction
# --- set default values for parameters ---
# 
$noplot=""; $backGround="square"; $grid="square10"; $mu=1.; $lambda=1.; $pv="nc"; $ts="me"; 
$Rg=8.314/27.; $yield=1.e10; $basePress=0.0; $c0=2.0; $cl=1.0; $hgFlag=2; $hgVisc=4.e-2; $rho=1.;
# turn off Q:
$Rg=8.314/27.; $yield=1.e10; $basePress=0.0; $c0=0.0; $cl=0.0; $hgFlag=0; $hgVisc=4.e-2;
$apr=0.0; $bpr=0.0; $cpr=0.0; $dpr=0.4;
$debug = 0;  $tPlot=.1; $bcn="sf"; $cons=0; $godunovOrder=2; $iw=2; 
$diss=.0; $dissOrder=2;  $filter=1; $filterOrder=6; $filterStages=2;
$tz = "poly"; $degreex=2; $degreet=2; $fx=2.; $fy=$fx; $fz=$fx; $ft=$fx;
$order = 2; $go="run"; 
$tFinal=5.; $cfl=.9; $dsf=.2; $p0=2.; $p1=1.; $modem=1; $moden=0; 
$en="max";
$ad=0.; $ad4=0.;  # art. diss for Godunov
$incompressible=1; # set to 1 for incompressible solids
$upwindSOS=0;      # set to 1 for upwind dissipation for second-order systems
$cdv=1.;           # divergence damping 
$known="strip"; $bc1="t"; $bc2="d"; $bc3="d"; $bc4="d"; $bc5="d"; $bc6="d"; 
$icase=1; $useCurlCurl=1; $projectLinearMode=0; 
#
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"known=s"=>\$known, \
 "tp=f"=>\$tPlot, "tz=s"=>\$tz, "show=s"=>\$show,"order=i"=>\$order,"debug=i"=>\$debug,"ad=f"=>\$ad,"ad4=f"=>\$ad4, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bc1=s"=>\$bc1,"bc2=s"=>\$bc2,"bc3=s"=>\$bc3,"bc4=s"=>\$bc4,"bc5=s"=>\$bc5,"bc6=s"=>\$bc6,\
 "go=s"=>\$go,"noplot=s"=>\$noplot,"iw=i"=>\$iw,\
  "mu=f"=>\$mu,"lambda=f"=>\$lambda,"rho=f"=>\$rho,"dtMax=f"=>\$dtMax, "cons=i"=>\$cons,"pv=s"=>\$pv,\
  "icase=f"=>\$icase,"p0=f"=>\$p0,"p1=f"=>\$p1,"modem=i"=>\$modem,"moden=i"=>\$moden,"ts=s"=>\$ts,"en=s"=>\$en,\
  "c0=f"=>\$c0,"cl=f"=>\$cl,"filter=i"=>\$filter,"filterOrder=i"=>\$filterOrder,"filterStages=i"=>\$filterStages,\
  "upwindSOS=i"=>\$upwindSOS,"cdv=f"=>\$cdv,"useCurlCurl=i"=>\$useCurlCurl,"projectLinearMode=i"=>\$projectLinearMode );
# -------------------------------------------------------------------------------------------------
if( $solver eq "best" ){ $solver="choose best iterative solver"; }
if( $tz eq "poly" ){ $tz="polynomial"; }else{ $tz="trigonometric"; }
if( $order eq "2" ){ $order = "second order accurate"; }else{ $order = "fourth order accurate"; }
if( $model eq "ins" ){ $model = "incompressible Navier Stokes"; }else\
                     { $model = "incompressible Navier Stokes\n Boussinesq model"; }
# 
if( $pv eq "nc" ){ $pv = "non-conservative"; }
if( $pv eq "c" ){ $pv = "conservative"; $cons=1; }
if( $pv eq "g" ){ $pv = "godunov"; }
if( $pv eq "h" ){ $pv = "hemp"; }
#
if( $ts eq "me" ){ $ts = "modifiedEquationTimeStepping"; }
if( $ts eq "fe" ){ $ts = "forwardEuler"; }
if( $ts eq "ie" ){ $ts = "improvedEuler"; }
if( $ts eq "ab" ){ $ts = "adamsBashforth2"; }
# 
if( $en eq "max" ){ $errorNorm="maximum norm"; }
if( $en eq "l1" ){ $errorNorm="l1 norm"; }
if( $en eq "l2" ){ $errorNorm="l2 norm"; }
# 
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
#
$grid
# -new: set-up stage: 
linear elasticity
$pv
if( $incompressible eq 0 ){ $cmd="compressible solid"; }else{ $cmd="incompressible solid"; }
$cmd
 continue
# 
# -- set the time-stepping method:
$ts
# 
#  NOTE: assign lambda,mu BEFORE setting IC's since stress depends on lambda
SMPDE:lambda $lambda
SMPDE:mu $mu 
#
SMPDE:upwind for SOS $upwindSOS
SMPDE:divergence damping $cdv
SMPDE:use curl-curl bc $useCurlCurl
SMPDE:project linear mode $projectLinearMode
#
if( $pv eq "godunov" && $iw eq 2 ){ $cmds = "reduce interpolation width\n $iw"; }else{ $cmds="#"; }
$cmds
# 
# ---- CHOOSE RANDOM INITIAL CONDITIONS ----
  OBIC:user defined...
    random
  exit
# 
OBTZ:twilight zone flow 0 
# error norm: 
$errorNorm
# 
final time $tFinal
times to plot $tPlot
dissipation $diss
# order of dissipation $dissOrder
cfl $cfl
# use conservative difference $cons
boundary conditions
  if( $bc1 eq "t" ){ $bc1 = "tractionBC"; }else{ $bc1="displacementBC"; }
  if( $bc2 eq "t" ){ $bc2 = "tractionBC"; }else{ $bc2="displacementBC"; }
  if( $bc3 eq "t" ){ $bc3 = "tractionBC"; }else{ $bc3="displacementBC"; }
  if( $bc4 eq "t" ){ $bc4 = "tractionBC"; }else{ $bc4="displacementBC"; }
  if( $bc5 eq "t" ){ $bc5 = "tractionBC"; }else{ $bc5="displacementBC"; }
  if( $bc6 eq "t" ){ $bc6 = "tractionBC"; }else{ $bc6="displacementBC"; }
  bcNumber1=$bc1
  bcNumber2=$bc2
  bcNumber3=$bc3
  bcNumber4=$bc4
  bcNumber5=$bc5
  bcNumber6=$bc6
  # bcNumber1=tractionBC
  # bcNumber2=displacementBC
done
#
displacement scale factor $dsf
debug $debug
check errors 0
plot errors 0
# plot vorticity 1 
plot divergence 1 
# For displacement solvers plot velocity and stress: 
# if( $pv eq "non-conservative" || $pv eq "conservative" ){ $plotCommands = "plot velocity 1\n plot stress 1"; }else{ $plotCommands="*"; }
# $plotCommands
# 
#*********************************
show file options...
  OBPSF:compressed
  OBPSF:open
    $show
 # OBPSF:frequency to save 
  OBPSF:frequency to flush 50
exit
#**********************************
continue
$go