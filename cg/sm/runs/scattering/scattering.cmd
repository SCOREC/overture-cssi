#
# cgsm -- scattering of a Gaussian plane wave 
#
# Usage:
#   
#  cgsm [-noplot] scattering -g=<name> -bc1=[d|t] -bc2=[d|t] -tf=<tFinal> -tp=<tPlot> ...
#                    -bcn=[d|sf|mixed] -diss=<> -order=<2/4> -debug=<num> -bg=<backGround> -cons=[0/1] ...
#                    -pv=[nc|c|g|h] -godunovOrder=[1|2] -mu=<> -rho=<> -go=[run/halt/og]
# 
# --- set default values for parameters ---
# 
$noplot=""; $backGround="square"; $grid="square10"; $mu=1.; $lambda=1.; $pv="nc"; $ts="me"; 
$debug = 0;  $tPlot=.1; $bcn="sf"; $cons=0; $godunovOrder=2; $iw=2; 
$kx=1.; $ky=0.; $kz=0.; $k0=0.; $beta=10; $x0=0.; $y0=0.; $z0=0.;  
$ax=0.; $ay=0; $az=0.; # (0,0,0) = choose automatically 
$order = 2; $go="run"; 
$tFinal=5.; $cfl=.9; $dsf=.2;
$en="max";
$ad=0.; $ad4=0.;  # art. diss for Godunov
$incompressible=1; # set to 1 for incompressible solids
$upwindSOS=0;      # set to 1 for upwind dissipation for second-order systems
$cdv=1.;           # divergence damping 
$known="strip"; $bc1="t"; $bc2="d"; $bc3="d"; $bc4="t"; $bc5="t"; $bc6="t"; $bc7="t";
$icase=1; $useCurlCurl=1; $skipLastPressureSolve=0; 
$n=1; $m=1; # Jn, lambda_nm
$orderInTime=-1; # -1 = use default
$psolver="yale"; $iluLevels=1; $ogesDebug=0; $rtolp=1.e-3; $atolp=1.e-4;  # For the pressure solve
#
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"known=s"=>\$known, \
 "tp=f"=>\$tPlot, "tz=s"=>\$tz, "show=s"=>\$show,"order=i"=>\$order,"debug=i"=>\$debug,"ad=f"=>\$ad,"ad4=f"=>\$ad4, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,\
 "bc1=s"=>\$bc1,"bc2=s"=>\$bc2,"bc3=s"=>\$bc3,"bc4=s"=>\$bc4,"bc5=s"=>\$bc5,"bc6=s"=>\$bc6,"bc7=s"=>\$bc7,\
 "go=s"=>\$go,"noplot=s"=>\$noplot,"iw=i"=>\$iw,\
  "mu=f"=>\$mu,"lambda=f"=>\$lambda,"rho=f"=>\$rho,"dtMax=f"=>\$dtMax, "cons=i"=>\$cons,"pv=s"=>\$pv,\
  "kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"k0=f"=>\$k0,"beta=f"=>\$beta,"ts=s"=>\$ts,"en=s"=>\$en,\
  "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"ax=f"=>\$ax,"ay=f"=>\$ay,"az=f"=>\$az,\
  "upwindSOS=i"=>\$upwindSOS,"cdv=f"=>\$cdv,"useCurlCurl=i"=>\$useCurlCurl,"orderInTime=i"=>\$orderInTime,\
  "skipLastPressureSolve=i"=>\$skipLastPressureSolve,"psolver=s"=>\$psolver,"rtolp=f"=>\$rtolp,"atolp=f"=>\$atolp,\
  "iluLevels=i"=>\$iluLevels,"n=i"=>\$n,"m=i"=>\$m );
# -------------------------------------------------------------------------------------------------
if( $psolver eq "best" ){ $psolver="choose best iterative solver"; }
if( $psolver eq "mg" ){ $psolver="multigrid"; }
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
if( $ts eq "pc2" ){ $ts = "adamsPredictorCorrector2"; }
if( $ts eq "pc4" ){ $ts = "adamsPredictorCorrector4"; }
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
if( $orderInTime>0 ){ $cmd="SMPDE:time order of accuracy $orderInTime"; }else{ $cmd="#"; }
$cmd
SMPDE:skip last pressure solve $skipLastPressureSolve
#
if( $pv eq "godunov" && $iw eq 2 ){ $cmds = "reduce interpolation width\n $iw"; }else{ $cmds="#"; }
$cmds
# 
OBTZ:user defined known solution 
  gaussian plane wave 
    $kx $ky $kz $k0 $beta $x0 $y0 $z0 $ax $ay $az
done
knownSolutionInitialCondition
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
#
  pressure solver options
   $ogesSolver=$psolver; $ogesRtol=$rtolp; $ogesAtol=$atolp; $ogesIluLevels=$iluLevels;
   include $ENV{CG}/ins/cmd/ogesOptions.h
  exit
#
# use conservative difference $cons
boundary conditions
  if( $bc1 eq "t" ){ $bc1 = "tractionBC"; }else{ $bc1="displacementBC"; }
  if( $bc2 eq "t" ){ $bc2 = "tractionBC"; }else{ $bc2="displacementBC"; }
  if( $bc3 eq "t" ){ $bc3 = "tractionBC"; }else{ $bc3="displacementBC"; }
  if( $bc4 eq "t" ){ $bc4 = "tractionBC"; }else{ $bc4="displacementBC"; }
  if( $bc5 eq "t" ){ $bc5 = "tractionBC"; }else{ $bc5="displacementBC"; }
  if( $bc6 eq "t" ){ $bc6 = "tractionBC"; }else{ $bc6="displacementBC"; }
  if( $bc7 eq "t" ){ $bc7 = "tractionBC"; }else{ $bc7="displacementBC"; }
  bcNumber1=$bc1
  bcNumber2=$bc2
  bcNumber3=$bc3
  bcNumber4=$bc4
  bcNumber5=$bc5
  bcNumber6=$bc6
  bcNumber7=$bc7
  # bcNumber1=tractionBC
  # bcNumber2=displacementBC
done
#
displacement scale factor $dsf
debug $debug
check errors 1
plot errors 1
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

contour
  if( $kx ne 0 ){ $cmd="plot:u2"; }else{ $cmd="plot:u1"; }
  $cmd
  plot contour lines (toggle)
  vertical scale factor 0.
exit
$go


  erase
  displacement
  exit this menu

