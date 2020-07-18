#================================================================================================
#  cgmx example: Test the radiation boundary conditions
#
# Usage:
#   
#  cgmx [-noplot] rbc -g=<name> -ic=[gs,gp] -tf=<tFinal> -tp=<tPlot> -diss=<> -debug=<num> -cons=[0/1] ...
#                     -rbc=[abcEM2|rbcNonLocal|abcPML|perfectElectricalConductor|symmetry] ...
#                     -pmlWidth=<> -pmlStrength=<> -pmlPower=<> -useICBB=[0|1] -go=[run/halt/og]
#
# -ic : gs=Gaussian source, gp=Gaussian plane wave
#
# -pmlWidth, -pmlStrength, -pmlPower: The pml damping function sigma(s) = (pmlStrength)*(s)^(pmlPower) where 0 <= s <= 1
# 
#  Examples: 
#  -- square: Generate this grid with: ogen -noplot squareArg -order=4 -nx=128
#     cgmx rbc -g=square128.order4.hdf -x0=.5 -y0=.5 -rbc=abcEM2 -go=halt [OK
#     cgmx rbc -g=square128.order4.hdf -x0=.5 -y0=.5 -rbc=perfectElectricalConductor -go=halt [For comparison
#     cgmx rbc -g=square128.order4.hdf -x0=.5 -y0=.5 -rbc=abcPML -pmlWidth=21 -pmlStrength=50. -go=halt 
# 
#     cgmx rbc -g=square128.order2.hdf -x0=.5 -y0=.5 -rbc=abcPML -pmlWidth=21 -pmlStrength=50. -go=halt [ second-order
# 
#  -- semi-periodic square: (Grid from squarepn.cmd)
#     cgmx rbc -g=square128npx2y4.order4.hdf -rbc=abcEM2 -go=halt                
#     cgmx rbc -g=square128npx2y4.order4.hdf -rbc=rbcNonLocal -go=halt
#     cgmx rbc -g=square128npx2y4.order4.hdf -rbc=abcPML -pmlWidth=11 -pmlStrength=50. -go=halt
#
# -- 3D box:
#   cgmx rbc -g=boxLx2Ly2Lz2Factor4.order4.hdf -rbc=abcEM2 -x0=0. -y0=0. -z0=0. -go=halt    
#   cgmx rbc -g=boxLx2Ly2Lz2Factor4.order4.hdf -rbc=abcPML -x0=0. -y0=0. -z0=0. -pmlWidth=11 -go=halt 
#   cgmx rbc -g=boxLx2Ly2Lz2Factor4.order4.hdf -rbc=perfectElectricalConductor -x0=0. -y0=0. -z0=0. -go=halt
#
$tFinal=10.; $tPlot=.1; $diss=.1; $cfl=.9; $x0=0.0; $y0=0.; $z0=0.; $kx=1; $ky=0; $kz=0.; $eps=1; $mu=1;
$ax=0; $ay=0; $az=0; # For PW and GPW - all zero = use default 
$grid="sib1.order4.hdf"; $ic="gs"; $ks="none"; $ic="gs"; $beta=50; $debug=0; $cons=0; 
$cons=0; $go="halt"; $rbc="abcEM2"; $bcn="debug $debug"; 
$pmlWidth=11; $pmlStrength=50.; $pmlPower=4.;
$useICBB=0; # 1=use initial condition bounding box for letting a wave enter from outside 
$compareToShowFile=""; $useSosupDissipation=0;
$omega=5.; $amp=1.; $rampTime=.5;
$numberOfPoles=-1; # use default: options 21 or 31 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"diss=f"=>\$diss,"tp=f"=>\$tPlot,"show=s"=>\$show,"debug=i"=>\$debug, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bcn=s"=>\$bcn,"go=s"=>\$go,"noplot=s"=>\$noplot,"ic=s"=>\$ic,"bc=s"=>\$bc,\
  "dtMax=f"=>\$dtMax, "cons=i"=>\$cons,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"kx=i"=>\$kx,"ky=i"=>\$ky,"kz=i"=>\$kz,\
   "ks=s"=>\$ks,"rbc=s"=>\$rbc,"pmlWidth=f"=>\$pmlWidth,"pmlStrength=f"=>\$pmlStrength,"pmlPower=f"=>\$pmlPower,\
   "eps=f"=>\$eps,"ic=s"=>\$ic,"beta=f"=>\$beta,"debug=i"=>\$debug,"useICBB=i"=>\$useICBB,"ax=f"=>\$ax,"ay=f"=>\$ay,"az=f"=>\$az,\
   "compareToShowFile=s"=>\$compareToShowFile,"useSosupDissipation=i"=>\$useSosupDissipation,\
   "amp=f"=>\$amp,"omega=f"=>\$omega,"rampTime=f"=>\$rampTime,"numberOfPoles=i"=>\$numberOfPoles );
# -------------------------------------------------------------------------------------------------
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
#
$grid
#
NFDTD
#
coefficients $eps 1 all (eps,mu,grid/domain name)
#
# All boundaries get the far field BC: 
bc: all=$rbc
# Number of poles in the non-local rad-BC approximation
number of poles $numberOfPoles
#
$xa=-100; $xb=.5; $ya=-100; $yb=100; $za=-100; $zb=100; 
# $useICBB=1; # use initial condition bounding box
if( $useICBB ){ $cmd="initial condition bounding box $xa $xb $ya $yb $za $zb"; }else{ $cmd="#"; }
$cmd
adjust boundaries for incident field $useICBB all
#
# TEMP: 
# bc: square(0,0)=dirichlet
use conservative difference $cons
#
#  sigma(s) = (pmlStrength)*(s)^(pmlPower)
pml width,strength,power $pmlWidth $pmlStrength $pmlPower
# 
# Gaussian source: 100. 5. $x0 $y0 $z0 
if( $ic eq "gpw" || $ic eq "gp" ){ $cmd="Gaussian plane wave: $beta $x0 $y0 0 (beta,x0,y0,z0)\n gaussianPlaneWave"; }
if( $ic eq "gs" ){ $cmd="gaussianSource\n Gaussian source: $beta $omega $x0 $y0 $z0 $amp $rampTime"; }
$cmd 
#
$epsPW=$eps; $muPW=$mu; # parameters for the incident plane wave
plane wave coefficients $ax $ay $az $epsPW $muPW
# 
tFinal $tFinal
tPlot $tPlot
dissipation $diss
# order of dissipation 4
cfl $cfl
if( $ic eq "gp" ){ $cmd="plot errors 1"; }else{ $cmd="#"; }
$cmd 
#
use sosup dissipation $useSosupDissipation
#
debug $debug 
#*********************************
show file options...
  MXSF:compressed
  MXSF:open
    $show
 # MXSF:frequency to save 
  MXSF:frequency to flush 10
exit
#**********************************
# We can compute errors compared to another solution (e.g. to test far-field BC's)
if( $compareToShowFile ne "" ){ $cmd ="compare to show file 1\n reference show file: $compareToShowFile\n plot errors 1\n check errors 1"; }else{ $cmd="#"; }
$cmd 
# 
continue
# plot:Hz
if( $az==0 ){ $cmd="plot:Ey"; }else{ $cmd="plot:Ez"; }
#
$cmd="#";
# for thinBoxGrid: 
if( $grid =~ /^thin/ ){ $cmd="plot:Ey\n contour\n  add contour plane  0.00000e+00  0.00000e+00  1.00000e+00 -1.43862e-02 -8.55536e-03  1.17979e-02\n exit";}
$cmd
$go 


