#================================================================================================
#  cgmx example: Test the radiation boundary conditions
#
# Usage:
#   
#  cgmx [-noplot] rbc -g=<name> -ic=[gs,gp] -tf=<tFinal> -tp=<tPlot> -diss=<> -debug=<num> -cons=[0/1] ...
#                     -rbc=[abcEM2|rbcNonLocal|abcPML|perfectElectricalConductor|symmetry] ...
#                     -pmlWidth=<> -pmlStrength=<> -pmlPower=<> -go=[run/halt/og]
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
$tFinal=10.; $tPlot=.1; $diss=.1; $cfl=.9; $x0=0.0; $y0=0.; $z0=0.; $kx=1; $ky=0; $kz=0.; $eps=1; 
$grid="sib1.order4.hdf"; $ic="gs"; $ks="none"; $ic="gs"; $beta=50; $debug=0; $cons=0; 
$cons=0; $go="halt"; $rbc="abcEM2"; $bcn="debug $debug"; 
$pmlWidth=11; 
$pmlStrength=50.;
$pmlPower=4.; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"diss=f"=>\$diss,"tp=f"=>\$tPlot,"show=s"=>\$show,"debug=i"=>\$debug, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bcn=s"=>\$bcn,"go=s"=>\$go,"noplot=s"=>\$noplot,"ic=s"=>\$ic,"bc=s"=>\$bc,\
  "dtMax=f"=>\$dtMax, "cons=i"=>\$cons,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"kx=i"=>\$kx,"ky=i"=>\$ky,"kz=i"=>\$kz,\
   "ks=s"=>\$ks,"rbc=s"=>\$rbc,"pmlWidth=f"=>\$pmlWidth,"pmlStrength=f"=>\$pmlStrength,"pmlPower=f"=>\$pmlPower,\
   "eps=f"=>\$eps,"ic=s"=>\$ic,"beta=f"=>\$beta,"debug=i"=>\$debug );
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
# TEMP: 
# bc: square(0,0)=dirichlet
use conservative difference $cons
#
#  sigma(s) = (pmlStrength)*(s)^(pmlPower)
pml width,strength,power $pmlWidth $pmlStrength $pmlPower
# 
# Gaussian source: 100. 5. $x0 $y0 $z0 
if( $ic eq "gp" ){ $cmd="Gaussian plane wave: $beta $x0 $y0 0 (beta,x0,y0,z0)\n gaussianPlaneWave"; }
if( $ic eq "gs" ){ $cmd="gaussianSource\n Gaussian source: 100. 5. $x0 $y0 $z0"; }
$cmd 
# 
tFinal $tFinal
tPlot $tPlot
dissipation $diss
# order of dissipation 4
cfl $cfl
if( $ic eq "gp" ){ $cmd="plot errors 1"; }else{ $cmd="#"; }
$cmd 
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
continue
# plot:Hz
$go

