#================================================================================================
# cgmx example:  Scattering of a plane wave from a sphere (compute the errors)
#
#   cgmx sib.planeWaveBC
# Usage:
#   
#  cgmx [-noplot] sib.planeWaveBC -g=<name> -tf=<tFinal> -tp=<tPlot> -diss=<> -debug=<num> -cons=[0/1] ...
#                     -ic=[exact|zero] -dm=[none|gdm] -go=[run/halt/og]
#
# Examples:
#   cgmx sib.planeWaveBC -g=sib1.order4.hdf
#   cgmx sib.planeWaveBC -g=sib2.order4.hdf
#   cgmx sib.planeWaveBC -g=sibe2.order4.hdf  ** trouble, inverted ghost cells?? -> Ok if built with -nrExtra=8
#   cgmx sib.planeWaveBC -g=sibe4.order4.hdf
# 
# -- sosup
#   cgmx sib.planeWaveBC -g=sibe4.order4.ng3.hdf -method=sosup
#
# -- sphere with 3 patches:
#   cgmx sib.planeWaveBC -g=sphereInABoxe1.order2 -bg=backGround -radius=1.
#   cgmx sib.planeWaveBC -g=sphereInABoxe2.order2 -bg=backGround -radius=1.
#
# parallel: 
#   mpirun -np 2 $cgmxp sib.planeWaveBC -g=sib2.order4.hdf
#   mpirun -np 2 $cgmxp sib.planeWaveBC -g=sib4.order4.hdf
#   srun -N2 -n8 -ppdebug $cgmxp noplot sib.planeWaveBC -g=sib2.order4.hdf -tf=.5 -tp=.5
#
# NOTES:
#   (1) In this example we solve for the scattered field. 
#   (2) The grid-spacing on the sphere needs to be small enough in the normal direction. I think there is trouble
#       if the ghost lines cross to form negative volumes (?)
#================================================================================================
# 
$tFinal=10.; $tPlot=.1; $diss=.1; $dissOrder=4; $cfl=.9; $method="NFDTD";
$grid="sib1.order4.hdf"; $bg="box"; $ic="exact"; 
$radius=.5; 
$cons=0; $go="halt"; 
#
$useSosupDissipation=0; $sosupParameter=1.;  $sosupDissipationOption=1; $sosupDissipationFrequency=1;
$selectiveDissipation=0;
$modeGDM=-1; 
$dm="none"; $alphaP=1.; $a0=1.; $a1=0.; $b0=0.; $b1=1.;  # GDM parameters
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"diss=f"=>\$diss,"tp=f"=>\$tPlot,"show=s"=>\$show,"debug=i"=>\$debug, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bcn=s"=>\$bcn,"go=s"=>\$go,"bg=s"=>\$bg,\
  "dtMax=f"=>\$dtMax, "cons=i"=>\$cons,"method=s"=>\$method,"diss=f"=>\$diss,"dissOrder=i"=>\$dissOrder,\
  "radius=f"=>\$radius,"ic=s"=>\$ic,\
  "dm=s"=>\$dm,"alphaP=f"=>\$alphaP,"a0=f"=>\$a0,"a1=f"=>\$a1,"b0=f"=>\$b0,"b1=f"=>\$b1,"modeGDM=i"=>\$modeGDM,\
  "useSosupDissipation=i"=>\$useSosupDissipation,"sosupParameter=f"=>\$sosupParameter,\
  "sosupDissipationOption=i"=>\$sosupDissipationOption,"sosupDissipationFrequency=i"=>\$sosupDissipationFrequency,\
  "selectiveDissipation=i"=>\$selectiveDissipation );
# -------------------------------------------------------------------------------------------------
#
if( $dm eq "none" ){ $dm="no dispersion"; }
if( $dm eq"gdm" ){ $dm="GDM"; }
#
if( $method eq "sosup" ){ $diss=0.; }
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
# 
# $grid="sibe2.order4.hdf"; 
# $grid="sibe4.order4.hdf";
# $grid="sib2.order4.hdf";
# 
#
$grid
# 
$method
# dispersion model:
$dm
GDM mode: $modeGDM
GDM params $a0 $a1 $b0 $b1 all (a0,a1,b0,b1,domain-name)
#
# -- When we solve for the scattered field directly the PEC boundary conditions are
#    changed to be inhomogeneous: 
planeWaveBoundaryForcing
# ====  Initial condition ====
# planeWaveInitialCondition
# planeWaveScatteredFieldInitialCondition
if( $ic eq "exact" ){ $cmd="planeWaveScatteredFieldInitialCondition"; }else{ $cmd="zeroInitialCondition"; }
$cmd
# zeroInitialCondition
# ====
scatteringFromASphereKnownSolution
#
kx,ky,kz 1 0 0
#
scattering radius $radius
# 
# Gaussian plane wave: 100. -1.5 0. 0.
# bc: box=abcEM2
bc: all=perfectElectricalConductor
# **** TESTING
# **** bc: all=dirichlet
#****
bc: $bg=dirichlet
# bc: box=abcEM2
# bc: box=abcPML
# pml width,strength,power 5 30. 4
#******
# =====================================
# radius for checking errors 1.5
# =====================================
#
tFinal $tFinal
tPlot  $tPlot
# 
order of dissipation $dissOrder
dissipation  $diss
#
if( $selectiveDissipation eq "1" ){ $cmd="selective dissipation...\n  turn off rectangular\n continue"; }else{ $cmd="#"; }
$cmd 
#
use sosup dissipation $useSosupDissipation
sosup parameter $sosupParameter
sosup dissipation option $sosupDissipationOption
sosup dissipation frequency $sosupDissipationFrequency
#***********************
slow start interval -1.
#***********************
# divergence damping  0.0025   
#***********************
use conservative difference $cons 
cfl $cfl
#*********************************
show file options...
  MXSF:compressed
  MXSF:open
    $show
  MXSF:frequency to flush 10
exit
#**********************************
#
# plot scattered field 1
debug $debug
plot errors 1
check errors 1
#
continue
# add contours on the surface of the sphere 
# contour
#   add coordinate surface 3 2 0 
#   add coordinate surface 1 2 1 
# exit
$go 


