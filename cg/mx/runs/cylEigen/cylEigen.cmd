#==================================================================
#  Compute eigenfunctions of a 3d cylinder (check errors) 
# 
# run this with 
#
#    cgmx [-noplot] cyleigen -g=<gridName> -tf=<> -tp=<> 
#    mpirun -np 2 $cgmxp cyleigen
#
#==================================================================
$tFinal=1.; $tPlot=.1; $cfl=.95; $diss=.5; $show=" "; $method="NFDTD";
$n=1; $m=1; $k=1;   # n=Jn, m=m*theta, k=k*Pi*z
$grid="tube1.order4.hdf"; $tPlot=.1; $cons=0; $errorNorm=0; 
$dm="none"; 
# $alphaP=1.; $a0=1.; $a1=0.; $b0=0.; $b1=1.;  # GDM parameters
$npv=1; $alphaP=1.; $modeGDM=-1; 
@a0 = (); @a1=(); @b0=(); @b1=(); # these must be null for GetOptions to work, defaults are given below  
# 
$useSosupDissipation=0;
#
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"diss=f"=>\$diss,"tp=f"=>\$tPlot,"show=s"=>\$show,"debug=i"=>\$debug, \
 "cfl=f"=>\$cfl,"go=s"=>\$go,"noplot=s"=>\$noplot,"dtMax=f"=>\$dtMax,"n=i"=>\$n,"m=i"=>\$m,"method=s"=>\$method,\
   "dm=s"=>\$dm,"alphaP=f"=>\$alphaP,"a0=f{1,}"=>\@a0,"a1=f{1,}"=>\@a1,"b0=f{1,}"=>\@b0,"b1=f{1,}"=>\@b1,\
   "npv=i"=>\$npv,"modeGDM=i"=>\$modeGDM,"cons=i"=>\$cons,"useSosupDissipation=i"=>\$useSosupDissipation,\
   "errorNorm=i"=>\$errorNorm   );
# -------------------------------------------------------------------------------------------------
if( $dm eq "none" ){ $dm="no dispersion"; }
if( $dm eq"gdm" ){ $dm="GDM"; }
# Give defaults here for array arguments: 
if( $a0[0] eq "" ){ @a0=(1,0,0,0); }
if( $a1[0] eq "" ){ @a1=(0,0,0,0); }
if( $b0[0] eq "" ){ @b0=(0,0,0,0); }
if( $b1[0] eq "" ){ @b1=(0,0,0,0); }
# 
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
if( $method eq "sosup" ){ $diss=0.; }
# 
# 2nd-order:
# tube4.hdf
# $grid="tube1.hdf"; $tPlot=.01; 
# $grid="tube2.hdf"; $tPlot=.01; 
#
# tube3e.hdf
# $grid="tube1.order4.hdf"; $tPlot=.1; 
# $grid="tube2.order4.hdf"; $tPlot=.1; 
# tube2.order4
# tube3.order4
# tube4.order4
# --- longer cylinder ---
# $grid="tube1a.order4.hdf"; $n=2; $m=3; $k=3; 
# $grid="tube2a.order4.hdf"; $n=2; $m=3; $k=3;
# $grid="tube3a.order4.hdf"; $n=2; $m=3; $k=3;
# $grid="tube4a.order4.hdf"; $n=2; $m=3; $k=3; $tPlot=.5; $show="cylEigen44a.show"; 
#
$grid
#
modifiedEquationTimeStepping
#
$method
# dispersion model:
$dm
# GDM params $a0 $a1 $b0 $b1 all (a0,a1,b0,b1,domain-name)
GDM mode: $modeGDM
$domain="all"; 
$cmd="#"; 
if( $npv == 1 ){ $cmd = "GDM params $a0[0] $a1[0] $b0[0] $b1[0] all (a0,a1,b0,b1,domain-name)"; }
if( $npv == 2 ){ \
   $cmd  = "GDM domain name: $domain\n"; \
   $cmd .= " number of polarization vectors: $npv\n"; \
   $cmd .= " GDM coeff: 0 $a0[0] $a1[0] $b0[0] $b1[0] (eqn, a0,a1,b0,b1)\n"; \
   $cmd .= " GDM coeff: 1 $a0[1] $a1[1] $b0[1] $b1[1] (eqn, a0,a1,b0,b1)"; \
      }
$cmd
#
#**
solve for magnetic field 0
#**
#twilightZone
#trigonometric 
#TZ omega: 1. 1. 1. 1.
#
#***
bc: all=perfectElectricalConductor
#bc: box(0,2)=dirichlet
#bc: box(1,2)=dirichlet
#bc: cylinder(0,1)=dirichlet
#bc: cylinder(1,1)=dirichlet
#**
#****
annulusEigenfunctionInitialCondition
 $n $m $k    * n=Jn, m=m*theta, k=k*Pi*z
annulusEigenfunctionKnownSolution
#**
# 
use conservative difference $cons
#
#
tFinal $tFinal
tPlot  $tPlot
dissipation $diss
#
use sosup dissipation $useSosupDissipation
# accuracy in space 4
# accuracy in time 4
# order of dissipation 4
cfl $cfl
#***
#
debug 0
#*********************************
show file options...
MXSF:compressed
MXSF:open
  $show
# MXSF:frequency to save 1
MXSF:frequency to flush 10
exit
#**********************************
plot errors 1
check errors 1
error norm $errorNorm
#
continue
#
$go


continue
# ----------------------------------------------


contour
plot the grid
  exit this menu
pick to add boundary surface
add boundary surface 1 1 2 
set view:0 0 0 0 1 -0.144686 -0.372044 0.916869 -0.342547 0.888155 0.306337 -0.928293 -0.269748 -0.255946
component 6
set view:0 0 0 0 1 0.81443 -0.403978 0.416541 0.0223302 0.739138 0.673184 -0.579832 -0.53896 0.610997
add boundary surface 1 1 1 
add boundary surface 0 1 2 
set view:0 0 0 0 1 0.454147 -0.868533 0.198498 -0.771917 -0.272344 0.574432 -0.444853 -0.414101 -0.79412
add boundary surface 1 0 1 
add boundary surface 0 0 2 



movie mode
finish


plot:Ez
erase
contour
pick to delete contour planes
delete contour plane 0
exit
y-r 60
x+r 10 

movie mode
finish


continue

finish


