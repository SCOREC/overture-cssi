#
#  cgcssi command file: compute errors for flow in the twilight zone
# 
# Usage:
#    cgcssi [-noplot] tz.cmd -g=<grid> -tf=<final time> -tp=<tPlot> -tz=[poly|trig] -degreex=<> -degreet=<> ...
#           -cssiVariation=<> -show=<show file> -axisym=[0|1] -bc[1234]=[noSlip|slip|d|outflow|inflow|axisym] ...
#           -move=[off|shift|oscillate|rotate]  -go=[go|halt|og]
#    -cssiVariation : jameson, godunov, nonconservative
#    -axisym : 1=axisymmetric flow 
#    -bc1=[noSlip|slip|d|outflow|inflow|axisym] : set boundaries with bc=1 to a given BC.
# 
# Examples:
#     cgcssi tz.cmd -g=square20 -tf=1. -tp=.1 
#     cgcssi tz.cmd -g=cice2.hdf -tf=1. -tp=.1 
#     cgcssi tz.cmd -g=box20 -tf=1. -tp=.1 
# -- check slip wall BC's 
#    cgcssi -noplot tz.cmd -g=square20.order2 -tf=.2 -tp=.1 -bc1=slip -bc3=slip -go=go 
#       -> 0.200 8.91e-04 3.55e-04 2.99e-04 2.67e-03 6.19e+00 1.42e-03 8.79e-01 (      18,      18)
#    cgcssi -noplot tz.cmd -g=square40.order2 -tf=.2 -tp=.1 -bc1=slip -bc3=slip -go=go 
#        -> 0.200 3.52e-04 1.14e-04 9.56e-05 1.38e-03 6.19e+00 1.38e-03 4.27e+00 (      19,      19)
#    cgcssi -noplot tz.cmd -g=nonSquare20.order2  -tf=.2 -tp=.1 -bc1=slip -bc3=slip -go=go
#        -> 0.200 8.91e-04 3.55e-04 2.99e-04 2.67e-03 6.19e+00 1.42e-03 8.94e-01 (      18,      18) 
#    cgcssi -noplot tz.cmd -g=nonSquare40.order2  -tf=.2 -tp=.1 -bc1=slip -bc3=slip -go=go
#        -> 0.200 3.52e-04 1.14e-04 9.56e-05 1.38e-03 6.19e+00 1.38e-03 3.99e+00 (      18,      18)
#    cgcssi -noplot tz.cmd -g=box20.order2  -tf=.1 -tp=.05 -bc1=slip -bc3=slip -bc5=slip -go=go
#        -> Max errors:  1.215e-03  &  4.518e-04  &  5.603e-04  &  5.017e-04  &  4.330e-03  &
#    cgcssi -noplot tz.cmd -g=box40.order2  -tf=.1 -tp=.05 -bc1=slip -bc3=slip -bc5=slip -go=go
#        -> Max errors:  5.110e-04  &  2.146e-04  &  1.762e-04  &  1.463e-04  &  2.487e-03  &
#
# -- nononservative:
#     cgcssi tz.cmd -g=square10 -cssiVariation=nonconservative -tf=1. -tp=.1 -degreex=2 -degreet=2 [exact]
# -- axisymmetric: TROUBLE: need to fix interior TZ forcing for Godunov
#    cgcssi -noplot tz.cmd -g=axiSquare2.order2 -tf=.2 -tp=.1 -axisym=1 -bc3=axisym -go=go
#     -- square off axis: 
#    cgcssi tz.cmd -g=axiSquare2a.order2 -tf=1. -tp=.001 -axisym=1 
#    cgcssi tz.cmd -g=axiSquare8a.order2 -tf=1. -tp=.001 -axisym=1 
#
# -- Moving:
#    cgcssi tz.cmd -g=nonSquare32.order2 -bcOption=4 -tp=.1 -tf=.5 -tz=trig -bc1=slip -slopeLimiter=0 -reduceInterpWidth=2 -ad=0. -move=shift -go=halt
#    cgcssi tz.cmd -g=nonPlug8 -bcOption=4 -tp=.1 -tf=.5 -tz=trig -bc1=slip -slopeLimiter=0 -reduceInterpWidth=2 -ad=0. -move=shift -gridToMove=plug -go=halt
# 
# --- set default values for parameters ---
$grid="square20.hdf"; $show = " "; $backGround="square"; $noplot=""; $cssiVariation="godunov"; 
$tFinal=1.; $tPlot=.1; $cfl=.9; $debug=1; $tol=.2; $x0=.5; $dtMax=1.e10; $nbz=2; $ad=0.; 
$axisym=0; $bc1=""; $bc2=""; $bc3=""; $bc4="";  $bc5="";  $bc6=""; 
$tz = "poly"; $degreex=2; $degreet=2; $fx=1.; $fy=1.; $fz=1.; $ft=1.;
$ratio=2;  $nrl=1;  # refinement ratio and number of refinement levels
$go="halt"; $reduceInterpWidth=2; $slopeLimiter=1; 
$move="off"; $gridToMove="square"; $fullGridGenFreq=10; $vShiftx=1.; 
$bcOption=0; # 4 = 2nd order BC's
$orderOfExtrapForOutflow=2; $orderOfExtrapForGhost2=2; $orderOfExtrapForInterpNeighbours=2;
$orderOfExtrapForOutflow=3; $orderOfExtrapForGhost2=3; $orderOfExtrapForInterpNeighbours=3;
# 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"l=i"=> \$nrl,"r=i"=> \$ratio, "tf=f"=>\$tFinal,"debug=i"=> \$debug, \
            "tp=f"=>\$tPlot, "bg=s"=>\$backGround,"show=s"=>\$show,"noplot=s"=>\$noplot,"tz=s"=>\$tz,\
            "degreex=i"=>\$degreex, "degreet=i"=>\$degreet,"fx=f"=>\$fx,"fy=f"=>\$fy,"fz=f"=>\$fz,"ft=f"=>\$ft,\
            "cssiVariation=s"=>\$cssiVariation,"go=s"=>\$go,"axisym=i"=>\$axisym,"bcOption=i"=>\$bcOption,\
            "bc1=s"=>\$bc1,"bc2=s"=>\$bc2,"bc3=s"=>\$bc3,"bc4=s"=>\$bc4,"bc5=s"=>\$bc5,"bc6=s"=>\$bc6,\
            "reduceInterpWidth=i"=> \$reduceInterpWidth,"slopeLimiter=i"=>\$slopeLimiter,"ad=f"=>\$ad,\
            "move=s"=>\$move,"vShiftx=f"=>\$vShiftx,"gridToMove=s"=>\$gridToMove );
# -------------------------------------------------------------------------------------------------
if( $tz eq "poly" ){ $tz="turn on polynomial"; }else{ $tz="turn on trigonometric"; }
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
if( $cssiVariation eq "godunov" ){ $cssiVariation="compressible Navier Stokes (Godunov)"; }
if( $cssiVariation eq "jameson" ){ $cssiVariation="compressible Navier Stokes (Jameson)"; }  
if( $cssiVariation eq "nonconservative" ){ $cssiVariation="compressible Navier Stokes (non-conservative)"; }  
# -- set faces with bc=1 : 
if( $bc1 eq "noSlip" ){ $bc1="bcNumber1=noSlipWall"; }\
 elsif( $bc1 eq "slip" ){ $bc1="bcNumber1=slipWall"; }\
 elsif( $bc1 eq "inflow" ){ $bc1="bcNumber1=superSonicInflow"; }\
 elsif( $bc1 eq "outflow"  ){ $bc1="bcNumber1=superSonicOutflow"; }\
 elsif( $bc1 eq "d"  ){ $bc1="bcNumber1=dirichletBoundaryCondition"; }\
 elsif( $bc1 eq "axisym"  ){ $bc1="bcNumber1=axisymmetric"; }\
 else{  $bc1="#"; }
# -- set faces with bc=2 : 
if( $bc2 eq "noSlip" ){ $bc2="bcNumber2=noSlipWall"; }\
 elsif( $bc2 eq "slip" ){ $bc2="bcNumber2=slipWall"; }\
 elsif( $bc2 eq "inflow" ){ $bc2="bcNumber2=superSonicInflow"; }\
 elsif( $bc2 eq "outflow"  ){ $bc2="bcNumber2=superSonicOutflow"; }\
 elsif( $bc2 eq "d"  ){ $bc2="bcNumber2=dirichletBoundaryCondition"; }\
 elsif( $bc2 eq "axisym"  ){ $bc2="bcNumber2=axisymmetric"; }\
 else{  $bc2="#"; }
# -- set faces with bc=3 : 
if( $bc3 eq "noSlip" ){ $bc3="bcNumber3=noSlipWall"; }\
 elsif( $bc3 eq "slip" ){ $bc3="bcNumber3=slipWall"; }\
 elsif( $bc3 eq "inflow" ){ $bc3="bcNumber3=superSonicInflow"; }\
 elsif( $bc3 eq "outflow"  ){ $bc3="bcNumber3=superSonicOutflow"; }\
 elsif( $bc3 eq "d"  ){ $bc3="bcNumber3=dirichletBoundaryCondition"; }\
 elsif( $bc3 eq "axisym"  ){ $bc3="bcNumber3=axisymmetric"; }\
 else{  $bc3="#"; }
# -- set faces with bc=4 : 
if( $bc4 eq "noSlip" ){ $bc4="bcNumber4=noSlipWall"; }\
 elsif( $bc4 eq "slip" ){ $bc4="bcNumber4=slipWall"; }\
 elsif( $bc4 eq "inflow" ){ $bc4="bcNumber4=superSonicInflow"; }\
 elsif( $bc4 eq "outflow"  ){ $bc4="bcNumber4=superSonicOutflow"; }\
 elsif( $bc4 eq "d"  ){ $bc4="bcNumber4=dirichletBoundaryCondition"; }\
 elsif( $bc4 eq "axisym"  ){ $bc4="bcNumber4=axisymmetric"; }\
 else{  $bc4="#"; }
# -- set faces with bc=5 : 
if( $bc5 eq "noSlip" ){ $bc5="bcNumber5=noSlipWall"; }\
 elsif( $bc5 eq "slip" ){ $bc5="bcNumber5=slipWall"; }\
 elsif( $bc5 eq "inflow" ){ $bc5="bcNumber5=superSonicInflow"; }\
 elsif( $bc5 eq "outflow"  ){ $bc5="bcNumber5=superSonicOutflow"; }\
 elsif( $bc5 eq "d"  ){ $bc5="bcNumber5=dirichletBoundaryCondition"; }\
 elsif( $bc5 eq "axisym"  ){ $bc5="bcNumber5=axisymmetric"; }\
 else{  $bc5="#"; }
# -- set faces with bc=6 : 
if( $bc6 eq "noSlip" ){ $bc6="bcNumber6=noSlipWall"; }\
 elsif( $bc6 eq "slip" ){ $bc6="bcNumber6=slipWall"; }\
 elsif( $bc6 eq "inflow" ){ $bc6="bcNumber6=superSonicInflow"; }\
 elsif( $bc6 eq "outflow"  ){ $bc6="bcNumber6=superSonicOutflow"; }\
 elsif( $bc6 eq "d"  ){ $bc6="bcNumber6=dirichletBoundaryCondition"; }\
 elsif( $bc6 eq "axisym"  ){ $bc6="bcNumber6=axisymmetric"; }\
 else{  $bc6="#"; }
#
# Here is the overlapping grid to use:
$grid 
#
  $cssiVariation
#    compressible Navier Stokes (Godunov)  
#*  compressible Navier Stokes (Jameson)
  define integer parameter SlopeLimiter $slopeLimiter
#
#   one step
#   add extra variables
#     1
  exit
# -- twilightzone options:
  $tz
  degree in space $degreex
  degree in time $degreet
  frequencies (x,y,z,t)   $fx $fy $fz $ft
# 
#
#  do not use iterative implicit interpolation
#
  final time $tFinal
  times to plot $tPlot
# no plotting
  plot and always wait
#
  show file options
    compressed
 # specify the max number of parallel hdf sub-files: 
    OBPSF:maximum number of parallel sub-files 2
    open
     $show
    frequency to flush
      20
    exit
# -- specify which variables will appear in the show file:
#     showfile options...
#     OBPSF:show variable: rho 1
#     OBPSF:show variable: u 0
#     OBPSF:show variable: v 0
#     OBPSF:show variable: w 0
#     OBPSF:show variable: T 0
#     OBPSF:show variable: Mach Number 0
#     OBPSF:show variable: p 1
#     close show file options
# 
 # no plotting
#
  if( $axisym eq "1" ){ $cmd="turn on axisymmetric flow"; }else{ $cmd="#"; }
  $cmd
#
  if( $reduceInterpWidth eq 2 ){ $cmds = "reduce interpolation width\n 2"; }else{ $cmds = "#"; }
  $cmds
#
#  --------------------------- MOVING ----------------------------
#
if( $move eq "off" ){ $cmd="turn off moving grids"; }else{ $cmd="turn on moving grids"; }
 $cmd
* ---
 frequency for full grid gen update $fullGridGenFreq
#
# (1,0,0) = tangent to motion
if( $move eq "shift" || $move eq "off" ){ $moveCmds="translate\n 1  0. 0.\n $vShiftx"; }\
elsif( $move eq "oscillate" ){ $moveCmds ="oscillate\n 1. .5 0\n  1.\n .25\n 0."; }\
else{ $moveCmds = "rotate\n 0. 0. 0.\n $freq $ramp"; }
if( $move ne "off" ){ $cmds = "specify grids to move\n $moveCmds\n done\n done"; }else{ $cmds="#"; }
$cmds 
#  ------------------------ END MOVING ---------------------------
#
  boundary conditions
    all=dirichletBoundaryCondition
    $bc1
    $bc2
    $bc3
    $bc4
    $bc5
    $bc6
 # all=noSlipWall uniform(T=.3572)
#    all=slipWall 
#     $backGround(0,0)=superSonicInflow uniform(r=2.6667,u=1.25,e=10.119,s=0)
#     $backGround(1,0)=superSonicOutflow
#     $backGround(0,1)=slipWall
#     $backGround(1,1)=slipWall
    done
#
# ---NOTE: set extrapolation conditions
  boundary conditions...
   order of extrap for outflow $orderOfExtrapForOutflow (-1=default)
   order of extrap for 2nd ghost line $orderOfExtrapForGhost2
   order of extrap for interp neighbours $orderOfExtrapForInterpNeighbours
  done
#
  OBPDE:slip wall boundary condition option $bcOption
#
  pde parameters
    mu
     0.0
    kThermal
     0.0
    heat release
      0.
    rate constant
      0.
   reciprocal activation energy
     1.
   OBPDE:artificial diffusion $ad $ad $ad $ad $ad $ad $ad
  done
  debug $debug
#
  turn off adaptive grids
#   save error function to the show file
#*  show amr error function
  order of AMR interpolation
      2
  error threshold
     $tol 
  regrid frequency
    $regrid=$nbz*$ratio;
    $regrid
  change error estimator parameters
    default number of smooths
      1
    set scale factors     
      2 1 1 1 1 
    done
    exit
  change adaptive grid parameters
    refinement ratio
      $ratio
    default number of refinement levels
      $nrl
    number of buffer zones
      $nbz
#    width of proper nesting 
#      1
    grid efficiency
      .7 .5 
  exit
#
#   initial conditions
# *    smooth step function
#      step function
#       $xStep
# *        5.
#       r=2.6667 u=1.25 e=10.119 
#       r=1. e=1.786 s=0.
#   continue
continue
#
$go
