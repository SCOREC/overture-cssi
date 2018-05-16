# Cgins: free surface over a submerged cylinder
#
# Usage:
#   
#  cgins [-noplot] submergedCylinder -g=<name> -pGrad=<f> -surfaceTension=<f> -tf=<tFinal> -tp=<tPlot> ...
#        -solver=<yale/best> -order=<2/4> -model=<ins/boussinesq> -ts=[pc|im|afs] -debug=<num> ..,
#        -ad2=<0|1> -project=<0/1> -iv=[viscous/adv/full] -imp=<val> -rf=<val> ...
#        -smoothSurface=[0|1] -numberOfSurfaceSmooths=<i>
#        -go=[run/halt/og]
# 
#  -surfaceTension : surface tension coefficient
#  -pAtmosphere : atmosphere pressure
#  -iv : implicit variation : viscous=viscous terms implicit, adv=viscous + advection, full=full linearized version
#  -imp : .5=CN, 1.=BE, 0.=FE
#  -rf : refactor frequency
#  -go : run, halt, og=open graphics
#  -ad2 : turn on or off the 2nd order artificial dissipation 
#  -dg : the name of thee grid to deform or -dg="share=<num>" to choose all grids with a given share value 
#  -df, -da : deformation frequency and amplitude. 
# 
# Examples: (Grid from freeSurfaceGrid2d.cmd)
# 
#  cgins submergedCylinder -g=submergedCylinderGride2.order2.ml1.hdf -dg="share=100" -nu=.05 -tf=2. -tp=.01 -model=ins -go=halt 
#  cgins submergedCylinder -g=freeSurfaceGrid2de4.order2.ml1 -dg="share=100" -nu=.05 -tf=2. -tp=.01 -model=ins -go=halt 
# -- turn on "gravity" : 
#  cgins submergedCylinder -g=submergedCylinderGride2.order2.ml1.hdf -dg="share=100" -nu=.01 -tf=2. -tp=.01 -model=ins -go=halt -surfaceTension=.001 -pGrad=-5. -ad2=1
# --- set default values for parameters ---
# 
$grid="halfCylinder.hdf"; $backGround="backGround"; $bcn="noSlipWall"; $pGrad=0.; 
$deformingGrid="ice"; $deformFrequency=2.; $deformAmplitude=1.; $deformationType="free surface"; 
$tFinal=1.; $tPlot=.1; $cfl=.9; $nu=.05; $Prandtl=.72; $thermalExpansivity=.1; 
$gravity = "0. 0. 0."; 
$model="ins"; $ts="pc"; $noplot=""; $implicitVariation="viscous"; $refactorFrequency=100; 
$debug = 0;   $maxIterations=100; $tol=1.e-16; $atol=1.e-16; 
$tz = "none"; $degreex=2; $degreet=2; $fx=1.; $fy=1.; $fz=1.; $ft=1.; $dtMax=.5; 
$order = 2; $fullSystem=0; $go="halt"; 
$ogesDebug=0; $project=0; $cdv=1.; $ad2=0; $ad22=2.; 
$psolver="yale"; $solver="yale"; 
$iluLevels=1; $ogesDebug=0; 
$rtolp=1.e-4; $atolp=1.e-5;  # tolerances for the pressure solve
$rtol=1.e-4; $atol=1.e-5;    # tolerances for the implicit solver
$bc="a"; 
$surfaceTension=.1; $pAtmosphere=0.;
$smoothSurface=1; $numberOfSurfaceSmooths=3;
$freeSurfaceOption="none"; 
$generatePastHistory=0;
# Decouple implicit BCs (e.g. free surface) so we can solve scalar velociity implicit equations
$decoupleImplicitBoundaryConditions=0;
$predictorOrder=0; # 0=use default 
$useNewTimeSteppingStartup=1;  # this will regenerate past time grids and solutions
$surfacePredictor="leap-frog";
# 
#
# use warnings;
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"degreex=i"=>\$degreex, "degreet=i"=>\$degreet, "model=s"=>\$model,\
 "tp=f"=>\$tPlot,"solver=s"=>\$solver,"psolver=s"=>\$psolver, "tz=s"=>\$tz, "show=s"=>\$show,\
 "order=i"=>\$order,"debug=i"=>\$debug, \
 "ts=s"=>\$ts,"nu=f"=>\$nu,"cfl=f"=>\$cfl, "bg=s"=>\$backGround,"fullSystem=i"=>\$fullSystem, "go=s"=>\$go,\
 "noplot=s"=>\$noplot,"dtMax=f"=>\$dtMax,"project=i"=>\$project,"rf=i"=> \$refactorFrequency,"bcn=s"=>\$bcn,\
 "iv=s"=>\$implicitVariation,"dtMax=f"=>\$dtMax,"ad2=i"=>\$ad2,"ad22=f"=>\$ad22,"imp=f"=>\$implicitFactor,\
  "bc=s"=>\$bc,"dg=s"=>\$deformingGrid,"dt=s"=>\$deformationType,"da=f"=>\$deformAmplitude,"df=f"=>\$deformFrequency,\
  "surfaceTension=f"=>\$surfaceTension,"pAtmosphere=f"=>\$pAtmosphere,"pGrad=f"=>\$pGrad,\
  "smoothSurface=i"=>\$smoothSurface,"numberOfSurfaceSmooths=i"=>\$numberOfSurfaceSmooths,\
  "rtol=f"=>\$rtol,"atol=f"=>\$atol,"rtolp=f"=>\$rtolp,"atolp=f"=>\$atolp );
# -------------------------------------------------------------------------------------------------
# 
$kThermal=$nu/$Prandtl; 
if( $solver eq "best" ){ $solver="choose best iterative solver"; }
if( $solver eq "mg" ){ $solver="multigrid"; }
if( $psolver eq "best" ){ $psolver="choose best iterative solver"; }
if( $psolver eq "mg" ){ $psolver="multigrid"; }
if( $model eq "ins" ){ $model = "incompressible Navier Stokes"; }else\
                     { $model = "incompressible Navier Stokes\n Boussinesq model"; }
#
if( $tz eq "none" ){ $tz="turn off twilight zone"; }
if( $tz eq "poly" ){ $tz="turn on twilight zone\n turn on polynomial"; $cdv=0.; }
if( $tz eq "trig" ){ $tz="turn on twilight zone\n turn on trigonometric"; $cdv=0.; }
if( $order eq "2" ){ $order = "second order accurate"; }else{ $order = "fourth order accurate"; }
# 
if( $ts eq "fe" ){ $ts="forward Euler";  }
if( $ts eq "be" ){ $ts="backward Euler"; }
if( $ts eq "im" ){ $ts="implicit";       }
if( $ts eq "pc" ){ $ts="adams PC";       }
if( $ts eq "mid"){ $ts="midpoint";       }  
if( $ts eq "pc4" ){ $ts="adams PC order 4"; $useNewImp=0; } # NOTE: turn off new implicit for fourth order
if( $ts eq "afs"){ $ts="approximate factorization"; $newts=1;  $useNewImp=0;}
# 
if( $implicitVariation eq "viscous" ){ $implicitVariation = "implicitViscous"; }\
elsif( $implicitVariation eq "adv" ){ $implicitVariation = "implicitAdvectionAndViscous"; }\
elsif( $implicitVariation eq "full" ){ $implicitVariation = "implicitFullLinearized"; }\
else{ $implicitVariation = "implicitFullLinearized"; }
if( $project eq "1" ){ $project = "project initial conditions"; }else{ $project = "do not project initial conditions"; }
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
#
# specify the overlapping grid to use:
$grid
# Specify the equations we solve:
  $model
  define real parameter kappa $kThermal
  define real parameter thermalExpansivity $thermalExpansivity
  define real parameter adcBoussinesq 0. 
  # Here is the surface tension
  define real parameter surfaceTension $surfaceTension
  define real parameter pAtmosphere $pAtmosphere
  done
# 
  show file options
    compressed
    open
      $show
    frequency to flush
      20 
    exit
# -- twilightzone options:
  $tz
  degree in space $degreex
  degree in time $degreet
  frequencies (x,y,z,t)   $fx $fy $fz $ft
# 
# choose the time stepping:
  $ts
use new time-stepping startup $useNewTimeSteppingStartup
# test:
if( $predictorOrder eq 1 ){ $predictorOrder = "first order predictor"; }else{ $predictorOrder="#"; }
$predictorOrder
# 
#****************************
 turn on moving grids
  specify grids to move
    deforming body
      user defined deforming body
        $deformationType
         debug
            $debug
        restrict to y direction
        velocity order of accuracy\n $gridEvolutionVelocityAccuracy
        acceleration order of accuracy\n $gridEvolutionAccelerationAccuracy
        generate past history $generatePastHistory
        # turn on surface smoothing:
        smooth surface $smoothSurface
	free surface predictor $surfacePredictor
        use known solution for initial conditions $useKnown
        number of surface smooths: $numberOfSurfaceSmooths
	past time dt: $tPlot
      done
      if( $deformingGrid =~ /^share=/ ){ $deformingGrid =~ s/^share=//; \
                 $deformingGrid="choose grids by share flag\n $deformingGrid"; };
      $deformingGrid
   done
 done
#**************************
##  useNewImplicitMethod
  $implicitVariation
  refactor frequency $refactorFrequency
  choose grids for implicit
    all=implicit
#     square=explicit
    done
  final time $tFinal
  times to plot $tPlot
  cfl $cfl
  dtMax $dtMax
# 
# 
  plot and always wait
 # no plotting
#
# Here is where we turn on gravity as a constant pressure gradient in the  negative y direction : 
if( $pGrad != 0 ){ $cmds ="user defined forcing\n constant forcing\n 2 $pGrad\n  done\n exit";}else{ $cmds="*"; }
$cmds
#
  pde parameters
    nu $nu
    kThermal $kThermal
    gravity
      $gravity
# 
    OBPDE:second-order artificial diffusion $ad2
    OBPDE:ad21,ad22  $ad22, $ad22
    OBPDE:divergence damping  $cdv 
    OBPDE:expect inflow at outflow
  done
  OBPDE:decouple implicit boundary conditions $decoupleImplicitBoundaryConditions
#
  maximum number of iterations for implicit interpolation
     10 
#***************************************************
#
  pressure solver options
   # $ogesDebug=$debug; 
   $ogesSolver=$psolver; $ogesRtol=$rtolp; $ogesAtol=$atolp; $ogesIluLevels=$iluLevels; $ogesDtol=1e20; 
   include $ENV{CG}/ins/cmd/ogesOptions.h
  exit
#
  implicit time step solver options
   $ogesSolver=$solver; $ogesRtol=$rtol; $ogesAtol=$atol; $ogesIluLevels=1; 
   include $ENV{CG}/ins/cmd/ogesOptions.h
  exit
#
#***************************************************
#-  pressure solver options
#-     $solver
#-     relative tolerance
#-       $rtol
#-     absolute tolerance
#-       $atol
#-    exit
#-  implicit time step solver options
#-     $solver
#-     relative tolerance
#-       $rtol
#-     absolute tolerance
#-       $atol 
#-    exit
# 
  boundary conditions
    $u=1.; $T=1.; 
    bcNumber5=noSlipWall, uniform(T=$T)
    all=$bcn, uniform(T=$T)
    bcNumber4=freeSurfaceBoundaryCondition
    bcNumber1=inflowWithVelocityGiven, uniform(u=$u,T=0.)
#    bcNumber2=outflow
    bcNumber3=slipWall
   bcNumber2=outflow ,  pressure(1.*p+0.1*p.n=0.), userDefinedBoundaryData
   # set the pressure at outflow to a linear profile.
    # p = p0*(y-y1)/(y0-y1) + p1*(y-y0)/(y1-y0). (linear function: p(y0)=p0, p(y1)=p1).
    $p0=-$pGrad; $p1=0.; 
    $deltay=0.;
    $y0=-1.+$deltay; $y1=0.+$deltay; 
    pressure profile
       $p0 $p1 $y0 $y1
    done
    pressure profile
       $p0 $p1 $y0 $y1
    done
    # bcNumber1=symmetry
    # bcNumber2=symmetry
    # bcNumber1=slipWall
    ## bcNumber3=inflowWithPressureAndTangentialVelocityGiven, uniform(p=0.)
  done
# 
  initial conditions
   if( $tz eq "turn off twilight zone" ){ $ic = "uniform flow\n p=1., u=$u, T=0. \n done";}else{ $ic = "done"; }
   $ic 
  done
  debug $debug
  $project
#
#-    output options...
#-    frequency to save probes 2
#-    create a probe
#-      file name probeFile.dat
#-      nearest grid point to -.5 .1 0.
#-      exit
#-    close output options
#
  exit
  $go
