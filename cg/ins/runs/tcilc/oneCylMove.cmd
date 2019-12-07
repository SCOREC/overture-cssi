# **********************************************************************************************************
# cgins command file: Moving cylinder in a square
#
#    cgins [-noplot] oneCylMove -g=<grid-name> -nu=<> -tp=<> -tf=<> -show=<> -debug=<> -project=[0|1] ...
#          -ts=[pc|im|afs|ss|bdf|imex] -solver=[best|mg|yale] -psolver=[best|mg|yale] -cpn=<> -inflow=[uniform|ramp] ...
#          -implicitVariation=[viscous/adv/full] -implicitFactor=<val> -refactorFrequency=<> ...
#          -plotResiduals=[0|1] -orderInTime=[] -ao=[centered|upwind|bweno] -upwindOrder=[1|2|...]
#
#  -implicitVariation : viscous=viscous terms implicit, full=full linearized version
#  -implicitFactor : .5=CN, 1.=BE, 0.=FE
#  -refactorFrequency : refactor frequency
#  -cpn : coefficient of p.n in the outflow BC. Normally increase for longer domains. 
#
# Examples:
#   cgins oneCylMove -g=oneCylGride2.order2.hdf -nu=.01 -tp=.1 -tf=20 -go=halt
#
# ********************************************************************************************************
#
$grid="oneCylGride2.order2"; $show = " "; $tFinal=5.; $tPlot=.1; $nu=.1; $cfl=.9; $inflow="uniform"; 
# 
$implicitVariation="viscous"; $impGrids="all=explicit"; $newts=0;
$debug = 1;  $debugp=0; $debugi=0; $opav=1; $ssr=0; $plotResiduals=0; $flushFrequency=4; 
$maxIterations=100; $tol=1.e-16; $atol=1.e-16; 
$tz = "none"; $degreex=2; $degreet=2; $fx=1.; $fy=1.; $fz=1.; $ft=1.; $dtMax=.05; 
$fullSystem=0; $go="halt"; $move=0;  $moveOnly=0; $freq=1.; 
$show=" "; $restart="";  $outflowOption="neumann"; 
$psolver="choose best iterative solver"; $solver="choose best iterative solver"; 
# $psolver="yale"; $solver="yale"; 
$iluLevels=1; $ogesDebug=0; $project=0; 
$ts="im"; 
$implicitFactor=.5;
$freqFullUpdate=10; # frequency for using full ogen update in moving grids 
$cdv=1.; $ad2=0; $ad21=1.; $ad22=1.;  $ad4=0; $ad41=1.; $ad42=1.; 
#
$orderInTime=-1;
#
$ao="centered"; $upwindOrder=-1;
#
$rtolp=1.e-3; $atolp=1.e-4;  # tolerances for the pressure solve
$rtol=1.e-4; $atol=1.e-5;    # tolerances for the implicit solver
$tolFactor=1.; # scale solver and psolver tolerances by this factor 
#
$refactorFrequency=10000; $recomputeDt=10000; 
$cp0=1.; $cpn=.1;  # For mixed pressure BC
$rate=.25;  #   number of rotations per second 
$rampOrder=3;   $ta=0.; $tb=1.;
# -- for Kyle's AF scheme:
$afit = 10;  # max iterations for AFS
$aftol=1e-2;
$filter=0; $filterFrequency=1; $filterOrder=6; $filterStages=2; 
$cdv=1;  $cDt=.25;
# -- for steady state solver
$maxIterations=200; $plotIterations=50; 
$ogmgAutoChoose=1;
# 
#
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"degreex=i"=>\$degreex, "degreet=i"=>\$degreet, "model=s"=>\$model,\
 "tp=f"=>\$tPlot, "solver=s"=>\$solver, "psolver=s"=>\$psolver, "tz=s"=>\$tz, "show=s"=>\$show,"debug=i"=>\$debug, \
 "ts=s"=>\$ts,"nu=f"=>\$nu,"cfl=f"=>\$cfl, "bg=s"=>\$backGround,"fullSystem=i"=>\$fullSystem, "go=s"=>\$go,\
 "noplot=s"=>\$noplot,"project=i"=>\$project,"rf=i"=> \$refactorFrequency,"impGrids=s"=>\$impGrids,\
 "implicitVariation=s"=>\$implicitVariation,"dtMax=f"=>\$dtMax,"implicitFactor=f"=>\$implicitFactor,\
 "rtol=f"=>\$rtol,"atol=f"=>\$atol,"rtolp=f"=>\$rtolp,"atolp=f"=>\$atolp,"restart=s"=>\$restart,\
 "freqFullUpdate=i"=>\$freqFullUpdate,"move=i"=>\$move,"moveOnly=i"=>\$moveOnly,\
  "freq=f"=>\$freq,"debugp=i"=>\$debugp,"debugi=i"=>\$debugi,"opav=i"=>\$opav,"ssr=i"=>\$ssr,"recomputeDt=i"=>\$recomputeDt,\
  "refactorFrequency=i"=>\$refactorFrequency,"ad2=i"=>\$ad2,"ad21=f"=>\$ad21,"ad22=f"=>\$ad22,\
  "ad4=i"=>\$ad4,"ad41=f"=>\$ad41,"ad42=f"=>\$ad42,"outflowOption=s"=>\$outflowOption,"cpn=f"=>\$cpn,\
  "ogmgAutoChoose=i"=>\$ogmgAutoChoose,"inflow=s"=>\$inflow,"iluLevels=i"=>\$iluLevels,\
  "maxIterations=i"=>\$maxIterations,"plotIterations=i"=>\$plotIterations,"plotResiduals=i"=>\$plotResiduals,\
  "orderInTime=i"=>\$orderInTime,"ao=s"=>\$ao,"upwindOrder=i"=>\$upwindOrder,"flushFrequency=i"=>\$flushFrequency,\
  "rate=f"=>\$rate,"rampOrder=i"=>\$rampOrder,"tb=f"=>\$tb,"tolFactor=f"=>\$tolFactor );
# -------------------------------------------------------------------------------------------------
if( $solver eq "best" ){ $solver="choose best iterative solver"; }
if( $solver eq "mg" ){ $solver="multigrid"; }
if( $psolver eq "best" ){ $psolver="choose best iterative solver"; }
if( $psolver eq "mg" ){ $psolver="multigrid"; }
if( $tz eq "none" ){ $tz="turn off twilight zone"; }
if( $tz eq "poly" ){ $tz="turn on twilight zone\n turn on polynomial"; $cdv=0.; }
if( $tz eq "trig" ){ $tz="turn on twilight zone\n turn on trigonometric"; $cdv=0.; }
if( $model eq "ins" ){ $model = "incompressible Navier Stokes"; }else\
                     { $model = "incompressible Navier Stokes\n Boussinesq model"; }
# 
if( $ts eq "fe" ){ $ts="forward Euler";  }
if( $ts eq "be" ){ $ts="backward Euler"; }
if( $ts eq "im" ){ $ts="implicit";       }
if( $ts eq "pc" ){ $ts="adams PC";       }
if( $ts eq "pc4" ){ $ts="adams PC order 4"; }
if( $ts eq "mid"){ $ts="midpoint";       }  
if( $ts eq "afs"){ $ts="approximate factorization"; $newts = 1;}
if( $ts eq "ss"){ $ts="steady state RK-line"; }
if( $ts eq "bdf" ){ $ts="implicit BDF"; }
if( $ts eq "imex" ){ $ts="implicit explicit multistep"; }
# 
if( $ao eq "centered" ){ $ao="centered advection"; }
if( $ao eq "upwind" ){ $ao="upwind advection"; }
if( $ao eq "bweno" ){ $ao="bweno advection"; }
#
if( $implicitVariation eq "viscous" ){ $implicitVariation = "implicitViscous"; }\
elsif( $implicitVariation eq "adv" ){ $implicitVariation = "implicitAdvectionAndViscous\n useNewImplicitMethod"; }\
elsif( $implicitVariation eq "full" ){ $implicitVariation = "implicitFullLinearized\n useNewImplicitMethod"; }\
else{ $implicitVariation = "implicitFullLinearized\n useNewImplicitMethod"; }
if( $newts eq "1" ){ $newts = "use new advanceSteps versions"; }else{ $newts = "#"; }
if( $project eq "1" ){ $project = "project initial conditions"; }else{ $project = "do not project initial conditions"; }
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
if( $opav eq "0" ){ $opav = "do not average coarse grid equations"; }\
  elsif( $opav eq "2" ){ $opav ="do not average coarse curvilinear grid equations"; }else{ $opav = "#"; }
if( $ssr eq 1 ){ $ssr="show smoothing rates"; }else{ $ssr="#"; }
#
$rtol = $rtol*$tolFactor; $atol=$atol*$tolFactor; $rtolp = $rtolp*$tolFactor; $atolp=$atolp*$tolFactor;
#
$grid
#
  incompressible Navier Stokes
  exit
#
  turn off twilight zone 
#
  refactor frequency $refactorFrequency
  implicit factor $implicitFactor 
  $implicitVariation
  dtMax $dtMax
#  -- choose time-stepping method:
  $ts
  $newts
  if( $orderInTime eq 4 ){ $cmd="fourth order accurate in time\n BDF order 4"; }else{ $cmd="#"; }
  $cmd
  # advection option:
  $ao
  upwind order: $upwindOrder
  #
  # -- for the AFS scheme:
  compact finite difference
  # -- convergence parameters for the af scheme
  max number of AF corrections $afit
  AF correction relative tol $aftol
  # optionally turn this on to improve stability of the high-order AF scheme by using 2nd-order dissipation at the boundary
  OBPDE:use boundary dissipation in AF scheme 1
  # -- for steady state solver:
  max iterations $maxIterations
  plot iterations $plotIterations
  plot residuals $plotResiduals
#
#------------------- MOVE GRIDS ----------------
$pi = 4.*atan2(1.,1.);
 turn on moving grids
  specify grids to move
#
    matrix motion 
      rotate around a line 
      point on line: 0 0 0 
      tangent to line: 0 0 1 
      #  Ramped rotation:
      #   Time function is the line  f(t) = 2*pi*t composed with a ramp function
      edit time function
        $a1 = 2.*$pi; 
        linear parameters: 0,$a1 (a0,a1)
        add composed function
          ramp function
          ramp end values: 0,$rate (start,end)
          ramp times: $ta,$tb (start,end)
          ramp order: $rampOrder
          printf("rampOrder=$rampOrder\n");
          exit
        exit
      exit
#-    rotate
#-      0 0 0
#-      $rate $rampTime
      annulus1
    done
 done
 #
##  OBPDE:use new fourth order boundary conditions 1
#
  choose grids for implicit
   all=implicit
   # square=explicit
  done
#
  final time $tFinal
  times to plot $tPlot
  cfl $cfl
#
  show file options
    compressed
    OBPSF:maximum number of parallel sub-files 8
    open
     $show
    frequency to flush
      $flushFrequency
    exit
#
  no plotting
  plot and always wait
#
    maximum number of iterations for implicit interpolation
      10
#
  recompute dt every $recomputeDt
  refactor frequency $refactorFrequency
# 
  pde parameters
    nu
     $nu 
    #  turn on 2nd-order AD here:
    OBPDE:second-order artificial diffusion $ad2
    OBPDE:ad21,ad22 $ad21 , $ad22
    OBPDE:fourth-order artificial diffusion $ad4
    OBPDE:ad41,ad42 $ad41 , $ad42
    if( $outflowOption eq "neumann" ){ $cmd = "use Neumann BC at outflow"; }else{ $cmd="use extrapolate BC at outflow"; }
    $cmd
  done
#
  pressure solver options
   $ogesSolver=$psolver; $ogesRtol=$rtolp; $ogesAtol=$atolp; $ogesIluLevels=$iluLevels; $ogesDebug=$debugp; $ogmgDebug=$debugp; $ogmgCoarseGridSolver="best"; 
   include $ENV{CG}/ins/cmd/ogesOptions.h
  exit
#
  implicit time step solver options
   $ogesSolver=$solver; $ogesRtol=$rtol; $ogesAtol=$atol; $ogesDebug=$debugi;
   include $ENV{CG}/ins/cmd/ogesOptions.h
  exit
#
  boundary conditions
#
    all=noSlipWall
    # square(1,1)=slipWall
#
    done
  initial conditions
  uniform flow
    p=0., u=0.
  done
  $project
  debug
   $debug
  continue
  $go



  movie mode
  finish
