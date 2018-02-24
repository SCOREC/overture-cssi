#
# Rotated grid for the elastic piston FSI example (cgmp)
#
#   yc  ---------------------------------------------
#       |                                           |
#       |                fluid                      |
#       |                                           |
#   yb  ---------------------------------------------
#       |                                           |
#       |                solid                      |
#       |                                           |
#       |                                           |
#   ya  ---------------------------------------------
#       xa                                         xb 
#
# Use the rotation: \xv = (\xv_0 - \cv) . [ \cos theta, - \sin theta; \sin theta, \cos theta]
# where \cv = [.5*(xb+xa), yb]
#
# usage: ogen [noplot] elasticPiston -factor=<num> -order=[2/4/6/8] -interp=[e/i] -per=[0|1]
# Options:
#   -per = 0 = no-periodic in x,  1=periodic in x
# 
# examples:
#     ogen -noplot elasticPistonGridRotated -factor=1
#     ogen -noplot elasticPistonGridRotated -factor=2
#     ogen -noplot elasticPistonGridRotated -factor=4
#
# -- periodic in x:
#    ogen -noplot elasticPistonGridRotated -factor=1 -per=1
#    ogen -noplot elasticPistonGridRotated -factor=2 -per=1
#    ogen -noplot elasticPistonGridRotated -factor=4 -per=1
#    ogen -noplot elasticPistonGridRotated -factor=8 -per=1
#
# --- reduced grid points in x:
#     ogen -noplot elasticPistonGridRotated -dsx=.2 -prefix=elasticPistonGridRotatedX5 -factor=4 
#
$prefix="elasticPistonRotated";
$order=2; $factor=1; $interp="e"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $t=0; $per=0; 
$xa=0.; $xb=1.; $ya=-.5; $yb=.0; $yc=1; 
$dsx=-1.;   # optionally set a different value for $ds in the x-direction
$thetad=30.; # give in degrees first
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,"yc=f"=> \$yc,\
            "t=f"=> \$t,"interp=s"=> \$interp,"name=s"=> \$name,"per=i"=>\$per,"dsx=f"=>\$dsx,"prefix=s"=>\$prefix,\
            "thetad=f"=>\$thetad);
# 
$cx=.5*($xa+$xb); $cy=$yb;
$pi = 4.*atan2(1.,1.);
$theta=$thetad*$pi/180.;
#
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $per eq 1 ){ $suffix .= "p"; }
if( $name eq "" ){$name = $prefix . $thetad . "DegreesGrid" . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
if( $dsx eq -1 ){ $dsx=$ds; }
# 
$dw = $order+1; $iw=$order+1; 
#
$bcInterface=100;  # bc for interfaces
$shareInterface=100;        # share value for interfaces
#
$nx = int( ($xb-$xa)/$dsx +1.5 ); 
$nyf = int( ($yc-$yb)/$ds +1.5 ); 
$nys = int( ($yb-$ya)/$ds +1.5 ); 
#
$degree=3;
$n=$nx;      # *** For now this must match the number of grid points
$hx=($xb-$xa)/($nx-1);
$hyf=($yc-$yb)/($nyf-1);
$hys=($yb-$ya)/($nys-1);
#
$cmd0="";
$cmd1="";
$cmd2="";
$cmd3="";
$cmd4="";
$cmd5="";
$cmd6="";
for($i=0; $i<$nx; $i++){$x0=$i*$hx+$xa; $y0=$ya; \
    $x=$cx+($x0-$cx)*cos($theta)-($y0-$cy)*sin($theta); $y=$cy+($x0-$cx)*sin($theta)+($y0-$cy)*cos($theta); \
    $cmd0=$cmd0 . "$x $y\n";}# lower curve 
for($i=0; $i<$nx; $i++){$x0=$i*$hx+$xa; $y0=$yb; \
    $x=$cx+($x0-$cx)*cos($theta)-($y0-$cy)*sin($theta); $y=$cy+($x0-$cx)*sin($theta)+($y0-$cy)*cos($theta); \
    $cmd1=$cmd1 . "$x $y\n";}# middle curve
for($i=0; $i<$nx; $i++){$x0=$i*$hx+$xa; $y0=$yc; \
    $x=$cx+($x0-$cx)*cos($theta)-($y0-$cy)*sin($theta); $y=$cy+($x0-$cx)*sin($theta)+($y0-$cy)*cos($theta); \
    $cmd2=$cmd2 . "$x $y\n";}# top curve
# wdh: flip orientation of solid left and right boundaries
for($i=0; $i<$nys; $i++){$x0=$xa; $y0=$yb-$i*$hys; \
    $x=$cx+($x0-$cx)*cos($theta)-($y0-$cy)*sin($theta); $y=$cy+($x0-$cx)*sin($theta)+($y0-$cy)*cos($theta); \
    $cmd3=$cmd3 . "$x $y\n";}# solid left curve
# wdh: flip orientation of solid left and right boundaries
for($i=0; $i<$nys; $i++){$x0=$xb; $y0=$yb-$i*$hys; \
    $x=$cx+($x0-$cx)*cos($theta)-($y0-$cy)*sin($theta); $y=$cy+($x0-$cx)*sin($theta)+($y0-$cy)*cos($theta); \
    $cmd4=$cmd4 . "$x $y\n";}# solid right curve
for($i=0; $i<$nyf; $i++){$x0=$xa; $y0=$yb+$i*$hyf; \
    $x=$cx+($x0-$cx)*cos($theta)-($y0-$cy)*sin($theta); $y=$cy+($x0-$cx)*sin($theta)+($y0-$cy)*cos($theta); \
    $cmd5=$cmd5 . "$x $y\n";}# fluid left curve
for($i=0; $i<$nyf; $i++){$x0=$xb; $y0=$yb+$i*$hyf; \
    $x=$cx+($x0-$cx)*cos($theta)-($y0-$cy)*sin($theta); $y=$cy+($x0-$cx)*sin($theta)+($y0-$cy)*cos($theta); \
    $cmd6=$cmd6 . "$x $y\n";}# solid right curve
#
create mappings
#
# Make the solid domains from a TFI so that we can later change the initial shape
#
  nurbs (curve)
    parameterize by index (uniform)
    enter points
    $nx $degree
    $cmd0
    lines
      $nx
    mappingName
      lowerBoundary
    exit
  nurbs (curve)
    parameterize by index (uniform)
    enter points
    $nx $degree
    $cmd1
    lines
      $nx
    mappingName
      middleBoundary
    exit
  nurbs (curve)
    parameterize by index (uniform)
    enter points
    $nx $degree
    $cmd2
    lines
      $nx
    mappingName
      topBoundary
    exit
  nurbs (curve)
    parameterize by index (uniform)
    enter points
    $nys $degree
    $cmd3
    lines
      $nys
    mappingName
      solidLeftBoundary
    exit
  nurbs (curve)
    parameterize by index (uniform)
    enter points
    $nys $degree
    $cmd4
    lines
      $nys
    mappingName
      solidRightBoundary
    exit
  nurbs (curve)
    parameterize by index (uniform)
    enter points
    $nyf $degree
    $cmd5
    lines
      $nyf
    mappingName
      fluidLeftBoundary
    exit
  nurbs (curve)
    parameterize by index (uniform)
    enter points
    $nyf $degree
    $cmd6
    lines
      $nyf
    mappingName
      fluidRightBoundary
    exit
#
# Grid for the elastic domain: 
#   
# *** WDH*** I reversed the orientation of the solidLeftBoundary and solidRightBoundary
#  since otherwise the TFI mapping automatically flipped the top and bottom curves 
#  to math the corners.
  tfi
    choose bottom curve (r_2=0)
      middleBoundary
    choose top curve    (r_2=1)
      lowerBoundary
    choose left curve   (r_1=0)
      solidLeftBoundary
    choose right curve  (r_1=1)
      solidRightBoundary
    lines
      $nx $nys
    mappingName
      elasticStrip
    boundary conditions
      if( $per eq 0 ){ $cmd="1 2 $bcInterface 3"; }else{ $cmd="-1 -1 $bcInterface 3"; }
      $cmd
    share 
       0 0 $shareInterface 0
   exit
# Interface grid for the fluid
  hyperbolic
    Start curve:middleBoundary
    $nr =5; 
    $dist= $ds*($nr-1);
    backward
    distance to march $dist
    lines to march $nr
    # $nxc = int( $nx*1.5 ); ## TEST 
    points on initial curve $nx
    BC: left outward splay
    BC: right outward splay
    outward splay 0.0, 0.0 (left, right for outward splay BC)
    #
    # -- periodicity: Note: setting hype BC's periodic --> will set derivativePeriodic 
    if( $per eq 0 ){ $cmd="#"; }else { $cmd ="BC: left periodic\n BC: right periodic"; }
    $cmd 
    boundary conditions
      if( $per eq 0 ){ $cmd="1 2 $bcInterface 0"; }else{ $cmd="-1 -1 $bcInterface 0"; }
      $cmd
    # 
    # BC: left match to a mapping
    # fluidLeftBoundary
    # BC: right match to a mapping
    # fluidRightBoundary
    # boundary condition options...
    #wdh -- try this : Feb 11, 2018
    march along normals 1
    generate
    share
     1 2 $shareInterface 0
    name fluidInterface
#
    mapping parameters
      robust inverse 1
    close mapping dialog
    # open graphics
  exit
# Background grid for the fluid
 tfi
  choose bottom curve (r_2=0)
    middleBoundary
  choose top curve    (r_2=1)
    topBoundary
  choose left curve   (r_1=0)
    fluidLeftBoundary
  choose right curve  (r_1=1)
    fluidRightBoundary
  lines
    $nx $nyf
  boundary conditions
   if( $per eq 0 ){ $cmd="1 2 0 4  "; }else{ $cmd="-1 -1 0 4"; }
   $cmd
  share 
    1 2 0 0 
  mappingName
    fluidBackGround
  exit
exit this menu
#
generate an overlapping grid
  elasticStrip
  fluidBackGround
  fluidInterface
  done choosing mappings
  change parameters
    specify a domain
     solidDomain
     elasticStrip
    done
    specify a domain
     fluidDomain
     fluidBackGround
     fluidInterface
    done
    # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      2 2 2 2 2 2
    exit
  # open graphics
  compute overlap
exit
#
# save an overlapping grid
save a grid (compressed)
$name
elasticPistonGridRotated
exit

