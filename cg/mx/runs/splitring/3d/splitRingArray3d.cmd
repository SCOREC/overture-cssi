$factor=1;  $ds0=.1;
$za=-1.; $zb=1.; $zc=0.;
$ng=1; 
$sharp=20; 
$order=2; $factor=1; $interp = "e";  $ml=0; # default values
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; } 
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
* get command line arguments
GetOptions("name=s"=> \$name,"zb=f"=>\$zb,"order=i"=>\$order,"factor=f"=> \$factor,"interp=s"=> \$interp);
*
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
*
$suffix = ".order$order";
if( $name eq "" ){$name = "splitRingArray3d" . "$interp$factor" . $suffix . ".hdf";}
*
* domain parameters:
$ds = $ds0/$factor; # target grid spacing
# nr = number of lines in normal directions to boundaries
$nr = max( 5 + $ng + 2*($order-2), 2**($ml+2) );
$nr = intmg( $nr );
#
$pi = 4.*atan2(1.,1.);
#
#
$wallBC=7;   # BC for walls of the box 
$wallShare=7; 
# 
$stretchFactor=.75; # scale nDist to account for stretching
$nDist= $stretchFactor*($nr-2)*$ds;
if( $order eq 4 ){ $nDist=$nDist+2.*$ds; } 
# 
*
$xa = -.8;
$xb =  .8;
$ya = -.8;
$yb =  .8;
$za = -.8;
$zb =  .8;
create mappings
*
  lofted surface
    smoothed polygon sections
    $tStretchaLB=0.05;  $tStretchbLB=.5*$sharp;  # stretching at corners in tangential directions
    flat tip profile 
    edit profile
      show parameters
      vertices
        4
        0.0       .5
        0.0525        .5
        0.0525       -.5
        ${ds}       -.5
      sharpness 
        $sharp
        $sharp
        $sharp 
        $sharp 
     t-stretch
      .0 $tStretchbLB
      $tStretchaLB $tStretchbLB
      $tStretchaLB $tStretchbLB
      .0 $tStretchbLB 
    exit 
    edit section
      vertices
       6 
       0.0  -0.5
       0.5  -0.5
       0.5  0.5
       -0.5  0.5
       -0.5  -0.5
       0.0  -0.5 
      * 
      n-dist
       fixed normal distance 
       $nDist 
      sharpness 
        $sharp
        $sharp
        $sharp
        $sharp
        $sharp
        $sharp
     t-stretch
      .0 $tStretchbLB
      $tsa = 2.*$tStretchaLB; # increase t-stretch weighting since there are twice as many corners
      $tsa $tStretchbLB
      $tsa $tStretchbLB
      $tsa $tStretchbLB
      $tsa $tStretchbLB
      .0 $tStretchbLB 
     periodicity
        2 
    exit 
    lines
      18 101
    mappingName
      cap 
    exit 
  box
    mappingName
      backGround
    set corners
     $xa $xb $ya $yb $za $zb
    lines
      $nx=int( ($xb-$xa)/$ds+1.5 );
      $ny=int( ($yb-$ya)/$ds+1.5 );
      $nz=int( ($zb-$za)/$ds+1.5 );
      $nx $ny $nz
    boundary conditions
      1 2 3 4 5 6
    exit 
      builder
        create surface grid...
          Start curve:cap
            plot options...
            plot boundary lines on reference surface 1
            close plot options
            picking:choose initial curve
            surface grid options...
            initial curve:points on surface 
            $delta=-1*$ds; # reduce cap width since ghost lines are added
            # increase cap width for fourth order and coarser grids -- ghost points are on surface anyway:
            if( $order eq 4 && $factor <4 ){ $delta=$delta-2.5*$ds; } 
            $capWidthFraction=.5;   # approx fraction of the end covered by the face grids 
            $xalb=-0.5;
            $xblb=0.5;
            $yalb=-0.5;
            $zblb=0.0525;
            $x1=$capWidthFraction*$xblb-$delta; $y1=$capWidthFraction*$yalb-$delta; $z1=$zblb; 
            choose point on surface 0  -$x1 $y1 $z1
            choose point on surface 0    0  $y1 $z1
            choose point on surface 0   $x1 $y1 $z1 
            done 
            $ns = intmg( $x1/$ds +1.5 ); 
            $ns = $ns + ($ns % 2);  # Make ns even so we avoid the singularity
            points on initial curve $ns
            $dist=-2*$y1 - $ds;   # ******************************* NOTE - $ds 
            distance to march $dist
              $nMarch = intmg( abs($dist)/$ds +1.5 );
            lines to march $nMarch 
            generate 
            fourth order
            $numGhost = int( $order/2 + .5 ); 
            boundary offset $numGhost $numGhost $numGhost $numGhost  (l r b t) 
            name capFace 
          exit
        exit
  reparameterize
    transform which mapping?
      cap
    restrict parameter space
      set corners
        .0 .6 0. 1. 
      exit
    $nTheta = intmg( ( ($xblb-$xalb) + ($yblb-$yalb) )/$ds +1.5 ); 
    $ns = intmg( ( ($zblb-$zalb) + ($yblb-$yalb) )/$ds +1.5 ); 
    lines 
       $ns $nTheta 
    mappingName
      capNoEnds
    exit 
* ----------------------  cross section -----------------------------------
  smoothedPolygon
      boundary conditions
      -1 -1 $wallBC 0
      vertices
       6 
       0.0  -0.0525 
       -0.0525  -0.0525
       -0.0525  0.0525
       0.0525  0.0525
       0.0525  -0.0525
       0.0  -0.0525 
    *
    * 
      sharpness
        * sharpen the corners as the grid is refined.  
        $sharp
        $sharp
        $sharp
        $sharp
        $sharp
        $sharp 
     n-dist
       fixed normal distance 
        $nDist
     t-stretch
      .0 $tStretchbLB
      $tsa = 2.*$tStretchaLB; # increase t-stretch weighting since there are twice as many corners
      $tsa $tStretchbLB
      $tsa $tStretchbLB
      $tsa $tStretchbLB
      $tsa $tStretchbLB
      .0 $tStretchbLB 
    mappingName
      crossSection1
    exit
  spline (3D)
    enter spline points
    32
      0.05 0.15 0.0
      0.05 0.154166666667 0.0
      0.05 0.158333333333 0.0
      0.05 0.1625 0.0
      -0.0625 0.275 0.0
      -0.0666666666667 0.275 0.0
      -0.0708333333333 0.275 0.0
      -0.075 0.275 0.0
      -0.375 0.275 0.0
      -0.379166666667 0.275 0.0
      -0.383333333333 0.275 0.0
      -0.3875 0.275 0.0
      -0.5 0.1625 0.0
      -0.5 0.158333333333 0.0
      -0.5 0.154166666667 0.0
      -0.5 0.15 0.0
      -0.5 -0.15 0.0
      -0.5 -0.154166666667 0.0
      -0.5 -0.158333333333 0.0
      -0.5 -0.1625 0.0
      -0.3875 -0.275 0.0
      -0.383333333333 -0.275 0.0
      -0.379166666667 -0.275 0.0
      -0.375 -0.275 0.0
      -0.075 -0.275 0.0
      -0.0708333333333 -0.275 0.0
      -0.0666666666667 -0.275 0.0
      -0.0625 -0.275 0.0
      0.05 -0.1625 0.0
      0.05 -0.158333333333 0.0
      0.05 -0.154166666667 0.0
      0.05 -0.15 0.0
    mappingName
      sweepCurve1
    exit
  sweep
    specify scaling factors
      32
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
    choose reference mapping
      crossSection1
    choose sweep curve
      sweepCurve1
    lines
      $nz=int( 1.9500000000000002/$ds+1.5 );
      81 7 $nz
    boundary conditions
      -1 -1  $wallBC 0 0  0
    share
      0 0 $wallShare 0 0 0 
    mappingName
      resonatorBody1
    exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEnds
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     90,0
      0.050000000000000044,0,0
  shift
    0, 0.15, 0.0
  mappingName
    capNoEndsTop1 
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFace
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     90,0
      0.050000000000000044,0,0
  shift
    0, 0.15, 0.0
  mappingName
    capFaceTop1 
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEnds
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     -90,0
      0.050000000000000044,0,0
  shift
    0, -0.15, 0.0
  mappingName
    capNoEndsBottom1
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFace
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     -90,0
      0.050000000000000044,0,0
  shift
    0, -0.15, 0.0
  mappingName
    capFaceBottom1
  exit
  hyperbolic
    Start curve:capFaceBottom1
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march -$nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( $capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
       0 0   0 0 $wallBC 0   
    share
      0  0 0  0 $wallShare 0   
    mappingName 
      capFaceBottom1Body
  exit
  hyperbolic
    Start curve:capNoEndsBottom1
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march  $nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( 3*$capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
      0 0 -1 -1   $wallBC 0
    share
      0  0 0  0 $wallShare 0
    mappingName 
      capNoEndsBottom1Body
  exit
  hyperbolic
    Start curve:capFaceTop1
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march -$nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( $capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
       0 0   0 0 $wallBC 0   
    share
      0  0 0  0 $wallShare 0   
    mappingName 
      capFaceTop1Body
  exit
  hyperbolic
    Start curve:capNoEndsTop1
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march  $nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( 3*$capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
      0 0 -1 -1   $wallBC 0
    share
      0  0 0  0 $wallShare 0
    mappingName 
      capNoEndsTop1Body
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEndsBottom1Body
  rotate
     0,0
     -0.5,0.0,0.0
  rotate
     0,1 
     -0.5,0.0,0.0
  mappingName
    capNoEndsBottom1BodyRotated
  lines
    15 51 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFaceBottom1Body
  rotate
     0,0
     -0.5,0.0,0.0
  rotate
     0,1 
     -0.5,0.0,0.0
  mappingName
    capFaceBottom1BodyRotated
  lines
    20 20 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEndsTop1Body
  rotate
     0,0
     -0.5,0.0,0.0
  rotate
     0,1 
     -0.5,0.0,0.0
  mappingName
    capNoEndsTop1BodyRotated
  lines
    15 51 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFaceTop1Body
  rotate
     0,0
     -0.5,0.0,0.0
  rotate
     0,1 
     -0.5,0.0,0.0
  mappingName
    capFaceTop1BodyRotated
  lines
    20 20 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      resonatorBody1
  rotate
     0,0
     -0.5,0.0,0.0
  rotate
     0,1 
     -0.5,0.0,0.0
  mappingName
    resonatorBody1Rotated
  exit 
  $wallBC = 8;
  $wallShare=8; 
  spline (3D)
    enter spline points
    32
      0.05 0.15 -0.5
      0.05 0.154166666667 -0.5
      0.05 0.158333333333 -0.5
      0.05 0.1625 -0.5
      -0.0625 0.275 -0.5
      -0.0666666666667 0.275 -0.5
      -0.0708333333333 0.275 -0.5
      -0.075 0.275 -0.5
      -0.375 0.275 -0.5
      -0.379166666667 0.275 -0.5
      -0.383333333333 0.275 -0.5
      -0.3875 0.275 -0.5
      -0.5 0.1625 -0.5
      -0.5 0.158333333333 -0.5
      -0.5 0.154166666667 -0.5
      -0.5 0.15 -0.5
      -0.5 -0.15 -0.5
      -0.5 -0.154166666667 -0.5
      -0.5 -0.158333333333 -0.5
      -0.5 -0.1625 -0.5
      -0.3875 -0.275 -0.5
      -0.383333333333 -0.275 -0.5
      -0.379166666667 -0.275 -0.5
      -0.375 -0.275 -0.5
      -0.075 -0.275 -0.5
      -0.0708333333333 -0.275 -0.5
      -0.0666666666667 -0.275 -0.5
      -0.0625 -0.275 -0.5
      0.05 -0.1625 -0.5
      0.05 -0.158333333333 -0.5
      0.05 -0.154166666667 -0.5
      0.05 -0.15 -0.5
    mappingName
      sweepCurve2
    exit
  sweep
    specify scaling factors
      32
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
    choose reference mapping
      crossSection1
    choose sweep curve
      sweepCurve2
    lines
      $nz=int( 1.9500000000000002/$ds+1.5 );
      81 7 $nz
    boundary conditions
      -1 -1  $wallBC 0 0  0
    share
      0 0 $wallShare 0 0 0 
    mappingName
      resonatorBody2
    exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEnds
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     90,0
      0.050000000000000044,0,0
  shift
    0, 0.15, -0.5
  mappingName
    capNoEndsTop2 
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFace
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     90,0
      0.050000000000000044,0,0
  shift
    0, 0.15, -0.5
  mappingName
    capFaceTop2 
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEnds
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     -90,0
      0.050000000000000044,0,0
  shift
    0, -0.15, -0.5
  mappingName
    capNoEndsBottom2
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFace
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     -90,0
      0.050000000000000044,0,0
  shift
    0, -0.15, -0.5
  mappingName
    capFaceBottom2
  exit
  hyperbolic
    Start curve:capFaceBottom2
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march -$nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( $capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
       0 0   0 0 $wallBC 0   
    share
      0  0 0  0 $wallShare 0   
    mappingName 
      capFaceBottom2Body
  exit
  hyperbolic
    Start curve:capNoEndsBottom2
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march  $nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( 3*$capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
      0 0 -1 -1   $wallBC 0
    share
      0  0 0  0 $wallShare 0
    mappingName 
      capNoEndsBottom2Body
  exit
  hyperbolic
    Start curve:capFaceTop2
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march -$nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( $capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
       0 0   0 0 $wallBC 0   
    share
      0  0 0  0 $wallShare 0   
    mappingName 
      capFaceTop2Body
  exit
  hyperbolic
    Start curve:capNoEndsTop2
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march  $nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( 3*$capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
      0 0 -1 -1   $wallBC 0
    share
      0  0 0  0 $wallShare 0
    mappingName 
      capNoEndsTop2Body
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEndsBottom2Body
  rotate
     -15,0
     -0.5,0.0,-0.5
  rotate
     0,1 
     -0.5,0.0,-0.5
  mappingName
    capNoEndsBottom2BodyRotated
  lines
    15 51 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFaceBottom2Body
  rotate
     -15,0
     -0.5,0.0,-0.5
  rotate
     0,1 
     -0.5,0.0,-0.5
  mappingName
    capFaceBottom2BodyRotated
  lines
    20 20 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEndsTop2Body
  rotate
     -15,0
     -0.5,0.0,-0.5
  rotate
     0,1 
     -0.5,0.0,-0.5
  mappingName
    capNoEndsTop2BodyRotated
  lines
    15 51 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFaceTop2Body
  rotate
     -15,0
     -0.5,0.0,-0.5
  rotate
     0,1 
     -0.5,0.0,-0.5
  mappingName
    capFaceTop2BodyRotated
  lines
    20 20 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      resonatorBody2
  rotate
     -15,0
     -0.5,0.0,-0.5
  rotate
     0,1 
     -0.5,0.0,-0.5
  mappingName
    resonatorBody2Rotated
  exit 
  $wallBC = 9;
  $wallShare=9; 
  spline (3D)
    enter spline points
    32
      0.05 0.15 0.5
      0.05 0.154166666667 0.5
      0.05 0.158333333333 0.5
      0.05 0.1625 0.5
      -0.0625 0.275 0.5
      -0.0666666666667 0.275 0.5
      -0.0708333333333 0.275 0.5
      -0.075 0.275 0.5
      -0.375 0.275 0.5
      -0.379166666667 0.275 0.5
      -0.383333333333 0.275 0.5
      -0.3875 0.275 0.5
      -0.5 0.1625 0.5
      -0.5 0.158333333333 0.5
      -0.5 0.154166666667 0.5
      -0.5 0.15 0.5
      -0.5 -0.15 0.5
      -0.5 -0.154166666667 0.5
      -0.5 -0.158333333333 0.5
      -0.5 -0.1625 0.5
      -0.3875 -0.275 0.5
      -0.383333333333 -0.275 0.5
      -0.379166666667 -0.275 0.5
      -0.375 -0.275 0.5
      -0.075 -0.275 0.5
      -0.0708333333333 -0.275 0.5
      -0.0666666666667 -0.275 0.5
      -0.0625 -0.275 0.5
      0.05 -0.1625 0.5
      0.05 -0.158333333333 0.5
      0.05 -0.154166666667 0.5
      0.05 -0.15 0.5
    mappingName
      sweepCurve3
    exit
  sweep
    specify scaling factors
      32
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
      1.
    choose reference mapping
      crossSection1
    choose sweep curve
      sweepCurve3
    lines
      $nz=int( 1.9500000000000002/$ds+1.5 );
      81 7 $nz
    boundary conditions
      -1 -1  $wallBC 0 0  0
    share
      0 0 $wallShare 0 0 0 
    mappingName
      resonatorBody3
    exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEnds
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     90,0
      0.050000000000000044,0,0
  shift
    0, 0.15, 0.5
  mappingName
    capNoEndsTop3 
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFace
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     90,0
      0.050000000000000044,0,0
  shift
    0, 0.15, 0.5
  mappingName
    capFaceTop3 
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEnds
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     -90,0
      0.050000000000000044,0,0
  shift
    0, -0.15, 0.5
  mappingName
    capNoEndsBottom3
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFace
  scale
    0.105, 0.105, 1.
  shift
      0.050000000000000044, 0., 0.
  rotate
     -90,0
      0.050000000000000044,0,0
  shift
    0, -0.15, 0.5
  mappingName
    capFaceBottom3
  exit
  hyperbolic
    Start curve:capFaceBottom3
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march -$nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( $capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
       0 0   0 0 $wallBC 0   
    share
      0  0 0  0 $wallShare 0   
    mappingName 
      capFaceBottom3Body
  exit
  hyperbolic
    Start curve:capNoEndsBottom3
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march  $nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( 3*$capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
      0 0 -1 -1   $wallBC 0
    share
      0  0 0  0 $wallShare 0
    mappingName 
      capNoEndsBottom3Body
  exit
  hyperbolic
    Start curve:capFaceTop3
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march -$nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( $capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
       0 0   0 0 $wallBC 0   
    share
      0  0 0  0 $wallShare 0   
    mappingName 
      capFaceTop3Body
  exit
  hyperbolic
    Start curve:capNoEndsTop3
    # target grid spacing .05 .05 (tang,normal, <0 : use default)
    backward 
    distance to march  $nDist 
    $linesToMarch = $nr-1; 
    lines to march $linesToMarch
    generate
    fourth order
    *boundary offset $numGhost $numGhost $numGhost $numGhost 0 0 (l r b t b f)
    lines
     $nx = intmg( $capWidthFraction*($xblb-$xalb)/$ds +1.5 ); 
     $ny = intmg( 3*$capWidthFraction*($yblb-$yalb)/$ds +1.5 ); 
     $nx $ny $nr
    boundary conditions
      0 0 -1 -1   $wallBC 0
    share
      0  0 0  0 $wallShare 0
    mappingName 
      capNoEndsTop3Body
  exit
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEndsBottom3Body
  rotate
     15,0
     -0.5,0.0,0.5
  rotate
     0,1 
     -0.5,0.0,0.5
  mappingName
    capNoEndsBottom3BodyRotated
  lines
    15 51 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFaceBottom3Body
  rotate
     15,0
     -0.5,0.0,0.5
  rotate
     0,1 
     -0.5,0.0,0.5
  mappingName
    capFaceBottom3BodyRotated
  lines
    20 20 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capNoEndsTop3Body
  rotate
     15,0
     -0.5,0.0,0.5
  rotate
     0,1 
     -0.5,0.0,0.5
  mappingName
    capNoEndsTop3BodyRotated
  lines
    15 51 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      capFaceTop3Body
  rotate
     15,0
     -0.5,0.0,0.5
  rotate
     0,1 
     -0.5,0.0,0.5
  mappingName
    capFaceTop3BodyRotated
  lines
    20 20 6
  exit 
  transform Mappings...
  rotate/scale/shift
  transform which mapping?
      resonatorBody3
  rotate
     15,0
     -0.5,0.0,0.5
  rotate
     0,1 
     -0.5,0.0,0.5
  mappingName
    resonatorBody3Rotated
  exit 
exit
generate an overlapping grid
  backGround
  capNoEndsBottom1BodyRotated
  capFaceBottom1BodyRotated
  capNoEndsTop1BodyRotated
  capFaceTop1BodyRotated 
  resonatorBody1Rotated
  capNoEndsBottom2BodyRotated
  capFaceBottom2BodyRotated
  capNoEndsTop2BodyRotated
  capFaceTop2BodyRotated 
  resonatorBody2Rotated
  capNoEndsBottom3BodyRotated
  capFaceBottom3BodyRotated
  capNoEndsTop3BodyRotated
  capFaceTop3BodyRotated 
  resonatorBody3Rotated
  done 
  compute overlap 
  exit
#
# save an overlapping grid
save a grid (compressed)
$name
SRR
exit

# 


generate an overlapping grid
  backGround
  resonator
  done

  compute overlap
