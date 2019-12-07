#
#   One Cylinder in a square
#
#
# usage: ogen [noplot] oneCylGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -blf=<num> -ml=<> -rgd=[fixed|var]
# 
#  -blf : boundary-layer-factor : blf>1 : make grid lines near boundary this many times smaller
#  -ml = number of (extra) multigrid levels to support
#  -rgd : var=variable : decrease radial grid distance as grids are refined. fixed=fix radial grid distance
# 
# Examples:
#
#  ogen -noplot oneCylGrid -order=2 -interp=e -factor=2
#  ogen -noplot oneCylGrid -order=2 -interp=e -factor=4
#  ogen -noplot oneCylGrid -order=2 -interp=e -factor=8  
#
#  ogen -noplot oneCylGrid -order=4 -interp=e -factor=2
#  ogen -noplot oneCylGrid -order=4 -interp=e -factor=4
#  ogen -noplot oneCylGrid -order=4 -interp=e -factor=8
#  ogen -noplot oneCylGrid -order=4 -interp=e -factor=16  
#  ogen -noplot oneCylGrid -order=4 -interp=e -factor=32
#
# Fixed radius grids 
#  ogen -noplot oneCylGrid -rgd=fixed -order=2 -interp=e -factor=2
#  ogen -noplot oneCylGrid -rgd=fixed -order=2 -interp=e -factor=4
#  ogen -noplot oneCylGrid -rgd=fixed -order=2 -interp=e -factor=8  
#
# multigrid:
#  ogen -noplot oneCylGrid -order=2 -interp=e -ml=2 -factor=2
#  ogen -noplot oneCylGrid -order=2 -interp=e -ml=3 -factor=4
#  ogen -noplot oneCylGrid -order=2 -interp=e -ml=3 -factor=8
#  ogen -noplot oneCylGrid -order=2 -interp=e -ml=3 -factor=16
#  ogen -noplot oneCylGrid -order=2 -interp=e -ml=4 -factor=32 
#  ogen -noplot oneCylGrid -order=2 -interp=e -ml=4 -factor=64
# 
#  ogen -noplot oneCylGrid -order=4 -interp=e -ml=2 -factor=2
#  ogen -noplot oneCylGrid -order=4 -interp=e -ml=3 -factor=4
#  ogen -noplot oneCylGrid -order=4 -interp=e -ml=3 -factor=8
#  ogen -noplot oneCylGrid -order=4 -interp=e -ml=3 -factor=16  
#  ogen -noplot oneCylGrid -order=4 -interp=e -ml=4 -factor=32  
#  ogen -noplot oneCylGrid -order=4 -interp=e -ml=4 -factor=64  
#  ogen -noplot oneCylGrid -order=4 -interp=e -ml=5 -factor=128 
#  ogen -noplot oneCylGrid -order=4 -interp=e -ml=5 -factor=256 
#
$prefix="oneCylGrid";  $rgd="var"; $numCyl=2; 
$order=2; $factor=1; $interp="i"; $ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-1.5; $xb=1.5; $ya=-1.5; $yb=1.5; 
$blf=1.;  # blf=1 : means no stretching
$innerRad=.5; $deltaRadius0=.35; # radius for rgd fixed
#
$x1=.5; $y1= 0.  # center of cylinder 1
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"ml=i"=>\$ml,"blf=f"=> \$blf, "prefix=s"=> \$prefix, \
	    "x1=f"=> \$x1,"y1=f"=> \$y1,"x2=f"=> \$x2,"y2=f"=> \$y2,"rgd=s"=> \$rgd,\
            "deltaRadius0=s"=> \$deltaRadius0,"numCyl=i"=>\$numCyl );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
if( $rgd eq "fixed" ){ $prefix = $prefix . "Fixed"; }
$suffix = ".order$order"; 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
# 
$dw = $order+1; $iw=$order+1; 
# parallel ghost lines: for ogen we need at least:
#       .5*( iw -1 )   : implicit interpolation 
#       .5*( iw+dw-2 ) : explicit interpolation
$parallelGhost=($iw-1)/2;
if( $interp eq "e" ){  $parallelGhost=($iw+$dw-2)/2; }
if( $parallelGhost<1 ){ $parallelGhost=1; } 
minimum number of distributed ghost lines
  $parallelGhost
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
#
$ds = .1/$factor; 
$pi=4.*atan2(1.,1.);
#
create mappings
  rectangle
    set corners
      $xa $xb $ya $yb
    lines
     $nx = intmg( ($xb-$xa)/$ds + 1.5 ); 
     $ny = intmg( ($yb-$ya)/$ds + 1.5 ); 
     $nx $ny
    boundary conditions
      1 2 3 4 
    mappingName
     square
    exit
#
  Annulus
    inner radius
      $innerRad
    # decrease the radius as the grids get finer -- but support number of requested MG levels 
    # $nr = max( 16, $ml2*4 )+1;
    # $nr = intmg( $nr ); 
    # MG levels are done below in stretching 
    $nr = 17;   # make large enough to keep away interpolation from the boundary layer
    $outerRad = $innerRad + $ds*$nr/1.5;   # divide by 1.5 to account for stretching
    $radialFactor=1.25; # add more points in radial direction
    if( $rgd eq "fixed" ){ $outerRad = $innerRad + $deltaRadius0; $nr=int( $radialFactor*$deltaRadius0/$ds + 1.5 ); }
    outer radius
      $outerRad
    centre for annulus
      $x1 $y1
    lines
      $nTheta = intmg( 2.*$pi*($outerRad+$innerRad)*.5/$ds+.5 );
      $nTheta $nr
    boundary conditions
      -1 -1 1 0
    mappingName
      # When there is no stretching just use the Annulus mapping
      if( $blf==1 ){ $unStretched1="annulus1"; $stretched1="annulus1Stretched"; }else{ $unStretched1="unstretched-annulus1"; $stretched1="annulus1"; }
      $unStretched1
    exit
 # stretch the annulus *********
 #
 # Stretch coordinates
  stretch coordinates
    transform which mapping?
      $unStretched1
    STRT:multigrid levels $ml
    $dxMin = $ds/$blf; 
    Stretch r2:exp to linear
    STP:stretch r2 expl: min dx, max dx $dxMin $ds
    stretch grid
#
    mappingName
      $stretched1
    exit
 #
  exit
  generate an overlapping grid
    square
    annulus1
    done
    change parameters
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      if( $order eq 4 ) { $ngp=3; }else{ $ngp=2; }
      $ng $ng $ng $ngp $ng $ng 
    exit
    compute overlap
# pause
    exit
  save a grid (compressed)
    $name
    oneCylGrid
  exit




