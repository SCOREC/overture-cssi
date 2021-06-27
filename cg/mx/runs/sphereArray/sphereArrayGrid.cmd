#
#  Ogen: Grid for an array pf spheres
#
# usage: ogen [-noplot] sphereArrayGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -nrExtra=<> -ml=<> -ns=<>
#
#  nrExtra: extra lines to add in the radial direction on the sphere grids 
#  ml = number of (extra) multigrid levels to support
#  ns : number of spheres, 1 or 2 
#  xa, xb, ya, yb, za, zb : bounds on the channel
# 
# examples:
#     ogen -noplot sphereArrayGrid -factor=1 -order=2
#     ogen -noplot sphereArrayGrid -factor=1 -interp=e -order=2
#     ogen -noplot sphereArrayGrid -factor=1 -order=4
# 
#     ogen -noplot sphereArrayGrid -factor=2 -interp=e -order=4
#
$prefix="sphereArrayGrid"; 
# Background: 
$xa=-2; $xb=2; $ya=-1.75; $yb=1.75; $za=-1.75; $zb=1.75; 
#
$nrExtra=3; $loadBalance=0; $ml=0; $ns=2; 
#
# -- Centers of the spheres:
# $x1=-.6; $y1=-.6; $z1=-.6; # center of sphere 1
# $x2=+.6; $y2=+.6; $z2=+.6; # center of sphere 2
$x1=-.5; $y1=-.5; $z1=+.5; # center of sphere 1
$x2=+.5; $y2=+.5; $z2=-.5; # center of sphere 2
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $dse=0.; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=i"=> \$factor,"nrExtra=i"=> \$nrExtra,"interp=s"=> \$interp,\
            "loadBalance=i"=>\$loadBalance,"xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "za=f"=>\$za,"zb=f"=>\$zb,"ml=i"=>\$ml,"ns=i"=>\$ns,"x1=f"=>\$x1,"y1=f"=>\$y1,"z1=f"=>\$z1,\
            "x2=f"=>\$x2,"y2=f"=>\$y2,"z2=f"=>\$z2 );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; $dse=1.; }
# 
$suffix = ".order$order"; 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $ns eq 1 ){ $prefix = "oneSphereInAChannel"; }
$name = $prefix . "$interp$factor" . $suffix . ".hdf";
# 
$ds=.1/$factor;
$pi=4.*atan2(1.,1.);
# 
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
#
create mappings
# 
# add extra points for stretching 
$nr=3+ $nrExtra +$order; if( $interp eq "e" ){ $nr=$nr+$order; } 
$nr = intmg( $nr );
#
$gridNames="*"; 
$numberOfSpheres=0;  $sphereShare=1; $sphereStretchb=10.; 
$sphereRadius=.5; 
$sphereBC=7; 
#
# ---- Loop to generate spheres
#
#  Array of sphere centres: 
@xv=($x1,$x2);
@yv=($y1,$y2);
@zv=($z1,$z2);
$nx=2; $ny=2; $nz=2; # number of spheres in x, y and z directions
$dist=1.5; # distance between spheres
$x0=-.75; $y0=-.75; $z0=-.75;  # offsets for first sphere
$ns=0;  # counts spheres
for( $iz=0; $iz<$nz; $iz++ ){\
for( $iy=0; $iy<$ny; $iy++ ){\
for( $ix=0; $ix<$nx; $ix++ ){\
    $xv[$ns]=$x0+$ix*$dist; $yv[$ns]=$y0+$iy*$dist; $zv[$ns]=$z0+$iz*$dist; $ns=$ns+1;\
   }}}
# add one more at center
$xv[$ns]=0.; $yv[$ns]=0.; $zv[$ns]=0.; $ns=$ns+1;
#
$cmd=""; 
for( $i=0; $i<$ns; $i++ ){\
  $cmd .= "include sphereArray.h\n"; \
 }
$cmd .="#";
$cmd
#- #
#- # sphere 1: 
#- # 
#- # $xSphere=-.5; $ySphere=-.4; $zSphere=-.4; 
#- # $xSphere=-.6; $ySphere=-.6; $zSphere=-.6; 
#- $xSphere=$x1; $ySphere=$y1; $zSphere=$z1; 
#- include $ENV{Overture}/sampleGrids/sphere.h
#- #
#- # sphere 2: 
#- # $xSphere= .5; $ySphere= .4; $zSphere= .4; 
#- # $xSphere= .6; $ySphere= .6; $zSphere= .6; 
#- $xSphere=$x2; $ySphere=$y2; $zSphere=$z2; 
#- if( $ns eq 2 ){ $cmd = "include $ENV{Overture}/sampleGrids/sphere.h"; }else{ $cmd="#"; }
#- $cmd
#- #include sphere.h
#
# Here is the back ground grid
#
Box
  set corners
    $xa $xb $ya $yb $za $zb
  lines
    $nx = intmg( ($xb-$xa)/$ds +1.5);
    $ny = intmg( ($yb-$ya)/$ds +1.5);
    $nz = intmg( ($zb-$za)/$ds +1.5);
    $nx $ny $nz
  boundary conditions
    1 2 3 4 5 6 
  mappingName
    backGround
  exit
#
exit
#
generate an overlapping grid
  backGround
  $gridNames
#  sphere1
#  northPole1
#  southPole1
  done
  change parameters
    * improve quality of interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
  exit
 # change the plot
 # open graphics
 compute overlap
exit
* save an overlapping grid
save a grid (compressed)
$name
sphereArrayGrid
exit



