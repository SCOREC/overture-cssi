#
#  Grid for an SOLID ellipsoid in a box (grids for inside and outside the ellipsoid)
#
# usage: ogen [noplot] solidEllipsoidInABoxGrid -a=<> -b=<> -c=<> -factor=<num> -order=[2/4/6/8] -interp=[e/i] -nrExtra=<>
#                 -angle1=<f> -rotationAxis1=[0|1|2] -angle2=<f> -rotationAxis2=[0|1|2]
#                 -xShift=<> -yShift=<> -zShift=<> -rgd=[fixed|var] -ml=<> ...
#                 -stretchFactor=<> -per=[nnn|ppp|npp] -overlapOption=[default|maximize]
#
#  a,b,c : semi-major axes : note: end patches are in the z-direction corresponding to "c"
#
#  -nrExtra: extra lines to add in the radial direction on the sphere grids 
#  -rgd : var=variable : decrease radial grid distance as grids are refined. fixed=fix radial grid distance
#  -ml = number of (extra) multigrid levels to support
#  - box :  if non zero on input then set -xa=xb=-ya=yb=-za=zb=box
#
# examples:
#     ogen -noplot solidEllipsoidInABoxGrid -order=2 -factor=4
#  rotated 
#     ogen -noplot solidEllipsoidInABoxGrid -order=2 -per=npp -a=.375 -b=.4 -c=.5 -angle1=45 -rotationAxis1=1 -factor=4
#     ogen -noplot solidEllipsoidInABoxGrid -order=2 -per=npp -a=.375 -b=.4 -c=.5 -angle1=45 -rotationAxis1=1 -angle2=45 -rotationAxis2=0 -factor=4
# 
# $a=1; $b=1.5; $c=2.; 
# Case 1 
# $a=1/4.; $b=1.5/4.; $c=2./4.; 
# Case 2: sphere 
$a=2/4.; $b=2/4.; $c=2./4.; 
# Case:
$a=.75/2.; $b=.75/2.; $c=1./2.;
# 
$xa=-.75; $xb=.75; $ya=-.75; $yb=.75; $za=-.75; $zb=.75; $nrMin=5; $nrExtra=0; $rgd="var"; $name=""; $per="nnn"; 
# specify up to two rotations and a shift 
$angle1=0; $rotationAxis1=0; $angle2=0; $rotationAxis2=1; $xShift=0; $yShift=0; $zShift=0; 
# 
$box=0.; # if non zero on input then set -xa=xb=-ya=yb=-za=zb=box
$order=2; $factor=1; $interp="e"; $ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $dse=0.;
$ds0=.1; # target grid spacing for $factor=1 
# parameter defining the width of the orthographic patch: 
# $orthographicPatchParameter=.55; 
$orthographicPatchParameter=.65;   # Dec 12, 2020 -- increased 
$stretchFactor=4.; # stretch grid lines by this factor at the sphere boundary
$deltaRadius0=.25; # do not make larger than .3 or troubles with cgmx
$rFactor=.75;   # Scale radial distance by this amount to make radial spacing smaller
$suffix=""; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
$prefix="solidEllipsoidInABoxGrid";
$overlapOption="default"; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=i"=> \$factor,"nrExtra=i"=>\$nrExtra,"nrMin=i"=>\$nrMin,\
            "interp=s"=> \$interp,"rgd=s"=> \$rgd,"deltaRadius0=f"=>\$deltaRadius0,"name=s"=>\$name,\
            "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,"za=f"=>\$za,"zb=f"=>\$zb,"ml=i"=>\$ml,\
            "stretchFactor=f"=>\$stretchFactor,"box=f"=>\$box,"suffix=s"=>\$suffix,"numGhost=i"=>\$numGhost,\
            "prefix=s"=> \$prefix,"a=f"=>\$a,"b=f"=>\$b,"c=f"=>\$c,"per=s"=>\$per,"ds0=f"=>\$ds0,\
	    "orthographicPatchParameter=f"=>\$orthographicPatchParameter,"angle1=f"=>\$angle1,"angle2=f"=>\$angle2, \
            "rotationAxis1=i"=>\$rotationAxis1,"rotationAxis2=i"=>\$rotationAxis2,"rFactor=f"=>\$rFactor,\
            "xShift=f"=>\$xShift,"yShift=f"=>\$yShift,"zShift=f"=>\$zShift,"overlapOption=s"=>\$overlapOption );
# 
if( $box ne 0 ){ $xa=-$box; $xb=$box; $ya=-$box; $yb=$box; $za=-$box; $zb=$box; }
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; $dse=1.; }
# 
if( $rgd eq "fixed" ){ $prefix = $prefix . "Fixed"; $sphereWidth=$deltaRadius0; }else{ $sphereWidth=-1.; }
$suffix .= ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){ $name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=$ds0/$factor;
# 
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
#
# ---------------------------------------
# turn off graphics
# ---------------------------------------
$dw = $order+1;
$iw = $dw;
$parallelGhost=($dw+1)/2;
if( $interp eq "e" ){ $parallelGhost=($dw+$iw-2)/2; }
minimum number of distributed ghost lines
  $parallelGhost
#
create mappings
$pi=4.*atan2(1.,1.);
# number of points to use in the radial direction : $nrExtra is used for stretching 
$nr=$nrMin + ($order-2); 
if( $interp eq "e" ){ $nr=$nr+$order+$nrExtra; } 
# the coarsest MG grid is 4 pts
$nr = max( $nr, 2**($ml+2) ); 
#
$gridNames="*"; 
#
# ------ Create Ellipsoid Surface Mappings for Hyperbolic Volume Grids ----
#
include $ENV{CG}/mx/runs/dielectricBodies/ellipsoidSurfaces.h
#
#
# printf("rFactor=$rFactor\n");
# pause 
$rDist = $rFactor*($nr-1)*$ds; 	# note $rFactor -- to make radial spacing a bit smaller 
$interfaceShare=100; 
#
#    ---- Construct INNER hyperbolic grids for the ellipsoid body and two patches on the poles ----
#
$nrSave=$nr;
$nr = $nrMin + ($order-2);
if( $interp eq "e" ){ $nr = $nrMin + $order +$nrExtra; }
$rDist = $rFactor*($nr-1)*$ds; 	# note $rFactor -- to make radial spacing a bit smaller 
$directionToMarch="backward"; 
#
$ellipsoidBodyName="ellipsoidInner";
$ellipsoidSouthPoleName="southPoleInner";
$ellipsoidNorthPoleName="northPoleInner";
# 
include $ENV{CG}/mx/runs/dielectricBodies/ellipsoidVolumeGrids.h
#
#    ---- Construct OUTER hyperbolic grids for the ellipsoid body and two patches on the poles ----
#
$nr=$nrSave; 
$rDist = $rFactor*($nr-1)*$ds; 	# note $rFactor -- to make radial spacing a bit smaller 
$directionToMarch="forward"; 
#
$ellipsoidBodyName="ellipsoidOuter";
$ellipsoidSouthPoleName="southPoleOuter";
$ellipsoidNorthPoleName="northPoleOuter";
# 
include $ENV{CG}/mx/runs/dielectricBodies/ellipsoidVolumeGrids.h 
#
# Here is the background box
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
    $bcCmd = "1 2 3 4 5 6"; 
    if( $per eq "ppp" ){ $bcCmd = "-1 -1 -1 -1 -1 -1"; }
    if( $per eq "npp" ){ $bcCmd = "1 2 -1 -1 -1 -1"; }
    $bcCmd
  mappingName
    backGround
  exit
# 
sub convertToNurbs\
{ local($old,$new,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift)=@_; \
  $cmds = "nurbs \n" . \
   "interpolate from mapping with options\n" . \
   " $old \n" . \
   " parameterize by index (uniform)\n" . \
   " number of ghost points to include\n $numGhost\n" . \
   " choose degree\n" . \
   "  3 \n" . \
   "done\n" . \
   "rotate \n" . \
   " $angle1 $rotationAxis1 \n" . \
   " 0. 0. 0.\n" . \
   "rotate \n" . \
   " $angle2 $rotationAxis2 \n" . \
   " 0. 0. 0.\n" . \
   "shift\n" . \
   " $xShift $yShift $zShift\n" . \
   "mappingName\n" . \
   " $new\n" . \
   "exit"; \
}
#
#  Convert to NURBS and optionally rotate and shift 
#
$numGhost=$ng+1; # N.B. to avoid negative volumes in the ghost points interpolate ghost too in Nurbs.
convertToNurbs(ellipsoidInner,ellipsoidInnerGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs(southPoleInner,southPoleInnerGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs(northPoleInner,northPoleInnerGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs(ellipsoidOuter,ellipsoidOuterGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs(southPoleOuter,southPoleOuterGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs(northPoleOuter,northPoleOuterGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
#
# Here is the INNER background box
#
Box
  set corners
    # Choose bounds on inner background grid 
    $delta = 2*$ds; 
    $xai=-$a+$delta; $xbi=$a-$delta;
    $yai=-$b+$delta; $ybi=$b-$delta;
    $zai=-$c+$delta; $zbi=$c-$delta; 
    $xai $xbi $yai $ybi $zai $zbi
  lines
    $nx = intmg( ($xbi-$xai)/$ds +1.5);
    $ny = intmg( ($ybi-$yai)/$ds +1.5);
    $nz = intmg( ($zbi-$zai)/$ds +1.5);
    $nx $ny $nz
  boundary conditions
    0 0 0 0 0 0
  mappingName
    backGroundInner
  exit
#**********************************
exit
#
generate an overlapping grid
   backGround
   ellipsoidOuterGrid
   southPoleOuterGrid
   northPoleOuterGrid
   # inside the ellipsoid: 
   backGroundInner
   ellipsoidInnerGrid
   southPoleInnerGrid
   northPoleInnerGrid 
 done
# 
  change parameters
# 
    specify a domain
      # --- outer domain name:
      outerDomain 
       # grids in the domain:
       backGround
       ellipsoidOuterGrid
       southPoleOuterGrid
       northPoleOuterGrid 
      done
    specify a domain
       # --- inner domain name:
      innerDomain 
       # grids in the domain:
       backGroundInner
       ellipsoidInnerGrid
       southPoleInnerGrid
       northPoleInnerGrid 
      done
   interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    # 
    # open graphics
    # Maximum overlap to ensure the mask's match at the interface *fix me*
    if( $overlapOption eq "maximize" ){ $cmd="maximize overlap"; }else{ $cmd="#"; }
      $cmd 
  exit
#
  # open graphics  
#
  compute overlap
 #  change the plot
 #    toggle grid 0 0
    # open graphics
 #    pause
  # exit
#
exit
# save an overlapping grid
save a grid (compressed)
$name
solidEllipsoidInABoxGrid
exit


