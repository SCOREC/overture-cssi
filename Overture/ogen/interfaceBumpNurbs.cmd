#**************************************************************************
#
#  Build grids for an interface with one or more bumps (for cgmx)
#   ** This version reads the interface grids generated by interfaceBumpHype.cmd ****
#      This version can then be run in parallel to generate very fine grids 
#
# usage: ogen [noplot] interfaceBumpNurbs -factor=<num> -order=[2/4/6/8] -interp=[e/i] -surfaceGrids=<name> ...
#    -factor0=<>
#
# -factor0 : surface grids were made with this factor value
#
# examples:
#
#  ogen -noplot interfaceBumpNurbs -interp=e -order=4 -curve=afm3d3 -surfaceGrids=afm3d3SurfaceGrids.hdf -factor0=4 -factor=4 -ds0=.4 -za=-14. -zb=6
#
#  ogen -noplot interfaceBumpNurbs -interp=e -order=4 -curve=afm3d3 -surfaceGrids=afm3d3SurfaceGrids.hdf -factor0=4 -factor=8 -ds0=.4 -za=-14. -zb=6
# OK: 140M pts
# srun -N2 -n8 -ppdebug $ogenp -noplot interfaceBumpNurbs -interp=e -order=4 -curve=afm3d3 -surfaceGrids=afm3d3SurfaceGrids.hdf -factor0=4 -factor=8 -ds0=.4 -za=-14. -zb=6
# 
# Finer grid: 1.1B pts
# submit.p -jobName="cgmx" -bank=windpowr -out="interfaceBumpNurbs16.out" -walltime=2:00 -submit=0 -cmd='srun -N16 -n64 -ppbatch $ogenp -noplot interfaceBumpNurbs -interp=e -order=4 -curve=afm3d3 -surfaceGrids=afm3d3SurfaceGrids.hdf -factor0=4 -factor=16 -ds0=.4 -za=-14. -zb=6'
#
# Higher domain: 
# srun -N2 -n8 -ppdebug $ogenp -noplot interfaceBumpNurbs -interp=e -order=4 -curve=afm3d3 -surfaceGrids=afm3d3SurfaceGrids.hdf -factor0=4 -factor=8 -ds0=.4 -za=-14. -zb=11 -name="interfaceBumpNurbsafm3d3z11e8.order4.hdf"
#
# Todo: (batch needed)
# srun -N4 -n16 -ppdebug $ogenp -noplot interfaceBumpNurbs -interp=e -order=4 -curve=afm3d3 -surfaceGrids=afm3d3SurfaceGrids.hdf -factor0=4 -factor=16 -ds0=.4 -za=-14. -zb=6
#
#
# -- Bigger central patch: (increase -zb=8 -> -zb=12)
#  ogen -noplot interfaceBumpNurbs -interp=e -order=4 -curve=m588MidPlus -surfaceGrids=m588MidPlusSurfaceGrids.hdf -factor0=4 -factor=4 -ds0=.4 -za=-18. -zb=8.
#  srun -N2 -n8 -ppdebug $ogenp -noplot interfaceBumpNurbs -interp=e -order=4 -curve=m588MidPlus -surfaceGrids=m588MidPlusSurfaceGrids.hdf -factor0=4 -factor=4 -ds0=.4 -za=-18. -zb=12.
# 
# finer: batch needed: cpu>30min  : 610M pts CPU=2.96e+03
#  srun -N8 -n16 -ppdebug $ogenp -noplot interfaceBumpNurbs -interp=e -order=4 -curve=m588MidPlus -surfaceGrids=m588MidPlusSurfaceGrids.hdf -factor0=4 -factor=8 -ds0=.4 -za=-18. -zb=12.
# 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $extra=0; $curve="1bump"; $surfaceGrids="afm3d3SurfaceGrids.hdf";
$xa=-1.; $xb=1.; $ya=-1.; $yb=1.; $za=-.5; $zb=.5; $ds0=.05; $factor0=4; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
        "za=f"=> \$za,"zb=f"=> \$zb,"interp=s"=> \$interp,"name=s"=> \$name,"curve=s"=> \$curve,\
        "extra=i"=>\$extra,"ds0=f"=> \$ds0,"surfaceGrids=s"=> \$surfaceGrids,"factor0=f"=> \$factor0 );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; $extra=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; $extra=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; $extra=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; $extra=$extra+$order; }
# *wdh* 090409 if( $factor > 4 ){ $extra=$extra+8; }  # make interface grids a bit wider for higher resolution cases
# 
$suffix = ".order$order"; 
$curveName = ""; 
if( $curve eq "afm1" ){ $curveName="One"; }elsif( $curve eq "flat" ){ $curveName="Flat"; }else{ $curveName=$curve; }
if( $name eq "" ){ $name = "interfaceBumpNurbs$curveName" . "$interp$factor" . $suffix . ".hdf";}
# 
# domain parameters:  
$ds = $ds0/$factor; # target grid spacing
# 
# parallel ghost lines: for ogen we need at least:
#       .5*( iw -1 )   : implicit interpolation 
#       .5*( iw+dw-2 ) : explicit interpolation
$dw = $order+1; $iw=$order+1;
$parallelGhost=($iw-1)/2;
if( $interp eq "e" ){  $parallelGhost=($iw+$dw-2)/2; }
if( $parallelGhost<1 ){ $parallelGhost=1; } 
minimum number of distributed ghost lines
  $parallelGhost
#
#
$bcInterface=100;  # bc for interfaces
$ishare=100;
#
create mappings 
if( $curve eq "afm3d1" ){ $amp=0.; $afm3dSurface="afm.smallMiddlePatch.dat"; }
# if( $curve eq "afm3d2" ){ $amp=0.; $afm3dSurface="/home/henshaw.0/cgDoc/nif/afm/afm.upperMiddleRight.dat"; }
if( $curve eq "afm3d2" ){ $amp=0.; $afm3dSurface="afm.upperMiddleRight.dat"; }
if( $curve eq "afm3d3" ){ $amp=0.; $afm3dSurface="afm.M588Mid.dat"; }
if( $curve eq "afm3d1" || $curve eq "afm3d2" || $curve eq "afm3d3" ){ $cmds = "include afm3d.cmd"; }
# Bigger middle section:
if( $curve eq "m588MidPlus" ){ $amp=0.; $afm3dSurface="afm.M588MidPlus.dat"; $cmds = "include afm3d.cmd"; }
$cmds 
#
#* -- nurbs --
# 
$pi =4.*atan2(1.,1.); $cmd="";
create mappings 
# 
  $curveFactor=1.1; # increase points from a flat surface
  $nx = int( $curveFactor*($xb-$xa)/$ds +1.5);
  $ny = int( $curveFactor*($yb-$ya)/$ds +1.5);
#
# -- fix me: this assumes the surface grids were saved with factor=4 and nr=10:
  #$nr0=10;    # number of lines in the normal direction in the file
  #$factor0=4; # file "factor"
  $nr0=5+$extra; 
  $nr = int( ($nr0-1)*$factor/$factor0 + 1.5);
#
open a data-base
  $surfaceGrids
  open an old file read-only
  get all mappings from the data-base
  close the data-base
#
#  -- set the number of lines for this resolution: 
  change a mapping
  lowerInterface
    lines
      $nx $ny $nr
  exit
#
  change a mapping
  upperInterface
    lines
      $nx $ny $nr
  exit
#
  $xar=$xa; $xbr=$xb; $yar=$ya; $ybr=$yb; $zar=$za; $zbr=$zMax; 
Box
    mappingName
      lower
  set corners
    $xar $xbr $yar $ybr $zar $zbr
  lines
    $nx = int( ($xbr-$xar)/$ds +1.5);
    $ny = int( ($ybr-$yar)/$ds +1.5);
    $nz = int( ($zbr-$zar)/$ds +1.5);
    $nx $ny $nz
    boundary conditions
      1 2 3 4 5 0 
    share
      1 2 3 4 0 0 
  exit
#
  $xar=$xa; $xbr=$xb; $yar=$ya; $ybr=$yb; $zar=$zMin; $zbr=$zb; 
Box
    mappingName
      upper
  set corners
    $xar $xbr $yar $ybr $zar $zbr
  lines
    $nx = int( ($xbr-$xar)/$ds +1.5);
    $ny = int( ($ybr-$yar)/$ds +1.5);
    $nz = int( ($zbr-$zar)/$ds +1.5);
    $nx $ny $nz
    boundary conditions
      1 2 3 4 0 6
    share
      1 2 3 4 0 0
  exit
#
exit
#
generate an overlapping grid 
  lower
  upper
  lowerInterface
  upperInterface
  done 
#
  change parameters
    specify a domain
 # domain name:
      lowerDomain 
 # grids in the domain:
      lower
      lowerInterface
      done
    specify a domain
 # domain name:
      upperDomain 
 # grids in the domain:
      upper
      upperInterface
      done
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    exit
# 
  # open graphics
  compute overlap
# 
  exit
#
maximum number of parallel sub-files
  4
save an overlapping grid
$name
interfaceBump
exit
