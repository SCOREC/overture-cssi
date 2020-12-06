#**************************************************************************
#
#  A grid for an interface calculation between two boxes
#
# Usage: 
#         ogen [noplot] twoBoxesInterface [options]
# where options are
#     -factor=<num>      : use this factor for all directions (if given)
#     -xFactor=<num>     : default grid spacing in x-direction is multiplied by this factor
#     -yFactor=<num>     : default grid spacing in y-direction is multiplied by this factor
#     -zFactor=<num>     : default grid spacing in z-direction is multiplied by this factor
#     -xFactorScale      : scalee grid spacing by this factro in x-direction (e.g. to make finer)
#     -order=[2/4/6/8]  : order of accuracy 
#     -interp=[e/i]     : implicit or explicit interpolation
#     -bc=[d|p]
# Examples: 
#    ogen -noplot twoBoxesInterface -xFactor=1 -order=2 -interp=e    (creates twoBoxesInterfacee111.order2.hdf)
#    ogen -noplot twoBoxesInterface -factor=1 -order=2 -interp=e
#    ogen -noplot twoBoxesInterface -factor=2 -order=2 -interp=e
#    ogen -noplot twoBoxesInterface -factor=4 -order=2 -interp=e
#    ogen -noplot twoBoxesInterface -factor=8 -order=2 -interp=e
# 
#    ogen -noplot twoBoxesInterface -order=4 -interp=e -factor=1
#    ogen -noplot twoBoxesInterface -order=4 -interp=e -factor=2
#    ogen -noplot twoBoxesInterface -order=4 -interp=e -factor=4
#    ogen -noplot twoBoxesInterface -order=4 -interp=e -factor=8
#    ogen -noplot twoBoxesInterface -order=4 -interp=e -factor=16
# -- periodic
#    ogen -noplot twoBoxesInterface -order=4 -interp=e -bc=p -factor=2
# -- rotated:
#    ogen -noplot twoBoxesInterface -factor=1 -order=2 -angle=45 -name="twoBoxesInterfaceRotated1.order2.hdf" 
#    ogen -noplot twoBoxesInterface -factor=1 -order=4 -angle=45 -name="twoBoxesInterfaceRotated1.order4.hdf" 
# 
#**************************************************************************
$prefix="twoBoxesInterface"; $name=""; 
$order=2; $bc="d"; 
$factor=-1; $xFactor=1; $yFactor=1; $zFactor=1; 
$xFactorScale=1.; # scale number of grid lines in x by this factor
$interp="i"; $interpType = "implicit for all grids";
$angle=0;  $rotAxis=1;   # rotation axis 0, 1, or 2 
$angle2=0; $rotaAxis2=2; # optional second rotation
#
$bcLeft = "1 100  1  1  1 1"; 
$bcRight= "100 1  1  1  1 1"; 
# 
$bclp = "1 100 -1 -1  -1 -1 ";  # periodic in the y/z-direction
$bcrp= "100 1 -1 -1  -1 -1 ";  # periodic in the y/z-direction
# 
$xa=-1;    $xb=1.; 
$ya= 0.;   $yb=.5; 
$za= 0.;   $zb=.5; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
#
# 
# get command line arguments
GetOptions("order=i"=>\$order,"factor=i"=> \$factor,"xFactor=i"=> \$xFactor,"yFactor=i"=> \$yFactor,"zFactor=i"=> \$zFactor,\
           "interp=s"=> \$interp,\
           "name=s"=>\$name,"bc=s"=>\$bc,"numGhost=i"=>\$numGhost,"prefix=s"=>\$prefix,"xFactorScale=f"=>\$xFactorScale,\
	   "angle=f"=> \$angle,"rotAxis=i"=> \$rotAxis,"angle2=f"=> \$angle2,"rotAxis2=i"=> \$rotAxis2,\
           "xa=f"=> \$xa,"xb=f"=> \$xb );
# 
if( $factor>0 ){ $xFactor=$factor; $yFactor=$factor; $zFactor=$factor; }
if( $order eq 2 ){ $orderOfAccuracy="second order"; $ng=2; }\
elsif( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
# 
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $bc eq "p" ){ $suffix .= "p"; } # periodic
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $factor>0 && $name eq "" ){ $name = $prefix . "$interp$factor" . $suffix . ".hdf"; }
if( $name eq "" ){ $name = $prefix . "$interp$xFactor$yFactor$zFactor" . $suffix . ".hdf"; }
# printf("name=[$name]\n"); 
if( $bc eq "p" ){ $bcLeft=$bclp; $bcRight=$bcrp; }
# 
# 
#
#
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
#
#
#**************************************************************************
#
#
# domain parameters:  
$dsx= $xFactorScale*(.1/$xFactor);       # target grid spacing in the x direction
$dsy= .1/$yFactor;                       # target grid spacing in the y direction
$dsz= .1/$zFactor;                       # target grid spacing in the y direction
#
create mappings
#
#  here is the left grid
#
  Box
    $xar=$xa;  $xbr=0.0; 
    set corners
     $xar $xbr $ya $yb $za $zb
    lines
      $nx=int( ($xbr-$xar)/$dsx+1.5 );
      $ny=int( ($yb-$ya)/$dsy+1.5 );
      $nz=int( ($zb-$za)/$dsz+1.5 );
      $nx $ny $nz
    boundary conditions
      $bcLeft
    share
 # for now interfaces are marked with share>=100 
      0 100 0 0 0 0
    mappingName
      leftBox0
  exit
#
#
  Box
    $xar= 0.0;  $xbr=$xb;
    set corners
     $xar $xbr $ya $yb $za $zb
    lines
      $nx=int( ($xbr-$xar)/$dsx+1.5 );
      $nx $ny $nz
    boundary conditions
      $bcRight
    share
 # for now interfaces are marked with share>=100 
      100 0 0 0 0 0 
    mappingName
      rightBox0
  exit
#
#
  rotate/scale/shift
    transform which mapping?
      leftBox0
    rotate
     # rotate about the x, y, or z-axis
     $angle $rotAxis 
     $yr = ($ya+$yb)*.5; 
     $zr = ($za+$zb)*.5; 
     0. $yr $zr
    # optionally rotate a second time
    if( $angle2 ne 0 ){ $rotateCommand ="rotate\n $angle2 $rotAxis2\n  0. $yr $zr"; }else{ $rotateCommand="#"; }
    $rotateCommand
    # pause
    mappingName
      leftBox
    exit
#
  rotate/scale/shift
    transform which mapping?
      rightBox0
    rotate
     $angle $rotAxis 
     0. $yr $zr
    # rotate a second time
   $rotateCommand
    mappingName
      rightBox
    exit
#
  exit this menu
#
generate an overlapping grid
  leftBox
  rightBox
  done 
  change parameters
 # define the domains -- these will behave like independent overlapping grids
    specify a domain
 # domain name:
      leftDomain 
 # grids in the domain:
      leftBox
      done
    specify a domain
 # domain name:
      rightDomain 
 # grids in the domain:
      rightBox
      done
    order of accuracy
     $orderOfAccuracy
 # choose implicit or explicit interpolation
    interpolation type
      $interpType
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
  exit
# pause
  compute overlap
# pause
  exit
#
save an overlapping grid
  $name
  twoBoxesInterface
exit
