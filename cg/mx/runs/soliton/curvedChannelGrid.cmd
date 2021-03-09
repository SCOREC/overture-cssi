#
# Grid for a curved channel
#
#
# usage: ogen [noplot] rectangleArg -factor=<num> -order=[2/4/6/8] -xa= -xb= -ya= -yb= -prefix=<> -name=<> -periodic=[p|np|pn]
# 
# examples:
#    ogen -noplot curvedChannelGrid -order=2 -factor=2
#    ogen -noplot curvedChannelGrid -order=2 -factor=4
#
#
$prefix="curvedChannel"; $xa=-1.; $xb=1.; $ya=-1.; $yb=1.;
$Nx=-1; # if set, use this many grid points in x 
$adjustCC=0; # apply cell-=centered adjustment for painted interface 
$order=2; $factor=1; $ds0=.1; # default values
$orderOfAccuracy = "second order"; $ng=2;  $periodic=""; $name=""; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
$extraLines=0; 
# 
# get command line arguments
# Getopt::Long::Configure("prefix_pattern=(--rectangleArg|--|-)");
GetOptions( "order=i"=>\$order,"factor=f"=>\$factor,"xa=f"=>\$xa,"xb=f"=>\$xb,"ds0=f"=>\$ds0,\
            "ya=f"=>\$ya,"yb=f"=>\$yb,"ybx=f"=>\$ybx,\
            "periodic=s"=>\$periodic,"name=s"=>\$name,"prefix=s"=>\$prefix,"numGhost=i"=> \$numGhost, \
	    "extraLines=i"=> \$extraLines,"Nx=i"=>\$Nx,"adjustCC=i"=>\$adjustCC );
# printf("rectangleArg: factor=$factor xa=$xa xb=$xb ya=$ya yb=$yb ybx=$ybx\n");
# pause
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix="";
if( $periodic eq "p" ){ $suffix = "p"; }
if( $periodic eq "np" ){ $suffix = "np"; }
if( $periodic eq "pn" ){ $suffix = "pn"; }
$suffix .= ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; }
if( $name eq "" ){$name = $prefix . $factor . $suffix . ".hdf";}
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
# 
$ds=$ds0/$factor;
# 
create mappings
#
# printf("ds=$ds\n");
include curvedChannel.h
#
 mapping from normals
  extend normals from which mapping?
  curvedChannelCurve
  $nDist=.1; 
  normal distance
    $nDist
  lines
    $ns =int($arcLength/$ds + 1.5 );
    $nr =int($nDist/$ds + 2.5 );
    $nr = max($nr,3+$ng); 
    $ns $nr
  boundary conditions
    1 2 3 4
  mappingName
    curvedChannel
  exit
#
exit
generate an overlapping grid
    curvedChannel
  done
  change parameters
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ngp=$ng+1; 
      $ng $ng $ng $ngp $ng $ng 
  exit
#  display intermediate results
  compute overlap
  # open graphics
# 
  # display computed geometry
  exit
#
save an overlapping grid
$name
curvedChannel
exit

