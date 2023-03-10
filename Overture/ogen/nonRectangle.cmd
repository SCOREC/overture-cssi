#
# Non-Rectangle 
#
#
# usage: ogen [noplot] nonRectangle -factor=<num> -order=[2/4/6/8] -xa= -xb= -ya= -yb= -prefix=<> -name=<>
# 
# examples:
#    ogen -noplot nonRectangle -prefix=nonRect6x2y -order=2 -xa=-3. -xb=3. -ya=-1. -yb=1. -factor=1
#
#
$prefix="nonRectangle"; $xa=-1.; $xb=1.; $ya=-1.; $yb=1.;
$order=2; $factor=1; # default values
$orderOfAccuracy = "second order"; $ng=2;  $periodic=""; $name=""; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "periodic=s"=>\$periodic,"name=s"=>\$name,"prefix=s"=>\$prefix );
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
if( $name eq "" ){$name = $prefix . $factor . $suffix . ".hdf";}
# 
$ds=.1/$factor;
# 
create mappings
#
rectangle
  set corners
    $xa $xb $ya $yb
  lines
    $nx = int( ($xb-$xa)/$ds +1.5 );
    $ny = int( ($yb-$ya)/$ds +1.5 );
    $nx $ny
  boundary conditions
     if( $periodic eq "p" ){ $bc ="-1 -1 -1 -1"; }\
     elsif( $periodic eq "np" ){ $bc ="1 2 -1 -1"; }\
     elsif( $periodic eq "pn" ){ $bc ="-1 -1 3 4"; }else{ $bc="1 2 3 4"; }
    $bc
  mappingName
    rectangleCartesian
  exit
  rotate/scale/shift
    mappingName
    rectangle
  exit
#
exit
generate an overlapping grid
    rectangle
  done
  change parameters
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
  exit
#  display intermediate results
  compute overlap
# 
  display computed geometry
  exit
#
save an overlapping grid
$name
rectangle
exit

