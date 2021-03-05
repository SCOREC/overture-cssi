#
# Rectangle (taking arguments)
#
#
# usage: ogen [noplot] rectangleArg -factor=<num> -order=[2/4/6/8] -xa= -xb= -ya= -yb= -prefix=<> -name=<> -periodic=[p|np|pn]
# 
# examples:
#    ogen -noplot rectangleArg -prefix=rect6x2y -order=2 -xa=-3. -xb=3. -ya=-1. -yb=1. -factor=1
#
#    ogen -noplot rectangleArg -factor=1 -order=2 -xa=-2. -xb=2. -ya=-2. -yb=2. -name="rect4x4y1.hdf"
#    ogen -noplot rectangleArg -factor=2 -order=2 -xa=-2. -xb=2. -ya=-2. -yb=2. -name="rect4x4y2.hdf"
#    ogen -noplot rectangleArg -factor=2 -order=2 -xa=-10. -xb=10. -ya=-10. -yb=10. -name="rect20x20y2.hdf"
# 
#    ogen -noplot rectangleArg -factor=4 -order=2 -xa=0. -xb=2. -ya=0. -yb=1. -name="rect2x1y4.hdf"
#  -- square with bottom at y=0 for axisymetric problems:
#    ogen -noplot rectangleArg -factor=2 -order=2 -xa=-0. -xb=1. -ya=0. -yb=1. -name="axiSquare2.order2.hdf"
#    ogen -noplot rectangleArg -factor=2 -order=2 -xa=-0. -xb=1. -ya=.2 -yb=1.2 -name="axiSquare2a.order2.hdf"
#    ogen -noplot rectangleArg -factor=4 -order=2 -xa=-0. -xb=1. -ya=.2 -yb=1.2 -name="axiSquare4a.order2.hdf"
#    ogen -noplot rectangleArg -factor=8 -order=2 -xa=-0. -xb=1. -ya=.2 -yb=1.2 -name="axiSquare8a.order2.hdf"
#
$prefix="rectangle"; $xa=-1.; $xb=1.; $ya=-1.; $yb=1.;
$Nx=-1; # if set, use this many grid points in x 
$adjustCC=0; # apply cell-=centered adjustment for painted interface 
$order=2; $factor=1; # default values
$orderOfAccuracy = "second order"; $ng=2;  $periodic=""; $name=""; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
$extraLines=0; 
# 
# get command line arguments
# Getopt::Long::Configure("prefix_pattern=(--rectangleArg|--|-)");
GetOptions( "order=i"=>\$order,"factor=f"=>\$factor,"xa=f"=>\$xa,"xb=f"=>\$xb,\
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
# 
$ds=.1/$factor;
# 
create mappings
#
rectangle
  set corners
    $xa $xb $ya $yb
  lines
    $nx = int( ($xb-$xa)/$ds +1.5 + $extraLines );
    $ny = int( ($yb-$ya)/$ds +1.5 + $extraLines );
    if( $Nx > 0 ){ $nx=$Nx; }
    # Adjust $nx so that points x=-.5 and x=.5 are at cell centers if $xb=-$xa and $xb = integer (for bamx interface code)
    if( $adjustCC ne 0 ){ $nx= int( ($xb-$xa)*( 10*$factor+1 ) + 1.5 ); }
    $nx $ny
  boundary conditions
     if( $periodic eq "p" ){ $bc ="-1 -1 -1 -1"; }\
     elsif( $periodic eq "np" ){ $bc ="1 2 -1 -1"; }\
     elsif( $periodic eq "pn" ){ $bc ="-1 -1 3 4"; }else{ $bc="1 2 3 4"; }
    $bc
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
      $ngp=$ng+1; 
      $ng $ng $ng $ngp $ng $ng 
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

