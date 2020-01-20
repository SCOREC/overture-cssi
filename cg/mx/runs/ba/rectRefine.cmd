#
# Rectangle with refinements 
#
#
# usage: ogen [noplot] rectRefine -factor=<num> -order=[2/4/6/8] -xa= -xb= -ya= -yb= -prefix=<> -name=<> -periodic=[p|np|pn]
# 
# examples:
#    ogen -noplot rectRefine -factor=1 -order=2 
#    ogen -noplot rectRefine -factor=2 -order=2 
#
# 
#
$prefix="rectRefine";
$xa=-1.; $xb=1.; $ya=-1.; $yb=1.;
$xar=-.5; $xbr=.5; $yar=-.5; $ybr=.5;
$order=2; $factor=1; # default values
$orderOfAccuracy = "second order"; $ng=2;  $periodic=""; $name=""; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,\
            "xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "xar=f"=> \$xar,"xbr=f"=> \$xbr,"yar=f"=> \$yar,"ybr=f"=> \$ybr,\
            "periodic=s"=>\$periodic,"name=s"=>\$name,"prefix=s"=>\$prefix,"numGhost=i"=> \$numGhost );
# == add an extra ghost 
if( $order eq 2 ){ $orderOfAccuracy="second order"; $ng=2; }\
elsif( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=5; }
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
    $nx = int( ($xb-$xa)/$ds +1.5 );
    $ny = int( ($yb-$ya)/$ds +1.5 );
    $nx $ny
  boundary conditions
     if( $periodic eq "p" ){ $bc ="-1 -1 -1 -1"; }\
     elsif( $periodic eq "np" ){ $bc ="1 2 -1 -1"; }\
     elsif( $periodic eq "pn" ){ $bc ="-1 -1 3 4"; }else{ $bc="1 2 3 4"; }
    $bc
  mappingName
    rectangle
exit
#     refinement 
rectangle
  set corners
    $xar $xbr $yar $ybr
  lines
    $nx = int( ($xbr-$xar)/$ds +1.5 );
    $ny = int( ($ybr-$yar)/$ds +1.5 );
    $nx $ny
  boundary conditions
    0 0 0 0 
  mappingName
    refinement1
exit
#
exit
generate an overlapping grid
    rectangle
    refinement1
  done
  change parameters
    # 
    #discretization width
    #  all
    #  5 5 5
    #interpolation width
    #  all
    #  all
    #  2 2 2 
#
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ngp=$ng+1; 
      $ng $ng $ng $ngp $ng $ng 
  exit
  # open graphics
  compute overlap
# 
#  display computed geometry
  exit
#
save an overlapping grid
$name
rectRefine
exit

