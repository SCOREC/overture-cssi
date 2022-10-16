#
# Smoothed polygon for a "washer"
#
#
# usage: ogen [-noplot] washer -factor=<num> -order=[2/4/6/8] -blf=<num> -ml=<>  -numGhost=<i> 
# 
#  -ml = number of (extra) multigrid levels to support
# 
# examples:
#     ogen -noplot washer -order=2 -factor=1
#     ogen -noplot washer -order=4 -factor=1
#     ogen -noplot washer -order=4 -factor=2
# 
#
$prefix="washer";  $rgd="var";
 $bcSquare="d"; # old way
 $periodic="";  # new way 
 $dw=""; $iw=""; 
$order=2; $factor=1; $interp="e"; $ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "explicit for all grids";
$name=""; 
$blf=1;  # this means no stretching
$numGhost=-1;  # if this value is set, then use this number of ghost points
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"ml=i"=>\$ml,"blf=f"=> \$blf, "prefix=s"=> \$prefix,\
            "cx=f"=>\$cx,"cy=f"=>\$cy,"rgd=s"=> \$rgd,"bcSquare=s"=>\$bcSquare,"numGhost=i"=>\$numGhost,\
            "iw=i"=>\$iw,"dw=i"=>\$dw,"periodic=s"=>\$periodic );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=3; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=4; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
if( $rgd eq "fixed" ){ $prefix = $prefix . "Fixed"; }
# if( $bcSquare eq "p" ){ $prefix = $prefix . "p"; }
$suffix=""; 
if( $periodic eq "p" ){ $suffix = "p"; }
if( $periodic eq "np" ){ $suffix = "np"; }
if( $periodic eq "pn" ){ $suffix = "pn"; }
if( $iw eq "" ){ $suffix .= ".order$order"; }else{ $suffix .= ".Iw$iw" . "Dw$dw" . "$bc"; }
# $suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $blf ne 1 ){ $suffix .= ".s$blf"; }
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
$pi = 4.*atan2(1.,1.);
# 
if( $dw eq "" ){ $dw = $order+1; $iw=$order+1; }
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
create mappings
  smoothedPolygon
    lines
    $nTheta = intmg( 7.5/$ds + 1 );
    $nr = intmg( 1./$ds + 1 );
    $nTheta $nr
  mappingName
    washer
  exit
#
exit
generate an overlapping grid
    washer
  done
  change parameters
    # choose implicit or explicit interpolation
    interpolation type
      $interpType
    #
    # order of accuracy 
    #   $orderOfAccuracy
    # -- set the discretization width and interpolation width --
    $cmd =" order of accuracy\n $orderOfAccuracy";
    if( $dw ne "" ){ $cmd="discretization width\n all\n $dw $dw $dw\n interpolation width\n all\n all\n $iw $iw $iw"; }
    $cmd
    #
    ghost points
      all
      # $ngp = $ng+1;
      $ngp = $ng;
      $ng $ng $ng $ngp $ng $ng
  exit
#   display intermediate results
  compute overlap
# plot
#   query a point 
#     interpolate point 1
#     check interpolation coords 1
#     pt: grid,i1,i2,i3: 1 5 6 0
# 
#*  display computed geometry
  exit
#
# save an overlapping grid
save a grid (compressed)
# printf(" name=$name\n");
$name
washer
exit

