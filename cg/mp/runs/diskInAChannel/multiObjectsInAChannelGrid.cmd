#
# Grid for multiple objects in a channel
#
#
# usage: ogen [noplot] multiObjectsInAChannelGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -gridOption=<num>
# 
# examples:
#   ogen multiObjectsInAChannelGrid -prefix=multiObjectsGrid -interp=e -factor=2 -gridOption=3
#
#
$prefix="multiObjectsInAChannelGrid";
$order=2; $factor=1; $interp="e"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $ml=0; 
$name=""; 
$xa=0; $xb=10; $ya=-2.0; $yb=2.0
$numberOfVolumeSmooths=0; 
$gridOption=3;
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"ra=f"=> \$ra,"rb=f"=> \$rb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"per=i"=>\$per,\
            "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "prefix=s"=>\$prefix,"numDisks=i"=>\$numDisks,\
            "gridOption=i"=>\$gridOption );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $per eq 1 ){ $suffix .= "p"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
$pi = 4.*atan2(1.,1.);
# 
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
#
# $bcInterface=100;  # bc for interfaces
# $shareInterface=100;        # share value for interfaces
#
$fluidDomain="";
$solidDomain="";
$specifyDomainCmd="";
create mappings
#
#
$count=0;  # counts number of objects
if ($gridOption eq 0) {$cmd = "include makeObjectsGrid0.h";}
if ($gridOption eq 1) {$cmd = "include makeObjectsGrid1.h";}
if ($gridOption eq 2) {$cmd = "include makeObjectsGrid2.h";}
if ($gridOption eq 3) {$cmd = "include makeObjectsGrid3.h";}
$cmd
#
#  --- fluid Channel Grid ---
#
 rectangle
  set corners
    $xa $xb $ya $yb
  lines
    $nx = int( ($xb-$xa)/$ds + .5 ); 
    $ny = int( ($yb-$ya)/$ds + .5 ); 
    $nx $ny
  boundary conditions
    1 2 3 4 
  share 
    0 0 0 0 
  mappingName
    fluidChannel
  exit
#
exit this menu
#
$fluidDomain=$fluidDomain . "fluidChannel" . "\n";
$specifyDomainCmd=$specifyDomainCmd \
  . "specify a domain\nfluidDomain\n" \
  . $fluidDomain \
  . "done\n";
$allDomains=$fluidDomain . $solidDomain;
generate an overlapping grid
  $allDomains
  done choosing mappings
  change parameters
    $specifyDomainCmd
    # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      2 2 2 2 2 2
    exit
  # open graphics
  compute overlap
exit
#
# save an overlapping grid
save a grid (compressed)
$name
exit
exit
