#
# Grid for two elastic disks in a channel
#
#
# usage: ogen [noplot] multiDisksInAChannelGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] 
# Options:
#    -ra,rb
# 
# examples:
#   ogen -noplot multiDisksInAChannelGrid -numDisks=2 -prefix=twoDisksGrid -interp=e -factor=2
#
#
$prefix="diskInAChannelGrid";
$order=2; $factor=1; $interp="e"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $ml=0; 
$name=""; 
$xa=-3; $xb=3; $ya=-1.75; $yb=1.75; 
$numDisks=2; # so far only 2 allowed 
# Disk 1:
$ra1=.4; $rb1=.8;   # disk is hollow 
$cx1=-1; $cy1=-.5; # centre of disk 1
# Disk 2:
$ra2=.4; $rb2=.8;   # disk is hollow 
$cx2= 1; $cy2= .5; # center of disk 2
$numberOfVolumeSmooths=0; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"ra=f"=> \$ra,"rb=f"=> \$rb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"per=i"=>\$per,\
            "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "prefix=s"=>\$prefix,"numDisks=i"=>\$numDisks );
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
#
# $bcInterface=100;  # bc for interfaces
# $shareInterface=100;        # share value for interfaces
#
create mappings
#
#
$count=0;  # counts number of disks
#
#  ------ Disk 1 ---------
#
$ra=$ra1; $rb=$rb1; $cx=$cx1; $cy=$cy1; 
include makeDisk.h
#
#  ------ Disk 2 ---------
#
$ra=$ra2; $rb=$rb2; $cx=$cx2; $cy=$cy2; 
include makeDisk.h
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
generate an overlapping grid
  solidDisk1
  solidDisk2
  fluidChannel
  fluidInterface1
  fluidInterface2
  done choosing mappings
  change parameters
    specify a domain
     solidDomain1
       solidDisk1
    done
    specify a domain
     solidDomain2
       solidDisk2
    done
    specify a domain
     fluidDomain
       fluidChannel
       fluidInterface1
       fluidInterface2
    done
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
multiDisks
exit

