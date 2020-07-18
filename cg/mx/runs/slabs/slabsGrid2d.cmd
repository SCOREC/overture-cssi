#
# Grid for multiple rectangular slabs in a 2d channel
#
#
# usage: ogen [-noplot] slabsGrid2d -factor=<num> -order=[2/4/6/8] -interp=[e/i] -ml=<> ...
#             -numSlabs=<i> -x1=<f> -x2=<f> -x3=<f> ...
#             -xa=<> -xb=<> -ya=<> -yb=<> -periodic=[n|p]
# 
#  -numSlabs : number of interior slabs  (current max =10) 
#  -x1, x2, x3, x4 : x locations of blocks.
#  -xa, -xb, -ya, -yb : bounds on the back ground grid
#  -angle : offset initial angles by this many degrees (future option)
#
#                                   
#           yb  +-----------------------------------+
#               |        |     |     |     |        |
#               |        |     |     |     |        |
#               |        |     |     |     |        |
#               |        |slab1|slab2|     |        |
#               |        |     |     |     |        |
#               |        |     |     |     |        |
#               |        |     |     |     |        |
#           ya  +-----------------------------------+
#               xa      x1     x2    x3    x4       xb
# 
# examples:
#  One interior slab:  
#   ogen -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=2 -interp=e -factor=2
#   ogen -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=4 -interp=e -factor=4
# 
#  Two interior slabs:  
#   ogen -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -order=2 -interp=e -factor=2
# 
#  Three interior slabs 
#   ogen -noplot slabsGrid2d -numSlabs=3 -x1=-1 -x2=-.5 -x3=0 -x4=.5 -order=2 -interp=e -factor=2
#
#
#
$prefix="slabsGrid2d";  $rgd="var"; $angle=0.; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
$periodic="n"; 
$order=2; $factor=1; $interp="i"; $ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-2.5; $xb=2.5; $ya=-.5; $yb=.5; 
$numSlabs=2; 
# Default block positions: 
$x1=-1; $x2=-.5; $x3=0; $x4=.5; $x5=1.; $x6=1.5; $x7=2.; $x8=2.5; $x9=3.; $x10=3.5; $x11=4.; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"ml=i"=>\$ml, "prefix=s"=> \$prefix,"numGhost=i"=> \$numGhost,\
            "numSlabs=i"=>\$numSlabs,"x1=f"=>\$x1,"x2=f"=>\$x2,"x3=f"=>\$x3,"x4=f"=>\$x4,"x5=f"=>\$x5,\
            "x6=f"=>\$x6,"x7=f"=>\$x7,"x8=f"=>\$x8,"x9=f"=>\$x9,"x10=f"=>\$x10,"x11=f"=>\$x11,"periodic=s"=> \$periodic );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$prefix = $prefix . "$numSlabs" . "Slabs";
if( $rgd eq "fixed" ){ $prefix = $prefix . "Fixed"; }
$suffix = ".order$order"; 
if( $periodic eq "p" ){ $suffix .= "p"; }
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
#
@xbv = ( $x1, $x2, $x3, $x4, $x5, $x6, $x7, $x8, $x9, $x10, $x11 );
# 
$ds=.1/$factor;
$pi = 4.*atan2(1.,1.);
# 
$dw = $order+1; $iw=$order+1; 
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
#
# ------- Left background grid -----
#
if( $periodic eq "p" ){ $bcTopBot = "-1 -1"; }else{ $bcTopBot = "3 4"; }
# 
rectangle
  set corners
    $xal=$xa; $xbl=$x1; 
    $xal $xbl $ya $yb
  lines
    $nx = intmg( ($xbl-$xal)/$ds +1.5 ); 
    $ny = intmg( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
  boundary conditions
    # 1 100 -1 -1 
    1 100 $bcTopBot 
  share 
    0  100 0 0
  mappingName
   leftBackGround
exit
#
# ------- Right background grid -----
#
rectangle
  set corners
    $xal=$xbv[$numSlabs]; $xbl= $xb;
    $xal $xbl $ya $yb
  lines
    $nx = intmg( ($xbl-$xal)/$ds +1.5 ); 
    $ny = intmg( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
  boundary conditions
    $share=100+$numSlabs; 
    # $share 2 -1 -1 
    $share 2 $bcTopBot
  share 
    $share 0 0 0
  mappingName
   rightBackGround
exit
##
# ------- Middle slabs  -----
#
$cmd="#"; $share=100; $grids="#"; 
for( $i=0; $i<$numSlabs; $i++ ){\
 $shareLeft=100+$i; $shareRight=$shareLeft+1; \
 $cmd .= "\n rectangle"; \
 $cmd .= "\n  set corners"; \
    $xal=$xbv[$i]; $xbl=$xbv[$i+1]; \
    $cmd .= "\n $xal $xbl $ya $yb"; \
  $cmd .= "\n lines"; \
    $nx = intmg( ($xbl-$xal)/$ds +1.5 );  \
    $ny = intmg( ($yb -$ya )/$ds +1.5 );  \
    $cmd .= "\n $nx $ny"; \
  $cmd .= "\n boundary conditions"; \
  $cmd .= "\n   $shareLeft $shareRight $bcTopBot"; \
  $cmd .= "\n share "; \
  $cmd .= "\n   $shareLeft $shareRight 0 0"; \
  $cmd .= "\n  mappingName"; \
  $cmd .= "\n  block$i"; \
  $cmd .= "\n exit"; \
  $grids .= "\n block$i"; \
}
#
# -- Now build the blocks: 
$cmd 
#
exit
#
#  --- generate the overlapping grid ---
#
generate an overlapping grid
    leftBackGround
    $grids
    rightBackGround
  done
  change parameters
    # -- Define domains from left to right ---
    specify a domain
      leftDomain
        leftBackGround
    done
    # create commands to define a separate domain for each slab
    $cmd="#"; 
    for( $i=0; $i<$numSlabs; $i++ ){\
      $cmd .="\n specify a domain"; \
      $cmd .="\n blockDomain$i "; \
      $cmd .="\n   block$i"; \
      $cmd .="\n done"; \
    }
    $cmd 
    # 
    specify a domain
      rightDomain
        rightBackGround
    done
    # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ngp = $ng+1;
      $ngp $ngp $ngp $ngp $ngp $ngp
  exit
  # open graphics
  compute overlap
#*  display computed geometry
  exit
#
# save an overlapping grid
save a grid (compressed)
$name
slabsGrid
exit


