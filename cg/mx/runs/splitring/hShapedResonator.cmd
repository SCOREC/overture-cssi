#
# usage: ogen [-noplot] hShapedResonator -factor=<num> -order=[2/4/6/8] -interp=[e/i] 
#
# draw an H shaped resonator 
#
$order=2; $factor=1; $interp="e";  $ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "explicit for all grids";
$prefix="hShapedResonator"; $name=""; 
$H=1.; $W=.5; $T=.1; $d=$T/3.;
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"interp=s"=> \$interp,"name=s"=> \$name,\
            "prefix=s"=> \$prefix,"H=f"=> \$H,"W=f"=>\$W,"T=f"=>\$T);
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $name eq "" ){$name = $prefix . "$factor" . $suffix . ".hdf";}
#
# 
$ds=.02/$factor;
$xa = -1.;
$xb =  1.;
$ya = -1.;
$yb =  1.; 
$x0 = 0.;
$y0 = 0.;
$pi = 4.*atan2(1.,1.);
$nx = ($xb-$xa)/$ds;
$ny = ($yb-$ya)/$ds; 
$nl = 100*int(1.5*$factor);
$sharp=100;
# 
# 
create mappings
#
  rectangle 
    set corners
      $xa $xb $ya $yb
    lines 
      $nx $ny
    boundary conditions
      -1 -1 -1 -1
    mappingName
      backGround 
  exit
  SmoothedPolygon
    vertices
      37
      $x=$x0+$W/2.       ;$y=$y0-$H/2.+$d ;
      $x $y
      $x=$x0+$W/2.       ;$y=$y0+$H/2.-$d ;
      $x $y
      $x=$x0+$W/2.       ;$y=$y0+$H/2.    ;
      $x $y
      $x=$x0+$W/2.-$d    ;$y=$y0+$H/2.    ;
      $x $y
      $x=$x0+$W/2.-$T+$d ;$y=$y0+$H/2.    ;
      $x $y
      $x=$x0+$W/2.-$T    ;$y=$y0+$H/2.    ;
      $x $y
      $x=$x0+$W/2.-$T    ;$y=$y0+$H/2.-$d ;
      $x $y
      $x=$x0+$W/2.-$T    ;$y=$y0+$T/2.+$d ;
      $x $y
      $x=$x0+$W/2.-$T    ;$y=$y0+$T/2.    ;
      $x $y
      $x=$x0+$W/2.-$T-$d ;$y=$y0+$T/2.    ;
      $x $y
      $x=$x0-$W/2.+$T+$d ;$y=$y0+$T/2.    ;
      $x $y
      $x=$x0-$W/2.+$T    ;$y=$y0+$T/2.    ;
      $x $y
      $x=$x0-$W/2.+$T    ;$y=$y0+$T/2.+$d ;
      $x $y
      $x=$x0-$W/2.+$T    ;$y=$y0+$H/2.-$d ;
      $x $y
      $x=$x0-$W/2.+$T    ;$y=$y0+$H/2.    ;
      $x $y
      $x=$x0-$W/2.+$T-$d ;$y=$y0+$H/2.    ;
      $x $y
      $x=$x0-$W/2.+$d    ;$y=$y0+$H/2.    ;
      $x $y
      $x=$x0-$W/2.       ;$y=$y0+$H/2.    ;
      $x $y
      $x=$x0-$W/2.       ;$y=$y0+$H/2.-$d ;
      $x $y
      $x=$x0-$W/2.       ;$y=$y0-$H/2.+$d ;
      $x $y
      $x=$x0-$W/2.       ;$y=$y0-$H/2.    ;
      $x $y
      $x=$x0-$W/2.+$d    ;$y=$y0-$H/2.    ;
      $x $y
      $x=$x0-$W/2.+$T-$d ;$y=$y0-$H/2.    ;
      $x $y
      $x=$x0-$W/2.+$T    ;$y=$y0-$H/2.    ;
      $x $y
      $x=$x0-$W/2.+$T    ;$y=$y0-$H/2.+$d ;
      $x $y
      $x=$x0-$W/2.+$T    ;$y=$y0-$T/2.-$d ;
      $x $y
      $x=$x0-$W/2.+$T    ;$y=$y0-$T/2.    ;
      $x $y
      $x=$x0-$W/2.+$T+$d ;$y=$y0-$T/2.    ;
      $x $y
      $x=$x0+$W/2.-$T-$d ;$y=$y0-$T/2.    ;
      $x $y
      $x=$x0+$W/2.-$T    ;$y=$y0-$T/2.    ;
      $x $y
      $x=$x0+$W/2.-$T    ;$y=$y0-$T/2.-$d ;
      $x $y
      $x=$x0+$W/2.-$T    ;$y=$y0-$H/2.+$d ;
      $x $y
      $x=$x0+$W/2.-$T    ;$y=$y0-$H/2.    ;
      $x $y
      $x=$x0+$W/2.-$T+$d ;$y=$y0-$H/2.    ;
      $x $y
      $x=$x0+$W/2.-$d    ;$y=$y0-$H/2.    ;
      $x $y
      $x=$x0+$W/2.       ;$y=$y0-$H/2.    ;
      $x $y
      $x=$x0+$W/2.       ;$y=$y0-$H/2.+$d ;
      $x $y 
    n-stretch   
      1. 2.0 0.  
    periodicity
      2
    sharpness   
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
      $sharp         
    t-stretch   
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
      0 10      
    n-dist      
      fixed normal distance
        -0.025
    lines 
      $nl 5
    boundary conditions
      -1 -1  2 0
    mappingName
      resonator
  exit 
exit
generate an overlapping grid
    backGround
    resonator
  done
  change parameters
    ghost points
      all
      $ng $ng $ng $ng $ng $ng
    order of accuracy
      $orderOfAccuracy
    interpolation type
      $interpType
    exit
  compute overlap
  exit
save a grid
$name
exit
exit
