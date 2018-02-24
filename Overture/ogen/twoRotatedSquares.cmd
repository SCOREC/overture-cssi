#
# Two adjacent squares ROTATED
#
#
# usage: ogen [noplot] twoRotatedSquares -factor=<num> -order=[2/4/6/8] -interp=[e/i] -per=[0|1]
# 
# examples:
#     ogen -noplot twoRotatedSquares -order=2 -interp=i -factor=1 
#     ogen -noplot twoRotatedSquares -order=2 -interp=i -factor=2 
#     ogen -noplot twoRotatedSquares -order=2 -interp=i -factor=4 
# 					       				       
#     ogen -noplot twoRotatedSquares -order=2 -interp=e -factor=2 
#     ogen -noplot twoRotatedSquares -order=4 -interp=e -factor=2 
# -- periodic 
#     ogen -noplot twoRotatedSquares -order=2 -interp=e -per=1 -factor=2
#
$xa=-1.; $xb=1.; $ya=-1.; $yb=1.; $per=0; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=i"=> \$factor,"interp=s"=> \$interp,"per=i"=> \$per);
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=5; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $per eq 1 ){ $suffix .= "p"; }
$name = "twoRotatedSquares" . "$interp$factor" . $suffix . ".hdf";
# 
$ds=.1/$factor;
$width = ($order-2)/2;
if( $interp eq "e" ){ $width=$width+1.; }
$overlap = $ds*$width + $ds*.125;
# 
create mappings
  rectangle
    set corners
     $xas=$xa; $xbs=$overlap; $yas=$ya; $ybs=$yb
     $xas $xbs $ya $yb 
    lines
      $nx=int( ($xbs-$xas)/$ds+1.5 );
      $ny=int( ($ybs-$yas)/$ds+1.5 );
      $nx $ny
    boundary conditions
      if( $per eq 0 ){ $cmd="1 0 2 2"; }else{ $cmd="1 0 -1 -1"; }
      $cmd
    share
      0 0 1 2
    mappingName
      leftSquare0
    exit
# -- rotate 
  rotate/scale/shift
   transform which mapping?
      leftSquare0
    rotate
      45
      0 0 
    mappingName
     leftSquare
 exit
# 
  rectangle
    set corners
     $xas=-$overlap; $xbs=$xb; $yas=$ya; $ybs=$yb
     $xas $xbs $ya $yb 
    lines
      $ds2 = $ds*1.12345;
      $nx=int( ($xbs-$xas)/$ds2+1.5 );
      $ny=int( ($ybs-$yas)/$ds2+1.5 );
      $nx $ny
    boundary conditions
      if( $per eq 0 ){ $cmd="0 1 2 2"; }else{ $cmd="0 1 -1 -1"; }
      $cmd
    share
      0 0 1 2
    mappingName
      rightSquare0
    exit
# -- rotate 
  rotate/scale/shift
   transform which mapping?
      rightSquare0
    rotate
      45
      0 0 
    mappingName
     rightSquare
 exit
#
  exit
#
generate an overlapping grid
  leftSquare
  rightSquare
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
save a grid (compressed)
$name
sis
exit
