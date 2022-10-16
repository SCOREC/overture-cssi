*
* Two adjacent squares
*
*
* usage: ogen [noplot] twoSquares -factor=<num> -order=[2/4/6/8] -interp=[e/i] -per=[0|1]
* 
* examples:
*     ogen noplot twoSquares -order=2 -interp=i -factor=1 
*     ogen noplot twoSquares -order=2 -interp=i -factor=2 
*     ogen noplot twoSquares -order=2 -interp=i -factor=4 
* 					       				       
*     ogen noplot twoSquares -order=2 -interp=e -factor=2 
*     ogen noplot twoSquares -order=4 -interp=e -factor=2 
*
$xa=-1.; $xb=1.; $ya=-1.; $yb=1.; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$per=0; 
* 
* get command line arguments
GetOptions( "order=i"=>\$order,"factor=i"=> \$factor,"interp=s"=> \$interp,"per=i"=>\$per);
* 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=5; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
* 
$suffix = ".order$order"; 
if( $per eq 1 ){ $suffix .= "p"; }
$name = "twoSquares" . "$interp$factor" . $suffix . ".hdf";
* 
$ds=.1/$factor;
$width = ($order-2)/2;
if( $interp eq "e" ){ $width=$width+1.; }
$overlap = $ds*$width + $ds*.125;
* 
create mappings
  rectangle
    set corners
     $delta=0.; # -1.e-3; # for testing ogen perturb the grid 
     $xas=$xa; $xbs=$overlap; $yas=$ya+$delta; $ybs=$yb+$delta
     $xas $xbs $yas $ybs 
    lines
      $nx=int( ($xbs-$xas)/$ds+1.5 );
      $ny=int( ($ybs-$yas)/$ds+1.5 );
      $nx $ny
    boundary conditions
      if( $per eq 0 ){ $cmd="1 0 2 2"; }else{ $cmd="1 0 -1 -1"; }
      $cmd
      # 1 0 2 2 
    share
      0 0 1 2
    mappingName
      leftSquare
    exit
* 
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
      # 0 1 2 2 
    share
      0 0 1 2
    mappingName
      rightSquare
    exit
*
  exit
*
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
