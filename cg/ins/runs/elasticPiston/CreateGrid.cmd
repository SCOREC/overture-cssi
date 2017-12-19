$factor=1; 
GetOptions("factor=f"=> \$factor);
$name="";  $ml=0; 
$ds = .1/$factor;
$xa = 0.; $xb = 1.; $ya = 0.; $yb = 1.;
$ml2 = 2**$ml; 
$order=2; $orderOfAccuracy = "second order"; $ng=2; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
#
create mappings
#
  rectangle
    $nx=intmg( ($xb-$xa)/$ds );
    $ny=intmg( ($yb-$ya)/$ds );
    set corners
      $xa $xb $ya $yb 
    lines
      $nx $ny 
    boundary conditions
      -1 -1 3 4
    share
      0 0 0 0
    mappingName
      square
    exit
  exit
generate an overlapping grid
  square
  done
  change parameters
    ghost points
      all
      2 2 2 2
  exit
  compute overlap
exit
save an overlapping grid
$name = "square" . "$factor" . ".hdf";
$name
exit
exit