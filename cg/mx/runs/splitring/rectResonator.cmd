#
# 
#
#
# usage: ogen [-noplot] rectResonator -factor=<num> -order=[2/4/6/8] -interp=[e/i] -ml=<>
# 
#  -ml = number of (extra) multigrid levels to support
# 
$order=2; $factor=4; $interp="e";  $ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "explicit for all grids";
$prefix="rectResonator"; $name=""; 
$nr=3;
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"interp=s"=> \$interp,"name=s"=> \$name,"ml=i"=>\$ml,\
            "prefix=s"=> \$prefix );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $name eq "" ){$name = $prefix . "$factor" . $suffix . ".hdf";}
#
sub corner {\
  my ($x1,$x2,$x3,$y1,$y2,$y3,$n1,$n2,$gName) = @_; \
  my $cmd = "";\
  $cmd .=  "SmoothedPolygon\n";\
  $cmd .=  "  vertices\n";\
  $cmd .=  "    3\n";\
  $cmd .=  "    $x1  $y1 \n";\
  $cmd .=  "    $x2  $y2 \n";\
  $cmd .=  "    $x3  $y3 \n";\
  $cmd .=  "  sharpness\n";\
  $cmd .=  "    10 \n";\
  $cmd .=  "    10 \n";\
  $cmd .=  "    10 \n";\
  $cmd .=  "  n-stretch\n";\
  $cmd .=  "    1. 2. 0.\n";\
  $cmd .=  "  t-stretch\n";\
  $cmd .=  "    0 10\n";\
  $cmd .=  "    0 10\n";\
  $cmd .=  "    0 10\n";\
  $cmd .=  "  \$nDist=-\$nr*\$ds;\n";\
  $cmd .=  "  n-dist\n";\
  $cmd .=  "    fixed normal distance\n";\
  $cmd .=  "    \$nDist\n";\
  $cmd .=  "  \$nl = 4*int(1.5*$factor);\n";\
  $cmd .=  "  lines \n";\
  $cmd .=  "    \$nl 7\n";\
  $cmd .=  "  boundary conditions\n";\
  $cmd .=  "    0 0 2 0\n";\
  $cmd .=  "  share\n";\
  $cmd .=  "    2 3 4 0\n";\
  $cmd .=  "  mappingName\n";\
  $cmd .=  " $gName\n";\
  $cmd .=  "exit\n";\
  return $cmd;\
} 
sub edge {\
  my ($x1,$x2,$y1,$y2,$n2,$gName) = @_; \
  my $dist;\
  if( abs($y2-$y1) > 0. ) {  $dist=abs($y2-$y1) }\
  else{  $dist = abs($x2-$x1) };\
  my $cmd = "";\
  $cmd .= "";\
  $cmd .=  "SmoothedPolygon\n";\
  $cmd .=  "  vertices\n";\
  $cmd .=  "    2\n";\
  $cmd .=  "    $x1  $y1 \n";\
  $cmd .=  "    $x2  $y2 \n";\
  $cmd .=  "  n-stretch\n";\
  $cmd .=  "    1. 2. 0.\n";\
  $cmd .=  "  \$nDist=-\$nr*\$ds;\n";\
  $cmd .=  "  n-dist\n";\
  $cmd .=  "    fixed normal distance\n";\
  $cmd .=  "    \$nDist\n";\
  $cmd .=  "  \$nl = \$dist/\$ds;\n";\
  $cmd .=  "  lines \n";\
  $cmd .=  "    \$nl 7\n";\
  $cmd .=  "  boundary conditions\n";\
  $cmd .=  "    0 0 2 0\n";\
  $cmd .=  "  share\n";\
  $cmd .=  "    2 3 4 0\n";\
  $cmd .=  "  mappingName\n";\
  $cmd .=  " $gName\n";\
  $cmd .=  "exit\n";\
  return $cmd;\
} 
#
# 
$ds=.01/$factor;
$xa = -1.;
$xb =  1.;
$ya = -1.;
$yb =  1.; 
$pi = 4.*atan2(1.,1.);
$nx = ($xb-$xa)/$ds;
$ny = ($yb-$ya)/$ds; 
$x0=-.3;
$y0 = 0.0;
$H = 0.8;
$H2 = .3;
$L = .1;
$W = .8;
$d1 = $L/3.; 
$d2 = $L/2.+$L/4. ;
$yt1 = $H/2.;
@v1 = ($x0-$L,$y0+$yt1);
@v2 = ($x0-$L,$y0-$yt1) ;
@v3 = ($x0-$L+$W,$y0-$yt1); 
@v4 = ($x0-$L+$W,$y0-$yt1+$H2) ;
@v5 = ($x0-2.*$L+$W,$y0-$yt1+$H2) ;
@v6 = ($x0-2.*$L+$W,$y0-$yt1+$L) ;
@v7 = ($x0,$y0-$yt1+$L) ;
@v8 = ($x0,$y0+$yt1-$L) ;
@v9 = ($x0+$W-2.*$L,$y0+$yt1-$L) ;
@v10 = ($x0+$W-2.*$L,$y0+$yt1-$H2) ;
@v11 = ($x0+$W-$L,$y0+$yt1-$H2) ;
@v12 = ($x0+$W-$L,$y0+$yt1) ; 
# 
create mappings
#
rectangle 
  set corners
    $xa $xb $ya $yb
  lines 
    $nx $ny
  boundary conditions
    1 1 1 1
  mappingName
    backGround 
exit
$cmd="";
$cmd .= corner($v1[0]+$d2,$v1[0],$v1[0],$v1[1],$v1[1],$v1[1]-$d2,10,5,'C1');
$cmd .= edge($v1[0],$v2[0],$v1[1]-$d1,$v2[1]+$d1,5,'C2'); 
#
$cmd .= corner($v2[0],$v2[0],$v2[0]+$d2,$v2[1]+$d2,$v2[1],$v2[1],30,5,'C3') ;
$cmd .= edge($v2[0]+$d1,$v3[0]-$d1,$v2[1],$v3[1],5,'C4') ;
#
$cmd .= corner($v3[0]-$d2,$v3[0],$v3[0],$v3[1],$v3[1],$v3[1]+$d2,30,5,'C5') ;
$cmd .= edge($v3[0],$v4[0],$v3[1]+$d1,$v4[1]-$d1,5,'C6') ;
# cap
$cmd .= corner($v4[0],$v4[0],$v4[0]-$d2,$v4[1]-$d2,$v4[1],$v4[1],30,5,'C7');
$cmd .= corner($v5[0]+$d2,$v5[0],$v5[0],$v5[1],$v5[1],$v5[1]-$d2,30,5,'C8');
#
$cmd .= edge($v5[0],$v6[0],$v5[1]-$d1,$v6[1]+$d1,5,'C9') ;
$cmd .= corner($v6[0],$v6[0],$v6[0]-$d2,$v6[1]+$d2,$v6[1],$v6[1],30,5,'C10');
#
$cmd .= edge($v6[0]-$d1,$v7[0]+$d1,$v6[1],$v7[1],5,'C11') ;
$cmd .= corner($v7[0]+$d2,$v7[0],$v7[0],$v7[1],$v7[1],$v7[1]+$d2,30,5,'C12');
#
$cmd .= edge($v7[0],$v8[0],$v7[1]+$d1,$v8[1]-$d1,5,'C13') ;
$cmd .= corner($v8[0],$v8[0],$v8[0]+$d2,$v8[1]-$d2,$v8[1],$v8[1],30,5,'C14');
#
$cmd .= edge($v8[0]+$d1,$v9[0]-$d1,$v8[1],$v9[1],5,'C15') ;
$cmd .= corner($v9[0]-$d2,$v9[0],$v9[0],$v9[1],$v9[1],$v9[1]-$d2,30,5,'C16');
$cmd .= edge($v9[0],$v10[0],$v9[1]-$d1,$v10[1]+$d1,5,'C17');
# cap 
$cmd .= corner($v10[0],$v10[0],$v10[0]+$d2,$v10[1]+$d2,$v10[1],$v10[1],30,5,'C18') ;
$cmd .= corner($v11[0]-$d2,$v11[0],$v11[0],$v11[1],$v11[1],$v11[1]+$d2,30,5,'C19') ;
#
$cmd .= edge($v11[0],$v12[0],$v11[1]+$d1,$v12[1]-$d1,5,'C20');
$cmd .= corner($v12[0],$v12[0],$v12[0]-$d2,$v12[1]-$d2,$v12[1],$v12[1],30,5,'C21') ;
# closeall
$cmd .= edge($v12[0]-$d1,$v1[0]+$d1,$v12[1],$v1[1],5,'C22') ;
$cmd 
exit
generate an overlapping grid
    backGround
    C1
    C2
    C3 
    C4 
    C5
    C6
    C7
    C8
    C9
    C10
    C11
    C12
    C13
    C14
    C15
    C16
    C17
    C18
    C19
    C20
    C21
    C22
  done
  change parameters
 # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      # add an extra ghost for singular problems: 
      $ng $ng $ng $ng $ng $ng 
  exit
#  display intermediate results
# change parameters
#   minimize overlap
#   exit 
  compute overlap 
#*  display computed geometry 
  exit
# 
# save an overlapping grid
save a grid (compressed)
$name
backGround
exit

