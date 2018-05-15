$order=2; $factor=1; $interp = "e"; 
$orderOfAccuracy="second order";
$ng=2; 
# -- convert a number so that it is a power of 2 plus 1 --
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
* get command line arguments
GetOptions("name=s"=> \$name,"zb=f"=>\$zb,"order=i"=>\$order,"factor=f"=> \$factor,"interp=s"=> \$interp);
*
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
*
$suffix = ".order$order";
if( $name eq "" ){$name = "splitRing3d" . "$interp$factor" . $suffix . ".hdf";}
*
* domain parameters:
$ds0=.05;
$za=-1.; $zb=1.; $zc=0.;
$sharp=30; 
# surrounding box
$xa = -1.; $xb = 1.;
$ya = -1.; $yb = 1.;
$za = -1.; $zb = 1.;
# resonator size
$H = 1.; $L=0.3; $W=0.9;
$w = $L/40;
$dw = $w;
$D = 0.125;
# reference point
$x0 = -0.5; $y0 = 0.; $z0 = 0.; 
$ds = $ds0/$factor; # target grid spacing
# nr = number of lines in normal directions to boundaries
$nr = 5 + $ng + ($order-2);
#
$pi = 4.*atan2(1.,1.);
#
#
$wallBC=7;   # BC for walls of the box 
$wallShare=7; 
# 
$stretchFactor=.725; # scale nDist to account for stretching
$nDist= $stretchFactor*($nr-2)*$ds;
if( $order eq 4 ){ $nDist=$nDist+2.*$ds; } 
sub sweepMap {\
  my ($numPts,$curve,$name,$section) = @_; \
  my $i;\
  $cmd = "sweep \n";\
  $cmd .= "  specify scaling factors\n";\
  $cmd .= "   $numPts \n";\
  for( $i = 0; $i < $numPts; $i++){ $cmd .= " 1. \n"; }\
  $cmd .= " choose reference mapping \n";\
  $cmd .= " $section \n";\
  $cmd .= " choose sweep curve \n";\
  $cmd .= "   $curve \n";\
  $cmd .= " lines \n";\
  $cmd .= " \$nz=int( (\$H+2*\$W+2*\$L)/\$ds+1.5 ); \n";\
  $cmd .= " \$nx=int( 4*\$D/\$ds+1.5 ); \n";\
  $cmd .= " \$nx \$nr \$nz \n";\
  $cmd .= " boundary conditions \n";\
  $cmd .= " -1 -1  \$wallBC 0 0  0 \n";\
  $cmd .= "  share \n";\
  $cmd .= " 0 0 \$wallShare 0 0 0 \n";\
  $cmd .= "  mappingName \n";\
  $cmd .= "   $name \n";\
  $cmd .= "  exit \n";\
  return $cmd;\
} 
sub splineCurve {\
  my $i;\
  my ($name,$numPts,$xref,$yref,$zref) = @_; \
  my @xdata = @{ $xref }; \
  my @ydata = @{ $yref }; \
  my @zdata = @{ $zref }; \
  my $cmd = "spline (3D) \n"; \
  $cmd .= "  enter spline points \n";\
  $cmd .= "    $numPts \n";\
  for($i = 0; $i < $numPts; $i++){ $cmd .= " $xdata[$i] $ydata[$i] $zdata[$i] \n";}\
  $cmd .= "mappingName \n"; \
  $cmd .= "  $name \n";\
  $cmd .= "exit\n";\
  return $cmd; \
} 
sub crossSection {\
  my $cmd = "smoothedPolygon \n"; \
  $cmd .= "boundary conditions\n";\
  $cmd .= "  -1 -1 \$wallBC 0 \n";\
  $cmd .= "vertices\n";\
  $cmd .= "  6\n";\
  $cmd .= "@{[    0.]}  @{[-$D/2.]} \n";\
  $cmd .= "@{[-$D/2.]}  @{[-$D/2.]} \n";\
  $cmd .= "@{[-$D/2.]}  @{[ $D/2.]} \n";\
  $cmd .= "@{[ $D/2.]}  @{[ $D/2.]} \n";\
  $cmd .= "@{[ $D/2.]}  @{[-$D/2.]} \n";\
  $cmd .= "@{[    0.]}  @{[-$D/2.]} \n";\
  $cmd .= "  sharpness \n"; \
  $cmd .= "    \$sharp \n"; \
  $cmd .= "    \$sharp \n"; \
  $cmd .= "    \$sharp \n"; \
  $cmd .= "    \$sharp \n"; \
  $cmd .= "    \$sharp \n"; \
  $cmd .= "    \$sharp \n"; \
  $cmd .= " n-dist \n"; \
  $cmd .= "   fixed normal distance \n";\
  $cmd .= "    \$nDist \n"; \
  $cmd .= " lines \n"; \
  $cmd .= " 51, 7 \n"; \
  $cmd .= " n-stretch \n"; \
  $cmd .= " 1., 1., 0. \n"; \
  $cmd .= " t-stretch \n"; \
  $cmd .= "  .0   \$tStretchbLB \n"; \
  $cmd .= "  \$tsa = 2.*\$tStretchaLB; \n";\
  $cmd .= "  \$tsa \$tStretchbLB  \n"; \
  $cmd .= "  \$tsa \$tStretchbLB  \n"; \
  $cmd .= "  \$tsa \$tStretchbLB  \n"; \
  $cmd .= "  \$tsa \$tStretchbLB  \n"; \
  $cmd .= "  .0 \$tStretchbLB  \n"; \
  $cmd .= "mappingName \n"; \
  $cmd .= "  crossSection \n";\
  $cmd .= "exit\n";\
  return $cmd;\
}
sub capBottom {\
  my $map = $_[0]; \
  my $cmd = "transform Mappings...\n";\
  $cmd .= "rotate/scale/shift \n";\
  $cmd .= "transform which mapping? \n";\
  $cmd .= "    $map \n";\
  $cmd .= "scale \n";\
  $cmd .= "  \$D, \$D, 1.0 \n";\
  $cmd .= "shift \n";\
  $cmd .= "    @{[$x0+$W]}, 0., @{[-$ds]} \n";\
  $cmd .= "rotate \n";\
  $cmd .= "   -90,0 \n";\
  $cmd .= "    @{[$x0+$W]},0,0 \n";\
  $cmd .= "shift \n";\
  $cmd .= "  \$yshift = \$y0-\$H/2.+(\$W-\$L)/2.-0*\$ng*\$ds;\n";\
  $cmd .= "  0, \$yshift, \$z0 \n";\
  $cmd .= "mappingName \n";\
  $cmd .= "  $map"."Bottom \n";\
  $cmd .= "exit \n";\
  return $cmd;\
}
sub capTop {\
  my $map = $_[0]; \
  my $cmd = "transform Mappings...\n";\
  $cmd .= "rotate/scale/shift \n";\
  $cmd .= "transform which mapping? \n";\
  $cmd .= "    $map \n";\
  $cmd .= "scale \n";\
  $cmd .= "  \$D, \$D, 1. \n";\
  $cmd .= "shift \n";\
  $cmd .= "    @{[$x0+$W]}, 0., @{[-$ds]} \n";\
  $cmd .= "rotate \n";\
  $cmd .= "   90,0 \n";\
  $cmd .= "    @{[$x0+$W]},0,0 \n";\
  $cmd .= "shift \n";\
  $cmd .= "  \$yshift = \$y0+\$H/2.-(\$W-\$L)/2.+0*\$ng*\$ds;\n";\
  $cmd .= "  0, \$yshift, \$z0 \n";\
  $cmd .= "mappingName \n";\
  $cmd .= "  $map"."Top \n";\
  $cmd .= "exit \n";\
  return $cmd;\
}
sub bodyGrid {\
  my ($map,$part) = @_;\
  my $cmd = "hyperbolic \n";\
  $cmd .="  Start curve:$map \n";\
  $cmd .="  forward  \n";\
  if( $part eq "face" ){ $cmd .="  distance to march \$nDist \n";}\
  else{ $cmd .= " distance to march   -\$nDist  \n"; }\
  if( $part eq "face" ){ $cmd .=" \$linesToMarch = \$nr-1; \n"; }\
  elsif( $part eq "end" ){ $cmd .=" \$linesToMarch = \$nr-1; \n"; }\
  $cmd .="  lines to march \$linesToMarch \n";\
  $cmd .="  generate \n";\
  $cmd .="  fourth order \n";\
  if( $part eq "end" ){\
    $cmd .= "boundary offset 0 0 0 0 0 0 (l r b t b f) \n";\
    $cmd .= "boundary conditions \n";\
    $cmd .= "  0 0 -1 -1   \$wallBC 0 \n";\
    $cmd .= "share \n";\
    $cmd .= "  0  0 0  0 \$wallShare 0 \n";}\
  elsif( $part eq "face" ){\
    $cmd .= "boundary offset \$ng \$ng \$ng \$ng 0 0 (l r b t b f) \n";\
    $cmd .= "boundary conditions \n";\
    $cmd .= "   0 0   0 0 \$wallBC 0    \n";\
    $cmd .= "share \n";\
    $cmd .= "  0  0 0  0 \$wallShare 0   \n";}\
   if( $part eq "face" ){ \
     $cmd .= " lines \n";\
     $cmd .= "\$nx = int(.25*(\$xblb-\$xalb)/\$ds +1.5 ); \n";\
     $cmd .= "\$ny = int(.25*(\$yblb-\$yalb)/\$ds +1.5 ); \n";\
     $cmd .= "\$nx \$ny \$nr \n";}\
  elsif( $part eq "end" ){ \
     $cmd .= "lines \n";\
     $cmd .= "\$nx = int(.25*(1.-\$capWidthFraction)*( (\$xblb-\$xalb) + (\$yblb-\$yalb) )/\$ds +1.5   ); \n";\
     $cmd .= "\$ny = int( 4*\$D/\$ds+1.5 ); \n";\
     $cmd .= "\$nx \$ny \$nr \n";}\
  $cmd .= "mappingName  \n";\
  $cmd .= "    $map"."Body\n";\
  $cmd .= "exit\n";\
  return $cmd;\
}
sub resonatorGrid {\
  my $cmd = "";\
  @xpts = ( $x0 + $W, $x0 + $W, $x0 + $W, $x0 + $W,\
            $x0 + ($W+$L)/2. + 3*$dw, $x0 + ($W+$L)/2. + 2*$dw, $x0 + ($W+$L)/2. + 1*$dw, $x0 + ($W+$L)/2. + 0*$dw,\
            $x0 + ($W-$L)/2. - 0*$dw, $x0 + ($W-$L)/2. - 1*$dw, $x0 + ($W-$L)/2. - 2*$dw, $x0 + ($W-$L)/2. - 3*$dw,\
            $x0,  $x0, $x0 ,$x0 ,\
            $x0, $x0,   $x0 ,  $x0,\
            $x0 + ($W-$L)/2. - 3*$dw, $x0 + ($W-$L)/2. - 2*$dw, $x0 + ($W-$L)/2. - 1*$dw, $x0 + ($W-$L)/2. - 0*$dw,\
            $x0 + ($W+$L)/2. + 0*$dw, $x0 + ($W+$L)/2. + 1*$dw, $x0 + ($W+$L)/2. + 2*$dw, $x0 + ($W+$L)/2. + 3*$dw,\
            $x0 + $W, $x0 + $W,   $x0 + $W, $x0 + $W);\
  @ypts = ( $y0 + $H/2. - ($W-$L)/2. + 0*$dw, $y0 + $H/2. - ($W-$L)/2. + 1*$dw, $y0 + $H/2. - ($W-$L)/2. + 2*$dw, $y0 + $H/2. - ($W-$L)/2. + 3*$dw,\
            $y0 + $H/2.,  $y0 + $H/2., $y0 + $H/2., $y0 + $H/2., \
            $y0 + $H/2., $y0 + $H/2., $y0 + $H/2., $y0 + $H/2.,\
            $y0 + $H/2. - ($W-$L)/2. + 3*$dw, $y0 + $H/2. - ($W-$L)/2. + 2*$dw, $y0 + $H/2. - ($W-$L)/2. + 1*$dw, $y0 + $H/2. - ($W-$L)/2. + 0*$dw,\
            $y0 - $H/2. + ($W-$L)/2. - 0*$dw, $y0 - $H/2. + ($W-$L)/2. - 1*$dw, $y0 - $H/2. + ($W-$L)/2. - 2*$dw, $y0 - $H/2. + ($W-$L)/2. - 3*$dw,\
            $y0 - $H/2., $y0 - $H/2., $y0 - $H/2., $y0 - $H/2.,\
            $y0 - $H/2., $y0 - $H/2., $y0 - $H/2., $y0 - $H/2.,\
            $y0 - $H/2. + ($W-$L)/2. - 3*$dw, $y0 - $H/2. + ($W-$L)/2. - 2*$dw, $y0 - $H/2. + ($W-$L)/2. - 1*$dw, $y0 - $H/2. + ($W-$L)/2. - 0*$dw);\
  @zpts = ( $z0, $z0, $z0, $z0, \
            $z0, $z0, $z0, $z0, \
            $z0, $z0, $z0, $z0, \
            $z0, $z0, $z0, $z0, \
            $z0, $z0, $z0, $z0, \
            $z0, $z0, $z0, $z0, \
            $z0, $z0, $z0, $z0, \
            $z0, $z0, $z0, $z0);\
  $cmd .= splineCurve("sweepCurve",32,\@xpts,\@ypts,\@zpts);\
  $cmd .= sweepMap(32,"sweepCurve","resonatorBody","crossSection");\
  $cmd .= capTop("capNoEnds");\
  $cmd .= capTop("capFace");\
  $cmd .= capBottom("capNoEnds");\
  $cmd .= capBottom("capFace");\
  $cmd .= bodyGrid("capNoEndsBottom","end");\
  $cmd .= "  stretch coordinates  \n";\
  $cmd .= "   Stretch r1:exp to linear \n";\
  $cmd .= "   STP:stretch r1 expl: parameters 0.01 1 0 (a,b,c) \n";\
  $cmd .= "   close r1 stretching parameters \n";\
  $cmd .= "   exit \n";\
  $cmd .= bodyGrid("capFaceBottom","face");\
  $cmd .= bodyGrid("capFaceTop","face");\
  $cmd .= bodyGrid("capNoEndsTop","end");\
  $cmd .= "  stretch coordinates  \n";\
  $cmd .= "   Stretch r1:exp to linear \n";\
  $cmd .= "   STP:stretch r1 expl: parameters 0.01 1 0 (a,b,c) \n";\
  $cmd .= "   close r1 stretching parameters \n";\
  $cmd .= "   exit \n";\
return $cmd;\
} 
# 
create mappings
*
  lofted surface
    smoothed polygon sections
    $tStretchaLB=0.05;  $tStretchbLB=.5*$sharp;  # stretching at corners in tangential directions
    flat tip profile 
    edit profile
      $s = $D/2.;
      $s2 = $D/4.;
      vertices
        6
        0.       .5
        $s2       .5
        $s        .5
        $s       -.5
        $s2      -.5
        0.       -.5
      sharpness 
        $sharp
        $sharp
        $sharp
        $sharp
        $sharp
        $sharp
      t-stretch
      .0 $tStretchbLB
      .0 $tStretchbLB
      $tStretchaLB $tStretchbLB
      $tStretchaLB $tStretchbLB
      .0 $tStretchbLB 
      .0 $tStretchbLB 
    exit 
    edit section 
      vertices
       6 
       $xc = 0.;
       $yc = 0.; 
       $cmd = "\
       @{[$xc - 0/2.]}  @{[$yc - 1/2.]} \n\
       @{[$xc + 1/2.]}  @{[$yc - 1/2.]} \n\
       @{[$xc + 1/2.]}  @{[$yc + 1/2.]} \n\
       @{[$xc - 1/2.]}  @{[$yc + 1/2.]} \n\
       @{[$xc - 1/2.]}  @{[$yc - 1/2.]} \n\
       @{[$xc - 0/2.]}  @{[$yc - 1/2.]} \n";
       $cmd 
      n-dist
       fixed normal distance 
       $nDist 
      sharpness 
        $sharp
        $sharp
        $sharp
        $sharp
        $sharp
        $sharp
     t-stretch
      .0 $tStretchbLB
      $tsa = 2.*$tStretchaLB; # increase t-stretch weighting since there are twice as many corners
      $tsa $tStretchbLB
      $tsa $tStretchbLB
      $tsa $tStretchbLB
      $tsa $tStretchbLB
      .0 $tStretchbLB 
     periodicity
        2 
     n-stretch
       1. 1. 0.  
    exit 
    mappingName
      cap 
    exit 
  box
    mappingName
      backGround
    set corners
     $xa $xb $ya $yb $za $zb
    lines
      $nx=int( ($xb-$xa)/(1.*$ds)+1.5 )+1;
      $ny=int( ($yb-$ya)/(1.*$ds)+1.5 )+1;
      $nz=int( ($zb-$za)/(1.*$ds)+1.5 )+1;
      $nx $ny $nz
    boundary conditions
      1 2 3 4 5 6
    exit 
   builder
     create surface grid...
       Start curve:cap 
         plot options...
         plot boundary lines on reference surface 1
         close plot options
         picking:choose initial curve
         surface grid options...
         initial curve:points on surface 
         $delta= ($ng+1)*$ds; # reduce cap width since ghost lines are added
         # increase cap width for fourth order and coarser grids -- ghost points are on surface anyway:
         if( $order eq 4 && $factor <4 ){ $delta=$delta-2.5*$ds; } 
         $capWidthFraction=.65;   # approx fraction of the end covered by the face grids 
         $xalb=-1/2;
         $xblb= 1/2;
         $yalb=-1/2;
         $yblb= 1/2;
         $zblb=$s;
         $x1=$capWidthFraction*$xblb+$delta; $y1=$capWidthFraction*$yalb-$delta; $z1=$zblb; 
         choose point on surface 0  -$x1 $y1 $z1
         choose point on surface 0    0  $y1 $z1
         choose point on surface 0   $x1 $y1 $z1 
         done 
         $ns = int( 2*$x1/$ds +1.5 );
         $ns = $ns + ($ns % 2);  # Make ns even so we avoid the singularity
         points on initial curve $ns
         *BC: left (forward) fix x, float y and z
         *BC: right (forward) fix x, float y and z 
         $dist=-2*$y1-($ng-2)*$ds;  
         distance to march $dist
         $nMarch = int( abs($dist)/$ds +1.5 );
         lines to march $nMarch 
         generate 
         fourth order
         # Ensure that the ghost points lie on the surface:
         ##$numGhost = $ng; 
         ### boundary offset $numGhost $numGhost $numGhost $numGhost  (l r b t) 
         name capFace 
       exit
     exit
  reparameterize
    transform which mapping?
      cap
    restrict parameter space
      set corners
        .0 .6 0. 1. 
      exit
    $nTheta = int( ( ($xblb-$xalb) + ($yblb-$yalb) )/$ds +1.5 ); 
    $ns =     int( ( ($zblb-$zalb) + ($yblb-$yalb) )/$ds +1.5 ); 
    lines 
       $ns $nTheta 
    mappingName
      capNoEnds 
    exit 
* ----------------------  cross section ----------------------------------- 
  $cmd = crossSection(); 
  $cmd
  $cmd = resonatorGrid(); 
  $cmd 
exit
generate an overlapping grid
  backGround 
  capFaceBottomBody
  capFaceTopBody 
  stretched-capNoEndsBottomBody
  stretched-capNoEndsTopBody 
  resonatorBody 
  done 
  change parameters
    # choose implicit or explicit interpolation
    interpolation type
      $interp
    # -- set the discretization width and interpolation width --
    $cmd =" order of accuracy\n $orderOfAccuracy";
    $cmd
    #
    ghost points
      all
      $ng $ng $ng $ng $ng $ng
  exit 
  project ghost points on shared sides
  make adjustments for nearby boundaries
  compute overlap 
  exit 
#
# save an overlapping grid
save a grid (compressed)
$name
SRR
exit

