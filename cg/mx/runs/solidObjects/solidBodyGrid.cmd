#
#   Grid for a solid body in a channel (for CGmx)
#
#
# usage: ogen [-noplot] solidBodyGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -numGhost=<i> ...
#                                     -body=[C|L|R|P|I]
# 
#  -body = body shape: "C"=C shape, L=L-shape, etc.
#
# examples:
# I 
# ogen -noplot solidBodyGrid -prefix=IBodyGrid -body=I -interp=e -order=2 -factor=4
# 
# P 
# ogen -noplot solidBodyGrid -prefix=PBodyGrid -body=P -interp=e -order=2 -factor=4
# 
# R 
# ogen -noplot solidBodyGrid -prefix=RBodyGrid -body=R -interp=e -order=2 -factor=4
# 
# SPLIT-RING 
# ogen -noplot solidBodyGrid -prefix=CBodyGrid -body=C -interp=e -order=2 -factor=4
#
$prefix="solidBodyGrid";
$body="C"; # body shape  
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-1.5; $xb=1.5; $ya=-1.5; $yb=1.5; $ml=0; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"numGhost=i"=> \$numGhost,"body=s"=> \$body,"prefix=s"=> \$prefix );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.05/$factor;
#
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
#
create mappings
#
  # ----- Background ----
  rectangle
    mappingName
     backGround
   set corners
    $xa $xb $ya $yb
   lines
    $nx = int( ($xb-$xa)/$ds +1.5 ); 
    $ny = int( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
    boundary conditions
      # 0 0 0 0
      1 2 -1 -1 
  exit
#
$gridNames="#"; 
#
# -----------------------------------------
# -------- build I-beam grid --------------
# -----------------------------------------
#
$iBeamNames="#"; 
$count=0;
#- $centerHeight=.75; $centerWidth=.25;
#- $edgeHeight=.25;   $edgeWidth=1.;
#- # -- we only build the hyperbolic grids once -- they are rotated/translated below 
#- include iBeamSolid.h
#- # 
#- # ---------------------- IBEAM 1 -------------------------------
#- $xr=.0; $yr=0; # center for rotation
#- $xShift=2; $yShift=0; $angle=0; $xMin=-$radX; $xMax=-$xMin; $yMin=-$radY; $yMax=-$yMin;
#- include buildIBeamGrids.h
# 
# open graphics
#   
# ---------------------------------------
# ------- build base split-ring ---------
# ---------------------------------------
#
$H=1.; $W=1.; $w=.4; $d=.4;
$numPointsPerSegment=2;   # number of extra points per segment (default=3) increase to make sharper corners
$arcLengthScaleFactor=.9; # scale arcLenght by this value for computing number of grid point default =1
$numRadial=5;             # number of point in radial direction is $numRadial + $order/2  (default=6)
$splitRingNames="#"; 
$count=0;
# include splitRingSolid.h
if( $body eq "C" ){ $cmd="include $ENV{'CG'}/mx/runs/solidObjects/splitRingSolid.h"; }else{ $cmd="#"; }
$cmd 
#
# ------------------ Split Ring ---------------
$xr=.0; $yr=0; # center for rotation
$xShift=0;  $yShift=0; $angle=0;
$xMin=-$radX; $xMax=-$xMin; $yMin=-$radY; $yMax=-$yMin;
if( $body eq "C" ){ $cmd="include $ENV{'CG'}/mx/runs/solidObjects/buildSplitRingGrids.h"; }else{ $cmd="#"; }
$cmd 
# include buildSplitRingGrids.h
#  
# ------ "R" shape -------------
#
$fluidGridNames="#";
$solidGridNames="#"; 
$solidDomains="#";
$shareInterface=99;
$bcInterface=99;
$count=0; 
$degree=3; # degree of nurbs
#
#   ---- "R" shape ----
#   The width of the R is about 1.25, the height is about 2.5 
$edgeWidth=0.5; $middleWidth=0.25; 
$slantLength=0.5; $bottomHeight=0.5; 
$topHeight=0.25;
$cx=-.25; $cy=$edgeWidth;
$coreRadius=.25; $xCore=.5*$middleWidth; $yCore=0.; # offset core from center
# --- define curve ---
if( $body eq "R" ){ $cmd="include $ENV{'CG'}/mx/runs/solidObjects/solidRBeamCurve.h"; }else{ $cmd="#"; }
$cmd
#  --- build grids ---
if( $body eq "R" ){ $cmd="include $ENV{'CG'}/mx/runs/solidObjects/solidBodyGrids.h"; }else{ $cmd="#"; }
$cmd 
#
#   ---- "P" shape ----
#
$edgeWidth=0.5; $middleWidth=0.25; 
$slantLength=0.5; $bottomHeight=0.5; 
$topHeight=0.25;
$cx=-.25; $cy=$edgeWidth;
$coreRadius=.25; $xCore=.5*$middleWidth; $yCore=0.; # offset core from center
# --- define curve ---
if( $body eq "P" ){ $cmd="include $ENV{'CG'}/mx/runs/solidObjects/solidPBeamCurve.h"; }else{ $cmd="#"; }
$cmd
# include solidPBeamCurve.h
#  --- build grids ---
if( $body eq "P" ){ $cmd="include $ENV{'CG'}/mx/runs/solidObjects/solidBodyGrids.h"; }else{ $cmd="#"; }
$cmd
#
#   ---- "I" shape ----
#
# printf("centerHeight=$centerHeight, centerWidth=$centerWidth, edgeWidth=$edgeWidth, edgeHeight=$edgeHeight\n");
#
$centerHeight=1.25; $centerWidth=0.5; $edgeWidth=1.75; $edgeHeight=0.5; 
#
$angle=0; 
$cx=0; $cy=0;
$coreRadius=.15; $xCore=0; $yCore=.5*$centerHeight+.25*$edgeHeight; # offset core from center
# --- define curve ---
if( $body eq "I" ){ $cmd="include $ENV{'CG'}/mx/runs/solidObjects/solidIBeamCurve.h"; }else{ $cmd="#"; }
$cmd
#  --- build grids ---
if( $body eq "I" ){ $cmd="include $ENV{'CG'}/mx/runs/solidObjects/solidBodyGrids.h"; }else{ $cmd="#"; }
$cmd
# 
#
exit 
#
generate an overlapping grid
  backGround
#  innerBackGround
  $gridNames
  $iBeamNames
  $splitRingNames
  $fluidGridNames
  $solidGridNames  
  done
  change parameters
    # ---- Specify domains -----
    specify a domain
      # domain name:
      outerDomain 
      # grids in the background domain:
      backGround
      $gridNames
      $fluidGridNames
     done
#
    specify a domain
      # domain name:
      innerDomain
      # grids in the domain:
      $splitRingNames
      $solidGridNames      
    done
#
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    exit
   # pause
   # open graphics
  compute overlap
exit
#
save an overlapping grid
$name
solidBodyGrid
exit
