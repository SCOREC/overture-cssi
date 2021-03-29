#
#  A lattice of LAYERED solid disks
#     
# usage: 
# ogen [-noplot] layeredDiskArrayGrid -prefix=<s> -factor=<num> -order=[2/4/6/8] -interp=[e/i] ...
#       -nCylx=<num> -nCyly=<num> -rad=<num> -dist=<num> -periodic=[p|np|pn] -xa=<f> -xb=<f> -ya=<f> -yb=<f> -deltaRadius=<>
#       -numLayers=<i> -layerWidth=<f,f>
# 
#   nCylx : number of cylinders in the x direction
#   nCyly : number of cylinders in the y direction
#   rad   : radius of the cylinder
#   dist  : distance between centers (x and y)
#   deltaX : horizontal dist between centers (over-rides dist of set)
#   deltaY : vertical dist between centers (over-rides dist of set)
#   xa,xb, ya,yb : bounds on outer box
#   deltaRadius : if set, use this as the radial width of the grids
#   numLayers : number of layers per disk
#   layerWidth : list of width of the layers
# 
# examples:
#  4 disks with 1 layer around each 
#  ogen -noplot layeredDiskArrayGrid -prefix=fourLayeredDisksGrid -layerWidth=.1 -xa=-1.75 -xb=1.75 -ya=-1.5 -yb=1.5 -nCylx=2 -nCyly=2 -deltaX=1.5 -deltaY=1.5 -order=2 -periodic=np -factor=2
#
#
$prefix = "layeredDiskArrayGrid"; 
$order=2; $factor=1; $interp="i";  $periodic=""; # default values
$orderOfAccuracy = "second order"; $ng=2;
$interp = "e";            # explicit interpolation by default 
$numGhost=-1;             # if >0 use this many ghost 
$nCylx=2; $nCyly=2;       # number of cylinders in each direction
$rad=.4;                  # radius of each cylinder
$numLayers=1;             # number of layers
@layerWidth = ();         # width of each layer 
$dist=1.;                 # spacing between cylinder centres
$deltaX=""; $deltaY="";   # x and y distance between disks (use $dist by default)
$xa=-1.5; $xb=1.5; $ya=-1.5; $yb=1.5;  # background grid
$deltaRadius=-1;          # if set, use this as the radial width of the grids 
# 
# get command line arguments
GetOptions( "prefix=s"=>\$prefix,"rad=f"=>\$rad,"dist=f"=>\$dist,"order=i"=>\$order,"numGhost=i"=>\$numGhost,\
            "factor=i"=> \$factor,"interp=s"=> \$interp,"periodic=s"=>\$periodic,"deltaX=f"=>\$deltaX,"deltaY=f"=>\$deltaY,\
            "nCylx=i"=> \$nCylx,"nCyly=i"=> \$nCyly,"xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "deltaRadius=f"=>\$deltaRadius,"numLayers=i"=> \$numLayers,"layerWidth=f{1,}"=>\@layerWidth );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
if( $interp eq "i" ){ $interpType = "implicit for all grids"; }
#
if( $layerWidth[0] eq "" ){ for( $i=0; $i<$numLayers; $i++ ){ $layerWidth[$i]=.2; } }
# 
$suffix="";
if( $deltaRadius > 0 ){ $prefix .= "Fixed"; }
if( $periodic eq "p"  ){ $prefix .= "p"; }
if( $periodic eq "np" ){ $prefix .= "np"; }
if( $periodic eq "pn" ){ $prefix .= "pn"; }
$suffix .= ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; }
if( $name eq "" ){$name = $prefix . $interp . $factor . $suffix . ".hdf";}
## $suffix = ".order$order"; 
## $name = $prefix . "$nCylx" ."x" . "$nCyly" ."y" . "$interp$factor" . $suffix . ".hdf";
# 
if( $deltaX eq "" ){ $deltaX=$dist; }  # x- spacing between cylinder centres
if( $deltaY eq "" ){ $deltaY=$dist; }  # y- spacing between cylinder centres
# 
#
# Define a subroutine to convert the number of grid points
sub getGridPoints\
{ local($n1,$n2)=@_; \
  $nx=int(($n1-1)*$factor+1.5); $ny=int(($n2-1)*$factor+1.5); \
}
#
#**************************************************************************
#
#
# domain parameters:  
$ds = .025/$factor; # target grid spacing
$pi=4*atan2(1.,1.);
$bcInterface=100;  # bc for interfaces
#
#
create mappings 
#
# Disks are ..
$disk=0;   # number the disks as 1,2,3...
# $layer=0;  # count layers for each disk as 1,2,3
$share=99; 
$domainCommands="#"; # define domains 
#
$radius=$rad; 
$minRadius=$radius;
$maxRadius=$radius; for( $i=0; $i<$numLayers; $i++ ){ $maxRadius=$maxRadius+$layerWidth[$i]; }
$aveRadius=.5*($minRadius+$maxRadius);
#
# ======================================================
# Define a function to build an inner-annulus, inner square, layers and outer-annulus.
# usage:
#   makeDisk(radius,xCenter,yCenter)
# =====================================================
sub makeDisk\
{ local($radius,$xc,$yc)=@_; \
  $count = $count + 1; $layer=1; \
  $innerAnnulus="innerAnnulus$count";     \
  $innerSquare="innerSquare$count";     \
  $outerAnnulus="outerAnnulus$count";     \
  $innerMappingNames .=  "   $innerSquare\n" . "   $innerAnnulus\n"; \
  $outerMappingNames .=  "   $outerAnnulus\n"; \
  $domainCommands  .=  "\n specify a domain\n innerDomain$count\n $innerSquare\n $innerAnnulus\n done"; \
  $nr = 3 + $order + $ng; \
  if( $deltaRadius<0 ){ $deltaRadius=$ds*$nr; } \
  $outerRadius=$radius; $innerRadius=$outerRadius-$deltaRadius; \
  $nTheta=int( (2.*$pi*$aveRadius)/$ds+1.5 ); \
  $nr=int( ($deltaRadius)/$ds+1.5 ); \
  $share= $share+1; \
  $commands = \
  "  annulus \n" . \
  "    mappingName \n" . \
  "      $innerAnnulus \n" . \
  "    centre\n" . \
  "     $xc $yc\n" .   \
  "    boundary conditions \n" . \
  "      -1 -1 0 $bcInterface \n" . \
  "    share\n" . \
  "     * material interfaces are marked by share>=100\n" . \
  "      0 0 0 $share \n" . \
  "    inner and outer radii \n" . \
  "      $innerRadius $outerRadius\n" . \
  "    lines \n" . \
  "      $nTheta $nr \n" . \
  "    exit \n"; \
  $dels = $ds*($order-2); \
  $xas=$xc-$innerRadius-$dels;  $xbs=$xc+$innerRadius+$dels; \
  $yas=$yc-$innerRadius-$dels;  $ybs=$yc+$innerRadius+$dels; \
  $nxs=int( ($xbs-$xas)/$ds+1.5 ); \
  $nys=int( ($ybs-$yas)/$ds+1.5 ); \
  $commands = $commands . \
  "  rectangle \n" . \
  "    mappingName\n" . \
  "      $innerSquare\n" . \
  "    set corners\n" . \
  "     $xas $xbs $yas $ybs \n" . \
  "    lines\n" . \
  "      $nxs $nys\n" . \
  "    boundary conditions\n" . \
  "      0 0 0 0\n" . \
  "    exit \n"; \
  $layer=1; \
  for( $layer=1; $layer<=$numLayers; $layer++ ) \
  { $layerName="disk$count" . "layer$layer";    \
    $innerMappingNames .= "$layerName\n"; \
    $layerDomainName="disk$count" . "Layer$layer" . "Domain";    \
    $domainCommands  .=  "\n specify a domain\n $layerDomainName\n $layerName\n done"; \
    $innerRadius=$outerRadius; \
    $outerRadius=$innerRadius+$layerWidth[$layer-1]; \
    $nr=int( ($layerWidth[$layer-1])/$ds+1.5 ); \
    $sharePrev= $share; $share=$share+1; \
    $commands = $commands . \
    "*\n" . \
    " # Define layer  \n" . \
    "  annulus \n" . \
    "    mappingName \n" . \
    "      $layerName \n" . \
    "    inner and outer radii \n" . \
    "      $innerRadius $outerRadius\n" . \
    "    centre\n" . \
    "      $xc $yc\n" .   \
    "    lines \n" . \
    "      $nTheta $nr\n" . \
    "    boundary conditions \n" . \
    "      -1 -1 $bcInterface $bcInterface \n" . \
    "    share\n" . \
    "     * material interfaces are marked by share>=100\n" . \
    "      0 0 $sharePrev $share   \n" . \
    "    exit \n"; \
    };\
  $innerRadius=$outerRadius; \
  $outerRadius=$innerRadius+$deltaRadius; \
  $nr=int( ($deltaRadius)/$ds+1.5 ); \
  $commands = $commands . \
  "*\n" . \
  "  annulus \n" . \
  "    mappingName \n" . \
  "      $outerAnnulus \n" . \
  "    inner and outer radii \n" . \
  "      $innerRadius $outerRadius\n" . \
  "    centre\n" . \
  "      $xc $yc\n" .   \
  "    lines \n" . \
  "      $nTheta $nr\n" . \
  "    boundary conditions \n" . \
  "      -1 -1 $bcInterface 0 \n" . \
  "    share\n" . \
  "     * material interfaces are marked by share>=100\n" . \
  "      0 0 $share 0   \n" . \
  "    exit \n"; \
}
# ========================================================
#
#
# ===================================================
#   Make an array of cylinders
#     makeDiskArray(radius,nCylx,nCyly,x0,dx0,y0,dy0)
#     
# Make cylinders at centers
#      (x0+i*dx0,y0+j*dy0)  i=0,..,nx0, j=0,..,ny0
# 
# Result: $commands
# =====================================================
sub makeDiskArray \
{ local($radius,$nCylx,$nCyly,$x0,$dx0,$y0,$dy0)=@_; \
  local $cmds; $cmds=""; \
  for( $j=0; $j<$nCyly; $j++ ){ \
  for( $i=0; $i<$nCylx; $i++ ){ \
    makeDisk($radius,$x0+$i*$dx0,$y0+$j*$dy0); \
    $cmds = $cmds . $commands; \
  }}\
  $commands=$cmds; \
}
#
# 
$x0=-($nCylx-1)*$deltaX*.5; $y0=-($nCyly-1)*$deltaY*.5;  
makeDiskArray($rad,$nCylx,$nCyly,$x0,$deltaX,$y0,$deltaY);
$commands
# 
# Background grid 
#
if( $periodic eq "p" ){ $bc ="-1 -1 -1 -1"; }\
elsif( $periodic eq "np" ){ $bc ="1 2 -1 -1"; }\
elsif( $periodic eq "pn" ){ $bc ="-1 -1 3 4"; }else{ $bc="1 2 3 4"; }
# 
  rectangle 
    mappingName
      backGround
    $outerMappingNames = "backGround\n" . $outerMappingNames;
    set corners
     $xa $xb $ya $yb 
    lines
      $nx=int( ($xb-$xa)/$ds+1.5 );
      $ny=int( ($yb-$ya)/$ds+1.5 );
      $nx $ny
    boundary conditions
      $bc
    exit 
#
  exit this menu 
#
generate an overlapping grid 
  backGround
  $mappingNames= $innerMappingNames . $outerMappingNames;
  $mappingNames
  done 
#
  change parameters 
    # -- define the outer domain ----
    specify a domain
      # domain name:
      outerDomain 
      # grids in the domain:
      $outerMappingNames
      done
    # --- define inner domains ---
      $domainCommands 
    #
    #specify a domain
    #  # domain name:
    #  innerDomain 
    #  # grids in the domain:
    #  $innerMappingNames
    #done
    #
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    exit 
#    display intermediate results
# pause
#
  # open graphics
  compute overlap
# pause
  exit
save a grid (compressed)
$name
layeredDiskArrayGrid
exit



  for( $layer=0; $layer<$numLayers; $layer++ ) \
  { \
    $layerName="disk$count" . "layer$layer";    \
    $innerMappingNames .= "$layerName\n"; \
    $layerDomainName="disk$count" . "Layer$layer" . "Domain";    \
    $domainCommands  .=  "\n specify a domain\n $layerDomainName\n $layerName\n done"; \
    $innerRadius=$outerRadius; \
    $outerRadius=$innerRadius+$layerWidth[$layer]; \
    $nr=int( ($layerWidth[$layer])/$ds+1.5 ); \
    $sharePrev= $share; $share=$share+1; \
    $commands = $commands . \
    "*\n" . \
    " # Define layer  \n" . \
    "  annulus \n" . \
    "    mappingName \n" . \
    "      $layerName \n" . \
    "    inner and outer radii \n" . \
    "      $innerRadius $outerRadius\n" . \
    "    centre\n" . \
    "      $xc $yc\n" .   \
    "    lines \n" . \
    "      $nTheta $nr\n" . \
    "    boundary conditions \n" . \
    "      -1 -1 $bcInterface $bcInterface \n" . \
    "    share\n" . \
    "     * material interfaces are marked by share>=100\n" . \
    "      0 0 $sharePrev $share   \n" . \
    "    exit \n"; \
  }\
