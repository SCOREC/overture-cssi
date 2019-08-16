#
#  An arrangement of material disks 
#
# usage: 
# ogen [noplot] multiDiskGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -nCylx=<num> -nCyly=<num> -rad=<num> -dist=<num>
# 
#   nCylx : number of cylinders in the x direction
#   nCyly : number of cylinders in the y direction
#   rad   : radius of the cylinder
#   dist  : distance between centers
# 
# examples:
#     ogen -noplot multiDiskGrid -interp=e -order=2 -factor=1
#     ogen -noplot multiDiskGrid -interp=e -order=2 -factor=2
#     ogen -noplot multiDiskGrid -interp=e -order=2 -factor=4
#     ogen -noplot multiDiskGrid -interp=e -order=2 -factor=8
#     ogen -noplot multiDiskGrid -interp=e -order=2 -factor=16
# 
#     ogen -noplot multiDiskGrid -interp=e -order=4 -numGhost=3 -factor=4 
#     ogen -noplot multiDiskGrid -interp=e -order=4 -numGhost=3 -factor=8
# 
#
$prefix="multiDiskGrid"; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
# $nCylx=2; $nCyly=2;       # number of cylinders in each direction
# $rad=.3;                  # radius of each cylinder
# $dist=1.;                 # spacing between cylinder centres
$bc="1 1 -1 -1"; # for the backGround grid
# -- back-ground rectangle: 
$xa=-4.; $xb=-$xa;
$ya=-3.; $yb=-$ya;
# 
# get command line arguments
GetOptions( "rad=f"=>\$rad,"dist=f"=>\$dist,"order=i"=>\$order,"factor=i"=> \$factor,"interp=s"=> \$interp,\
             "nCylx=i"=> \$nCylx,"nCyly=i"=> \$nCyly,"prefix=s"=>\$prefix,\
	     "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,"numGhost=i"=> \$numGhost );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; }
$name = $prefix . "$interp$factor" . $suffix . ".hdf";
# 
$deltaY=$dist; $deltaX=$dist;   # spacing between cylinder centres
# 
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
$pi=4.*atan2(1.,1.); # 3.141592653;
$bcInterface=100;  # bc for interfaces
#
#
create mappings 
#
#
# ======================================================
# Define a function to build an inner-annulus, outer-annulus and inner square.
# usage:
#   makeDisk(radius,xCenter,yCenter)
# =====================================================
sub makeDisk\
{ local($radius,$xc,$yc)=@_; \
  $count = $count + 1; \
  $innerAnnulus=$innerDomain . "InnerAnnulus$count";     \
  $innerSquare =$innerDomain . "InnerSquare$count";     \
  $outerAnnulus=$outerDomain . "OuterAnnulus$count";     \
  $innerMappingNames = $innerMappingNames . "   $innerSquare\n". "   $innerAnnulus\n"; \
  $outerMappingNames = $outerMappingNames . "   $outerAnnulus\n"; \
  $nr = 2*($order-2) + 5; \
  $deltaRadius=$ds*$nr; \
  $outerRadius=$radius; $innerRadius=$outerRadius-$deltaRadius; \
  $nx=int( (2.*$pi*$outerRadius)/$ds+1.5 ); \
  $ny=int( ($deltaRadius)/$ds+1.5 ); \
  $nTheta=$nx; \
  $share=100+$count; \
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
  "      $nx $ny \n" . \
  "    exit \n"; \
  $extraWidth = ($order-2)*$ds;  \
  $xar=$xc-$innerRadius-$extraWidth;  $xbr=$xc+$innerRadius+$extraWidth; \
  $yar=$yc-$innerRadius-$extraWidth;  $ybr=$yc+$innerRadius+$extraWidth; \
  $nx=int( ($xbr-$xar)/$ds+1.5 ); \
  $ny=int( ($ybr-$yar)/$ds+1.5 ); \
  $commands = $commands . \
  "  rectangle \n" . \
  "    mappingName\n" . \
  "      $innerSquare\n" . \
  "    set corners\n" . \
  "     $xar $xbr $yar $ybr \n" . \
  "    lines\n" . \
  "      $nx $ny\n" . \
  "    boundary conditions\n" . \
  "      0 0 0 0\n" . \
  "    exit \n"; \
  $innerRadius=$outerRadius; \
  $outerRadius=$innerRadius+$deltaRadius; \
  $nx=$nTheta; \
  $ny=int( ($deltaRadius)/$ds+1.5 ); \
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
  "      $nx $ny\n" . \
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
# ===================================================
#   Make an circle of disks
#     makeDiskArray(radius,numCyl,circleRadiu)
#     
# Result: $commands
# =====================================================
sub makeDiskCircle \
{ local($radius,$numCyl,$circleRadius)=@_; \
  local $cmds; $cmds=""; \
  for( $j=0; $j<$numCyl; $j++ ){ \
    $theta = 2.*$pi*($j)/$numCyl; \
    $x0 = $circleRadius*cos($theta); $y0=$circleRadius*sin($theta); \
    makeDisk($radius,$x0+$i*$dx0,$y0+$j*$dy0); \
    $cmds = $cmds . $commands; \
  }\
  $commands=$cmds; \
}
#
#
# 
$x0=-($nCylx-1)*$deltaX*.5; $y0=-($nCyly-1)*$deltaY*.5;  
# makeDiskArray($rad,$nCylx,$nCyly,$x0,$deltaX,$y0,$deltaY);
#
# ========================OUTER DISK AND OUTER CIRCLE ===================================================================
# -- make a big disk to hold the smaller disks
#
$count=0; $innerMappingNames=""; $outerMappingNames=""; $innerDomain="bigDisk"; $outerDomain="backGround"; 
$radius=2.; $x0=0.; $y0=0.; 
makeDisk($radius,$x0,$y0);
$commands
$bigDiskNames     = $innerMappingNames; 
$backGroundNames  = $outerMappingNames; 
#
# --- make a circle of smaller disks
#
$innerMappingNames=""; $outerMappingNames="";  $innerDomain="bigDiskCircle"; $outerDomain="bigDisk"; 
$rad=.3; $circleRadius=1.5; $numCyl=12; 
makeDiskCircle($rad,$numCyl,$circleRadius);
$commands
$bigDiskCircleNames = $innerMappingNames; 
$bigDiskNames       = $bigDiskNames . $outerMappingNames; 
#
# ============================INNER DISK AND INNER CIRCLE ==========================================================
# -- make a medium size disk to hold the smaller disks
#
$innerMappingNames=""; $outerMappingNames=""; $innerDomain="mediumDisk"; $outerDomain="bigDisk"; 
$radius=1.; $x0=0.; $y0=0.; 
makeDisk($radius,$x0,$y0);
$commands
$mediumDiskNames = $innerMappingNames; 
$bigDiskNames    = $bigDiskNames . $outerMappingNames; 
#
# --- make a circle of smaller disks
#
$innerMappingNames=""; $outerMappingNames="";  $innerDomain="mediumDiskCircle"; $outerDomain="mediumDisk"; 
$rad=.2; $circleRadius=.7; $numCyl=8; 
makeDiskCircle($rad,$numCyl,$circleRadius);
$commands
$mediumDiskCircleNames = $innerMappingNames; 
$mediumDiskNames       = $mediumDiskNames . $outerMappingNames; 
#
# 
#
  rectangle 
    mappingName
      backGround
      $backGroundNames = "backGround" . "\n" . $backGroundNames;
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
printf("backGroundNames:\n$backGroundNames\n");
printf("bigDiskNames:\n$bigDiskNames\n");
printf("bigDiskCircleNames:\n$bigDiskCircleNames\n");
printf("mediumDiskNames:\n$mediumDiskNames\n");
printf("mediumiskCircleNames:\n$mediumDiskCircleNames\n");
# pause
#
generate an overlapping grid 
  # Here is a list of all grids: 
  $mappingNames= $backGroundNames . $bigDiskNames . $bigDiskCircleNames . $mediumDiskNames . $mediumDiskCircleNames; 
  $mappingNames
  done 
#
  change parameters 
    # ---- Specify domains -----
    specify a domain
      # domain name:
      outerDomain 
      # grids in the background domain:
      $backGroundNames 
     done
#
    specify a domain
      # domain name:
      bigDiskDomain 
      # grids in the domain:
      $bigDiskNames
    done
#
    specify a domain
      # domain name:
      bigDiskCircleDomain 
      # grids in the domain:
      $bigDiskCircleNames
    done
#
    specify a domain
      # domain name:
      mediumDiskDomain 
      # grids in the domain:
      $mediumDiskNames
    done
#
    specify a domain
      # domain name:
      mediumDiskCircleDomain 
      # grids in the domain:
      $mediumDiskCircleNames
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
#  open graphics 
#    display intermediate results
# pause
#
# 
   compute overlap
# pause
  exit
save a grid (compressed)
$name
multiDiskGrid
exit
