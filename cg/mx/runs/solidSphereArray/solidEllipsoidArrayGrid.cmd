#
#  Grid for an ARRAY of SOLID ELLIPSOIDs in a box 
#
echo to terminal 0
# usage: ogen [noplot] solidEllipsoidArrayGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -nrExtra=<>
#                 -configFile=>s>-rgd=[fixed|var] -nsx=<i> -nsy=<i> -nsz=<i> -ml=<> -stretchFactor=<> -per=[nnn|ppp|npp]
#
#  -configFile=>s> : configuration file with ellipsoid parameters (see examples below)
#  -nsx, -nsy, -nsz : build an array of solid ellipsoids, this many in each direction
#  -nrExtra: extra lines to add in the radial direction on the sphere grids 
#  -rgd : var=variable : decrease radial grid distance as grids are refined. fixed=fix radial grid distance
#  -ml = number of (extra) multigrid levels to support
#
# EXAMPLES:
# 
# E18 
#    ogen -noplot solidEllipsoidArrayGrid -prefix=18SolidEllipsoidsGrid -nsx=3 -nsy=2 -nsz=3 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -factor=4
# 
# E12 - OK 
# ogen -noplot solidEllipsoidArrayGrid -prefix=12SolidEllipsoidsGrid -nsx=2 -nsy=2 -nsz=3 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -factor=4
# 
# E4 
# ogen -noplot solidEllipsoidArrayGrid -prefix=4SolidEllipsoidsGrid -nsx=2 -nsy=2 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -factor=4
#
# New way: use a config file
$configFile = "solidEllipsoidArrayConfig3d.h";  # default configuration file 
# 
$nsx=1; $nsy=1; $nsz=1; # number of solids along each axis 
# sphere: 
$a=1/2.; $b=1/2.; $c=1./2.;
# 
$xa=-2.5; $xb=2.5; $ya=-.75; $yb=.75; $za=-.75; $zb=.75; $nrMin=5; $nrExtra=0; $rgd="var"; $name=""; $per="nnn"; 
# specify up to two rotations and a shift 
$angle1=0; $rotationAxis1=0; $angle2=0; $rotationAxis2=1; $xShift=0; $yShift=0; $zShift=0; 
# 
$box=0.; # if non zero on input then set -xa=xb=-ya=yb=-za=zb=box
$order=2; $factor=1; $interp="e"; $interpType = "explicit for all grids";
$ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $dse=0.;
$ds0=.1; # target grid spacing for $factor=1 
$orthographicPatchParameter=.65; # parameter defining the width of the orthographic patch
$stretchFactor=4.; # stretch grid lines by this factor at the sphere boundary
$deltaRadius0=.25; # do not make larger than .3 or troubles with cgmx
$suffix=""; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
# maximize overlap so the mask matches on the interface -- need to fix ogen to do this automatically
$overlapOption="maximize"; # "default"; 
$prefix="solidEllipsoidArrayGrid";
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=i"=>\$factor,"nrExtra=i"=>\$nrExtra,"nrMin=i"=>\$nrMin,\
            "interp=s"=>\$interp,"rgd=s"=>\$rgd,"deltaRadius0=f"=>\$deltaRadius0,"name=s"=>\$name,\
            "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,"za=f"=>\$za,"zb=f"=>\$zb,"ml=i"=>\$ml,\
            "stretchFactor=f"=>\$stretchFactor,"box=f"=>\$box,"suffix=s"=>\$suffix,"numGhost=i"=>\$numGhost,\
            "prefix=s"=> \$prefix,"a=f"=>\$a,"b=f"=>\$b,"c=f"=>\$c,"per=s"=>\$per,"ds0=f"=>\$ds0,\
	    "orthographicPatchParameter=f"=>\$orthographicPatchParameter,"angle1=f"=>\$angle1,"angle2=f"=>\$angle2, \
            "rotationAxis1=i"=>\$rotationAxis1,"rotationAxis2=i"=>\$rotationAxis2,\
            "xShift=f"=>\$xShift,"yShift=f"=>\$yShift,"zShift=f"=>\$zShift,"configFile=s"=>\$configFile,\
            "overlapOption=s"=>\$overlapOption,"nsx=i"=>\$nsx,"nsy=i"=>\$nsy,"nsz=i"=>\$nsz );
# 
if( $box ne 0 ){ $xa=-$box; $xb=$box; $ya=-$box; $yb=$box; $za=-$box; $zb=$box; }
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "i" ){ $interpType = "implicit for all grids"; $dse=1.; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; $dse=1.; }
# 
if( $rgd eq "fixed" ){ $prefix = $prefix . "Fixed"; $sphereWidth=$deltaRadius0; }else{ $sphereWidth=-1.; }
$suffix .= ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){ $name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=$ds0/$factor;
# 
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
#
#
# ---------------------------------------
# turn off graphics
# ---------------------------------------
$dw = $order+1;
$iw = $dw;
$parallelGhost=($dw+1)/2;
if( $interp eq "e" ){ $parallelGhost=($dw+$iw-2)/2; }
minimum number of distributed ghost lines
  $parallelGhost
#
$pi=4.*atan2(1.,1.);
# 
# Include the configuration file
#
include $configFile
#
# 
# ------------------------ CREATE MAPPINGS ----------------------
# 
create mappings
#
# number of points to use in the radial direction : $nrExtra is used for stretching 
$nr=$nrMin + ($order-2); 
if( $interp eq "e" ){ $nr=$nr+$order+$nrExtra; } 
# the coarsest MG grid is 4 pts
$nr = max( $nr, 2**($ml+2) ); 
#
## $gridNames="*"; # IS THIS NEEED ??
#
# Here is the background box
#
Box
  set corners
    $xa $xb $ya $yb $za $zb
  lines
    $nx = intmg( ($xb-$xa)/$ds +1.5);
    $ny = intmg( ($yb-$ya)/$ds +1.5);
    $nz = intmg( ($zb-$za)/$ds +1.5);
    $nx $ny $nz
  boundary conditions
    $bcCmd = "1 2 3 4 5 6"; 
    if( $per eq "ppp" ){ $bcCmd = "-1 -1 -1 -1 -1 -1"; }
    if( $per eq "npp" ){ $bcCmd = "1 2 -1 -1 -1 -1"; }
    $bcCmd
  mappingName
    backGround
  exit
# 
#    Perl routine to convert a mapping into a Nurbs mapping and optionally rotate and shift
#
sub convertToNurbs\
{ local($old,$new,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift)=@_; \
  $cmds = "nurbs \n" . \
   "interpolate from mapping with options\n" . \
   " $old \n" . \
   " parameterize by index (uniform)\n" . \
   " number of ghost points to include\n $numGhost\n" . \
   " choose degree\n" . \
   "  3 \n" . \
   "done\n" . \
   "rotate \n" . \
   " $angle1 $rotationAxis1 \n" . \
   " 0. 0. 0.\n" . \
   "rotate \n" . \
   " $angle2 $rotationAxis2 \n" . \
   " 0. 0. 0.\n" . \
   "shift\n" . \
   " $xShift $yShift $zShift\n" . \
   "mappingName\n" . \
   " $new\n" . \
   "exit"; \
}
# ==========================================================================================
# domain 0: all outer grids
# domain 1: ellipsoid 1 inner grids 
# domain 2: ellipsoid 2 inner grids
# domain n: ellipsoid n inner grids
# ==========================================================================================
#
$numDomains=1; 
@domainName = (); # list of domains names
@domainGrids = (); # list of grids in a domain
#
$domainName[0] = "outerDomain";
$domainGrids[0] = "backGround";
#
# 
#
# ---- NOTES on include files ----
#   
#  buildMasterEllipsoidGrids.h
#    include ellipsoidSurfacePatch.h
#         - set $rDist 
#         - rest $nr for inner grids  based on $nrMin
#    include ellipsoidVolumeGrids.h
#        - generate hyperbolic grids (interior or exterior)
# 
#  constructEllipsoidGrids.h
#    - set shifts and rotations 
#    include buildEllipsoidGridsFromMaster.h
#              - convert to Nurbs
#              - build inner background grid 
# 
# --- Build a string containing the commands to shift and rotate the master grids ---
#
# NOTE: includes are processed in REVERSE ORDER 
$cmd=""; 
for( $i=0; $i<$numEllipsoids; $i++){ \
  $cmd .= "include \$ENV{CG}/mx/runs/solidSphereArray/constructEllipsoidGrids.h\n"; \
  $cmd .= "include \$ENV{CG}/mx/runs/solidSphereArray/buildMasterEllipsoidGrids.h\n"; \
  }
$cmd .="#";  
#
$masterCount=-1; # counts master ellipsoids 
$count=0; # counts ellipsoids 
$cmd
# 
# open graphics
## include $ENV{CG}/mx/runs/solidSphereArray/buildMasterEllipsoidGrids.h
# 
## include $ENV{CG}/mx/runs/solidSphereArray/constructEllipsoidGrids.h
#- #
#- #
#- #  --------ellipsoid 1 ---------
#- $angle1=0; $rotationAxis10; $angle2=0; $rotationAxis2=1; $xShift=.75; $yShift=0; $zShift=0; 
#- # 
#- include $ENV{CG}/mx/runs/solidSphereArray/buildEllipsoidGridsFromMaster.h
#- # -----------------------------
#- #  --------ellipsoid 2 ---------
#- $angle1=0; $rotationAxis10; $angle2=0; $rotationAxis2=1; $xShift=-.75; $yShift=0; $zShift=0; 
#- # 
#- include $ENV{CG}/mx/runs/solidSphereArray/buildEllipsoidGridsFromMaster.h
#- # -----------------------------
#- # 
exit
#
generate an overlapping grid
  $domainGrids[0]
  # make a list of all ellipsoid grids 
  $gridNames=""; for( $d=1; $d<$numDomains; $d++){ $gridNames .= $domainGrids[$d] . "\n"; } $gridNames .= "#"; 
  $gridNames 
#-  $domainGrids[1]
#-  $domainGrids[2]
#-   backGround
#-   ellipsoidOuterGrid
#-   southPoleOuterGrid
#-   northPoleOuterGrid
#-   # inside the ellipsoid: 
#-   backGroundInner
#-   ellipsoidInnerGrid
#-   southPoleInnerGrid
#-   northPoleInnerGrid 
 done
# 
  change parameters
# 
    specify a domain
      # --- outer domain name:
      $domainName[0]
      $domainGrids[0]
      done
    # make a list of all domain commands:
    $cmd ="";
    for( $d=1; $d<$numDomains; $d++){\
      $cmd .= "specify a domain\n $domainName[$d]\n $domainGrids[$d]\n done\n"; }
    $cmd .="#";
    #
    # -- define the domains ---
    #
    $cmd 
#-#
#-     specify a domain
#-        # --- inner domain name:
#-       $domainName[1]
#-       $domainGrids[1]
#-       done
#-# 
#-     specify a domain
#-       $domainName[2]
#-       $domainGrids[2]
#-      done
# 
   interpolation type
      $interpType
      printf("interpType=$interpType\n");
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    # 
    # Maximum overlap to ensure the mask's match at the interface *fix me*
    if( $overlapOption eq "maximize" ){ $cmd="maximize overlap"; }else{ $cmd="#"; }
      $cmd 
  exit
#
echo to terminal 1
  compute overlap
 #  change the plot
 #    toggle grid 0 0
    # open graphics
 #    pause
  # exit
#
exit
# save an overlapping grid
save a grid (compressed)
$name
solidEllipsoidArrayGrid
exit


