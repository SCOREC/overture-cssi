#
#  Array of solid sphere in a box 
#
# usage: ogen [-noplot] solidSphereArrayGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -nrExtra=<> ...
#         -rgd=[fixed|var] -outerBC=[dirichlet|yzPeriodic]
#
#  nrExtra: extra lines to add in the radial direction on the sphere grids 
#  -rgd : var=variable : decrease radial grid distance as grids are refined. fixed=fix radial grid distance
#  -outerBC : 
# 
# examples:
#     ogen -noplot solidSphereArrayGrid -order=2 -factor=1 
# 
# One sphere: 
#  ogen -noplot solidSphereArrayGrid -prefix=solidSphereGrid -xa=-4 -xb=4 -ya=-3 -yb=3 -za=-3 -zb=3 - -nsx=1 -nsy=1 -nsz=1 -interp=e -order=2 -factor=1
#  ogen -noplot solidSphereArrayGrid -prefix=solidSphereGrid -xa=-4 -xb=4 -ya=-3 -yb=3 -za=-3 -zb=3 - -nsx=1 -nsy=1 -nsz=1 -interp=e -order=2 -factor=2
# 
#   2 x 2 x 2 = 8 spheres: 
#  ogen -noplot solidSphereArrayGrid -nsx=2 -nsy=2 -nsz=2 -order=2 -factor=1 -go=og  
# 
#  3 x 3 x 3 = 27 spheres: 
#  ogen -noplot solidSphereArrayGrid -nsx=3 -nsy=3 -nsz=3 -order=2 -factor=1 -go=og  
# 
#  ogen -noplot solidSphereArrayGrid -nsx=3 -nsy=3 -nsz=3 -order=4 -factor=2 -interp=e -go=og  
#
# OLD: 
# --implicit:
#     ogen -noplot solidSphereArrayGrid -order=4 -factor=8 -nrMin=11 
# 
$prefix="solidSphereArrayGrid";
$xa=-4.; $xb=4.; $ya=-4; $yb=4; $za=-4.; $zb=4.; $nrMin=3; $nrExtra=0; $rgd="var"; $name=""; 
$order=2; $factor=1; $interp="e"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $dse=0.; 
$deltaRadius0=.3; # do not make larger than .3 or troubles with cgmx
$numGhost=-1;  # if this value is set, then use this number of ghost points
$outerBC="yzPeriodic"; # "dirichlet"; 
# 
$nsx=3; $nsy=3; $nsz=3; # number of spheres in x, y and z directions
$dist=2.5; # distance between spheres
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=i"=> \$factor,"nrExtra=i"=>\$nrExtra,"nrMin=i"=>\$nrMin,\
            "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,"za=f"=>\$za,"zb=f"=>\$zb,\
            "interp=s"=> \$interp,"rgd=s"=> \$rgd,"deltaRadius0=f"=>\$deltaRadius0,"name=s"=>\$name,\
            "numGhost=i"=>\$numGhost,"outerBC=s"=>\$outerBC,"prefix=s"=>\$prefix,\
	    "nsx=i"=>\$nsx,"nsy=i"=>\$nsy,"nsz=i"=>\$nsz );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; $dse=1.; }
# 
if( $rgd eq "fixed" ){ $prefix = $prefix . "Fixed"; }
$suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $name eq "" ){ $name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
#
# Here is the radial width of the spherical grids -- this will be fixed if rgd=fixed
# matching interface grids should be given distinct share values for now for cgmx -- fix me -- cgmp is different
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
create mappings
$pi=4.*atan2(1.,1.);
#
#
#
# Here is the backGround grid 
#
Box
  set corners
    $xa $xb $ya $yb $za $zb
  lines
    $nx = int( ($xb-$xa)/$ds +1.5);
    $ny = int( ($yb-$ya)/$ds +1.5);
    $nz = int( ($zb-$za)/$ds +1.5);
    $nx $ny $nz
  boundary conditions
    $cmd="1 2 3 4 5 6";
    if( $outerBC eq "yzPeriodic"){ $cmd="1 2 -1 -1 -1 -1"; } # periodic BC in y and z 
    $cmd
  mappingName
    outerBox
  exit
#
#
# number of points to use in the radial direction : $nrExtra is used for stretching 
$nr=$nrMin + $order; if( $interp eq "e" ){ $nr=$nr+$order+$nrExtra; } 
# fixed radial distance
if( $rgd eq "fixed" ){ $nr = $nr*($factor-1); }
#
$count=0;  # counts spheres and/or individual domains
# $innerGridNames="#"; 
@domainGridNames = ();
# domain=0 is the outer domain 
# 
$domainName[0] = "outerDomain"; 
$domainGridNames[0] = "outerBox"; 
#
$sphereRadius=1.; $radiusDir=1; 
$xSphere=0.; $ySphere=-1.25; $zSphere=0.; 
#
# ---- Loop to generate spheres
#
#  Array of sphere centres: 
@xv=();
@yv=();
@zv=();
$x0=0; $y0=0; $z0=0;  # offsets for first sphere
$ns=0;  # counts spheres
for( $iz=0; $iz<$nsz; $iz++ ){\
for( $iy=0; $iy<$nsy; $iy++ ){\
for( $ix=0; $ix<$nsx; $ix++ ){\
    $xv[$ns]=$x0+($ix-.5*($nsx-1))*$dist; $yv[$ns]=$y0+($iy-.5*($nsy-1))*$dist; $zv[$ns]=$z0+($iz-.5*($nsz-1))*$dist; $ns=$ns+1;\
   }}}
# add one more at center
### $xv[$ns]=0.; $yv[$ns]=0.; $zv[$ns]=0.; $ns=$ns+1;
#
$cmd=""; 
for( $i=0; $i<$ns; $i++ ){\
  $cmd .= "include constructSolidSphere.h\n"; \
 }
$cmd .="#";
# Execute commands to build sphere grids:
$cmd
# 
#- #  -------------- SPHERE 1 ---------------
#- $count = $count+1; 
#- $sphereRadius=1.; $radiusDir=1; 
#- $xSphere=0.; $ySphere=-1.25; $zSphere=0.; 
#- # 
#- include buildSolidSphere.h
#- # 
#- $domainName[$count] = "sphereDomain$count"; 
#- $domainGridNames[$count] = $innerNames; 
#- $domainGridNames[0] .= "\n" . $outerNames; 
#- # 
#- #  -------------- SPHERE 2 ---------------
#- $count = $count+1; 
#- $sphereRadius=1.; $radiusDir=1; 
#- $xSphere=2.5; $ySphere=-1.25; $zSphere=0.; 
#- # 
#- include buildSolidSphere.h
#- # 
#- $domainName[$count] = "sphereDomain$count"; 
#- $domainGridNames[$count] = $innerNames; 
#- $domainGridNames[0] .= "\n" . $outerNames; 
#- #
#- #  -------------- SPHERE 3 ---------------
#- $count = $count+1; 
#- $sphereRadius=1.; $radiusDir=1; 
#- $xSphere=0; $ySphere=1.25; $zSphere=0.; 
#- # 
#- include buildSolidSphere.h
#- # 
#- $domainName[$count] = "sphereDomain$count"; 
#- $domainGridNames[$count] = $innerNames; 
#- $domainGridNames[0] .= "\n" . $outerNames; 
#- #
exit
#
$numDomains = $count+1; 
generate an overlapping grid
  $allGridNames="#"; 
  for( $d=0; $d<$numDomains; $d++ ){ $allGridNames .= "\n" . $domainGridNames[$d]; }
  $allGridNames
#  $domainGridNames[0]
#  $domainGridNames[1]
#  $domainGridNames[2]
  done
# 
  change parameters
    # Make a string to define all domains: 
    $cmd=""; 
    for( $d=0; $d<$numDomains; $d++ ){\
      $cmd .= "specify a domain\n"; \
      $cmd .= "  $domainName[$d]\n"; \
      $cmd .= "  $domainGridNames[$d]\n"; \
      $cmd .= "done\n"; }
    $cmd .= "#";
    # Now execute the string to define all domains 
    $cmd 
    # 
#-     specify a domain
#-       # domain name:
#-       $domainName[1]
#-        # grids in the domain:
#-        $domainGridNames[1]
#-     done
#-     # 
#-     specify a domain
#-       # domain name:
#-       $domainName[2]
#-        # grids in the domain:
#-        $domainGridNames[2]
#-     done
# 
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
  exit
#
  # open graphics 
  compute overlap
#
exit
# save an overlapping grid
save a grid (compressed)
$name
solidSphereArrayGrid
exit
