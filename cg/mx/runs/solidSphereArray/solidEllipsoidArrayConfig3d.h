#
# Configuration file for solidEllipsoidArrayGrid.cmd
#
#    Arbitrary 3D array of ellipsoids:
#         $nsx by $nsy by $nsz 
#     
# $nsx=2; $nsy=2; $nsz=2;
$numEllipsoids=$nsx*$nsy*$nsz;
$edist = 1.3;  # distance between centres
#
# ------ first make a list of difference ellipsoid shapes --------
#   Note: probably best to have c=longest semi-axis so end patches are smaller
#
$i=0;
# flat ellipsoid 
$aee[$i]= .55;  $bee[$i]= .5;   $cee[$i]=.40; $rot1e[$i]=  0; $rot2e[$i]= 90; $i=$i+1; 
$aee[$i]= .375; $bee[$i]= .5;   $cee[$i]=.5;  $rot1e[$i]= 20; $rot2e[$i]=  0; $i=$i+1; 
$aee[$i]= .55;  $bee[$i]= .4;   $cee[$i]=.4;  $rot1e[$i]= 45; $rot2e[$i]= 45; $i=$i+1;  
$aee[$i]= .55;  $bee[$i]= .4;   $cee[$i]=.4;  $rot1e[$i]=  0; $rot2e[$i]=  0; $i=$i+1;  
#
# thin ellipsoid
$aee[$i]= .375; $bee[$i]= .375; $cee[$i]=.5;  $rot1e[$i]=  0; $rot2e[$i]=  0; $i=$i+1; 
$aee[$i]= .375; $bee[$i]= .5;   $cee[$i]=.5;  $rot1e[$i]= 20; $rot2e[$i]=  0; $i=$i+1; 
# flat ellipsoid 
$aee[$i]= .50;  $bee[$i]= .5;   $cee[$i]=.35; $rot1e[$i]=  0; $rot2e[$i]=  0; $i=$i+1;  
$aee[$i]= .55;  $bee[$i]= .4;   $cee[$i]=.4;  $rot1e[$i]=  0; $rot2e[$i]=  0; $i=$i+1;  
#
# thin ellipsoid
$aee[$i]= .55;  $bee[$i]= .5;   $cee[$i]=.40; $rot1e[$i]=  0; $rot2e[$i]= 90; $i=$i+1; 
$aee[$i]= .375; $bee[$i]= .5;   $cee[$i]=.5;  $rot1e[$i]= 20; $rot2e[$i]=  0; $i=$i+1; 
$aee[$i]= .55;  $bee[$i]= .4;   $cee[$i]=.4;  $rot1e[$i]= 45; $rot2e[$i]= 45; $i=$i+1;  
$aee[$i]= .55;  $bee[$i]= .4;   $cee[$i]=.4;  $rot1e[$i]=  0; $rot2e[$i]=  0; $i=$i+1;  
$numEllipsoidTypes=$i;
#
# --- Cycle through the different ellipsoid types ---
# 
$x0=0; $y0=0; $z0=0;  # offsets for first ellipsoid
$xeMin=1e12; $xeMax=-$xeMin; $yeMin=$xeMin; $yeMax=-$yeMin; $zeMin=$xeMin; $zeMax=-$zeMin; 
$i=0;  # counts ellipsoids
for( $iz=0; $iz<$nsz; $iz++ ){\
for( $iy=0; $iy<$nsy; $iy++ ){\
for( $ix=0; $ix<$nsx; $ix++ ){\
    $xe[$i]=$x0+($ix-.5*($nsx-1))*$edist; $ye[$i]=$y0+($iy-.5*($nsy-1))*$edist; $ze[$i]=$z0+($iz-.5*($nsz-1))*$edist; \
    $xeMin = min($xeMin,$xe[$i]);  $xeMax = max($xeMax,$xe[$i]);  \
    $yeMin = min($yeMin,$ye[$i]);  $yeMax = max($yeMax,$ye[$i]);  \
    $zeMin = min($zeMin,$ze[$i]);  $zeMax = max($zeMax,$ze[$i]);  \
    $j = $i % $numEllipsoidTypes; \
    $ae[$i]=$aee[$j]; $be[$i]=$bee[$j]; $ce[$i]=$cee[$j]; $rot1[$i]=$rot1e[$j]; $rot2[$i]=$rot2e[$j]; \
    $i=$i+1;  \
   }}}
#
printf(" Bounds on ellipsoids: xeMin=$xeMin, xeMax=$xeMax, yeMin=$yeMin, yeMax=$yeMax, zeMin=$zeMin, zeMax=$zeMax\n");
# 
# Box bounds: 
# $xa=-1.5; $xb=1.5; $ya=-1.5; $yb=1.5; $za=-1.5; $zb=1.5;
$de = .5*$edist; # half distance between ellipsoids
$xa=$xeMin-$de-.25; $xb=$xeMax+$de+.25; $ya=$yeMin-$de-.125; $yb=$yeMax+$de+.125; $za=$zeMin-$de-.125; $zb=$zeMax+$de+.125;
