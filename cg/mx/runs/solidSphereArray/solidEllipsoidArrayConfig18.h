#
# Configuration file for solidEllipsoidArrayGrid.cmd
#
#    Arbitrary 3D array of elliposids:
#         $nex by $ney by $nez 
#     
$nex=2; $ney=2; $nez=2;
$numEllipsoids=$nex*$ney*$nez;
$edist = 1.3;  # distance between centres
# Box bounds: 
$xa=-1.5; $xb=1.5; $ya=-1.5; $yb=1.5; $za=-1.5; $zb=1.5;
#
# ------ list of difference ellipsoid parameters: --------
#   Note: probably best to have c=longest semi-axis so end patches are smaller
#
$i=0;
$de=.5*$edist; # half distance between ellipsoids 
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
$i=0;  # counts ellipsoids
for( $iz=0; $iz<$nez; $iz++ ){\
for( $iy=0; $iy<$ney; $iy++ ){\
for( $ix=0; $ix<$nex; $ix++ ){\
    $xe[$i]=$x0+($ix-.5*($nex-1))*$edist; $ye[$i]=$y0+($iy-.5*($ney-1))*$edist; $ze[$i]=$z0+($iz-.5*($nez-1))*$edist; \
    $j = $i % $numEllipsoidTypes; \
    $ae[$i]=$aee[$j]; $be[$i]=$bee[$j]; $ce[$i]=$cee[$j]; $rot1[$i]=$rot1e[$j]; $rot2[$i]=$rot2e[$j]; \
    $i=$i+1;  \
   }}}
