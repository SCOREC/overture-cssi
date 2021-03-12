# Define a few spheres and ellipsoids 
forcing options...
define material region...
material file: $matFile2
$xc1=-1; $yc1=-1; $zc1=1; $ae1=.5; $be1=.5; $ce1=.5;
ellipsoid: $xc1 $yc1 $zc1 $ae1 $be1 $ce1 (x0,y0,z0, a,b,c)
# 
material file: $matFile2
$xc2=0; $yc2=0; $zc2=0; $ae2=1.; $be2=.5; $ce2=.75;
ellipsoid: $xc2 $yc2 $zc2 $ae2 $be2 $ce2 (x0,y0,z0, a,b,c)
# 
material file: $matFile2
$xc3=1; $yc3=.5; $zc3=.5; $ae3=.5; $be3=1; $ce3=1;
ellipsoid: $xc3 $yc3 $zc3 $ae3 $be3 $ce3 (x0,y0,z0, a,b,c)
continue
