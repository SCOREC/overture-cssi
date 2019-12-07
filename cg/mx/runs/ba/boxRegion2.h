# Define a 2 box material region for the BA MX solver
forcing options...
define material region...
material file: $matFile2
  box: .4 .8 .4 .8 0 1 (xa,xb,ya,yb,za,zb)
# 
material file: $matFile3
  box: .5 .9 .1 .6 0 1 (xa,xb,ya,yb,za,zb)
continue
