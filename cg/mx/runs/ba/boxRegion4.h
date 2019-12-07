#
# Define a 4 box material region for the BA MX solver
# 
forcing options...
define material region...
material file: $matFile2
  box: -.3 .0 -.3 .0 0 1 (xa,xb,ya,yb,za,zb)
# 
material file: $matFile2
  box: -.3 .0 .0 .3 0 1 (xa,xb,ya,yb,za,zb)
# 
material file: $matFile2
  box: .0 .3 -.3 .0 0 1 (xa,xb,ya,yb,za,zb)
# 
material file: $matFile2
  box: .0 .3 .0 .3 0 1 (xa,xb,ya,yb,za,zb)
continue
