#
# Define a 4 box material region for the BA MX solver
# 
forcing options...
define material region...
material file: $matFile2
  box: -.45 -.05 -.35 -.05 0 1 (xa,xb,ya,yb,za,zb)
  box:  .05  .35  .05  .35 0 1 (xa,xb,ya,yb,za,zb)
  box:  .45   .8  -.35 -.05  0 1 (xa,xb,ya,yb,za,zb)
# 
material file: $matFile3
  box: -.45 -.05  .05 .35 0 1 (xa,xb,ya,yb,za,zb)
  box:  .05  .35 -.35 -.05 0 1 (xa,xb,ya,yb,za,zb)
  box:  .45   .8   .05 .35 0 1 (xa,xb,ya,yb,za,zb)
continue
