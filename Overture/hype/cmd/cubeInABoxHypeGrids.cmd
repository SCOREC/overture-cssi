* ===================================================================
*
*   Grids for a cube in a box 
* 
* ==================================================================
*
create mappings 
* 
*
*    Build body fitted grids around the exterior of a cube
*
* Make planes to form the side of a cube
*
  plane or rhombus
    specify plane or rhombus by three points
      0. 0. 0. 1. 0. 0. 0. 1. 0
    lines
      5 5 
    mappingName
      frontPlane
    exit
*
  plane or rhombus
    specify plane or rhombus by three points
      1. 0. 0. 1. 0. -1 1. 1. 0
    lines
      5 5 
    mappingName
      rightPlane
    exit
*
  plane or rhombus
    specify plane or rhombus by three points
      0. 1. 0. 1. 1. 0. 0. 1. -1
    lines
      5 5 
    mappingName
      topPlane
    exit
*
*
  plane or rhombus
    specify plane or rhombus by three points
      0. 0. -1.  0. 0. 0.  0. 1. -1.
    lines
      5 5 
    mappingName
      leftPlane
    exit
*
    plane or rhombus
      specify plane or rhombus by three points
        0. 0. 0.  0. 0. -1. 1. 0. 0. 
      lines
        5 5 
      mappingName
        bottomPlane
      exit
*
   plane or rhombus
     specify plane or rhombus by three points
       0. 0. -1.  0. 1. -1.  1. 0. -1. 
     lines
       5 5 
     mappingName
       backPlane
     exit
*
*  Join the planes into a single CAD surface
*
  composite surface
    add all mappings
    mappingName
    corner
* pause
    determine topology
      deltaS .1
      build edge curves
      merge edge curves
      triangulate
* pause
      exit
    exit
  builder
*
* 4-edge surface grid around front
*
    create surface grid... 
      choose edge curve 1 
      choose edge curve 2 
      choose edge curve 3 
      choose edge curve 0 
      done 
* pause
      points on initial curve 41
* 
      picking:choose interior matching curve
      choose edge curve 8 0.000000e+00 1.000000e+00 -5.555556e-01 
      done
      choose edge curve 9 0.000000e+00 0.000000e+00 -5.555556e-01 
      done
      choose edge curve 4 1.000000e+00 0.000000e+00 -4.444444e-01 
      done
      choose edge curve 6 1.000000e+00 1.000000e+00 -5.555556e-01 
      done
*
      forward and backward
      distance to march .15 .15  (forward,backward)
*      points on initial curve 41
      marching options...
      uniform dissipation 0.2
      spacing: geometric
      geometric stretch factor 1.1
      lines to march 9, 9 
      generate
   pause
      exit
*
* 4-edge volume grid around front
*
    create volume grid...
      backward
      distance to march .2
      * target grid spacing $ds $ds (tang,normal, <0 : use default)
      lines to march 9
*       decrease uniform dissipation near very sharp corners!
      uniform dissipation 0.01
      spacing: geometric
      geometric stretch factor 1.1
      volume smooths 5
      generate
      Share Value: back    1
      use robust inverse
  pause
      exit
*
* 4-edge surface grid around back
*
    create surface grid...
      choose edge curve 7 5.555556e-01 1.000000e+00 -1.000000e+00 
      choose edge curve 10 0.000000e+00 5.555556e-01 -1.000000e+00 
      choose edge curve 11 4.444444e-01 0.000000e+00 -1.000000e+00 
      choose edge curve 5 1.000000e+00 4.444444e-01 -1.000000e+00 
      done
* 
      points on initial curve 41
*
      points on initial curve 41
      picking:choose interior matching curve
      choose edge curve 6 1.000000e+00 1.000000e+00 -5.555556e-01 
      done
      choose edge curve 4 1.000000e+00 0.000000e+00 -4.444444e-01 
      done
      choose edge curve 8 0.000000e+00 1.000000e+00 -5.555556e-01 
      done
      choose edge curve 9 0.000000e+00 0.000000e+00 -5.555556e-01 
      done
      forward and backward
      distance to march .15 .15  (forward,backward) 
      lines to march 9 9
      marching options...
      uniform dissipation 0.2
      spacing: geometric
      geometric stretch factor 1.1
      close marching options
      generate
*  pause
      exit
*
* 4-edge volume grid around back
*
    create volume grid...
      backward
      distance to march .2
      lines to march 9
*       decrease uniform dissipation near very sharp corners!
      uniform dissipation 0.01
      spacing: geometric
      geometric stretch factor 1.1
      volume smooths 5
      generate
      Share Value: back    1
      use robust inverse
*  pause
      exit
*
*  upper right edge surface grid
*
    create surface grid...
      choose edge curve 6
      done
      forward and backward
      BC: left (forward) fix z, float x and y
      BC: right (forward) fix z, float x and y
      BC: left (backward) fix z, float x and y
      BC: right (backward) fix z, float x and y
      spacing: geometric
      geometric stretch factor 1.1 
      distance to march 0.15 .15  (forward,backward) 
      lines to march 9 9
      points on initial curve 19
      Start curve parameter bounds .05 .95
      generate
*  pause
      exit
*
* upper right edge volume grid
*
    create volume grid...
      distance to march .2
      lines to march 7
      uniform dissipation 0.05
      backward
      generate
      Share Value: back    1
*  pause
      exit
*
*  lower right edge surface grid
*
    create surface grid...
      choose edge curve 4
      done
      forward and backward
      BC: left (forward) fix z, float x and y
      BC: right (forward) fix z, float x and y
      BC: left (backward) fix z, float x and y
      BC: right (backward) fix z, float x and y
      spacing: geometric
      geometric stretch factor 1.1 
      distance to march 0.15 .15  (forward,backward) 
      lines to march 9 9
      points on initial curve 19
      Start curve parameter bounds .05 .95
      generate
*  pause
      exit
*
* lower right edge volume grid
*
    create volume grid...
      distance to march .2
      lines to march 7
      uniform dissipation 0.05
      backward
      generate
*  pause
      Share Value: back    1
      exit
*
*  lower left edge surface grid
*
    create surface grid...
      choose edge curve 9 0.000000e+00 0.000000e+00 -5.555556e-01
      done
      forward and backward
      BC: left (forward) fix z, float x and y
      BC: right (forward) fix z, float x and y
      BC: left (backward) fix z, float x and y
      BC: right (backward) fix z, float x and y
      spacing: geometric
      geometric stretch factor 1.1 
      distance to march 0.15 .15  (forward,backward) 
      lines to march 9 9
      points on initial curve 19
      Start curve parameter bounds .05 .95
      generate
*  pause
      exit
*
* lower left edge volume grid
*
    create volume grid...
      distance to march .2
      lines to march 7
      uniform dissipation 0.05
      backward
      generate
      Share Value: back    1
*  pause
      exit
*
*  upper left edge surface grid
*
    create surface grid...
      choose edge curve 8 0.000000e+00 1.000000e+00 -5.555556e-01
      done
      forward and backward
      BC: left (forward) fix z, float x and y
      BC: right (forward) fix z, float x and y
      BC: left (backward) fix z, float x and y
      BC: right (backward) fix z, float x and y
      spacing: geometric
      geometric stretch factor 1.1 
      distance to march 0.15 .15  (forward,backward) 
      lines to march 9 9
      points on initial curve 19
      Start curve parameter bounds .05 .95
      generate
*  pause
      exit
*
* lupper left edge volume grid
*
    create volume grid...
      distance to march .2
      lines to march 7
      uniform dissipation 0.05
      backward
      generate
      Share Value: back    1
*  pause
      exit
*
exit
exit    


