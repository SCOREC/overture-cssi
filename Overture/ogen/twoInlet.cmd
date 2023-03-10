*
* create a grid to demonstrate various features
*
create mappings
  * make a back ground grid
  rectangle
    specify corners
      0 0 2. 1.
    lines
      61 31  
    mappingName
     backGroundGrid
    share
      1 2 3 4
    pause
  exit
  * make an annulus
  Annulus
    centre for annulus
      1. .5
    inner radius
     .2
    outer radius
     .4
    lines
      41 9
    mappingName
      annulus
    boundary conditions
      -1 -1 1 0
    pause
  exit
  * the inlet (on the right) will consist of two 
  * smoothed polygons
  SmoothedPolygon
    mappingName
      inlet-top
    vertices
     3
     2. .85
     2. .65 
     2.25 .65 
   n-dist
     fixed normal distance
     -.175  .2
   sharpness
     10.
     10.
     10.
   t-stretch
     0. 10.
     1. 10.
     0. 10.
   lines
     25 11 
   boundary conditions
     0 1 1 0
   * One boundary here should match one boundary of 
   * the backGroundGrid, while another boundary 
   * should match a boundary on the inlet-bottom.
   * Set share flag to match corresponding share values
   share
     0 5 2 0
    pause
   exit
* 
  SmoothedPolygon
    mappingName
      inlet-bottom
    vertices
     3
     2. .15 
     2. .35 
     2.25 .35 
   lines
     25 11 
   n-dist
     fixed normal distance
      .175  .2
   sharpness
     10.
     10.
     10.
   t-stretch
     0. 10.
     1. 10.
     0. 10.
   boundary conditions
     0 1 1 0
   * One boundary here should match one boundary 
   * of the backGroundGrid, while another boundary 
   * should match a bounbdary on the inlet-bottom.
   * Set share flag to match corresponding share values
   share
     0 5 2 0
   exit
  * here is an outlet grid made in the poor man's way
  rectangle
    specify corners
      -.35 .3  .05 .7
    lines
      15 15
    mappingName
     outlet
    boundary conditions
      1 0 1 1
    pause
  exit
  * here is another inlet grid made in the poor man's way
  rectangle
    specify corners
       .85  .925  1.15 1.35
    lines
      11 11
    mappingName
     inlet2
    boundary conditions
      1 1 0 1
    pause
  exit
  * now look at the mappings
  view mappings
    backGroundGrid
    annulus
    inlet-top
    inlet-bottom
    outlet
    inlet2
    *
    * The grid is plotted with boundaries coloured
    *  by the boundary condition number. Here we 
    * should check that all interpolation boundaries 
    * are 0 (blue), all physical boundaries are positive 
    * and periodic boundaries are black
    * pause
    *
    * now we plot the boundaries by share value
    * The sides that correspond to the same boundary 
    * should be the same colour
    colour boundaries by share value 
    pause
    erase and exit
  exit
generate an overlapping grid
* put the nonconforming grid first to be a lower 
* priority than the back-ground
    outlet
    inlet2
    backGroundGrid
    annulus
    inlet-top
    inlet-bottom
  done
  change parameters
    prevent hole cutting
      backGroundGrid
        all
      outlet
        all
      inlet2
        all
    done
    ghost points
      all
      2 2 2 2 2 2
  exit
  compute overlap
  pause
  exit
*
save an overlapping grid
  twoInlet.hdf
  twoInlet
exit
