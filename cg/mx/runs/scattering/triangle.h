#
#  ----- Build a triangle ------
#
  SmoothedPolygon
    vertices
      5
       .866025 0.
       .866025 .5
      0. 0 .
       .866025 -.5
       .866025 0.
    boundary conditions
      -1 -1 1 0
    n-dist
    fixed normal distance
      $nr=9+$order;
      $nDist=-($nr-3)*$ds;
      $nDist
    n-stretch
      1. 2.  0
    t-stretch
      0. 0.
      .15 15
      .15 15
      .15 15
      .15 15
    lines
      $nTheta = int( 3./$ds + 1.5 );
      $nTheta $nr
    mappingName
      triangleBase
    exit
