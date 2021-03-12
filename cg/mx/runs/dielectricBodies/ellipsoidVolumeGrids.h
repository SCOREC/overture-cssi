# ---- Construct three volume grids for the ellipsoid ---
# 
#   This file should follow ellipsoidSurfaces.h 
# 
#    ---- Hyperbolic grid for the ellipsoid body ----
#
#  
  hyperbolic
    # Start curve:ellipsoid
    Start curve:ellipsoidSurfaceNurbs
    $directionToMarch
    distance to march $rDist
    $nrm = $nr-1; 
    lines to march $nrm  
    BC: bottom outward splay
    BC: top outward splay
    generate
    share
      0 0 0 0 $interfaceShare 0
    fourth order
    # Keep ghost points on the surface: 
    # bc: -1 -1 0 0 2 0 
    $numGhostBO = int( $order/2 + .5 ); 
    boundary offset 0 0 $numGhostBO $numGhostBO 0 0 (l r b t b f)
    # boundary offset 0 0 1 1 0 0 (l r b t b f)
    name $ellipsoidBodyName
    #
    # open graphics
  exit
  #
  #  ---north-pole hyperbolic grid ----
  # 
  hyperbolic
    Start curve:northPoleSurface
    $directionToMarch
    distance to march $rDist
    $nrm = $nr-1; 
    lines to march $nrm  
    BC: bottom outward splay
    BC: top outward splay
    generate
    share
      0 0 0 0 $interfaceShare 0
    fourth order
    # Keep ghost points on the surface: 
    boundary offset $numGhostBO $numGhostBO $numGhostBO $numGhostBO 0 0  (l r b t b f)
    # boundary offset 1 1 1 1 0 0 (l r b t b f)
    name $ellipsoidNorthPoleName
    # open graphics 
 exit
  #
  #  ---south-pole hyperbolic grid ----
  # 
  hyperbolic
    Start curve:southPoleSurface
    $directionToMarch
    distance to march $rDist
    $nrm = $nr-1; 
    lines to march $nrm  
    BC: bottom outward splay
    BC: top outward splay
    generate
    fourth order
    # Keep ghost points on the surface: 
    boundary offset $numGhostBO $numGhostBO $numGhostBO $numGhostBO 0 0  (l r b t b f)
    share
      0 0 0 0 $interfaceShare 0
    name $ellipsoidSouthPoleName
    #open graphics 
 exit
