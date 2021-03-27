# 
# ------ Create Ellipsoid Surface Mappings for Hyperbolic Volume Grids ---- 
#
$crossSectionName         = "ellipsoidSurface$masterCount";
$ellipsoidSurfaceNurbs    = "ellipsoidSurfaceNurbs$masterCount";
$ellipsoidSingularSurface = "ellipsoidSingularSurface$masterCount";
$northPoleSurface         = "northPoleSurface$masterCount";
$southPoleSurface         = "southPoleSurface$masterCount";
#
  CrossSection 
    mappingName 
      $crossSectionName
      # ellipsoidSurface
    # 
    ellipse 
    a,b,c for ellipse 
    #
    $a $b $c 
    surface or volume (toggle)
    # 
    # remove singularities:
    # $baseAxialOffset=.02;
    $baseAxialOffset=.025; # March 19, 2021
    $deltaAxial=$baseAxialOffset*($factor**.5);  # Dec 12, 2020  -- try this to decrease small cells near poles
    $axialStart=$deltaAxial; $axialEnd=1.0-$deltaAxial;
    ##    $axialStart=.025; $axialEnd=.975;
    ##    $axialStart=.01; $axialEnd=.99;
    $axialFactor=$axialEnd-$axialStart; 
    start axial
      $axialStart
    end axial
      $axialEnd
    show parameters
    no singularity at start
    no singularity at end
    # cubic interpolate between cross sections:
    cubic
    lines
      # Fix me: 
      $arc1 = $pi*( 3*($a+$b) - sqrt( ($a+3*$b)*(3*$a+$b) ) );   # approx arc-length of an ellipse from Ramanujan
      $nTheta = intmg( $arc1/$ds +1.5);
      $ab = .5*($a+$b); # max($a,$b); 
      # $arc2 = $pi*( 3*($c+$ab) - sqrt( ($c+3*$ab)*(3*$c+$ab) ) ); # approx arc-length from Ramanujan
      # $nAxial = intmg( $arc2/$ds +1.5);
      # $nAxial = intmg( $axialFactor*(2*$c)/$ds +1.5);
      # try this: 
      $nAxial = intmg( $axialFactor*(2*$c+ .25*($a+$b))/$ds +1.5);
      $nTheta $nAxial
    exit
#
#  Turn cross-section surface into a NURBS and reparameterize by arc-length
#
  nurbs (surface)
    interpolate from mapping with options
      $crossSectionName
      # ellipsoidSurface
    parameterize by chord length
    done
    mappingName
      $ellipsoidSurfaceNurbs
      # ellipsoidSurfaceNurbs
    exit  
#
# ------ Ellipsoid with singular ends for cosntructing north and south pole patches ----
#
  CrossSection 
    mappingName 
      $ellipsoidSingularSurface
      # ellipsoidSingularSurface
    # 
    ellipse 
    a,b,c for ellipse 
      $a $b $c 
    surface or volume (toggle)
    # cubic interpolate between cross sections:
    cubic
    lines
      $nTheta $nAxial
  exit
  #
  #  ---- north-pole surface patch ---
  #
  reparameterize
    transform which mapping?
      $ellipsoidSingularSurface
    mappingName
      $northPoleSurface
    orthographic
      choose north or south pole
        +1
      specify sa,sb
        # Make sa bigger if $a is smaller so we get a bigger wrap around 
        # $sa=.95; $sb=.75;
        $abRatio = $a/$b-1; 
        if( $orthographicPatchParameter eq "" ){ $orthographicPatchParameter=.55; } #
        $sab0=$orthographicPatchParameter + ($order-2)*.1;   # value of sa, sb for a sphere 
        # $sab0=.75;
        $sa = $sab0*(1 + .1*$abRatio); $sb = $sab0*(1 + .1*$abRatio); 
        $sa $sb
      exit
    lines
      # $nTheta1 =intmg( .8*$sa*$pi*$a/$ds +1.5 );     # use $a
      # $nTheta2 =intmg( .8*$sb*$pi*$b/$ds +1.5 );     # use $b 
      $nTheta1 =intmg( .75*$sa*$pi*$a/$ds +1.5 );     # use $a
      $nTheta2 =intmg( .75*$sb*$pi*$b/$ds +1.5 );     # use $b 
      $nTheta1 $nTheta2 
    share
      0 0 0 0 1 0
    exit
  #
  #  ---- south-pole surface patch ---
  #
  reparameterize
    transform which mapping?
      $ellipsoidSingularSurface
    mappingName
      $southPoleSurface
    orthographic
      choose north or south pole
        -1
      specify sa,sb
        $sa $sb
      exit
    lines
      $nTheta1 $nTheta2 
    share
      0 0 0 0 1 0
  exit
