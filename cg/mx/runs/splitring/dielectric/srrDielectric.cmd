$ds0=0.0125; 
$H=1.; $W=1.; $L=.4;
$D=.2; 
$nr=8;
$factor=1;
$ds=$ds0/$factor;
$xa=-2.; $xb=2;
# some default parameters: 
$numberOfVolumeSmooths=50;
$ya=-2.; $yb=2.;
$ng=2;
$interp="implicit for all grids";
$name="dielectricRes" ."$factor".".hdf";
sub resonator {\
  my ($x0,$y0,$num) = @_; \
  my $cmd = "";\
  $cmd .= "create mappings \n";\
  $cmd .= "  nurbs (curve) \n";\
  $cmd .= "  enter control points \n";\
  $cmd .= "  2 \n";\
  $cmd .= "  25  \n";\
  $cmd .= "  0.025 \n";\
  $cmd .= "  0.05 \n";\
  $cmd .= "  0.1 \n";\
  $cmd .= "  0.15 \n";\
  $cmd .= "  0.2 \n";\
  $cmd .= "  0.25 \n";\
  $cmd .= "  0.3 \n";\
  $cmd .= "  0.35 \n";\
  $cmd .= "  0.4 \n";\
  $cmd .= "  0.45 \n";\
  $cmd .= "  0.5 \n";\
  $cmd .= "  0.5 \n";\
  $cmd .= "  0.55 \n";\
  $cmd .= "  0.6 \n";\
  $cmd .= "  0.65 \n";\
  $cmd .= "  0.7 \n";\
  $cmd .= "  0.75 \n";\
  $cmd .= "  0.8 \n";\
  $cmd .= "  0.85 \n";\
  $cmd .= "  0.9  \n";\
  $cmd .= "  0.95 \n";\
  $cmd .= "  0.975 \n";\
  $cmd .= "  $x0          ,$y0             ,1.   \n";\
  $cmd .= "  $x0           ,@{[ $y0-$H/4    ]},1.   \n";\
  $cmd .= "  $x0           ,@{[ $y0-$H/2    ]},1.   \n";\
  $cmd .= "  @{[$x0+$W/2   ]},@{[$y0-$H/2     ]},1. \n";\
  $cmd .= "  @{[$x0+$W     ]},@{[$y0-$H/2     ]},1. \n";\
  $cmd .= "  @{[$x0+$W     ]},@{[$y0-($H-$L)/2]},1. \n";\
  $cmd .= "  @{[$x0+$W     ]},@{[$y0-$H/2+$L  ]},1. \n";\
  $cmd .= "  @{[$x0+$W-$D/2]},@{[$y0-$H/2+$L  ]},1. \n";\
  $cmd .= "  @{[$x0+$W-$D  ]},@{[$y0-$H/2+$L  ]},1. \n";\
  $cmd .= "  @{[$x0+$W-$D  ]},@{[$y0-$H/2+$D  ]},1. \n";\
  $cmd .= "  @{[$x0+$W/2   ]},@{[$y0-$H/2+$D  ]},1. \n";\
  $cmd .= "  @{[$x0+$D     ]},@{[$y0-$H/2+$D  ]},1. \n";\
  $cmd .= "  @{[$x0+$D     ]},@{[$y0        ]},1. \n";\
  $cmd .= "  @{[$x0+$D     ]},@{[ $y0+$H/2.-$D ]},1. \n";\
  $cmd .= "  @{[$x0+$W/2   ]},@{[ $y0+$H/2.-$D ]},1. \n";\
  $cmd .= "  @{[$x0+$W-$D  ]},@{[ $y0+$H/2.-$D ]},1. \n";\
  $cmd .= "  @{[$x0+$W-$D  ]},@{[ $y0+$H/2.-$L ]},1. \n";\
  $cmd .= "  @{[$x0+$W-$D/2]},@{[ $y0+$H/2.-$L ]},1. \n";\
  $cmd .= "  @{[$x0+$W     ]},@{[ $y0+$H/2.-$L ]},1. \n";\
  $cmd .= "  @{[$x0+$W     ]},@{[ $y0+($H-$L)/2]},1. \n";\
  $cmd .= "  @{[$x0+$W     ]},@{[ $y0+$H/2.    ]},1. \n";\
  $cmd .= "  @{[$x0+$W/2   ]},@{[ $y0+$H/2.    ]},1. \n";\
  $cmd .= "  $x0             ,@{[ $y0+$H/2.    ]},1. \n";\
  $cmd .= "  $x0             ,@{[ $y0+$H/4.    ]},1. \n";\
  $cmd .= "  $x0           ,$y0              ,1. \n";\
  $cmd .= "  \$nLines = int(1*(\$H+2*(\$W+\$L)-3*\$D)/\$ds); \n";\
  $cmd .= "  lines \n";\
  $cmd .= "    \$nLines  \n";\
  $cmd .= "  exit \n";\
  $cmd .= "  hyperbolic \n";\
  $cmd .= "    forward \n";\
  $cmd .= "    \$nDist=(\$nr-5)*\$ds; \n";\
  $cmd .= "    distance to march \$nDist \n";\
  $cmd .= "    \$nrm=\$nr-3;  \n";\
  $cmd .= "    lines to march \$nrm \n";\
  $cmd .= "    points on initial curve \$nLines \n";\
  $cmd .= "    uniform dissipation 0.05 \n";\
  $cmd .= "    volume smooths \$numberOfVolumeSmooths \n";\
  $cmd .= "    equidistribution 0 (in [0,1]) \n";\
  $cmd .= "    # \n";\
  $cmd .= "    spacing: geometric \n";\
  $cmd .= "    geometric stretch factor 1.05  \n";\
  $cmd .= "    # \n";\
  $cmd .= "    generate \n";\
  $cmd .= "    # open graphics \n";\
  $cmd .= "    boundary conditions \n";\
  $cmd .= "      -1 -1 101 0 0 0 \n";\
  $cmd .= "    share  \n";\
  $cmd .= "       0 0 101 0 0 0 \n";\
  $cmd .= "    name outer$num  \n";\
  $cmd .= "  exit \n";\
  $cmd .= "  hyperbolic \n";\
  $cmd .= "    backward \n";\
  $cmd .= "    \$nDist=(\$nr-5)*\$ds; \n";\
  $cmd .= "    distance to march \$nDist \n";\
  $cmd .= "    \$nrm=\$nr-3;  \n";\
  $cmd .= "    lines to march \$nrm \n";\
  $cmd .= "    points on initial curve \$nLines \n";\
  $cmd .= "    uniform dissipation 0.05 \n";\
  $cmd .= "    volume smooths \$numberOfVolumeSmooths \n";\
  $cmd .= "    equidistribution 0 (in [0,1]) \n";\
  $cmd .= "    # \n";\
  $cmd .= "    spacing: geometric \n";\
  $cmd .= "    geometric stretch factor 1.05  \n";\
  $cmd .= "    # \n";\
  $cmd .= "    generate \n";\
  $cmd .= "    # open graphics \n";\
  $cmd .= "    boundary conditions \n";\
  $cmd .= "      -1 -1 101 0 0 0 \n";\
  $cmd .= "    share  \n";\
  $cmd .= "       0 0 101 0 0 0 \n";\
  $cmd .= "    name inner$num  \n";\
  $cmd .= "  exit  \n";\
  $cmd .= "# ------- inner background grid ----- \n";\
  $cmd .= "#  \n";\
  $cmd .= "  rectangle \n";\
  $cmd .= "    \$xai = @{[$x0       ]}; \n";\
  $cmd .= "    \$xbi = @{[$x0+$W   ]}; \n";\
  $cmd .= "    \$yai = @{[$y0-$H/2.]}; \n";\
  $cmd .= "    \$ybi = @{[$y0+$H/2.]}; \n";\
  $cmd .= "    set corners  \n";\
  $cmd .= "      \$xai \$xbi \$yai \$ybi \n";\
  $cmd .= "    lines \n";\
  $cmd .= "    \$nx = int( \$W/\$ds +1.5 );  \n";\
  $cmd .= "    \$ny = int( \$H/\$ds +1.5 );  \n";\
  $cmd .= "      \$nx \$ny \n";\
  $cmd .= "    boundary conditions \n";\
  $cmd .= "      0 0 0 0  \n";\
  $cmd .= "    mappingName \n";\
  $cmd .= "      innerBackGround$num  \n";\
  $cmd .= "  exit  \n";\
  return $cmd;\
} 
 $cmd = resonator(0.,0.,0);
 $cmd 
 $cmd = resonator(0.,1.25,1);
 $cmd 
 $cmd = resonator(0.,-1.25,2);
 $cmd 
 $cmd = resonator(-1.25,0.,3);
 $cmd 
 $cmd = resonator(-1.25,1.25,4);
 $cmd 
 $cmd = resonator(-1.25,-1.25,5);
 $cmd 
 rectangle
    set corners: $xa,$xb,$ya,$yb
    $nx=int(($xb-$xa)/$ds+1.5);
    $ny=int(($yb-$ya)/$ds+1.5);
    lines
      $nx, $ny 
    boundary conditions
     1 2 -1 -1
    mappingName
      backGround
  exit
  exit
generate an overlapping grid
  backGround
  outer0 
  innerBackGround0
  inner0
  outer1 
  innerBackGround1
  inner1
  outer2 
  innerBackGround2
  inner2
  outer3 
  innerBackGround3
  inner3
  outer4 
  innerBackGround4
  inner4
  outer5 
  innerBackGround5
  inner5
  done 
  change parameters
    specify a domain
      innerDomain 
      innerBackGround0
      inner0
      innerBackGround1
      inner1
      innerBackGround2
      inner2
      innerBackGround3
      inner3
      innerBackGround4
      inner4
      innerBackGround5
      inner5
    done
    specify a domain
      outerDomain
      backGround
      outer0
      outer1
      outer2
      outer3
      outer4
      outer5
    done 
    # choose implicit or explicit interpolation
    interpolation type
      $interp
    # -- set the discretization width and interpolation width --
    #$cmd =" order of accuracy\n $orderOfAccuracy";
    #$cmd
    #
    ghost points
      all
      $ng $ng $ng $ng $ng $ng
  exit 
  compute overlap 
  exit
save a grid (compressed)
$name
SRR
exit

