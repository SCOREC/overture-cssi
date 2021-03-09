#
#   Define the table title for regression tests by conv.p
#
if( $numberOfComponents eq 5 && $numberOfDimensions eq 2 )\
  { $title= "grid  \& N \&  \$E_x\$ \&  \$E_y\$   \& \$H_z\$  \&  \$\\Ev\$ \& \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; }
#
if( $numberOfComponents eq 8  && $numberOfDimensions eq 2 )\
  { $title= "grid  \& N \&  \$E_x\$ \&  \$E_y\$   \& \$H_z\$   \&  \$\\Ev\$ \&  \$Ex_t\$ \&  \$Ey_t\$ \& \$Hz_t\$ \& \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; }
#
if( $numberOfComponents eq 5 && $numberOfDimensions eq 3 )\
  { $title= "grid  \& N \&  \$E_x\$ \&  \$E_y\$ \& \$E_z\$    \&  \$\\Ev\$ \& \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; }
#
if( $numberOfComponents eq 8  && $numberOfDimensions eq 3 )\
  { $title= "grid  \& N \&  \$E_x\$ \&  \$E_y\$ \& \$E_z\$ \&    \&   \$\\Ev\$ \$Ex_t\$ \&  \$Ey_t\$ \& \$Ez_t\$ \& \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; }
#
#      ----- DISPERSIVE MODEL TITLE -----
# 
#  -- Dispersive model 2D ---
if( $numberOfComponents eq 5 && $numberOfDimensions eq 2 && $dm ne "none" && $nm eq "none" )\
  { $numberOfComponents=6; $title= "grid  \&  N  \&    \$E_x\$   \&    \$E_y\$    \&    \$H_z\$    \&   \$\\Ev\$   \&    \$|\\Pv|\$   \&    \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; }
#  -- Dispersive model 3D ---
if( $numberOfComponents eq 5 && $numberOfDimensions eq 3 && $dm ne "none" && $nm eq "none"  )\
  { $numberOfComponents=6; $title= "grid  \&  N  \&    \$E_x\$   \&    \$E_y\$    \&    \$E_z\$    \&   \$\\Ev\$   \&    \$|\\Pv|\$   \&    \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; }
#
#      ---- NONLINEAR MODEL TITLE -------
#
#  -- Nonlinear model 2D ---
if( $numberOfDimensions eq 2 && $nm ne "none" )\
  { $numberOfComponents=7; $title= "grid  \&  N  \&    \$E_x\$   \&    \$E_y\$    \&    \$H_z\$    \&   \$\\Ev\$   \&    \$|\\Pv|\$  \&   \$|\\Qv|\$   \&    \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; }
#  -- Nonlinear model 3D ---
if( $numberOfDimensions eq 3 && $nm ne "none" )\
  { $numberOfComponents=8; $title= "grid  \&  N  \&    \$E_x\$   \&    \$E_y\$    \&    \$E_z\$    \&   \$\\Ev\$   \&    \$|\\Pv|\$  \&    \$|\\Qv|\$   \&    \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; }
#   
#
#    ----------- BA-MAXWELL TITLE ------
if( $method eq "bamx" && $solveForAllFields eq 1  &&  $dm eq "none" )\
  { $title= "grid  \& N \&  \$E_x\$ \&  \$E_y\$ \& \$E_z\$ \&  \$H_x\$ \&  \$H_y\$ \& \$H_z\$  \&   \$\\Ev\$ \&   \$\\Hv\$  \&  \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; $numberOfComponents=9; }
#
if( $method eq "bamx" && $solveForAllFields eq 1  &&  $dm eq "gdm" )\
  { $title= "grid  \& N \&  \$E_x\$   \&  \$E_y\$     \& \$E_z\$    \&  \$H_x\$    \&  \$H_y\$    \& \$H_z\$   \&   \$\\Ev\$   \&   \$\\Hv\$   \&  \$\\Pv\$   \&  \$\\grad\\cdot\\Ev\/\\grad\\Ev\$"; $numberOfComponents=10; }
#
# BA-MAXWELL -- only display vector errors 
if( $method eq "bamx" && $solveForAllFields eq 1 &&  $dm eq "none" && $saveVectorErrors eq 1 )\
  { $title= "grid  \& N \&    \$\\Ev\$ \&   \$\\Hv\$ "; $numberOfComponents=9;  \
foreach $i (0,1,2,3,4,5,8) { $ignoreComponent[$i]=1; }\
}
# BA-GDM save all errors 
if( $method eq "bamx" && $solveForAllFields eq 1 && $dm eq "gdm" )\
  { $title= "grid  \& N \&  \$E_x\$ \&  \$E_y\$     \& \$E_z\$     \&  \$H_x\$       \&  \$H_y\$      \& \$H_z\$      \&   \$\\Ev\$    \&   \$\\Hv\$    \&  \$\\Pv\$"; $numberOfComponents=9; \
 } \
  }
# BA-GDM : Save vector errors only 
if( $method eq "bamx" && $solveForAllFields eq 1 && $dm eq "gdm" && $saveVectorErrors eq 1 )\
  { $title= "grid  \& N \&  \$\\Ev\$    \&   \$\\Hv\$    \&  \$\\Pv\$"; $numberOfComponents=9;  \
  foreach $i (0,1,2,3,4,5) { $ignoreComponent[$i]=1; } \
  }
# 
$labelMethod=$method;
if( $method eq "NFDTD" ){ $labelMethod="CGFD$order"; }
if( $method eq "bamx" ){ $labelMethod="BAMX$order"; }
$caption= "\\caption\{Cgmx, $name, method=$labelMethod, $normName norm, method=$method, ts=$ts, order=\$$order\$, \$t=$tf\$, cfl=\$$cfl\$, diss=\$$diss\$, dissOrder=\$$dissOrder\$, kx=$kx, ky=$ky, kz=$kz, cons=$cons, nm=$nm, dm=$dm, alphaP=$alphaP, npv=$npv, modeGDM=$modeGDM a0=($a0[0],$a0[1],$a0[2]), a1=($a1[0],$a1[1],$a1[2]), b0=($b0[0],$b0[1],$b0[2]), b1=($b1[0],$b1[1],$b1[2]) tz=$tz degreex=$degreex degreet=$degreet fx=$fx fy=$fy fz=$fz -ft=$ft $date}\\label\{table:$name" . $method . "Order$order$normName}"; 
if( $method eq "bamx" ){\
$caption= "\\caption\{Cgmx, $name, method=$labelMethod, $normName norm, method=$method, ts=$ts, order=\$$order\$, \$t=$tf\$, cfl=\$$cfl\$, diss=\$$diss\$, kx=$kx, ky=$ky, kz=$kz, dm=$dm, -matFile=$matFile -drOption=$drOption -solveForAllFields=$solveForAllFields $date}\\label\{table:$name" . $method . "Order$order$normName}"; \
}