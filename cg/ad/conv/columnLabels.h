#
#   Define the table column labels for regression tests by conv.p
#      These are usually common between different cases
#
# $numberOfComponents=1;
if( $numberOfComponents eq 1  )\
  { $title= "grid  \& N \&  \$ u \$ "; }
#
# BA-GDM : Save vector errors only 
if( $method eq "bamx" && $solveForAllFields eq 1 && $dm eq "gdm" && $saveVectorErrors eq 1 )\
  { $title= "grid  \& N \&  \$\\Ev\$    \&   \$\\Hv\$    \&  \$\\Pv\$"; $numberOfComponents=9;  \
  foreach $i (0,1,2,3,4,5) { $ignoreComponent[$i]=1; } \
  }
# 
$labelMethod=$method;
if( $method eq "NFDTD" ){ $labelMethod="CGFD$order"; }
if( $method eq "bamx" ){ $labelMethod="BAMX$order"; }
