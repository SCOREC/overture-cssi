#
#  Specify material files for the R P I grids 
#
#
$outerDomain    = "fluidDomain"; 
$baseDomainName = "solidDomain";
#  
$i=0; # (outer domain is domain 0)
# Outer-domain:
$domainName[$i]="$outerDomain"; $npv[$i]=0; $dmFile[$i]="baIsoEps1.txt"; $i=$i+1;
# 
for( $j=0; $j<3; $j++ ){\
  $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="mlaMat4levels.txt"; $i=$i+1; \
}
