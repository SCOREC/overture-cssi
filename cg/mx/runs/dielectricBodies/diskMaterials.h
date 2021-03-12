#
#  Specify material files for the array of solid disks
#
if( $outerDomain eq "" ){ $outerDomain="outerDomain"; } # default name of outer-domain
$baseDomainName = "innerDomain";
#  
$i=0; # (outer domain is domain 0)
# Outer-domain:
$domainName[$i]="$outerDomain"; $npv[$i]=0; $dmFile[$i]="baIsoEps1.txt"; $i=$i+1;
# 
# non-dispersive materials: 
# $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps2.txt"; $i=$i+1;
# $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps3.txt"; $i=$i+1;
#
# dispersive materials
$domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
#
$domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
#
$domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
#
$domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
#
for( $j=0; $j<200; $j++ )\
{\
$domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1; \
$domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1; \
}
