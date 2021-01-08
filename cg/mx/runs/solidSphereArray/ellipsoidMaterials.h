#
#  Specify material files for the array of solid ellipsoids 
#
$i=1; # start at domain 1 (outer domain is domain 0)
$domainName[$i]="ellipsoid$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
#
$domainName[$i]="ellipsoid$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
#
$domainName[$i]="ellipsoid$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
#
$domainName[$i]="ellipsoid$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1;
$domainName[$i]="ellipsoid$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1;
#
for( $j=0; $j<200; $j++ )\
{\
$domainName[$i]="ellipsoid$i"; $npv[$i]=1; $dmFile[$i]="gdmMaterial1.txt"; $i=$i+1; \
$domainName[$i]="ellipsoid$i"; $npv[$i]=2; $dmFile[$i]="gdmMaterial2.txt"; $i=$i+1; \
}

   
