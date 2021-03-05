eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
#
# perl program to PERFORM CGMX SCATTERING RUNS FOR BA HOMOGENIZATION
# 
use Getopt::Long; use Getopt::Std;

printf("\n");
printf("=================================================================================\n");
printf("runScat.p : perl program to PERFORM CGMX SCATTERING RUNS FOR BA HOMOGENIZATION\n");
printf("=================================================================================\n\n");

$verbose=1;
$gridFactor=-1;         # if set, run one resolution with this grid factor 
$resStart=1;        # start at this resolution
$numResolutions=2;
$order=4;   # order od accuracy 

$cmdFile = "dielectricBodies"; # "baScat" 

$numSlabs=1;
$case = "Block";           # single block case 
$matFile = "matIsoEps2";   # material file (without .txt)
$matFile1 = "baMatIsoEps2";   # material file (without .txt)
$matFile2 = "baMatIsoEps3";   # material file (without .txt)
$matFile3 = "baMatIsoEps4";   # material file (without .txt)
$matFile4 = "baMatIsoEps8";   # material file (without .txt)

# number of bodies to set materials for:
$numBodies=1;
# ------ material properities -------
#  File containing info on dispersive and non-dispersive materials: 
if( $materialFile  eq "" ){ $materialFile="diskMaterials.h"; }
# optionally supply material files on the command line to over-ride those in $materialFile : 
@matFileArray = ();  

$dm="none"; # dispersion model 

$tf=1; 
$polarization=1; # set to -1 to do both
$direction=1;    # set to -1 to do both
$kx=1;
$ic="gp"; # ic=gp or pw
$grid = "blockGrid"; # "blockGrid3d"; 
$rbc = "rbcNonLocal"; 

$xProbe=1;
$xGaussianPlaneWaveStartPoint = 2;   # becomes -2 for right-going and +2 for left-going 

$cgmxCmd = "";   # supply full cgmx command using TFINAL, FACTOR, AX, AY, AZ and XO, XPROBE, PROBEFILENAME 
# -- old way ---
foreach $arg ( @ARGV )
{
  if( $arg =~ /-solver=(.*)/ ){ $solver=$1; }
#  elsif( $arg =~ /-cgmxCmd=(.*)/ ){$cgmxCmd=$1; }
#  elsif( $arg =~ /-resStart=(.*)/ ){$resStart=$1; }
#  elsif( $arg =~ /-numResolutions=(.*)/ ){$numResolutions=$1; }
#  elsif( $arg =~ /-verbose=(.*)/ ){ $verbose=$1; }
#  elsif( $arg =~ /-case=(.*)/ ){ $case=$1;  }
#  elsif( $arg =~ /-order=(.*)/ ){ $order=$1;  }
#  elsif( $arg =~ /-matFile=(.*)/ ){ $matFile=$1; }
#  elsif( $arg =~ /-matFile1=(.*)/ ){ $matFile1=$1; }
#  elsif( $arg =~ /-matFile2=(.*)/ ){ $matFile2=$1; }
#  elsif( $arg =~ /-matFile3=(.*)/ ){ $matFile3=$1; }
#  elsif( $arg =~ /-matFile4=(.*)/ ){ $matFile4=$1; }
#  elsif( $arg =~ /-direction=(.*)/ ){ $direction=$1; }
#  elsif( $arg =~ /-polarization=(.*)/ ){ $polarization=$1; }
#  elsif( $arg =~ /-tf=(.*)/   ){ $tf=$1;  }
#  elsif( $arg =~ /-ic=(.*)/   ){ $ic=$1;  }
#  elsif( $arg =~ /-grid=(.*)/ ){ $grid=$1; }
#  elsif( $arg =~ /-rbc=(.*)/  ){$rbc=$1; }
#  elsif( $arg =~ /-dm=(.*)/   ) {  $dm=$1;  }
#  elsif( $arg =~ /-xProbe=(.*)/    ) {  $xProbe=$1;  }
#  elsif( $arg =~ /-numSlabs=(.*)/  ) {  $numSlabs=$1;  }
#  elsif( $arg =~ /-xGaussianPlaneWaveStartPoint=(.*)/    ) {  $xGaussianPlaneWaveStartPoint=$1;  }
  else
  {
    $testName=$arg;
  }
}

# *new* Read options using get GetOptions:
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "solver=s"=>\$solver,
            "cmdFile=s"=>\$cmdFile,
            "cgmxCmd=s"=>\$cgmxCmd,
            "resStart=s"=>\$resStart,
            "numResolutions=i"=>\$numResolutions,
            "gridFactor=i"=>\$gridFactor,
            "g=s"=>\$grid,"grid=s"=>\$grid,
            "ic=s"=>\$ic,
            "dm=s"=>\$dm,
            "rbc=s"=>\$rbc,
            "case=s"=>\$case,
            "order=i"=>\$order,
            "verbose=i"=>\$verbose,
            "tf=f"=>\$tf,
            "direction=i"=>\$direction,
            "polarization=i"=>\$polarization,
            "xProbe=f"=>\$xProbe,
            "xGaussianPlaneWaveStartPoint=f"=>\$xGaussianPlaneWaveStartPoint,
            "numSlabs=i"=>\$numSlabs,
            "matFile=s"=>\$matFile,
            "matFile1=s"=>\$matFile1,
            "matFile2=s"=>\$matFile2,
            "matFile3=s"=>\$matFile3,
            "matFile4=s"=>\$matFile4,
            "materialFile=s"=>\$materialFile,
            "numBodies=i"=>\$numBodies,
            "matFileArray=s{1,}"=>\@matFileArray );
# -------------------------------------------------------------------------------------------------


printf("runScat.p: direction=$direction, polarization=$polarization, numResolutions=$numResolutions, gridFactor=$gridFactor\n");
# printf("gridFactor=$gridFactor\n");
printf("cgmxCmd=[$cgmxCmd]\n");
# exit;

#
#  Scattering runs: single block
#
#  Variables:
#   $factor : grid resolution
#   $polarization : 1 or 2
#   $direction : 1 or 2 for forward of backward 
#   $kx, $ky, $kz : wave vector 



$program = "/home/henshw/cg.g/mx/bin/cgmx"; 

 # starting point for plane wave, forward and backward: 
@planeWaveStartPoint = ( -$xGaussianPlaneWaveStartPoint, $xGaussianPlaneWaveStartPoint);
@kxv = ($kx,-$kx); 


if( $ic eq "pw" ){ $case .= "PW"; }
$MatFile = $matFile; 
$MatFile =~ s/^([a-z])/\u$1/; # change first char to upper case 

if( $numSLabs eq "1" ) 
{
  $caseName = $case . $MatFile;  # **OLD WAY**
}
else
{
  $caseName = $case; # new way
}
if( $direction>0 ){ $dirStart=$direction; $dirEnd=$direction; }else{ $dirStart=1; $dirEnd=2; }
if( $polarization>0 ){ $polarStart=$polarization; $polarEnd=$polarization; }else{ $polarStart=1; $polarEnd=2; }

if( $order eq "4" )
{
   $gridOrder=".order4.ng3"; 
}
elsif( $order eq "2" )
{
   $gridOrder=".order2"; 
}
else
{
  printf("runScat.p:ERROR: unexpected order=$order\n");
  exit 1;
}  

# If $gridFactor is set then run one resolution with is factor 
if( $gridFactor>0 ){ $resStart=1; $numResolutions=1; }

# -------------------- LOOP OVER GRID RESOLUTIONS -----------------------
for( $i=$resStart; $i <= $numResolutions; $i++ )
{  
  if( $gridFactor>0 ){ $factor=$gridFactor; }else{ $factor = 2**($i); }

  # ----- Loop over directions: forward/backward 
  for( $dir=$dirStart; $dir<=$dirEnd; $dir++ )
  {
    # ----- Loop over polarizations
    for( $polar=$polarStart; $polar<=$polarEnd; $polar++ )
    {

      # Choose the polarization:
      #   Ev = [ax,ay,az]^T exp( ...  )
      if( $polar==1 ){ $ax=0; $ay=1; $az=0; }else{ $ax=0; $ay=0; $az=1; } 
      
      $x0 = $planeWaveStartPoint[$dir-1]; 
      $kxd = $kxv[$dir-1];
    
      if( $ic eq "pw" ){ $x0Cmd = "-xb=$x0"; }else{ $x0Cmd = "-x0=$x0"; } 
    
      $probeFileName = "$caseName" . "Polar$polar" . "Dir$dir" . "G$factor"; 
    
      if( $cgmxCmd ne "" )
      {
        # **new way
        $cmd = $cgmxCmd;
        # TFINAL, FACTOR, AX, AY, AZ and XO, XPROBE, PROBEFILENAME are replaced         
        $cmd =~ s/FACTOR/$factor/g;
        $cmd =~ s/TFINAL/$tf/g;
        $cmd =~ s/KX/$kxd/g;
        $cmd =~ s/AX/$ax/g;
        $cmd =~ s/AY/$ay/g;
        $cmd =~ s/AZ/$az/g;
        $cmd =~ s/X0/$x0Cmd/g;
        $cmd =~ s/XPROBE/$xProbe/g;
        $cmd =~ s/PROBEFILENAME/$probeFileName/;
      }
      elsif( $cmdFile eq "dielectricBodies" )
      {
        $cmd = "-noplot dielectricBodies -g=$grid" . "e$factor$gridOrder.hdf -tf=$tf -tp=1 -useSosupDissipation=1 -ic=$ic -kx=$kxd -ay=$ay -az=$az $x0Cmd -rbc=$rbc -background=OuterSquare -xLeftProbe=-$xProbe -xRightProbe=$xProbe -dm=$dm -probeFrequency=1 -probeFileName=$probeFileName -matFile=$matFile.txt -go=go";
        # $cmd = "-noplot dielectricBodies -g=$grid" . "e$factor.order4.ng3.hdf -tf=$tf -tp=1 -useSosupDissipation=1 -ic=$ic -kx=$kxd -ay=$ay -az=$az $x0Cmd -rbc=$rbc -background=OuterSquare -xLeftProbe=-$xProbe -xRightProbe=$xProbe -dm=$dm -probeFrequency=1 -probeFileName=$probeFileName -matFile=$matFile.txt -go=go";
      }
      elsif( $cmdFile eq "baScat" )
      {
        # bamx: 
        $cmd = "-noplot baScat.cmd -g=$grid$factor" . "np.order4.ng3.hdf -method=bamx -dm=$dm -numMatRegions=2 -regionFile=matBoxRegion.h -matFile=baIsoEps1.txt -matFile2=$matFile.txt -ts=rk4 -ic=gp -kx=$kxd -ay=$ay -az=$az $x0Cmd -beta=50 -solveForAllFields=1 -tp=1 -tf=$tf -diss=1 -xLeftProbe=-$xProbe -xRightProbe=$xProbe -probeFileName=$probeFileName -rbc=absorbing -useSuperGrid=1 -superGridWidth=.5 -go=go"; 
    
	if( $numSlabs eq "2" )
        {
           $cmd = "-noplot baScat.cmd -g=$grid$factor" . "np.order4.ng3.hdf -method=bamx -dm=$dm -numMatRegions=3 -regionFile=matTwoBoxRegion.h -matFile=baIsoEps1.txt -matFile2=$matFile2.txt -matFile3=$matFile3.txt -ts=rk4 -ic=gp -kx=$kxd -ay=$ay -az=$az $x0Cmd -beta=50 -solveForAllFields=1 -tp=1 -tf=$tf -diss=1 -xLeftProbe=-$xProbe -xRightProbe=$xProbe -probeFileName=$probeFileName -rbc=absorbing -useSuperGrid=1 -superGridWidth=.5 -go=go";          
        }
        # Order=2 - super-grid
        # $cmd = "-noplot baScat.cmd -g=$grid$factor" . "np.order2.hdf -method=bamx -numMatRegions=2 -regionFile=matBoxRegion.h -matFile=baIsoEps1.txt -matFile2=$matFile.txt -ts=rk4 -ic=gp -x0=-2 -beta=50 -solveForAllFields=1 -tp=1 -tf=$tf -diss=1 -xLeftProbe=-$xProbe -xRightProbe=$xProbe -probeFileName=$probeFileName -rbc=absorbing -useSuperGrid=1 -superGridWidth=.5 -go=go"; 
    
        # $cmd = "-noplot baScat.cmd -g=$grid$factor" . "np.order2.hdf -method=bamx -numMatRegions=2 -regionFile=matBoxRegion.h -matFile=baIsoEps1.txt -matFile2=$matFile.txt -ts=rk4 -ic=gp -x0=-2 -beta=50 -solveForAllFields=1 -tp=1 -tf=$tf -diss=1 -xLeftProbe=-$xProbe -xRightProbe=$xProbe -probeFileName=$probeFileName -go=go"; 
    
    
      }
      elsif( $cmdFile eq "slabs" )
      {
        $cmd = "-noplot slabs -g=$grid" . "e$factor.order4p.hdf -numSlabs=$numSlabs -tf=$tf -tp=1 -useSosupDissipation=1 -ic=$ic -kx=$kxd -ay=$ay -az=$az $x0Cmd -rbc=$rbc -xLeftProbe=-$xProbe -xRightProbe=$xProbe -dm=$dm -matFile1=$matFile1.txt -matFile2=$matFile2.txt -matFile3=$matFile3.txt -matFile4=$matFile4.txt -probeFrequency=1 -probeFileName=$probeFileName -matFile=$matFile.txt -go=go";
      }
      elsif( $cmdFile eq "mxScript" )
      {
        # --- use new master script ---
        ## if( $#matFileArray > 5 ){ printf("runScat.p: ERROR: fix me for more entries in matFileArray\n"); exit; }
        $matFileList=""; 
        for( $mf=0; $mf <= $#matFileArray; $mf++ )
        {
          $matFileList .= $matFileArray[$mf] . " "; # list of mat files
        }
        $cmd = "-noplot mxScript -g=$grid" . "e$factor$gridOrder.hdf -tf=$tf -tp=1 -useSosupDissipation=1 -ic=$ic -kx=$kxd -ay=$ay -az=$az $x0Cmd -rbc=$rbc -background=backGround -xLeftProbe=-$xProbe -xRightProbe=$xProbe -dm=$dm -numBodies=$numBodies -matFileArray=$matFileList -materialFile=$materialFile -probeFrequency=1 -probeFileName=$probeFileName -go=go";

      }
      else
      {
        printf("runscat:ERROR: Unknown cmd=[$cmd]\n");
        exit; 
      }    
    
      $runTimeOutput="runScat$probeFileName"; # pipe output from cgmx to this file 
    
      if( $verbose ){ printf("running:  $program -noplot $cmd >! $runTimeOutput.out\n"); }
      $returnValue = system("csh -f -c '$program -noplot $cmd >! $runTimeOutput.out'");

    } # end for polar
  } # end for dir
}

exit;

