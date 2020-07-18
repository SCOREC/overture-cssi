eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#!/usr/bin/perl
#
# perl program to PERFORM CGMX SCATTERING RUNS FOR BA HOMOGENIZATION
# 

printf("\n");
printf("================================================================================\n");
printf("perl program to PERFORM CGMX SCATTERING RUNS FOR BA HOMOGENIZATION\n");
printf("==============================================================================\n\n");

$verbose=1;
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

$dm="none"; # dispersion model 

$tf=1; 
$polarization=1; # set to -1 to do both
$direction=1;    # set to -1 to do both
$kx=1;
$ic="gp"; # ic=gp or pw
$grid = "blockGrid"; # "blockGrid3d"; 
$rbc = "rbcNonLocal"; 



foreach $arg ( @ARGV )
{
  if( $arg =~ /-solver=(.*)/ ){ $solver=$1; }
  elsif( $arg =~ /-cmdFile=(.*)/ ){$cmdFile=$1; }
  elsif( $arg =~ /-resStart=(.*)/ ){$resStart=$1; }
  elsif( $arg =~ /-numResolutions=(.*)/ ){$numResolutions=$1; }
  elsif( $arg =~ /-verbose=(.*)/ ){ $verbose=$1; }
  elsif( $arg =~ /-case=(.*)/ ){ $case=$1;  }
  elsif( $arg =~ /-order=(.*)/ ){ $order=$1;  }
  elsif( $arg =~ /-matFile=(.*)/ ){ $matFile=$1; }
  elsif( $arg =~ /-matFile1=(.*)/ ){ $matFile1=$1; }
  elsif( $arg =~ /-matFile2=(.*)/ ){ $matFile2=$1; }
  elsif( $arg =~ /-matFile3=(.*)/ ){ $matFile3=$1; }
  elsif( $arg =~ /-matFile4=(.*)/ ){ $matFile4=$1; }
  elsif( $arg =~ /-direction=(.*)/ ){ $direction=$1; }
  elsif( $arg =~ /-polarization=(.*)/ ){ $polarization=$1; }
  elsif( $arg =~ /-tf=(.*)/   ){ $tf=$1;  }
  elsif( $arg =~ /-ic=(.*)/   ){ $ic=$1;  }
  elsif( $arg =~ /-grid=(.*)/ ){ $grid=$1; }
  elsif( $arg =~ /-rbc=(.*)/  ){$rbc=$1; }
  elsif( $arg =~ /-dm=(.*)/   ) {  $dm=$1;  }
  elsif( $arg =~ /-numSlabs=(.*)/   ) {  $numSlabs=$1;  }
  else
  {
    $testName=$arg;
  }
}

#
#  Scattering runs: single block
#
#  Variables:
#   $factor : grid resolution
#   $polarization : 1 or 2
#   $direction : 1 or 2 for forward of backward 
#   $kx, $ky, $kz : wave vector 



$program = "/home/henshw/cg.g/mx/bin/cgmx"; 


@planeWaveStartPoint = ( -2, 2); # starting point for plane wave, forward and backward
@kxv = ($kx,-$kx); 
$xProbe=1; 


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

for( $i=$resStart; $i <= $numResolutions; $i++ )
{  
  $factor = 2**($i); 

  # ----- Loop over directions: forward/backward 
  for( $dir=$dirStart; $dir<=$dirEnd; $dir++ )
  {
    # ----- Loop over polarizations
    for( $polar=$polarStart; $polar<=$polarEnd; $polar++ )
    {

      # Choose the polarization:
      #   Ev = [ax,ay,az]^T exp( ...  )
      if( $polar==1 ){ $ay=1; $az=0; }else{ $ay=0; $az=1; } 
      
      $x0 = $planeWaveStartPoint[$dir-1]; 
      $kxd = $kxv[$dir-1];
    
      if( $ic eq "pw" ){ $x0Cmd = "-xb=$x0"; }else{ $x0Cmd = "-x0=$x0"; } 
    
      $probeFileName = "$caseName" . "Polar$polar" . "Dir$dir" . "G$factor"; 
    
      if( $cmdFile eq "dielectricBodies" )
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

