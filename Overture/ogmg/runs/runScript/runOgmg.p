eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#
# -------------------------------------------------------------------
# runOgmg.p : perl program to run ogmgt and vary parameters 
# -------------------------------------------------------------------
# 

use Getopt::Long; use Getopt::Std;

$numberOfParameters = @ARGV;
if ($numberOfParameters eq 0)
{
  
  printf("\n");
  printf("================================================================================\n");
  printf("This perl script will run ogmgt and vary parameters\n");
  printf("  Usage: \n");
  printf("    runOgmg.p -cmd=<file.cmd> -g=<gridName> -numOmega=<i> -omegaStart=<>f -omegaEnd=<> -matlabFile=<s> ... \n");
  printf("              -sm=<i> -preSmooth=<i> -postSmooth=<i> -orderCoarse=<i> -options=\"... ogmg options ....\" \n");
  printf("  where \n");
  printf("    -cmd=file.cmd  : ogmgt command file \n");
  printf("==============================================================================\n\n");
  exit;
  
}

$OVERTURECHECKOUT=$ENV{OvertureCheckout};
$ogmgt = "$OVERTURECHECKOUT/ogmg/ogmgt";  # command for ogmgt


$cmdFile="tz.cmd";
$grid="square1024.order2.hdf"; 
$check="ogmg.check";
$ic=2; #random
$levels=2;
$sm=rb;
$numOmega=2;
$omegaStart=1.;
$omegaEnd  =1.2; 
$matlabFile="ogmgResults.m"; # name of matlab output file
$cycle="V"; $preSmooth=1; $postSmooth=1; 
$opav=0; # 0=noGalerkin, 1=Galerkin coarse grid
$orderCoarse=-1; # can be 2 or 4 for fourth-order, -1=default
$options=""; 
foreach $arg ( @ARGV )
{
  if( $arg =~ /-cmd=(.*)/ )
  {
    $cmdFile = $1;
    printf("Using command file [%s]\n",$cmdFile);
  }
  elsif( $arg =~ /-levels=(.*)/)
{
  $levels = $1;
  printf("Setting levels=$levels\n");
}
 elsif( $arg =~ /-ic=(.*)/)
{
  $ic = $1;
  printf("Setting ic=$ic\n");
}
  elsif(  $arg =~ /-numOmega=(.*)/ )
  {
    $numOmega = $1;
    printf("Setting numOmega=$numOmega\n");
  }
 elsif(  $arg =~ /-sm=(.*)/ )
  {
    $sm = $1;
    printf("Setting sm=$sm\n");
  }
  elsif(  $arg =~ /-preSmooth=(.*)/ )
  {
    $preSmooth = $1;
    printf("Setting preSmooth=$preSmooth\n");
  }
  elsif(  $arg =~ /-postSmooth=(.*)/ )
  {
    $postSmooth = $1;
    printf("Setting postSmooth=$postSmooth\n");
  }
  elsif(  $arg =~ /-opav=(.*)/ )
  {
    $opav = $1;
    printf("Setting opav=$opav\n");
  }
  elsif(  $arg =~ /-g=(.*)/ )
  {
    $grid = $1;
    printf("Setting grid=$grid\n");
  }
  elsif(  $arg =~ /-cycle=(.*)/ )
  {
    $cycle = $1;
    printf("Setting cycle=$cycle\n");
  }
  elsif(  $arg =~ /-orderCoarse=(.*)/ )
  {
    $orderCoarse = $1;
    printf("Setting orderCoarse=$orderCoarse\n");
  }
  elsif(  $arg =~ /-matlabFile=(.*)/ )
  {
    $matlabFile = $1;
    printf("Setting matlabFile=$matlabFile\n");
  }
  elsif(  $arg =~ /-options=(.*)/ )
  {
    $options = $1;
    printf("Setting options=$options\n");
  }
}

$date= localtime();   # here is the date
printf("date=[$date]\n");

for( $i=0; $i<$numOmega; $i++ )
{
  $omega = $omegaStart + ($omegaEnd-$omegaStart)*$i/($numOmega-1); 

  # --- Here is the ogmgt command and command line arguments -------
  $cmd = "$ogmgt -noplot $cmdFile -ic=$ic -g=$grid -levels=$levels -cycle=$cycle -sm=$sm -autoSmooth=0 -omega=$omega -maxIt=8 -opav=$opav -opavCoarseGrid=$opav -orderCoarse=$orderCoarse -nsm=\"$preSmooth $postSmooth\""; 

  printf("run [%s]\n",$cmd);
  $startTime=time();
  $returnValue = system("csh -f -c 'nohup $cmd >! junk.out'");
  $cpuTime=time()-$startTime;
  printf("...returnValue=%g, cpu=%g (s)\n",$returnValue,$cpuTime);

  # --------------------------------------------------------------------------------
  # --------------------- PROCESS CURRENT CHECK FILE -------------------------------
  # --------------------------------------------------------------------------------
  open(CHECKFILE,"$check") || die "cannot open file $check!" ;

  # read the check file and look for the convergence rate
  printf("--- checkfile:\n");
  $iStart=5; $iEnd=8; # average these cycles
  $j=-1; 
  $crAve=0; $ecrAve=0; 
  while( <CHECKFILE> )
  {
    $line = $_;
    chop($line);
    # printf("%s\n",$line);   
    @token =  split(' ',$line); 
    $t=$token[0]; $wu=$token[4]; $cr=$token[6]; $ecr=$token[7]; 
    # printf(" t=%f, CR=%f, WU=%f, ECR=%f\n",$t,$cr,$wu,$ecr);

    if( $j >=$iStart && $j <= $iEnd )
    { 
      # printf("average this:  t=%f, CR=%f, WU=%f, ECR=%f\n",$t,$cr,$wu,$ecr);
      $crAve = $crAve + $cr; $ecrAve = $ecrAve + $ecr;
    }
    $j++; 
  }
  close(CHECKFILE);

  $crAve = $crAve/($iEnd-$iStart+1);
  $ecrAve = $ecrAve/($iEnd-$iStart+1);
  printf("omega=%f, Average the last few iters: CR=%f, ECR=%f\n",$omega,$crAve,$ecrAve);

  $omegav[$i]=$omega;
  $crv[$i]=$crAve;
  $ecrv[$i]=$ecrAve;

}

# ---- Write a matlab file -----
open(MFILE,">$matlabFile") || die print "unable to open $matlabFile\n";
print MFILE "%\n% File created by runOgmg.p on $date \n%\n"; 
print MFILE "fontSize=16; lineWidth=2;\n"; 
# --- omega ----
print MFILE "w = [";
for( $ii=0; $ii<$numOmega; $ii++ )
{
  print MFILE "$omegav[$ii] ";
}
print MFILE "];\n";
# -- convergence rate ----
print MFILE "cr = [";
for( $ii=0; $ii<$numOmega; $ii++ )
{
  print MFILE "$crv[$ii] ";
}
print MFILE "];\n";
# -- effective convergence rate ----
print MFILE "ecr = [";
for( $ii=0; $ii<$numOmega; $ii++ )
{
  print MFILE "$ecrv[$ii] ";
}
print MFILE "];\n";
print MFILE "figure;\n";
print MFILE "plot(w,cr,'b-x','LineWidth',lineWidth);\n"; 
print MFILE "set(gca,'fontSize',fontSize);\n"; 
print MFILE "title('MG-CR: $levels-grid $cycle\[$preSmooth,$postSmooth\] sm=$sm opav=$opav $grid');\n";
print MFILE "xlabel('\\omega');\n";
print MFILE "grid on;\n";
$plotName = $matlabFile; 
$plotName =~ s/\.m$//; # remove trailing ".m"
print MFILE "print('-dpdf','$plotName.pdf');\n"; 
close( MFILE );
print "Matlab results written to file [$matlabFile].\n";


exit;
# ===============================================================================================================
