eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#
# ------------------------------------------------------------------------
# runCgmx.p : perl program to run cgmx in parallel and compute speedups
# ------------------------------------------------------------------------
# 

use Getopt::Long; use Getopt::Std;

sub parseLogFile # ($i)
# Parse the cgmx  log file
{
  # --------------------------------------------------------------------------------
  # --------------------- PROCESS CURRENT LOG FILE -------------------------------
  # --------------------------------------------------------------------------------

  local($i);

  $i = $_[0];

  open(LOGFILE,"$logFile") || die "cannot open file $logFile!" ;

  # read the log file and look for CPU times for the "advance" stage
  # printf("--- logFile ($logFile) :\n");

  while( <LOGFILE> )
  {
    $line = $_;
    chop($line);
    # printf("%s\n",$line);   

    if( $line =~/ Grid: (.*$)/ )
    { 
      @stuff = split('/',$1);  # separate grid name by "/" so we can remove the prefix 
      $gridName[$i] = $stuff[$#stuff];  

      printf("gridName[$i]=$gridName[$i]\n"); 
      # printf("gridName=$gridName[$i] (line=$line -> $stuff[$#stuff])\n"); 
    }

    if( $line =~/ total number of grid points =([^,]*)/ )
    { 
      $gridPoints[$i] = $1;  
      printf("gridPoints[$i]=$gridPoints[$i]\n"); 
    }

    
    if( $line =~/max errors:\[([^ ]*),\]/ )
    { 
      $maxErrors[$i]=$1;  
      # printf("max errors: $maxErrors[$i]\n"); 
    }

    if( $line =~/==== final time=([^,]*)/ )
    { 
      $tFinal[$i] = $1;  
      # printf("tFinal=$tFinal[$i]\n"); 
    }
    if( $line =~ /^advance/ )
    {
      # printf("%s\n",$line);   
      @token =  split(' ',$line); 
      # printf("Tokens: [$token[0]] [$token[1]] [$token[2]]\n");

      $cpuTotal[$i]  =$token[1];
      $cpuPerStep[$i]=$token[2];
      $cpuPerStepPerPt[$i] = $token[3];
      $cpuPercent[$i]=$token[4];
      printf("Case $i np=$np: advance: cpu=$cpuTotal[$i], cpu/step=$cpuPerStep[$i], percent=$cpuPercent[$i] \n");

      last;
    }

  }
  close(LOGFILE);
}



$numberOfParameters = @ARGV;
if ($numberOfParameters eq 0)
{
  
  printf("\n");
  printf("================================================================================\n");
  printf("This perl script will run cgmx in parallel and compute speedups\n");
  printf("  Usage: \n");
  printf("    runCgmx.p -cmd=<file.cmd> -g=<gridName> -tf=<f> -matlabFile=<s> -title=<s>  \n");
  printf("              -np=[num. processors] -order=[2|4] \n");
  printf("              -scheme=[FD,FDA,UPC|SOSUP] -options=\"... ogmg options ....\" \n");
  printf("  where \n");
  printf("    -cmd=file.cmd  : cgmx command file \n");
  printf("==============================================================================\n\n");
  exit;
  
}

$CGBUILDPREFIX=$ENV{CGBUILDPREFIX};
$cgmx = "$CGBUILDPREFIX/mx/bin/cgmx";  # command for cgmx


# TEST -- read a configuration  file
#- 
#- $return = do './myConfig.pl';
#- 
#- print "numProcs=$numProcs\n";
#- print "myCmd=[$myCmd]\n";
#- 
#- exit;


# The configuration file holds the commands for running in parallel
$configFile = "boxEigenConfig.pl"; 

$runCase="singleRun";   # run 1 case only 

$caseName="square256UPC";  # case name - used int matlab file name 

$computerInfo = "CG6, 48 procs = 2 X 12 core (X 2 threads), Intel Xeon CPU E5-2697 v2 @ 2.70GHz, 132 Gb"; 

$cmdFile="boxEigen.cmd";
$grid="square128.order2.hdf"; 
$check="mx.check";
$logFile="mx.log";
$matlabFile ="cgmxTiming";
$title=""; 
$option=""; 
$order=2;
$scheme="FD";
$gridType="rectangular"; # rectangular or curvilinear 

$factor=3;

$numRuns=2;    # number of parallel runs 

$nd=2; # number of space dimensions
$scaling = "strong";

$np=0;   #   0 : means run in serial 

$tf=.5; 

foreach $arg ( @ARGV )
{
  if( $arg =~ /-cmd=(.*)/ )
  {
    $cmdFile = $1;
    printf("Using command file [%s]\n",$cmdFile);
  }
  elsif( $arg =~ /-tf=(.*)/)
  {
    $tf = $1;
    printf("Setting tf=$tf\n");
  }
  elsif(  $arg =~ /-g=(.*)/ )
  {
    $grid = $1;
    printf("Setting grid=$grid\n");
  }
  elsif(  $arg =~ /-matlabFile=(.*)/ )
  {
    $matlabFile = $1;
    printf("Setting matlabFile=$matlabFile\n");
  }
  elsif( $arg =~ /-title=(.*)/)
  {
    $title = $1;
    printf("Setting title=$title\n");
  }
  elsif( $arg =~ /-gridType=(.*)/)
  {
    $gridType = $1;
    printf("Setting gridType=$gridType\n");
  }
  elsif(  $arg =~ /-runCase=(.*)/ )
  {
    $runCase = $1;
    printf("Setting runCase=$runCase\n");
  }
  elsif(  $arg =~ /-order=(.*)/ )
  {
    $order = $1;
    printf("Setting order=$order\n");
  }
  elsif(  $arg =~ /-caseName=(.*)/ )
  {
    $caseName = $1;
    printf("Setting caseName=$caseName\n");
  }
  elsif(  $arg =~ /-configFile=(.*)/ )
  {
    $configFile = $1;
    printf("Setting configFile=$configFile\n");
  }
  elsif(  $arg =~ /-scheme=(.*)/ )
  {
    $scheme = $1;
    printf("Setting scheme=$scheme\n");
  }
  elsif(  $arg =~ /-scaling=(.*)/ )
  {
    $scaling = $1;
    printf("Setting scaling=$scaling\n");
  }
  elsif(  $arg =~ /-numRuns=(.*)/ )
  {
    $numRuns = $1;
    printf("Setting numRuns=$numRuns\n");
  }
  elsif(  $arg =~ /-factor=(.*)/ )
  {
    $factor = $1;
    printf("Setting factor=$factor\n");
  }
  elsif(  $arg =~ /-np=(.*)/ )
  {
    $np = $1;
    printf("Setting np=$np\n");
  }
  elsif(  $arg =~ /-nd=(.*)/ )
  {
    $nd = $1;
    printf("Setting nd=$nd\n");
  }


}

$date= localtime();   # here is the date
printf("date=[$date]\n");


# ---- READ THE CONFIGURATION FILE ----
$myConfigFile = "./$configFile"; 
$return = do $myConfigFile;
# print "return=[$return]\n";
print "caseName=[$caseName]\n";
print "title=[$title]\n";
print "myCmd=[$myCmd]\n";

if( $myCmd eq "" )
{
  print "runParallel: ERROR: the configuration file is invalid. myCmd has not been set.\n";
  exit; 
}

$matlabFile = $caseName . "ParPerf"; 
$matlabFileName = $matlabFile . "$scheme"; 


# if( $runCase eq "singleRun" )
if( $runCase eq "performance" )
{
  # ================== Run in parallel and save CPU times ===============

#-  if( $scheme eq "FD" ){
#-    $option = "-diss=0.0"; 
#-  }
#-  elsif( $scheme eq "FDA" ){
#-    $option = "-diss=.5"; 
#-  }
#-  elsif( $scheme eq "UPC" ){
#-    $option = "-useSosupDissipation=1";
#-  }
#-  elsif( $scheme eq "SOSUP" ){
#-    $option = "-method=sosup";
#-  }
#-  else{ printf("Unknown scheme=[$scheme]\n"); exit; }
#-
#-  #  $matlabFile = "$grid" . "ParPerf"; 


  for( $irun=0; $irun<$numRuns; $irun++ )
  {
    
    $np = $numProcs[$irun];
    
    # ---- READ THE CONFIGURATION FILE AGAIN (file may use $irun and $np) ----
    $myConfigFile = "./$configFile"; 
    $return = do $myConfigFile;    # $return = do './myConfig.pl';

    # print "configFile=[$configFile]\n";
    # $return = do './boxEigenConfig.pl';
    # print "return=[$return]\n";
    # print "myCmd=[$myCmd]\n";
    # $cmd = $myCmd; 

    $cmd = $myCmd; 

    # exit;
    
    if( $np <= 0 )
    {
      print "ERROR: irun=$irun: np=$np is INVALID -- check yuour config file for setting the numProcs array\n";
      exit;
    }

    $cmd = "mpirun -np $np " . $cmd;

    printf("\n+++++++++++++++++++++++++++ scheme=$scheme, order=$order, irun=$irun scaling=$scaling +++++++++++++++++++++++++++++++++++++++++++\n");
    printf("run [%s]\n",$cmd);
    $startTime=time();
    $returnValue = system("csh -f -c 'nohup $cmd'");
    $cpuTime=time()-$startTime;
    printf("...returnValue=%g, cpu=%g (s)\n",$returnValue,$cpuTime);
    if( $returnValue != 0 )
    {
      print "runParallel.p: ERROR return from running code. Check the command\n";
      exit 1; 
    }

    # --------------------------------------------------------------------------------
    # --------------------- PROCESS CURRENT LOG FILE -------------------------------
    # --------------------------------------------------------------------------------
    parseLogFile($irun);
 
    # $line = $maxErrors[$res];
    # @errv =  split(',',$line);   #extract errors in Ex, Ey, [Ez or Hz]

    # $maxErr =0; 
    # for( $j=0; $j<3; $j++){ $err=$errv[$j]; if( $err > $maxErr ){ $maxErr=$err; } }
    # printf(" irun=$irun: maxErr=$maxErr\n");
    # $err[$res]= $maxErr; 
  }

  # -----------------------------
  # ---- Write a matlab file -----
  # -----------------------------
  open(MFILE,">$matlabFileName.m") || die print "unable to open $matlabFile.m\n";
  print MFILE "%\n% File created by runCgmx.p on $date \n%\n"; 
  print MFILE "% grid=$gridName[0]\n";
  print MFILE "% grid points=$gridPoints[0]\n";
  print MFILE "% computerInfo: $computerInfo\n";
  print MFILE "% \n";

  # print MFILE "clf; clear;\n"; 
  # print MFILE "fontSize=18; lineWidth=2;\n"; 
  
  print MFILE "caseName ='$caseName';\n"; 
  print MFILE "scaling ='$scaling';\n"; 
  print MFILE "myTitle ='$title';\n"; 
  print MFILE "tFinal=$tFinal[0];\n";
  
  print MFILE "numProcs = [";  for( $irun=0; $irun<$numRuns; $irun++ ){  print MFILE "$numProcs[$irun] "; } print MFILE "];\n";

  ## print MFILE "h = [";  for( $res=0; $res<=$numRes; $res++ ){  print MFILE "$h[$res] "; } print MFILE "];\n";
  ## print MFILE "err = [";  for( $res=0; $res<=$numRes; $res++ ){  print MFILE "$err[$res] "; } print MFILE "];\n";

  print MFILE "cpuTotal = [";  for( $irun=0; $irun<$numRuns; $irun++ ){ print MFILE "$cpuTotal[$irun] "; } print MFILE "];\n";
  
  
  print MFILE "if( plotOption == 1 )\n";  # plot results if plotOption==0 

  # log-log plot
  print MFILE "  plot( numProcs,cpuTotal,'b-x','linewidth',lineWidth);\n";
  
  # print MFILE "  set(gca,'fontSize',fontSize);\n"; 
  if( $title eq "" ){ $plotTitle="CPU versus np, $scheme$order, $grid"; }else{ $plotTitle="$title"; }
  print MFILE "  title('$plotTitle');\n";
  print MFILE "  grid on;\n";
  
  $plotName = $matlabFileName; 
  $plotName =~ s/\.m$//; # remove trailing ".m"
  print MFILE "  %% print('-depsc','$plotName.eps');\n"; 
  
  print MFILE "end % end if plotOption\n"; 

  close( MFILE );
  print "Matlab results written to file [$matlabFileName.m].\n";


  # -----------------------------
  # ---- Write a LaTeX file -----
  # -----------------------------
  open(MFILE,">$matlabFileName.tex") || die print "unable to open $matlabFile.tex\n";

  print MFILE "%\n% File created by runCgmx.p on $date \n%\n";
  print MFILE "% caseName = $caseName\n";
  print MFILE "% title = $title\n";
  print MFILE "% \n";
  # print MFILE "{ % START table scope \n";

  print MFILE "\\begin{table}[hbt]\n";
  print MFILE "\\begin{center}\\tableFont\n";
  print MFILE "\\begin{tabular}{|c|c|c|} \\hline\n"; 
  print MFILE "  NP  & sec/step   &  speedup  \\\\   \\hline\\hline \n"; 
  for( $irun=0; $irun<$numRuns; $irun++ )
  {  
    $speedup = $cpuPerStep[0]/$cpuPerStep[$irun];
    # print MFILE "$numProcs[$irun]      &  $cpuPerStep[$irun]   & $speedup  \\\\\n";
    printf MFILE "%4d  & %8.1e   & %7.1f  \\\\\n",$numProcs[$irun],$cpuPerStep[$irun],$speedup;
  }
  print MFILE " \\hline \n"; 
  print MFILE "\\end{tabular}\n";   
  print MFILE "\\end{center}\n"; 		
  print MFILE "\\caption{$caseName, $gridName[0], $gridPoints[0] grid points, Computer: $computerInfo, $date.}\n"; 
  print MFILE "\\label{tab:$caseName} \n"; 
  print MFILE "\\end{table}\n"; 
  #   print MFILE "} % END table scope \n";

  close( MFILE );
  print "LaTeX results written to file [$matlabFileName.tex].\n";

}


#- if( $runCase eq "compare" )
#- {
#-   # ============ COMPARE SCHEMES =================
#-   $numCases=4;
#-   for( $i=0; $i<$numCases; $i++ )
#-   {
#-   
#-     if( $i == 0 ){
#-       $option = "-diss=0."; $scheme="FD"; 
#-     }
#-     elsif( $i==1 ){
#-       $option = "-diss=.5";  $scheme="FDA";
#-     }
#-     elsif( $i==2 ){
#-       $option = "-useSosupDissipation=1"; $scheme="UPC";
#-     }
#-     elsif( $i==3 ){
#-       $option = "-method=sosup"; $scheme="SOSUP";
#-     }
#-   
#-     $gridName = "$grid.order$order"; 
#-     if( $order eq 4 && ( $scheme eq "UPC" || $scheme eq "SOSUP" ) ){ $gridName .= ".ng3"; } # append ".ng3" 
#- 
#-     # --- Here is the cgmx command and command line arguments -------
#-     $cmd = "$cgmx -noplot $cmdFile -g=$gridName -mx=3 -my=2 -m=3 -n=2 -tf=$tf -tp=$tf $option -go=go >! $matlabFile$i.out"; 
#-   
#-     printf("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#-     printf("run [%s]\n",$cmd);
#-     $startTime=time();
#-     $returnValue = system("csh -f -c 'nohup $cmd'");
#-     $cpuTime=time()-$startTime;
#-     printf("...returnValue=%g, cpu=%g (s)\n",$returnValue,$cpuTime);
#-   
#-     # --------------------------------------------------------------------------------
#-     # --------------------- PROCESS CURRENT LOG FILE -------------------------------
#-     # --------------------------------------------------------------------------------
#-     parseLogFile($i);
#-   
#-   }
#-   
#-   # ---- Write a matlab file -----
#-   open(MFILE,">$matlabFile.m") || die print "unable to open $matlabFile.m\n";
#-   print MFILE "%\n% File created by runCgmx.p on $date \n%\n"; 
#-   # print MFILE "clf; clear;\n"; 
#-   # print MFILE "fontSize=18; lineWidth=2;\n"; 
#-   
#-   print MFILE "tFinal=$tFinal[0];\n";
#-   
#-   for( $ii=0; $ii<$numCases; $ii++ ){ $jj=$ii+1; print MFILE "maxErrors$jj=[$maxErrors[$ii]];\n"; }
#-   
#-   print MFILE "cpuTotal = [";
#-   for( $ii=0; $ii<$numCases; $ii++ ){ print MFILE "$cpuTotal[$ii] "; } print MFILE "];\n";
#-   
#-   print MFILE "cpuPerStep = [";
#-   for( $ii=0; $ii<$numCases; $ii++ ){ print MFILE "$cpuPerStep[$ii] ";} print MFILE "];\n";
#-   
#-   
#-   print MFILE "relCPU = [";
#-   for( $ii=0; $ii<$numCases; $ii++ ){ $relCPU[$ii]=$cpuTotal[$ii]/$cpuTotal[0]; print MFILE "$relCPU[$ii] ";} print MFILE "];\n";
#-   
#-   # Bar chart
#-   print MFILE "bar(diag(relCPU),'stacked');\n";
#-   
#-   print MFILE "set(gca,'fontSize',fontSize);\n"; 
#-   if( $title eq "" ){ $plotTitle="Cgmx: relative CPU times"; }else{ $plotTitle="Rel-CPU, $title"; }
#-   print MFILE "title('$plotTitle');\n";
#-   
#-   
#-   print MFILE "legend('FD','FDA','UPC','SOSUP','Location','northwest');\n";
#-   print MFILE "grid on;\n";
#-   
#-   
#-   #-print MFILE "xlabel('\\omega');\n";
#-   #-print MFILE "grid on;\n";
#-   
#-   $plotName = $matlabFile; 
#-   $plotName =~ s/\.m$//; # remove trailing ".m"
#-   print MFILE "print('-depsc','$plotName.eps');\n"; 
#-   
#-   close( MFILE );
#-   print "Matlab results written to file [$matlabFile.m].\n";
#- 
#- }
#- 
#- exit;
#- # ===============================================================================================================
#- 
