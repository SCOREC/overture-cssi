eval 'exec perl -S $0 ${1+"$@"}'
if 0;
#
# -------------------------------------------------------------------
# runCgmx.p : perl program to run cgmx and compare CPU times 
# -------------------------------------------------------------------
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
  printf("--- logFile:\n");

  while( <LOGFILE> )
  {
    $line = $_;
    chop($line);
    # printf("%s\n",$line);   

    if( $line =~/max errors:\[([^ ]*),\]/ ){ $maxErrors[$i]=$1;  printf("max errors: $maxErrors[$i]\n"); }

    if( $line =~/==== final time=([^,]*)/ ){ $tFinal[$i] = $1;  printf("tFinal=$tFinal[$i]\n"); }
    if( $line =~ /^advance/ )
    {
      # printf("%s\n",$line);   
      @token =  split(' ',$line); 
      # printf("Tokens: [$token[0]] [$token[1]] [$token[2]]\n");

      $cpuTotal[$i]  =$token[1];
      $cpuPerStep[$i]=$token[2];
      $cpuPerStepPerPt[$i] = $token[3];
      $cpuPercent[$i]=$token[4];
      printf("Case $i: advance: cpu=$cpuTotal[$i], cpu/step=$cpuPerStep[$i], percent=$cpuPercent[$i] \n");

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
  printf("This perl script will run cgmx and compare timings\n");
  printf("  Usage: \n");
  printf("    runCgmx.p -cmd=<file.cmd> -g=<gridName> -tf=<f> -matlabFile=<s> -title=<s>  \n");
  printf("              -runCase=[compare|performance] -order=[2|4] \n");
  printf("              -scheme=[FD,FDA,UPC|SOSUP] -options=\"... ogmg options ....\" \n");
  printf("  where \n");
  printf("    -cmd=file.cmd  : cgmx command file \n");
  printf("==============================================================================\n\n");
  exit;
  
}

$CGBUILDPREFIX=$ENV{CGBUILDPREFIX};
$cgmx = "$CGBUILDPREFIX/mx/bin/cgmx";  # command for cgmx

$runCase="compare";   # compare/performance

$cmdFile="boxEigen.cmd";
$grid="square128.order2.hdf"; 
$check="mx.check";
$logFile="mx.log";
$matlabFile ="cgmxTiming";
$title=""; 
$option=""; 
$order=2;
$scheme="FD"; 
$numRes=3;

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
  elsif(  $arg =~ /-scheme=(.*)/ )
  {
    $scheme = $1;
    printf("Setting scheme=$scheme\n");
  }
  elsif(  $arg =~ /-numRes=(.*)/ )
  {
    $numRes = $1;
    printf("Setting numRes=$numRes\n");
  }


}

$date= localtime();   # here is the date
printf("date=[$date]\n");

if( $runCase eq "performance" )
{
  # ================== Compute Errors and CPU versus grid spacing ===============

  if( $scheme eq "FD" ){
    $option = "-diss=0."; 
  }
  elsif( $scheme eq "FDA" ){
    $option = "-diss=.5"; 
  }
  elsif( $scheme eq "UPC" ){
    $option = "-useSosupDissipation=1";
  }
  elsif( $scheme eq "SOSUP" ){
    $option = "-method=sosup";
  }
  else{ printf("Unknown scheme=[$scheme]\n"); exit; }

  $matlabFile = "$grid" . "Perf"; 
  $matlabFileName = $matlabFile . "$scheme$order"; 

  for( $res=0; $res<$numRes; $res++ )
  {
    if( ($grid eq "square") || ($grid eq "nonSquare") ){
      $n = 2**(3+$res+1); 
      $h[$res]= 1.0/($n);
    }
    elsif( ($grid eq "box") || ($grid eq "nonBox") ){
      $n = 2**($res); 
      $h[$res]= .1/($n);
    }
    else{ printf("Unknown grid=[$grid]\n"); exit; }

    if( $scheme eq "FDA" ){
      $diss = .01/$h[$res]; # scale AD by 1/h
      $option = "-diss=$diss"; 
    }

    $gridName = "$grid$n.order$order"; 
    if( $order eq 4 && ( $scheme eq "UPC" || $scheme eq "SOSUP" ) ){ $gridName .= ".ng3"; } # append ".ng3" 

    # --- Here is the cgmx command and command line arguments -------
    $cmd = "$cgmx -noplot $cmdFile -g=$gridName -mx=4 -my=4 -m=4 -n=4 -tf=$tf -tp=$tf $option -go=go >! $matlabFileName$i.perf.out"; 

    printf("\n+++++++++++++++++++++++++++ scheme=$scheme, order=$order, res=$res +++++++++++++++++++++++++++++++++++++++++++\n");
    printf("run [%s]\n",$cmd);
    $startTime=time();
    $returnValue = system("csh -f -c 'nohup $cmd'");
    $cpuTime=time()-$startTime;
    printf("...returnValue=%g, cpu=%g (s)\n",$returnValue,$cpuTime);

    # --------------------------------------------------------------------------------
    # --------------------- PROCESS CURRENT LOG FILE -------------------------------
    # --------------------------------------------------------------------------------
    parseLogFile($res);
 
    $line = $maxErrors[$res];
    @errv =  split(',',$line);   #extract errors in Ex, Ey, [Ez or Hz]

    $maxErr =0; 
    for( $j=0; $j<3; $j++){ $err=$errv[$j]; if( $err > $maxErr ){ $maxErr=$err; } }
    printf(" res=$res: maxErr=$maxErr\n");
    $err[$res]= $maxErr; 
  }

  # ---- Write a matlab file -----
  open(MFILE,">$matlabFileName.m") || die print "unable to open $matlabFile.m\n";
  print MFILE "%\n% File created by runCgmx.p on $date \n%\n"; 
  # print MFILE "clf; clear;\n"; 
  # print MFILE "fontSize=18; lineWidth=2;\n"; 
  
  print MFILE "tFinal=$tFinal[0];\n";
  
  print MFILE "h = [";  for( $res=0; $res<=$numRes; $res++ ){  print MFILE "$h[$res] "; } print MFILE "];\n";
  print MFILE "err = [";  for( $res=0; $res<=$numRes; $res++ ){  print MFILE "$err[$res] "; } print MFILE "];\n";

  print MFILE "cpuTotal = [";  for( $res=0; $res<=$numRes; $res++ ){ print MFILE "$cpuTotal[$res] "; } print MFILE "];\n";
  
  
  print MFILE "if( plotOption == 1 )\n";  # plot results if plotOption==0 

  # log-log plot
  print MFILE "  loglog(h,err,'b-x','linewidth',lineWidth);\n";
  
  print MFILE "  set(gca,'fontSize',fontSize);\n"; 
  if( $title eq "" ){ $plotTitle="Err versus h, $scheme$order, $grid"; }else{ $plotTitle="$title"; }
  print MFILE "  title('$plotTitle');\n";
  print MFILE "  grid on;\n";
  
  $plotName = $matlabFileName; 
  $plotName =~ s/\.m$//; # remove trailing ".m"
  print MFILE "  %% print('-depsc','$plotName.eps');\n"; 
  
  print MFILE "end % end if plotOption\n"; 

  close( MFILE );
  print "Matlab results written to file [$matlabFileName.m].\n";
}


if( $runCase eq "compare" )
{
  # ============ COMPARE SCHEMES =================
  $numCases=4;
  for( $i=0; $i<$numCases; $i++ )
  {
  
    if( $i == 0 ){
      $option = "-diss=0."; $scheme="FD"; 
    }
    elsif( $i==1 ){
      $option = "-diss=.5";  $scheme="FDA";
    }
    elsif( $i==2 ){
      $option = "-useSosupDissipation=1"; $scheme="UPC";
    }
    elsif( $i==3 ){
      $option = "-method=sosup"; $scheme="SOSUP";
    }
  
    $gridName = "$grid.order$order"; 
    if( $order eq 4 && ( $scheme eq "UPC" || $scheme eq "SOSUP" ) ){ $gridName .= ".ng3"; } # append ".ng3" 

    # --- Here is the cgmx command and command line arguments -------
    $cmd = "$cgmx -noplot $cmdFile -g=$gridName -mx=3 -my=2 -m=3 -n=2 -tf=$tf -tp=$tf $option -go=go >! $matlabFile$i.out"; 
  
    printf("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("run [%s]\n",$cmd);
    $startTime=time();
    $returnValue = system("csh -f -c 'nohup $cmd'");
    $cpuTime=time()-$startTime;
    printf("...returnValue=%g, cpu=%g (s)\n",$returnValue,$cpuTime);
  
    # --------------------------------------------------------------------------------
    # --------------------- PROCESS CURRENT LOG FILE -------------------------------
    # --------------------------------------------------------------------------------
    parseLogFile($i);
  
  }
  
  # ---- Write a matlab file -----
  open(MFILE,">$matlabFile.m") || die print "unable to open $matlabFile.m\n";
  print MFILE "%\n% File created by runCgmx.p on $date \n%\n"; 
  # print MFILE "clf; clear;\n"; 
  # print MFILE "fontSize=18; lineWidth=2;\n"; 
  
  print MFILE "tFinal=$tFinal[0];\n";
  
  for( $ii=0; $ii<$numCases; $ii++ ){ $jj=$ii+1; print MFILE "maxErrors$jj=[$maxErrors[$ii]];\n"; }
  
  print MFILE "cpuTotal = [";
  for( $ii=0; $ii<$numCases; $ii++ ){ print MFILE "$cpuTotal[$ii] "; } print MFILE "];\n";
  
  print MFILE "cpuPerStep = [";
  for( $ii=0; $ii<$numCases; $ii++ ){ print MFILE "$cpuPerStep[$ii] ";} print MFILE "];\n";
  
  
  print MFILE "relCPU = [";
  for( $ii=0; $ii<$numCases; $ii++ ){ $relCPU[$ii]=$cpuTotal[$ii]/$cpuTotal[0]; print MFILE "$relCPU[$ii] ";} print MFILE "];\n";
  
  # Bar chart
  print MFILE "bar(diag(relCPU),'stacked');\n";
  
  print MFILE "set(gca,'fontSize',fontSize);\n"; 
  if( $title eq "" ){ $plotTitle="Cgmx: relative CPU times"; }else{ $plotTitle="Rel-CPU, $title"; }
  print MFILE "title('$plotTitle');\n";
  
  
  print MFILE "legend('FD','FDA','UPC','SOSUP','Location','northwest');\n";
  print MFILE "grid on;\n";
  
  
  #-print MFILE "xlabel('\\omega');\n";
  #-print MFILE "grid on;\n";
  
  $plotName = $matlabFile; 
  $plotName =~ s/\.m$//; # remove trailing ".m"
  print MFILE "print('-depsc','$plotName.eps');\n"; 
  
  close( MFILE );
  print "Matlab results written to file [$matlabFile.m].\n";

}

exit;
# ===============================================================================================================
