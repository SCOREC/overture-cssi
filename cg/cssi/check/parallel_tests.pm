# Parallel regression tests:
$CGBUILDPREFIX = $ENV{"CGBUILDPREFIX"};
$program = "$CGBUILDPREFIX/cssi/bin/cgcssi"; 

@cmdFiles=(
 	     "squarej.cssi",    # compressible, Jameson with TZ
 	     "squarej.cssi", 
	     "cicej.cssi",
	     "cicej.cssi",
	     "squareg.cssi",    # compressible, Godunov with TZ
	     "squareg.cssi",
	     "cicej.cssi",
	     "cicej.cssi"  
	     ); 
  # specify the number of processors to use in each of the above cases 
@numProc=(
             "1",
             "2",
             "1",
             "2",
             "1",
             "2",
             "1",
             "2",
             "1",
             "2",
             "1",
             "2"
	    );
