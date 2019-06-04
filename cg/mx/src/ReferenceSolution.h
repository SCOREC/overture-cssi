#ifndef REFERENCE_SOLUTION_H
#define REFERENCE_SOLUTION_H

//
//  Class to hold a "reference solution" (e.g. a solution from a show file)
//  for comparison to other solutions (on plossible different grids)
//

#include "Overture.h"

// forward declarations
class ShowFileReader;
class InterpolatePointsOnAGrid;

class ReferenceSolution
{
public:

ReferenceSolution();

~ReferenceSolution();

int setShowFileName( const aString & name );

realCompositeGridFunction & 
getSolution( real t, CompositeGrid & cgTarget, Range & C );


private:

aString nameOfShowFile;
ShowFileReader *showFileReader;
realCompositeGridFunction *uTarget;

int currentSolution;
real currentTime, tPlot;

InterpolatePointsOnAGrid *interpolatePointsOnAGrid;

};

  
#endif
