% convert memory (Mb) and processors and grid points to reals-per-grid-point
function mem( x,np,gp )


rpg = (x*1024.*1024.)/(gp*1e6*8.);

fprintf(' mem=%e MB, np=%d, grid-points=%d, reals/grid-point = rpg=%7.1f \n',x,np,gp,rpg);


