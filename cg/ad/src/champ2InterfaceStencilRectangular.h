// build the coefficient matrix for 9 point stencil.
// axis=0 if left and right. axis=1 if top and bottom.
// side=0 for left. side=1 for right.
// Real coefB[3][3]={0,0,0,0,0,0,0,0,0}; // initialize
// #define coefA(i,j) coefB[i+1][j+1]
Real dxl,dyl;
Real a00,a01,a02,a10,a11,a12,a20;
int is=1;


printF(" Running Sijia's Cartesian script\n ");

dxl = dx[0];
dyl = dx[1];

/*
if (side==0){
	is=1;
}
else{
	is=-1;
}
*/
if (side==0){
	is=-1;
}
else{
	is=1;
}

printF("Sl=%g, xl=%g, dyl=%g, side1=%i, axis1=%i, axis2=%i theta=%g, beta=%g\n", Sl,dxl, dyl,side,axis,axis2, theta,beta);


// MAPLE input 
if (axis==0){
	a00  = Sl;
 
	a01  = 0;
	 
	a02  = (double) (dxl * (beta - 1) * (Sl * dxl + 2 )) / 0.2e1;
	 
	a10  = theta * (Sl * dxl  + 1);
	 
	a11  = 0;
	 
	a12  = 0;
	 
	a20  = (double) (dxl * beta * (Sl * dxl + 2)) / 0.2e1;
}	 
else{
	a00  = Sl;
	a01  = 0;
	a02  = (dyl * (beta - 1) * (Sl * dyl + 2)) / 0.2e1;
	a10  = theta * (Sl * dyl + 1);
	a11  = 0;
	a20  = (beta * dyl * (Sl * dyl + 2)) / 0.2e1;
}


printF("a00=%g, a01=%g, a02=%g, a10=%g, a11=%g, a20=%g,\n", a00,a01,a02,a10,a11,a20);


if (axis==0){
	coefA( 0, 0) = a00 - 2/SQR(dxl) *a20 -2/SQR(dyl)*a02;
	coefA(-1, 0) =-is*a10/(2*dxl) + a20/SQR(dxl);
	coefA( 1, 0) = is*a10/(2*dxl) + a20/SQR(dxl);
	coefA( 0, 1) = a01/(2*dyl) + a02/SQR(dyl);
	coefA( 0,-1) =-a01/(2*dyl) + a02/SQR(dyl);
}
else{
	coefA( 0, 0) = a00 - 2/SQR(dyl) *a20 -2/SQR(dxl)*a02;
	coefA( 0,-1) =-is*a10/(2*dyl) + a20/SQR(dyl);
	coefA( 0, 1) = is*a10/(2*dyl) + a20/SQR(dyl);
	coefA( 1, 0) = a01/(2*dxl) + a02/SQR(dxl);
	coefA(-1, 0) =-a01/(2*dxl) + a02/SQR(dxl);


}