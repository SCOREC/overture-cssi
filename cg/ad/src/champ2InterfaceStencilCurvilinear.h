// fill in the matrix coefficients
//
if( debug & 8 )
  printF("========== Running Sijia's Mapped script ==========\n");
// Real coefB[3][3]={0,0,0,0,0,0,0,0,0}; // initialize
// #define coefA(i,j) coefB[i+1][j+1]

Real drl,dsl,drr,dsr;
Real a00,a01,a02,a10,a11,a12,a20,a21,a22;
//Real hI,SlI,b1LI=abs(b1L),b1RI=abs(b1R);
Real hI,SlI,b1LI=b1L,b1RI=b1R;

int is=1;

drl = drL[0];
dsl = drL[1];

drr = drR[0];
dsr = drR[1];

hI = drr;

if (side==0){
	is=-1;
	b1LI = is*b1L;
	b1RI = is*b1R;
	if (Sl<0){
		SlI=-1.*Sl/b1RI;
	}
	else{
		SlI=Sl/b1RI;
	}
}
else{
	is=1;
	if (Sl<0){
		SlI=-1.*Sl/b1RI;
	}
	else{
		SlI=Sl/b1RI;
	}
}

/*
if (side==0){
	is=-1;
	if (Sl<0){
		SlI=-1.*Sl/b1RI;
	}
	else{
		SlI=Sl/b1RI;
	}
}
else{
	is=1;
	if (Sl<0){
		SlI=-1.*Sl/b1RI;
	}
	else{
		SlI=Sl/b1RI;
	}
}
*/

if( debug & 8 )
{
	printF("    side1=%i, side2=%i, axis1=%i, axis2=%i, is=%i \n", side,side2,axis1,axis2,is );
	printF("    theta=%g, beta=%g, Sl=%g, SlI=%g, hI=%g \n", theta,beta,Sl,SlI,hI);
	printF("    drl=%g, dsl=%g, drr=%g, dsr=%g \n", drl, dsl,drr,dsr);
}

a00  = b1RI * SlI;
 
a01  = (double) (b1RI * (2 * hI * c12R * (b2L * theta - b2R) * (SlI * hI + 2) * b1Rr2 + b1RI * (-2 * hI * c12R * theta * (SlI * hI + 2) * b2Lr2 + 2 * hI * c12R * (SlI * hI + 2) * b2Rr2 + hI * (c2L * beta - c2R) * (SlI * hI + 2) * b1RI - ((-2 * SlI * hI * is - 2 * is) * c11R + hI * c1R * (SlI * hI + 2)) * (b2L * theta - b2R))) *  pow((double) b1RI, (double) (-2)) / c11R) / 0.2e1;
 
a02  = -b1RI * (double) hI * (double) (SlI * hI + 2) / b1RI * ((-c22L * beta / 0.2e1 + c22R / 0.2e1) * b1RI + c12R * (b2L * theta - b2R)) / c11R;
 
a10  = (double) (b1RI * (-2 * b1RI * hI * c12R * theta * (SlI * hI + 2) * b1Lr2 + 2 * b1LI * hI * c12R * theta * (SlI * hI + 2) * b1Rr2 + b1RI * (c1L * hI * beta * (SlI * hI + 2) * b1RI - b1LI * ((-2 * SlI * hI * is - 2 * is) * c11R + hI * c1R * (SlI * hI + 2)) * theta)) * pow((double) b1RI, (double) (-2)) / c11R) / 0.2e1;
 
a11  = (double) (b1RI * hI * (SlI * hI + 2) * (2 * c12L * beta * b1RI - 2 * b1LI * c12R * theta) / b1RI / c11R) / 0.2e1;
 
a12  = 0;
 
a20  = (double) (b1RI * hI * (SlI * hI + 2) * c11L * beta / c11R) / 0.2e1;
 
a21  = 0;
 
a22  = 0;
 
if( debug & 8 )
{
  printF("    a00=%g, a01=%g, a02=%g, a10=%g, a11=%g, a12=%g, a20=%g \n", a00,a01,a02,a10,a11,a12,a20);
}

if (axis==0){
	coefA( 0, 0) = a00 - 2.0/SQR(drl) *a20 -2.0/SQR(dsl)*a02;
	coefA(-1, 0) =-a10/(2.0*drl) + a20/SQR(drl);
	coefA( 1, 0) = a10/(2.0*drl) + a20/SQR(drl);
	coefA( 0, 1) = a01/(2.0*dsl) + a02/SQR(dsl);
	coefA( 0,-1) =-a01/(2.0*dsl) + a02/SQR(dsl);
	coefA( 1, 1) = a11/(4*drl*dsl);
	coefA( 1,-1) =-a11/(4*drl*dsl);
	coefA(-1, 1) =-a11/(4*drl*dsl);
	coefA(-1,-1) = a11/(4*drl*dsl);
}
else{
	coefA( 0, 0) = a00 - 2.0/SQR(dsl) *a20 -2.0/SQR(drl)*a02;
	coefA( 0,-1) =-a10/(2.0*dsl) + a20/SQR(dsl);
	coefA( 0, 1) = a10/(2.0*dsl) + a20/SQR(dsl);
	coefA( 1, 0) = a01/(2.0*drl) + a02/SQR(drl);
	coefA(-1, 0) =-a01/(2.0*drl) + a02/SQR(drl);
	coefA( 1, 1) = a11/(4.*drl*drl);
	coefA(-1, 1) =-a11/(4.*drl*drl);
	coefA( 1,-1) =-a11/(4.*drl*drl);
	coefA(-1,-1) = a11/(4.*drl*drl);


}


 
