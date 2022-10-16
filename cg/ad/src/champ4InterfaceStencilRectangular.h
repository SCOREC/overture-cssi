// build the coefficient matrix for 9 point stencil.
// axis=0 if left and right. axis=1 if top and bottom.
// side=0 for left. side=1 for right.
// Real coefB[3][3]={0,0,0,0,0,0,0,0,0}; // initialize
// #define coefA(i,j) coefB[i+1][j+1]
Real dxl,dyl,hI;
Real a00,a01,a02,a03,a04,a10,a11,a12,a13,a14,a20,a21,a22,a23,a24,a30,a31,a32,a33,a34,a40,a41,a42,a43,a44;
int is=1;
int shftI4=2;
printF("========== Running Sijia's 4th-order Rectangular script ==========\n");

dxl = dx[0];
dyl = dx[1];
hI = dxs;

if (side==0){
    is=-1;
}
else{
    is=1;
}

//printF("Sl=%g, xl=%g, dyl=%g, side1=%i, axis1=%i, axis2=%i theta=%g, beta=%g\n", Sl,dxl, dyl,side,axis,axis2, theta,beta);


// MAPLE input 
a00  = Sl;
 
a01  = 0;
 
a02  =  (hI * (beta - 1) * (Sl * hI + 2)) / 0.2e1;
 
a03  = 0;
 
a04  = pow(hI, 0.3e1) *   pow( (beta - 1),  2) * (Sl * hI + 0.4e1) / 0.24e2;
 
a10  = is * theta * (Sl * hI + 1);
 
a11  = 0;
 
a12  =  (hI * hI * is * theta * (beta - 1) * (Sl * hI + 3)) / 0.6e1;
 
a13  = 0;
 
a14  = 0;
 
a20  =  (beta * (Sl * hI + 2) * hI) / 0.2e1;
 
a21  = 0;
 
a22  =  beta * pow(hI, 0.3e1) *  (beta - 1) * (Sl * hI + 0.4e1) / 0.12e2;
 
a23  = 0;
 
a24  = 0;
 
a30  =  (beta * is * theta * (Sl * hI + 3) * hI * hI) / 0.6e1;
 
a31  = 0;
 
a32  = 0;
 
a33  = 0;
 
a34  = 0;
 
a40  = pow(hI, 0.3e1) * beta * beta * (Sl * hI + 0.4e1) / 0.24e2;
 
a41  = 0;
 
a42  = 0;
 
a43  = 0;
 
a44  = 0;


printF("     side1=%i, side2=%i, axis1=%i, axis2=%i, is=%i\n", side,side2,axis,axis2,is );
printF("     theta=%g, beta=%g, Sl=%g\n", theta,beta,Sl);
printF("     dxl=%g, dyl=%g, dxs=%g\n", dxl, dyl,dxs);
 
/* 
printF("    a00=%g, a01=%g, a02=%g, a03=%g, a04=%g \n", a00,a01,a02,a03,a04);
printF("    a10=%g, a11=%g, a12=%g, a13=%g, a14=%g \n", a10,a11,a12,a13,a14);
printF("    a20=%g, a21=%g, a22=%g, a23=%g, a24=%g \n", a20,a21,a22,a23,a24);
printF("    a30=%g, a31=%g, a32=%g, a33=%g, a34=%g \n", a30,a31,a32,a33,a34);
printF("    a40=%g, a41=%g, a42=%g, a43=%g, a44=%g \n", a40,a41,a42,a43,a44);
*/

Real Dcc[5][5] = {0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.,
                  0.,0.,1.,0.,0.,
                  0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.}; // initialize

Real Dx4cc[5][5]={0.,0., 1./(12.*dxl),0.,0.,
                  0.,0.,-8./(12.*dxl),0.,0.,
                  0.,0.,0.,0.,0.,
                  0.,0., 8./(12.*dxl),0.,0.,
                  0.,0.,-1./(12.*dxl),0.,0.}; 

Real Dy4cc[5][5]={0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.,
                  1./(12.*dyl),-8./(12.*dyl),0.,8./(12.*dyl),-1./(12.*dyl),
                  0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.};

Real Dxx4cc[5][5]={0.,0., -1./(12.*SQR(dxl)),0.,0.,
                   0.,0., 16./(12.*SQR(dxl)),0.,0.,
                   0.,0.,-30./(12.*SQR(dxl)),0.,0.,
                   0.,0., 16./(12.*SQR(dxl)),0.,0.,
                   0.,0., -1./(12.*SQR(dxl)),0.,0.};
Real Dyy4cc[5][5]={0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.,
                  -1./(12.*SQR(dyl)),16./(12.*SQR(dyl)),-30./(12.*SQR(dyl)),16./(12.*SQR(dyl)),-1./(12.*SQR(dyl)),
                  0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.};

Real Dxxx2cc[5][5]={0.,0.,(-1./2.)/(pow(dxl,3)),0.,0.,
                   0.,0.,       1./(pow(dxl,3)),0.,0.,
                   0.,0.,0.,0.,0.,
                   0.,0.,      -1./(pow(dxl,3)),0.,0.,
                   0.,0.,  (1./2.)/(pow(dxl,3)),0.,0.};
Real Dyyy2cc[5][5]={0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.,
                  (-1./2.)/(pow(dyl,3)),1./(pow(dyl,3)),0.,-1./(pow(dyl,3)),(1./2.)/(pow(dyl,3)),
                  0.,0.,0.,0.,0.,
                  0.,0.,0.,0.,0.}; 
Real Dxxxx2cc[5][5]={0.,0., 1./(pow(dxl,4)),0.,0.,
                     0.,0.,-4./(pow(dxl,4)),0.,0.,
                     0.,0., 6./(pow(dxl,4)),0.,0.,
                     0.,0.,-4./(pow(dxl,4)),0.,0.,
                     0.,0., 1./(pow(dxl,4)),0.,0.};
Real Dyyyy2cc[5][5]= {0.,0.,0.,0.,0.,
                      0.,0.,0.,0.,0.,
                      1./(pow(dyl,4)),-4./(pow(dyl,4)),6./(pow(dyl,4)) ,-4./(pow(dyl,4)),1./(pow(dyl,4)),
                      0.,0.,0.,0.,0.,
                      0.,0.,0.,0.,0.};

Real Dxy4cc[5][5]={ 1./(144.*dxl*dyl), -8./(144.*dxl*dyl),0.,  8./(144.*dxl*dyl),-1./(144.*dxl*dyl),
                   -8./(144.*dxl*dyl), 64./(144.*dxl*dyl),0.,-64./(144.*dxl*dyl), 8./(144.*dxl*dyl),
                    0.,0.,0.,0.,0.,
                    8./(144.*dxl*dyl),-64./(144.*dxl*dyl),0., 64./(144.*dxl*dyl),-8./(144.*dxl*dyl),
                   -1./(144.*dxl*dyl),  8./(144.*dxl*dyl),0., -8./(144.*dxl*dyl), 1./(144.*dxl*dyl)};

Real Dxxy4cc[5][5]={ -1./(144.*SQR(dxl)*dyl),   8./(144.*SQR(dxl)*dyl),0.,  -8./(144.*SQR(dxl)*dyl),  1./(144.*SQR(dxl)*dyl),
                     16./(144.*SQR(dxl)*dyl),-128./(144.*SQR(dxl)*dyl),0., 128./(144.*SQR(dxl)*dyl),-16./(144.*SQR(dxl)*dyl),
                    -30./(144.*SQR(dxl)*dyl), 240./(144.*SQR(dxl)*dyl),0.,-240./(144.*SQR(dxl)*dyl), 30./(144.*SQR(dxl)*dyl),
                     16./(144.*SQR(dxl)*dyl),-128./(144.*SQR(dxl)*dyl),0., 128./(144.*SQR(dxl)*dyl),-16./(144.*SQR(dxl)*dyl),
                     -1./(144.*SQR(dxl)*dyl),   8./(144.*SQR(dxl)*dyl),0.,  -8./(144.*SQR(dxl)*dyl),  1./(144.*SQR(dxl)*dyl)};

Real Dxyy4cc[5][5]={ -1./(144.*dxl*SQR(dyl)),  16./(144.*dxl*SQR(dyl)), -30./(144.*dxl*SQR(dyl)),  16./(144.*dxl*SQR(dyl)),-1./(144.*dxl*SQR(dyl)),
                      8./(144.*dxl*SQR(dyl)),-128./(144.*dxl*SQR(dyl)), 240./(144.*dxl*SQR(dyl)),-128./(144.*dxl*SQR(dyl)), 8./(144.*dxl*SQR(dyl)),
                      0.,0.,0.,0.,0.,
                     -8./(144.*dxl*SQR(dyl)), 128./(144.*dxl*SQR(dyl)),-240./(144.*dxl*SQR(dyl)), 128./(144.*dxl*SQR(dyl)),-8./(144.*dxl*SQR(dyl)),
                      1./(144.*dxl*SQR(dyl)), -16./(144.*dxl*SQR(dyl)),  30./(144.*dxl*SQR(dyl)), -16./(144.*dxl*SQR(dyl)), 1./(144.*dxl*SQR(dyl))};


Real Dxxyy4cc[5][5]={  1./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  30./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  1./(144.*SQR(dxl)*SQR(dyl)),
                     -16./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-16./(144.*SQR(dxl)*SQR(dyl)),
                      30./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 900./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 30./(144.*SQR(dxl)*SQR(dyl)),
                     -16./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-16./(144.*SQR(dxl)*SQR(dyl)),
                       1./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  30./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  1./(144.*SQR(dxl)*SQR(dyl))};

Real Dxyyy2cc[5][5]={-1./(24.*dxl*pow(dyl,3)),  2./(24.*dxl*pow(dyl,3)),0., -2./(24.*dxl*pow(dyl,3)), 1./(24.*dxl*pow(dyl,3)),
                      8./(24.*dxl*pow(dyl,3)),-16./(24.*dxl*pow(dyl,3)),0., 16./(24.*dxl*pow(dyl,3)),-8./(24.*dxl*pow(dyl,3)),
                      0.,0.,0.,0.,0.,
                     -8./(24.*dxl*pow(dyl,3)), 16./(24.*dxl*pow(dyl,3)),0.,-16./(24.*dxl*pow(dyl,3)), 8./(24.*dxl*pow(dyl,3)),
                      1./(24.*dxl*pow(dyl,3)), -2./(24.*dxl*pow(dyl,3)),0.,  2./(24.*dxl*pow(dyl,3)),-1./(24.*dxl*pow(dyl,3))};

Real Dxxxy2cc[5][5]={-1./(24.*pow(dxl,3)*dyl),  8./(24.*pow(dxl,3)*dyl),0., -8./(24.*pow(dxl,3)*dyl), 1./(24.*pow(dxl,3)*dyl),
                      2./(24.*pow(dxl,3)*dyl),-16./(24.*pow(dxl,3)*dyl),0., 16./(24.*pow(dxl,3)*dyl),-2./(24.*pow(dxl,3)*dyl),
                      0.,0.,0.,0.,0.,
                     -2./(24.*pow(dxl,3)*dyl), 16./(24.*pow(dxl,3)*dyl),0.,-16./(24.*pow(dxl,3)*dyl), 2./(24.*pow(dxl,3)*dyl),
                      1./(24.*pow(dxl,3)*dyl), -8./(24.*pow(dxl,3)*dyl),0.,  8./(24.*pow(dxl,3)*dyl),-1./(24.*pow(dxl,3)*dyl)};

//printF("hfwd1=%i,hfwd1=%i\n",shftI4,shftI4);

#define Dc(i,j) Dcc[i+shftI4][j+shftI4]
#define Dx4c(i,j) Dx4cc[i+shftI4][j+shftI4]
#define Dy4c(i,j) Dy4cc[i+shftI4][j+shftI4]
#define Dxx4c(i,j) Dxx4cc[i+shftI4][j+shftI4]
#define Dyy4c(i,j) Dyy4cc[i+shftI4][j+shftI4]
#define Dxxx2c(i,j) Dxxx2cc[i+shftI4][j+shftI4]
#define Dyyy2c(i,j) Dyyy2cc[i+shftI4][j+shftI4]
#define Dxxxx2c(i,j) Dxxxx2cc[i+shftI4][j+shftI4]
#define Dyyyy2c(i,j) Dyyyy2cc[i+shftI4][j+shftI4]
#define Dxy4c(i,j) Dxy4cc[i+shftI4][j+shftI4]
#define Dxxy4c(i,j) Dxxy4cc[i+shftI4][j+shftI4]    
#define Dxyy4c(i,j) Dxyy4cc[i+shftI4][j+shftI4]    
#define Dxxyy4c(i,j) Dxxyy4cc[i+shftI4][j+shftI4]    
#define Dxxxy2c(i,j) Dxxxy2cc[i+shftI4][j+shftI4]    
#define Dxyyy2c(i,j) Dxyyy2cc[i+shftI4][j+shftI4]        
   

if (axis==0){
    for (int i = -2; i < 3; i++) {
        for (int j = -2; j < 3; j++) {
        
            /*
            printF("    a00=%g, a01=%g, a02=%g, a03=%g, a04=%g \n", Dc(i,j),Dy4c(i,j),Dyy4c(i,j),Dyyy2c(i,j),Dyyyy2c(i,j));
            printF("    a10=%g, a11=%g, a12=%g, a13=%g\n", Dx4c(i,j),Dxy4c(i,j),Dxyy4c(i,j),Dxyyy2c(i,j));
            printF("    a20=%g, a21=%g, a22=%g\n", Dxx4c(i,j),Dxxy4c(i,j),Dxxyy4c(i,j));
            printF("    a30=%g, a31=%g \n", Dxxx2c(i,j),Dxxxy2c(i,j));
            printF("    a40=%g \n", Dxxxx2c(i,j));
            */

            coefA(i,j) = Dc(i,j)*a00+Dx4c(i,j)*a10+Dxx4c(i,j)*a20+Dxxx2c(i,j)*a30+Dxxxx2c(i,j)*a40\
                                    +Dy4c(i,j)*a01+Dyy4c(i,j)*a02+Dyyy2c(i,j)*a03+Dyyyy2c(i,j)*a04\
                                    +Dxxy4c(i,j)*a21+Dxyy4c(i,j)*a12\
                                    +Dxy4c(i,j)*a11+Dxxyy4c(i,j)*a22\
                                    +Dxxxy2c(i,j)*a31+Dxyyy2c(i,j)*a13;
        }
    } 
}
 else{
    for (int i = -2; i < 3; i++) {
        for (int j = -2; j < 3; j++) {
            coefA(j,i) = Dc(j,i)*a00+Dy4c(j,i)*a10+Dyy4c(j,i)*a20+Dyyy2c(j,i)*a30+Dyyyy2c(j,i)*a40\
                                    +Dx4c(j,i)*a01+Dxx4c(j,i)*a02+Dxxx2c(j,i)*a03+Dxxxx2c(j,i)*a04\
                                    +Dxyy4c(j,i)*a21+Dxxy4c(j,i)*a12\
                                    +Dxy4c(j,i)*a11+Dxxyy4c(j,i)*a22\
                                    +Dxyyy2c(j,i)*a31+Dxxxy2c(j,i)*a13;
        }
    } 

 }

/*
if (axis==0){ 
    coefA( 2, 0) = Dx4c(2, 0)*a10+Dxx4c(2, 0)*a20+Dxxx2c(2, 0)*a30+Dxxxx2c(2, 0)*a40+Dxyy4c( 2, 0)*a12+Dxxyy4c( 2, 0)*a22;
    coefA( 1, 0) = Dx4c(1, 0)*a10+Dxx4c(1, 0)*a20+Dxxx2c(1, 0)*a30+Dxxxx2c(1, 0)*a40+Dxyy4c( 1, 0)*a12+Dxxyy4c( 1, 0)*a22;
    coefA( 0, 0) = a00+Dxx4c(0,0)*a20+Dyy4c(0,0)*a02+Dxxxx2c(0,0)*a40+Dyyyy2c(0,0)*a04+Dxxyy4c(0,0)*a22;
    coefA(-1, 0) = Dx4c(-1, 0)*a10+Dxx4c(-1, 0)*a20+Dxxx2c(-1, 0)*a30+Dxxxx2c(-1, 0)*a40+Dxyy4c(-1, 0)*a12+Dxxyy4c(-1, 0)*a22;
    coefA(-2, 0) = Dx4c(-2, 0)*a10+Dxx4c(-2, 0)*a20+Dxxx2c(-2, 0)*a30+Dxxxx2c(-2, 0)*a40+Dxyy4c(-2, 0)*a12+Dxxyy4c(-2, 0)*a22;
    coefA( 0, 2) = Dy4c(0, 2)*a01+Dyy4c(0, 2)*a02+Dyyy2c(0, 2)*a03+Dyyyy2c(0, 2)*a04+Dxyy4c( 0, 2)*a21+Dxxyy4c( 0, 2)*a22;   
    coefA( 0, 1) = Dy4c(0,1)*a01+Dyy4c(0,1)*a02+Dyyy2c(0,1)*a03+Dyyyy2c(0,1)*a04+Dxyy4c( 0, 1)*a21+Dxxyy4c( 0, 1)*a22;  
    coefA( 0,-1) = Dy4c(0,-1)*a01+Dyy4c(0,-1)*a02+Dyyy2c(0,-1)*a03+Dyyyy2c(0,-1)*a04+Dxyy4c( 0,-1)*a21+Dxxyy4c( 0,-1)*a22;  
    coefA( 0,-2) = Dy4c(0,-2)*a01+Dyy4c(0,-2)*a02+Dyyy2c(0,-2)*a03+Dyyyy2c(0,-2)*a04+Dxyy4c( 0,-2)*a21+Dxxyy4c( 0,-2)*a22;  

    coefA( 1, 1) = Dxy4c( 1, 1)*a11+Dxxy4c( 1, 1)*a21+Dxyy4c( 1, 1)*a12+Dxxyy4c( 1, 1)*a22+Dxxxy2c( 1, 1)*a31+Dxyyy2c( 1, 1)*a13;
    coefA( 1,-1) = Dxy4c( 1,-1)*a11+Dxxy4c( 1,-1)*a21+Dxyy4c( 1,-1)*a12+Dxxyy4c( 1,-1)*a22+Dxxxy2c( 1,-1)*a31+Dxyyy2c( 1,-1)*a13;
    coefA(-1, 1) = Dxy4c(-1, 1)*a11+Dxxy4c(-1, 1)*a21+Dxyy4c(-1, 1)*a12+Dxxyy4c(-1, 1)*a22+Dxxxy2c(-1, 1)*a31+Dxyyy2c(-1, 1)*a13;
    coefA(-1,-1) = Dxy4c(-1,-1)*a11+Dxxy4c(-1,-1)*a21+Dxyy4c(-1,-1)*a12+Dxxyy4c(-1,-1)*a22+Dxxxy2c(-1,-1)*a31+Dxyyy2c(-1,-1)*a13;
    coefA( 1,-2) = Dxy4c( 1,-2)*a11+Dxxy4c( 1,-2)*a21+Dxyy4c( 1,-2)*a12+Dxxyy4c( 1,-2)*a22+Dxxxy2c( 1,-2)*a31+Dxyyy2c( 1,-2)*a13;
    coefA( 1, 2) = Dxy4c( 1, 2)*a11+Dxxy4c( 1, 2)*a21+Dxyy4c( 1, 2)*a12+Dxxyy4c( 1, 2)*a22+Dxxxy2c( 1, 2)*a31+Dxyyy2c( 1, 2)*a13;
    coefA(-1,-2) = Dxy4c(-1,-2)*a11+Dxxy4c(-1,-2)*a21+Dxyy4c(-1,-2)*a12+Dxxyy4c(-1,-2)*a22+Dxxxy2c(-1,-2)*a31+Dxyyy2c(-1,-2)*a13;
    coefA(-1, 2) = Dxy4c(-1, 2)*a11+Dxxy4c(-1, 2)*a21+Dxyy4c(-1, 2)*a12+Dxxyy4c(-1, 2)*a22+Dxxxy2c(-1, 2)*a31+Dxyyy2c(-1, 2)*a13;
    coefA( 2,-2) = Dxy4c( 2,-2)*a11+Dxxy4c( 2,-2)*a21+Dxyy4c( 2,-2)*a12+Dxxyy4c( 2,-2)*a22+Dxxxy2c( 2,-2)*a31+Dxyyy2c( 2,-2)*a13;
    coefA( 2,-1) = Dxy4c( 2,-1)*a11+Dxxy4c( 2,-1)*a21+Dxyy4c( 2,-1)*a12+Dxxyy4c( 2,-1)*a22+Dxxxy2c( 2,-1)*a31+Dxyyy2c( 2,-1)*a13;
    coefA( 2, 1) = Dxy4c( 2, 1)*a11+Dxxy4c( 2, 1)*a21+Dxyy4c( 2, 1)*a12+Dxxyy4c( 2, 1)*a22+Dxxxy2c( 2, 1)*a31+Dxyyy2c( 2, 1)*a13;
    coefA( 2, 2) = Dxy4c( 2, 2)*a11+Dxxy4c( 2, 2)*a21+Dxyy4c( 2, 2)*a12+Dxxyy4c( 2, 2)*a22+Dxxxy2c( 2, 2)*a31+Dxyyy2c( 2, 2)*a13;
    coefA(-2,-2) = Dxy4c(-2,-2)*a11+Dxxy4c(-2,-2)*a21+Dxyy4c(-2,-2)*a12+Dxxyy4c(-2,-2)*a22+Dxxxy2c(-2,-2)*a31+Dxyyy2c(-2,-2)*a13;
    coefA(-2,-1) = Dxy4c(-2,-1)*a11+Dxxy4c(-2,-1)*a21+Dxyy4c(-2,-1)*a12+Dxxyy4c(-2,-1)*a22+Dxxxy2c(-2,-1)*a31+Dxyyy2c(-2,-1)*a13;
    coefA(-2, 1) = Dxy4c(-2, 1)*a11+Dxxy4c(-2, 1)*a21+Dxyy4c(-2, 1)*a12+Dxxyy4c(-2, 1)*a22+Dxxxy2c(-2, 1)*a31+Dxyyy2c(-2, 1)*a13;
    coefA(-2, 2) = Dxy4c(-2, 2)*a11+Dxxy4c(-2, 2)*a21+Dxyy4c(-2, 2)*a12+Dxxyy4c(-2, 2)*a22+Dxxxy2c(-2, 2)*a31+Dxyyy2c(-2, 2)*a13;
    printF("    a00=%g, a01=%g, a02=%g, a03=%g, a04=%g \n", Dc(0,-2),Dy4c(0,-2),Dyy4c(0,-2),Dyyy2c(0,-2),Dyyyy2c(0,-2));
    printF("    a10=%g, a11=%g, a12=%g, a13=%g\n", Dx4c(0,-2),Dxy4c(0,-2),Dxyy4c(0,-2),Dxyy4c(0,-2));
    printF("    a20=%g, a21=%g, a22=%g\n", Dxx4c(0,-2),Dxxy4c(0,-2),Dxxyy4c(0,-2));
    printF("    a30=%g, a31=%g \n", Dxxx2c(0,-2),Dxxxy2c(0,-2));
    printF("    a40=%g \n", Dxxxx2c(0,-2));     
}
*/