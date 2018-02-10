// ============================================================
// ============================================================
#beginMacro initializeAnnulusEigenfunction()

 real sint = sin(omega*t), cost = cos(omega*t);
 real sintp = omega*cost, costp = -omega*sint;

 real sintm = sin(omega*(t-dt)), costm = cos(omega*(t-dt));

 real sr,si,psir[10],psii[10], ct,st,expt, ctm,stm,exptm;
 real ampH, ampE, ampHm, ampEm, ampHp, ampEp, ampHmp, ampEmp;
 real ampP[10], ampPm[10];
 if( dispersionModel==noDispersion )
 {
   if( numberOfDimensions==2 )
   {
     ampH  = cost;   ampHp  =-omega*sint;
     ampE  = sint;   ampEp  = omega*cost;

     ampHm = costm;  ampHmp =-omega*sintm;
     ampEm = sintm;  ampEmp = omega*costm;
   }
   else
   {
     ampE  = cost;   ampEp  = -omega*sint; 
     ampEm = costm;  ampEmp = -omega*sintm;
  }
   
 }
 else 
 {
   // --- DISPERSIVE ----
   DispersiveMaterialParameters & dmp = getDispersiveMaterialParameters(grid);

   // Evaluate the dispersion relation for "s"
   assert( dmp.numberOfPolarizationVectors<10 );
   
   const real kk = omega/c; //  *CHECK ME* 
   dmp.evaluateDispersionRelation( c,kk, sr, si, psir,psii ); 

   if( t<3.*dt )
     printF("--DISK-EIGEN-- (dispersive) t=%10.3e, sr=%g, si=%g psir[0]=%g psii[0]=%g\n",t,sr,si,psir[0],psii[0] );

   expt =exp(sr*t);
   st=sin(si*t)*expt; ct=cos(si*t)*expt;
   // const real stp= si*ct+sr*st , ctp=-si*st+sr*ct;

   const real tm=t-dt;
   exptm =exp(sr*tm);
   stm=sin(si*tm)*exptm; ctm=cos(si*tm)*exptm;
   // const real stmp= si*ctm+sr*stm , ctmp=-si*stm+sr*ctm;

   const real sNormSq = sr*sr+si*si;

   ampH = ct;   
   // ampHp = -si*st + sr*ct;

   // eps Ev_t = curl( Hv ) - alphaP*eps* SUM (Pv_j).t 
   //   Pv_j = psi_j * Ev   
   // eps*( 1 + alphaP*Sum psi_j) \Ev_t = curl ( Hv ) 

   // E = Re( (1/(eps*s) * 1/( 1+alphaP*psi) * ( ct + i sint ) )
   //   = Re( (phir+i*phii)*( ct + i sint )
   const real alphaP = dmp.alphaP;

   real psirSum=0., psiiSum=0.;
   for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
   {
     psirSum += psir[iv]; 
     psiiSum += psii[iv];
   }
   
   real chiNormSq = SQR(1.+alphaP*psirSum)+SQR(alphaP*psiiSum); //   | 1+alphaP*psi|^2 

   //  phi = (1/(eps*s) * 1/( 1+alphaP*psi)
   //      = (sr-i*si)*( 1+alphaP*psir - i*alphaP*psii)/(eps* sNormSq*chiNormSq )
   //      = phir +i*phii 

   real phir = ( sr*(1.+alphaP*psirSum)-si*alphaP*psiiSum)/( eps*sNormSq*chiNormSq );
   real phii = (-si*(1.+alphaP*psirSum)-sr*alphaP*psiiSum)/( eps*sNormSq*chiNormSq );
   
   ampE = phir*ct - phii*st;

   // P = Re( (psir+i*psii)*(phir+i*phii)*( ct + i sint ) )
   //   = Re( (psir+i*psii)*( phir*ct-phii*st + i*( phir*st +phii*ct )
   //   = psir*( phir*ct-phii*st) -psii*(  phir*st +phii*ct )
   for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
   {
     ampP[iv] = psir[iv]*(phir*ct-phii*st ) - psii[iv]*( phir*st +phii*ct);
   }
   
   // tm = t-dt 
   ampHm = ctm;  
   // ampHp = -si*stm + sr*ctm;
   ampEm = phir*ctm - phii*stm;
   for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
   {
      ampPm[iv] = psir[iv]*(phir*ctm-phii*stm ) - psii[iv]*( phir*stm +phii*ctm);
   }
   
 }

#endMacro 


//==================================================================================================
// Evaluate the annulus eigenfunction or it's error
// 
// OPTION: solution, error
//==================================================================================================
#beginMacro annulusEigenFunction(OPTION,J1,J2,J3)

//   printF(" I1.getBase(),uLocal.getBase(0),I1.getBound(),uLocal.getBound(0)=%i %i %i %i \n",
// 	 I1.getBase(),uLocal.getBase(0),I1.getBound(),uLocal.getBound(0));



 if( numberOfDimensions==2 )
 {
#include "besselPrimeZeros.h"

   const int n = int(initialConditionParameters[0]+.5);  // angular number, n=0,1,... --> Jn(omega*r)
   const int m = int(initialConditionParameters[1]+.5);  // radial number m=0,... 
   
   assert( m<mdbpz && n<ndbpz );
   
   real omega = besselPrimeZeros[n][m];  // m'th zero of Jn' (excluding r=0 for J0)
   
// printF("Annulus: Bessel function solution: n=%i, m=%i, omega=%e (c=%8.2e)\n",n,m,omega,c);

   const real epsilon=sqrt(REAL_EPSILON);
   
   real np1Factorial=1.;
   for( int k=2; k<=n+1; k++ )
     np1Factorial*=k;              //  (n+1)!

   int i1,i2,i3;
   real r,gr,xd,yd,zd,bj,bjp,rx,ry,theta,thetax,thetay;
   real cosTheta,sinTheta,bjThetax,bjThetay,uex,uey,cosn,sinn;

   // real sint = sin(omega*t), cost = cos(omega*t);
   // real sintp = omega*cost, costp = -omega*sint;

   // real sintm = sin(omega*(t-dt)), costm = cos(omega*(t-dt));

   initializeAnnulusEigenfunction();

   FOR_3D(i1,i2,i3,J1,J2,J3)
   {
     xd=X(i1,i2,i3,0);
     yd=X(i1,i2,i3,1);
     #If #DIM == "3" 
       zd=X(i1,i2,i3,2)
       sinPhi=sin(Pi*zd)
       cosPhi=cos(Pi*zd)
     #End
     r = sqrt(xd*xd+yd*yd);
     theta=atan2(yd,xd);
     // if( theta<0. ) theta+=2.*Pi;
     cosTheta=cos(theta);
     sinTheta=sin(theta);
     
     cosn=cos(n*theta);
     sinn=sin(n*theta);
     
     gr=omega*r;
     
     rx = cosTheta;  // x/r
     ry = sinTheta;  // y/r
     
     bj=jn(n,gr);  // Bessel function J of order n
     
     
     if( gr>epsilon )  // need asymptotic expansion for small gr ??
     {
       bjp = -jn(n+1,gr) + n*bj/gr;  // from the recursion relation for Jn'
       thetay= cosTheta/r;
       thetax=-sinTheta/r;
       if( dispersionModel==noDispersion )
       {
         uex =  (1./omega)*(omega*ry*bjp*cosn -n*bj*thetay*sinn); // Ex.t = Hz.y
         uey = -(1./omega)*(omega*rx*bjp*cosn -n*bj*thetax*sinn); // Ey.t = - Hz.x
       }
       else
       {
         uex =  (omega*ry*bjp*cosn -n*bj*thetay*sinn); // Ex.t = Hz.y
         uey = -(omega*rx*bjp*cosn -n*bj*thetax*sinn); // Ey.t = - Hz.x
       }
       
     
     }
     else
     {
       // Jn(z) = (.5*z)^n *( 1 - (z*z/4)/(n+1)! + .. 
     
     
       // At r=0 all the Jn'(0) are zero except for n=1
       // bjp = n==1 ? 1./2. : 0.;
       bjp = n==0 ? 0. : pow(.5,double(n))*pow(gr,n-1.)*( 1. - (gr*gr)/(4.*np1Factorial) );
     
       // bj/r = omega*bjp at r=0
       bjThetay= omega*bjp*cosTheta;
       bjThetax=-omega*bjp*sinTheta;
       if( dispersionModel==noDispersion )
       {
         uex =  (1./omega)*(omega*ry*bjp*cosn -n*bjThetay*sinn);  // Ex.t = Hz.y
         uey = -(1./omega)*(omega*rx*bjp*cosn -n*bjThetax*sinn);  // Ey.t = - Hz.x
       }
       else
       {
         uex =  (omega*ry*bjp*cosn -n*bjThetay*sinn);  // Ex.t = Hz.y
         uey = -(omega*rx*bjp*cosn -n*bjThetax*sinn);  // Ey.t = - Hz.x
       }
       
       
     
     }

     #If #OPTION == "solution"
       UHZ(i1,i2,i3)  = bj*cosn*ampH;
       UEX(i1,i2,i3) = uex*ampE; 
       UEY(i1,i2,i3) = uey*ampE; 
     
       if( method==nfdtd )
       {
         UMHZ(i1,i2,i3) = bj*cosn*ampHm;
         UMEX(i1,i2,i3) = uex*ampEm;
         UMEY(i1,i2,i3) = uey*ampEm;
       }
       else if( method==sosup )
       {
         uLocal(i1,i2,i3,hzt) = bj*cosn*ampHp;
         uLocal(i1,i2,i3,ext) = uex*ampEp;
         uLocal(i1,i2,i3,eyt) = uey*ampEp;
       }
       
       if( dispersionModel!=noDispersion )
       { // -- dispersive ---

         for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
         {
           const int pc= iv*numberOfDimensions;
       
           // Do this for now -- set all vectors to be the same: 
           pLocal(i1,i2,i3,pc  ) = uex*ampP[iv];
           pLocal(i1,i2,i3,pc+1) = uey*ampP[iv];
           if( method==nfdtd )
           {
             pmLocal(i1,i2,i3,pc  ) = uex*ampPm[iv];
             pmLocal(i1,i2,i3,pc+1) = uey*ampPm[iv];
           }
         }
         
       }
       

     #Elif #OPTION == "boundaryCondition"
       // *check me*
       uLocal(i1,i2,i3,hz) = bj*cosn*ampH;
       uLocal(i1,i2,i3,ex) = uex*ampE; 
       uLocal(i1,i2,i3,ey) = uey*ampE; 
       if( dispersionModel!=noDispersion )
       { // -- dispersive ---
         for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
         {
           const int pc= iv*numberOfDimensions;
           // Do this for now -- set all vectors to be the same: 
           pLocal(i1,i2,i3,pc  ) = uex*ampP[iv];
           pLocal(i1,i2,i3,pc+1) = uey*ampP[iv];
         }
         // uLocal(i1,i2,i3,pxc) = uex*ampP;
         // uLocal(i1,i2,i3,pyc) = uey*ampP;
       }
     
       if( method==sosup )
       {
         uLocal(i1,i2,i3,hzt) = bj*cosn*ampHp;
         uLocal(i1,i2,i3,ext) = uex*ampEp;
         uLocal(i1,i2,i3,eyt) = uey*ampEp;
       }

     #Elif #OPTION == "error"
       ERRHZ(i1,i2,i3) = UHZ(i1,i2,i3) - bj*cosn*ampH;
       ERREX(i1,i2,i3) = UEX(i1,i2,i3) - uex*ampE;  // Ex.t = Hz.y
       ERREY(i1,i2,i3) = UEY(i1,i2,i3) - uey*ampE;  // Ey.t = - Hz.x
       if( dispersionModel!=noDispersion )
       { // -- dispersive ---
         for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
         {
           const int pc= iv*numberOfDimensions;
           // Do this for now -- set all vectors to be the same: 
           errPolarization(i1,i2,i3,pc  ) = pLocal(i1,i2,i3,pc  ) - uex*ampP[iv];
           errPolarization(i1,i2,i3,pc+1) = pLocal(i1,i2,i3,pc+1) - uey*ampP[iv];
         }
       }

       if( method==sosup )
       {
         errLocal(i1,i2,i3,hzt) = uLocal(i1,i2,i3,hzt) - bj*cosn*ampHp;
         errLocal(i1,i2,i3,ext) = uLocal(i1,i2,i3,ext) - uex*ampEp;
         errLocal(i1,i2,i3,eyt) = uLocal(i1,i2,i3,eyt) - uey*ampEp;
       }
     #Else
       Overture::abort("error");
     #End
  
  }
 }
 else /* 3D */
 {
#include "besselZeros.h"

   const real cylinderLength=cylinderAxisEnd-cylinderAxisStart;

   const int n = int(initialConditionParameters[0]+.5);  // angular number, n=0,1,... --> Jn(omega*r)
   const int m = int(initialConditionParameters[1]+.5);  // radial number m=0,... 
   const int k = int(initialConditionParameters[2]+.5);  // axial number k=1,2,3,...
   
   assert( m<mdbz && n<ndbz );
   
   real lambda = besselZeros[n][m];  // m'th zero of Jn (excluding r=0 for J0)
   real omega = sqrt( SQR(k*Pi/cylinderLength) + lambda*lambda );
   
   printF("***Cylinder: Bessel function soln: n=%i, m=%i, k=%i, lambda=%e, omega=%e (c=%8.2e) [za,zb]=[%4.2f,%4.2f]\n",
          n,m,k,lambda,omega,c,cylinderAxisStart,cylinderAxisEnd);

   const real epsilon=sqrt(REAL_EPSILON);
   
   real np1Factorial=1.;
   for( int k=2; k<=n+1; k++ )
     np1Factorial*=k;              //  (n+1)!

   int i1,i2,i3;
   real r,gr,xd,yd,zd,bj,bjp,rx,ry,theta,thetax,thetay;
   real cosTheta,sinTheta,bjThetax,bjThetay,uex,uey,cosn,sinn,sinkz,coskz;
   
   initializeAnnulusEigenfunction();

   FOR_3D(i1,i2,i3,J1,J2,J3)
   {
     xd=X(i1,i2,i3,0);
     yd=X(i1,i2,i3,1);
     zd=(X(i1,i2,i3,2)-cylinderAxisStart)/cylinderLength; // *wdh* 040626 -- allow for any length
     
     sinkz=sin(Pi*k*zd);   
     coskz=cos(Pi*k*zd); 

     r = sqrt(xd*xd+yd*yd);
     theta=atan2(yd,xd);

     cosTheta=cos(theta);
     sinTheta=sin(theta);
     
     cosn=cos(n*theta);
     sinn=sin(n*theta);

     //  cost=cos(omega*t);
     
     gr=lambda*r;
     
     rx = cosTheta;  // x/r
     ry = sinTheta;  // y/r
     
     bj=jn(n,gr);  // Bessel function J of order n
     
     
     if( gr>epsilon )  // need asymptotic expansion for small gr ??
     {
       bjp = -jn(n+1,gr) + n*bj/gr;  // from the recursion relation for Jn'
       thetay= cosTheta/r;
       thetax=-sinTheta/r;
     
       uex = -(k*Pi/(cylinderLength*lambda*lambda))*( lambda*rx*bjp*cosn - n*bj*thetax*sinn );
       uey = -(k*Pi/(cylinderLength*lambda*lambda))*( lambda*ry*bjp*cosn - n*bj*thetay*sinn );
     
     
     }
     else
     {
       // Jn(z) = (.5*z)^n *( 1 - (z*z/4)/(n+1)! + .. 
     
     
       // At r=0 all the Jn'(0) are zero except for n=1
       // bjp = n==1 ? 1./2. : 0.;
       bjp = n==0 ? 0. : pow(.5,double(n))*pow(gr,n-1.)*( 1. - (gr*gr)/(4.*np1Factorial) );
     
       // bj/r = lambda*bjp at r=0
       bjThetay= lambda*bjp*cosTheta;
       bjThetax=-lambda*bjp*sinTheta;
     
       uex = -(k*Pi/(cylinderLength*lambda*lambda))*( lambda*rx*bjp*cosn -n*bjThetax*sinn);  // Ex.t = Hz.y
       uey = -(k*Pi/(cylinderLength*lambda*lambda))*( lambda*ry*bjp*cosn -n*bjThetay*sinn);  // Ey.t = - Hz.x
     
     }

     uex = uex*sinkz;
     uey = uey*sinkz;
     real uez = bj*cosn*coskz;

     #If #OPTION == "solution"

       UEX(i1,i2,i3) = uex*ampE;
       UEY(i1,i2,i3) = uey*ampE;
       UEZ(i1,i2,i3) = uez*ampE;
       // UEX(i1,i2,i3) = uex*sinkz*cost;
       // UEY(i1,i2,i3) = uey*sinkz*cost;
       // UEZ(i1,i2,i3) = bj*cosn*coskz*cost;
     
       UMEX(i1,i2,i3) = uex*ampEm;
       UMEY(i1,i2,i3) = uey*ampEm;
       UMEZ(i1,i2,i3) = uez*ampEm;
       // cost=cos(omega*(t-dt)); 
       // UMEX(i1,i2,i3) = uex*sinkz*cost;
       // UMEY(i1,i2,i3) = uey*sinkz*cost;
       // UMEZ(i1,i2,i3) = bj*cosn*coskz*cost;
     
       if( method==sosup )
       {
         assert( dispersionModel==noDispersion );
         
         sint=sin(omega*t); 
         uLocal(i1,i2,i3,ext) = -omega*uex*sint;
         uLocal(i1,i2,i3,eyt) = -omega*uey*sint;
         uLocal(i1,i2,i3,ezt) = -omega*uez*sint;
       }

       if( dispersionModel!=noDispersion )
       { // -- dispersive ---

         for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
         {
           const int pc= iv*numberOfDimensions;
       
           // Do this for now -- set all vectors to be the same: 
           pLocal(i1,i2,i3,pc  ) = uex*ampP[iv];
           pLocal(i1,i2,i3,pc+1) = uey*ampP[iv];
           pLocal(i1,i2,i3,pc+2) = uez*ampP[iv];
           if( method==nfdtd )
           {
             pmLocal(i1,i2,i3,pc  ) = uex*ampPm[iv];
             pmLocal(i1,i2,i3,pc+1) = uey*ampPm[iv];
             pmLocal(i1,i2,i3,pc+2) = uez*ampPm[iv];
           }
         }
       }

     #Elif #OPTION == "boundaryCondition"

       // *check me*
       uLocal(i1,i2,i3,ex) = uex*ampE;
       uLocal(i1,i2,i3,ey) = uey*ampE;
       uLocal(i1,i2,i3,ez) = uez*ampE;
       // uLocal(i1,i2,i3,ex) = uex*sinkz*cost;
       // uLocal(i1,i2,i3,ey) = uey*sinkz*cost;
       // uLocal(i1,i2,i3,ez) = bj*cosn*coskz*cost;
     
       if( method==sosup )
       {
         assert( dispersionModel==noDispersion );
         sint=sin(omega*t); 
         uLocal(i1,i2,i3,ext) = -omega*uex*sint;
         uLocal(i1,i2,i3,eyt) = -omega*uey*sint;
         uLocal(i1,i2,i3,ezt) = -omega*uez*sint;
       }

       if( dispersionModel!=noDispersion )
       { // -- dispersive ---

         for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
         {
           const int pc= iv*numberOfDimensions;
       
           // Do this for now -- set all vectors to be the same: 
           pLocal(i1,i2,i3,pc  ) = uex*ampP[iv];
           pLocal(i1,i2,i3,pc+1) = uey*ampP[iv];
           pLocal(i1,i2,i3,pc+2) = uez*ampP[iv];
         }
       }
     #Elif #OPTION == "error"


       ERREX(i1,i2,i3) = UEX(i1,i2,i3) - uex*ampE;
       ERREY(i1,i2,i3) = UEY(i1,i2,i3) - uey*ampE;
       ERREZ(i1,i2,i3) = UEZ(i1,i2,i3) - uez*ampE;

       // ERREX(i1,i2,i3) = UEX(i1,i2,i3) - uex*sinkz*cost;
       // ERREY(i1,i2,i3) = UEY(i1,i2,i3) - uey*sinkz*cost;
       // ERREZ(i1,i2,i3) = UEZ(i1,i2,i3) - bj*cosn*coskz*cost;

       if( dispersionModel!=noDispersion )
       { // -- dispersive ---
         for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
         {
           const int pc= iv*numberOfDimensions;
           // Do this for now -- set all vectors to be the same: 
           errPolarization(i1,i2,i3,pc  ) = pLocal(i1,i2,i3,pc  ) - uex*ampP[iv];
           errPolarization(i1,i2,i3,pc+1) = pLocal(i1,i2,i3,pc+1) - uey*ampP[iv];
           errPolarization(i1,i2,i3,pc+2) = pLocal(i1,i2,i3,pc+2) - uez*ampP[iv];
         }
       }
       if( method==sosup )
       {
         assert( dispersionModel==noDispersion );
         sint=sin(omega*t); 
         errLocal(i1,i2,i3,ext) = uLocal(i1,i2,i3,ext) + omega*uex*sint;
         errLocal(i1,i2,i3,eyt) = uLocal(i1,i2,i3,eyt) + omega*uey*sint;
         errLocal(i1,i2,i3,ezt) = uLocal(i1,i2,i3,ezt) + omega*uez*sint;
       }
     #Else
       Overture::abort("error");
     #End
  
  }
 }

#endMacro
