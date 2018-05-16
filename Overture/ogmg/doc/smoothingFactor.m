function [lambdaMax,lambdaNum,z1Max,z2Max,l1Max,l2Max,l3Max,l1z1Max,l1z2Max]=smoothingFactor(omega,cMin,nx,ny,order,nu)
%
%
% This version uses a brute force approach to find the max eigenvalue
% of the omega-RB-GS smoothing operator 
%
% Input:
%  w = omega
%  cMin
%
% Output:
%   lambdaMax : max(|lambda1|,|lambda2|,|lambda3|)
%  lambdaNum: 1,2 or 3 gives the lambda that achieved lambdaMax
%  z1Max,z2Max : position where the max occured.

xa=0.;  xb=1.;
ya=-1.; yb=0.;
hx=(xb-xa)/(nx-1);
hy=(yb-ya)/(ny-1);

if order==4 
  c41=16/15;
  c42=1/15;
else
  c41=1.;
  c42=0.;
end

c1=cMin;
cmax=1-c1;
c2=1-c1;

w=omega;

lambdaMax=0.;
l1Max=0;
l2Max=0;
l3Max=0;

z1Max=0;
z2Max=0;
lambdaNum=1;
l1z1Max=0;
l1z2Max=0;

for j=1:ny
for i=1:nx

  z1=xa+(i-1)*hx;
  z2=ya+(j-1)*hy;

  % lambda is maximized for z1 and z2s in [0,1]
  z2s=-z2;
  alpha   =c1*(  c41*z1 - 2*c42*z1^2 )+c2*(  c41*z2s - 2*c42*z2s^2 ) + c42;
  alphaBar=c1*( -c41*z1 - 2*c42*z1^2 )+c2*( -c41*z2s - 2*c42*z2s^2 ) + c42;  % z1-> -z1,  z2 -> -z2 
  s =1-w*(1-alpha);
  sb=1-w*(1-alphaBar);
  s22=( (sb+1)^2+(sb-1)*(1-s) )/4;
  if nu==1 
    lambda1=abs( s22 );
  else
    s11=( (s +1)^2+(sb-1)*(1-s) )/4;
    s12=( (sb+1)*(1-s )+s^2 -1  )/4;
    s21=( (s +1)*(1-sb)+sb^2-1  )/4;

    if nu==2
      lambda1=sqrt(abs( s12*s21+s22^2 ));
    else
      lambda1=(abs( s11*s12*s21+2*s12*s21*s22+s22^3 ))^(1/3);
    end
  end
%  lambda1=abs( 1 - w*(1-alphaBar)*( 1-w*(alpha-alphaBar)/4 ));


  alpha   =c1*(  c41*z1 - 2*c42*z1^2 )+c2*(  c41*z2 - 2*c42*z2^2 ) + c42;
  alphaBar=c1*( -c41*z1 - 2*c42*z1^2 )+c2*( -c41*z2 - 2*c42*z2^2 ) + c42;  % z1-> -z1,  z2 -> -z2 

  s =1-w*(1-alpha);
  sb=1-w*(1-alphaBar);

  s11=( (s +1)^2+(sb-1)*(1-s) )/4;
  s22=( (sb+1)^2+(sb-1)*(1-s) )/4;
  s12=( (sb+1)*(1-s )+s^2 -1  )/4;
  s21=( (s +1)*(1-sb)+sb^2-1  )/4;

  ssum=(s11+s22);
  sdiff=(s11-s22);
  sprod=(s12*s21);
  des=sqrt(sdiff^2+4*sprod);

  lambda2=abs( (1/2)*( ssum + des ) );

  lambda3=abs( (1/2)*( ssum - des ) );

  lm=max(lambda1,max(lambda2,lambda3));
  
  if lm>lambdaMax 
    lambdaMax=lm;
    z1Max=z1;
    z2Max=z2;
    if lambda1>=lambdaMax
      lambdaNum=1;
      z2Max=z2s;
    elseif lambda2>lambda3 
      lambdaNum=2;
    else
     lambdaNum=3;
    end
  end

  % keep track of the max of each lamba
  if lambda1>l1Max
    l1Max=lambda1;
    l1z1Max=z1;
    l1z2Max=z2s;
  end
  if lambda2>l2Max
    l2Max=lambda2;
  end
  if lambda3>l3Max
    l3Max=lambda3;
  end


end
end
