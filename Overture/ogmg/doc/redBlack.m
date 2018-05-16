% determine optimal relaxation parameters for GS-RB
%
set(gca,'FontSize',16);


nx=501; % 51; % 501; % 501
ny=501; % 51; % 501;
xa=1.;  xb=2.;
ya=0;   yb=1.;
hx=(xb-xa)/(nx-1);
hy=(yb-ya)/(ny-1);


order=4;  % 2 or 4



z=zeros(nx,ny);

[omega,c] = meshgrid(xa:hx:xb,ya:hy:yb);

for j=1:ny
for i=1:nx

  x=xa+(i-1)*hx;
  y=ya+(j-1)*hy;

%   z(i,j) = sin(3.14*(x+y));
% omega=x
% C=y
t0= x*x*y*y - 2.*(x-1);
t1=x*y*sqrt( x*x*y*y-4.*(x-1) );

z(i,j) = abs( .5*( t0+t1 ) );
z(i,j) = max( z(i,j), abs( .5*( t0-t1 ) ));

%  -- more general way

if order==4
  s  = 1-x*(1-(16/15)*y+(1/15)*(2*y*y-1));
  sb = 1-x*(1+(16/15)*y+(1/15)*(2*y*y-1));
else
  s  = 1-x*(1-y);
  sb = 1-x*(1+y);
end

s11=( (s+1)*(s+1) + (sb-1)*(1-s) )/4;
s12=( (s+1.)*(1-sb)+sb*sb-1 )/4;
s21=( (sb+1)*(1-s)+s*s-1 )/4;
s22=( (sb+1)*(sb+1)+(sb-1)*(1-s) )/4;

t0 = s11+s22;
t1 = sqrt( t0*t0-4*(s11*s22-s12*s21) );

z(i,j) = abs( .5*( t0+t1 ) );
z(i,j) = max( z(i,j), abs( .5*( t0-t1 ) ));

% muSmooth= .125*( x*x+4.*x-4. );
muSmooth=abs(s22);

z(i,j) = max( z(i,j), muSmooth );


end
end

% surf(z);
mesh(omega,c,z.');
zlabel('\mu');
ylabel('C');
xlabel('omega');

% print -depsc2 redBlackSurf.order2.eps
% print -depsc2 redBlackSurf.order4.eps

pause

% Now find the max(z) as a function of omega

zMax=zeros(nx,1);
xx=zeros(nx,1);




nz=51;
hz=.5/(nz-1);
cc=0:hz:.5;
omegaMin=zeros(nz,1);
muMin=zeros(nz,1);
mu1=zeros(nz,1);  % mu(omega=1)

for k=1:nz

cmin=cc(k);
cmax=1.-cmin;

muMin(k)=1;

for i=1:nx

  x=xa+(i-1)*hx;
  xx(i)=x;

  zm=0.;
  for j=1:ny
    y=ya+(j-1)*hy;
    if y<= cmax 
      zm=max(zm,z(i,j));
    end
  end
  if zm< muMin(k)
    omegaMin(k)=x;
    muMin(k)=zm;
  end
  zMax(i)=zm;
  if i==1
    mu1(k)=muMin(k);
  end 
end

% fprintf('k=%i cMin=%4.2f omegaMin=%5.3f muMin=%5.3f\n',k,cmin,omegaMin(k),muMin(k));

if mod(k,10)==0 
 plot(xx,zMax,'r');
 title(sprintf(' cMin=%4.2f \\omega_{min}=%5.3f \\mu_{min}=%5.3f order=%i',cmin,omegaMin(k),muMin(k)),order);
 ylabel('\mu');
 xlabel('\omega');
 grid on
 pause
% print -depsc2 redBlack.cut.order2.eps
% print -depsc2 redBlack.cut.order4.eps
end


end % for k

plot(cc,omegaMin,'r-o',cc,muMin,'b-x',cc,mu1,'g-+');
title(sprintf('Red-Black GS with relaxation order=%i',order));
xlabel('cMin');
legend('omegaMin','muMin','mu(1)');
grid on
% grid minor

% print -depsc2 redBlack.order2.eps
% print -depsc2 redBlack.order4.eps






