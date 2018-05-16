%
% Determine the eigenvalues and eigenvectors of the 4th order approximation
% to the 1D eigenvalue problem ** NEUMANN BC's ***
%
%           - u_xx = lambda u   x in (0,1)
%            u_x=0 at x=0,1
%         
%
%     u(1)                  u(n)
%      1  2                 n
%      +--+--+--+--+--+--+--+                
%     x=xa               xc xd              
%

clear;
clf;
set(gca,'FontSize',16);

n = 81;  % number of grid points
n= 41;

xa=0.;
xb=1.;

h=(xb-xa)/(n-1);
hsq=h*h;
hsq12=h*h*12.;

% Include 2 ghost points on each side
ng=n+4;
ia=3; 
ib=ng-2;
x=zeros(ng,1);
for i=1:ng 
  x(i)=xa+(i-ia)*h;
end

%  bcOption=0 : use BC D0 D+D-u(0) = 0 
%  bcOption=1 : use BC D+^5 u(-2) = 0 
bcOption=1;  
if bcOption==0 
  identifier='eqn';
else
  identifier='extrap';
end

adim=n;
a=zeros(adim,adim);

i=1;
% a(i,i-2)=  1./hsq12;  % u(-1) = 2*u(0)-u(1)
% a(i,i-1)=-16./hsq12;  % u(0)=0
a(i,i  )= 30/hsq12;
a(i,i+1)=-16./hsq12;
a(i,i+2)=  1./hsq12;
if bcOption==0 
 a(i,i+2)=a(i,i+2)+1./hsq12; % BC's in this case are equiv to u(-2)=u(2) u(-1)=u(1)
 a(i,i+1)=a(i,i+1)-16./hsq12;
elseif bcOption==1
 % u(i-2) = -80/3u0 + 40u1 -15 u2 + 8/3 u3
 % u(i-1) = -10/3u0 +  6u1  -2 u2+(1/3)*u3
 a(i,i  )=a(i,i  )+1./hsq12*(-80/3.) -16./hsq12*(-10./3.);
 a(i,i+1)=a(i,i+1)+1./hsq12*( 40.  ) -16./hsq12*(  6.   );
 a(i,i+2)=a(i,i+2)+1./hsq12*(-15.  ) -16./hsq12*( -2.   );
 a(i,i+3)=        +1./hsq12*( 8./3.) -16./hsq12*( 1./3. );
end

i=2;
% a(i,i-2)=  1./hsq12;  % u(0)=0
a(i,i-1)=-16./hsq12;
a(i,i  )= 30./hsq12;
a(i,i+1)=-16./hsq12;
a(i,i+2)=  1./hsq12;
if bcOption==0 
 a(i,i)=a(i,i)+1./hsq12; % BC's in this case are equiv to u(-2)=u(2) u(-1)=u(1)
elseif bcOption==1
 a(i,i-1)=a(i,i-1)+1./hsq12*(-10./3.);
 a(i,i  )=a(i,i  )+1./hsq12*(  6.   );
 a(i,i+1)=a(i,i+1)+1./hsq12*( -2.   );
 a(i,i+2)=a(i,i+2)+1./hsq12*( 1./3. );
end

for i=3:n-2
  a(i,i-2)=  1./hsq12;
  a(i,i-1)=-16./hsq12;
  a(i,i  )= 30./hsq12;
  a(i,i+1)=-16./hsq12;
  a(i,i+2)=  1./hsq12;
end 

i=n;
a(i,i-2)=  1./hsq12; 
a(i,i-1)=-16./hsq12; 
a(i,i  )= 30./hsq12;
% a(i,i+1)=-16./hsq12;
% a(i,i+2)=  1./hsq12;
if bcOption==0 
 a(i,i-2)=a(i,i-2)+1./hsq12;
 a(i,i-1)=a(i,i-1)-16./hsq12;
elseif bcOption==1
 % u(i-2) = -80/3u0 + 40u1 -15 u2 + 8/3 u3
 % u(i-1) = -10/3u0 +  6u1  -2 u2+(1/3)*u3
 a(i,i  )=a(i,i  )+1./hsq12*(-80/3.) -16./hsq12*(-10./3.);
 a(i,i-1)=a(i,i-1)+1./hsq12*(40.   ) -16./hsq12*(  6.   );
 a(i,i-2)=a(i,i-2)+1./hsq12*(-15.  ) -16./hsq12*( -2.   );
 a(i,i-3)=        +1./hsq12*( 8./3.) -16./hsq12*( 1./3. );
end

i=n-1;
a(i,i-2)=  1./hsq12;  % u(0)=0
a(i,i-1)=-16./hsq12;
a(i,i  )= 30./hsq12;
a(i,i+1)=-16./hsq12;
% a(i,i+2)=  1./hsq12;
if bcOption==0 
 a(i,i)=a(i,i)+1./hsq12;
elseif bcOption==1
 a(i,i+1)=a(i,i+1)+1./hsq12*(-10./3.);
 a(i,i  )=a(i,i  )+1./hsq12*(  6.   );
 a(i,i-1)=a(i,i-1)+1./hsq12*( -2.   );
 a(i,i-2)=a(i,i-2)+1./hsq12*( 1./3. );
end

%  a
%  pause

[v,L] = eig(a);

% Now plot the results

for i=1:adim
  lambda(i)=L(i,i);
end

% sort from smallest to largest in absolute value
[zz,iOrder] = sort(abs(lambda));

e=zeros(ng,n);
u=zeros(ng,n);
uTrue=zeros(ng,1);


% ------  Build eigenfunctions for plotting -------------
for j=1:n  % loop over eigenfunctions
  ie=iOrder(j);

  lambdaTrue= (pi*(j-1))^2;

  relErr = abs( lambda(ie)-lambdaTrue )/max(1,lambdaTrue);

  fprintf(' j=%i ie=%i lambda(ie)=%8.2e  true=%8.2e  relErr=%8.2e err/(m^6 h^4)=%8.2e \n',...
                j,ie,lambda(ie),lambdaTrue,relErr,abs(lambda(ie)-lambdaTrue )/max(h^4,((j-1)^6*h^4)));

 for i=1:n
   u(i+ia-1,j)=v(i,ie);
 end 
 % end points from BC's
 if bcOption==0
   u(ia-2,j)=u(ia+2,j);
   u(ia-1,j)=u(ia+1,j);
   u(ib+1,j)=u(ib-1,j);
   u(ib+2,j)=u(ib-2,j);
 else
   i=ia-1;
   u(i,j) = -(10./3.)*u(i+1,j) + 6.*u(i+2,j) -2.*u(i+3,j)+(1./3.)*u(i+4,j);
   i=1;
   u(i,j) = -(80./3.)*u(i+2,j) +40.*u(i+3,j)-15.*u(i+4,j)+(8./3.)*u(i+5,j);
   i=ng-1;
   u(i,j) = -(10./3.)*u(i-1,j) + 6.*u(i-2,j) -2.*u(i-3,j)+(1./3.)*u(i-4,j);
   i=ng;
   u(i,j) = -(80./3.)*u(i-2,j) +40.*u(i-3,j)-15.*u(i-4,j)+(8./3.)*u(i-5,j);
 end 
 % check:
 i=ia;
 resa= (u(i-2,j)-16.*u(i-1,j)+30.*u(i,j)-16.*u(i+1,j)+u(i+2,j))/(12.*h*h)-lambda(ie)*u(i,j);
 i=ib;
 resb= (u(i-2,j)-16.*u(i-1,j)+30.*u(i,j)-16.*u(i+1,j)+u(i+2,j))/(12.*h*h)-lambda(ie)*u(i,j);
 fprintf('LU-lambda*U at end points: %8.2e %8.2e \n',resa,resb);
 i=ia;
 resa= (u(i-2,j)-8.*u(i-1,j)+8.*u(i+1,j)-u(i+2,j))/(12.*h);
 i=ib;
 resb= (u(i-2,j)-8.*u(i-1,j)+8.*u(i+1,j)-u(i+2,j))/(12.*h);
 fprintf('Residual in Neumman BC at end points: %8.2e %8.2e \n',resa,resb);
 if bcOption==0
   i=ia;
   resa= -u(i-2,j)+2.*u(i-1,j)-2.*u(i+1,j)+u(i+2,j);
   i=ib;
   resb= -u(i-2,j)+2.*u(i-1,j)-2.*u(i+1,j)+u(i+2,j);
   fprintf('Residual in u_xxx BC at end points: %8.2e %8.2e \n',resa,resb);
 elseif bcOption==1
   i=ia;
   resa= u(i-2,j)-5.*u(i-1,j)+10.*u(i,j)-10.*u(i+1,j)+5.*u(i+2,j)-u(i+3,j);
   i=ib-1;
   resb= u(i-2,j)-5.*u(i-1,j)+10.*u(i,j)-10.*u(i+1,j)+5.*u(i+2,j)-u(i+3,j);
   fprintf('Residual in Extrap BC at end points: %8.2e %8.2e \n',resa,resb);
 end

 % arrange the eigenfunctions so they always are initially positive (easier then to compare)
 if u(2,j) < 0 
   u(:,j)=u(:,j)*(-1);
 end 

 uTrue=cos(pi*(j-1)*x);

 alpha=u(:,j).'*uTrue/(u(:,j).'*u(:,j));

 u(:,j)=u(:,j)*alpha; % scale to match uTrue: min || uTrue - alpha*v || 

 e(:,j)=u(:,j)-uTrue;  % error in the eigenfunction

 u(:,j)=u(:,j)*.5+j;


% plot(xu,vu(:,j),'r-o',xv,vv(:,j),'b-x');
% title(sprintf('Eigenfunction %i, \\lambda=%8.2e, m=%i, n=%i, d=%3.2fh_1=%3.2fh_2',j,lambda(ie),m,n,d/hu,d/hv));
% pause

end


k=0;
md=10;
md=5;
for j=1:md:n
  k=k+1;
  plot(x,u(:,j  ),'r-o', ...
       x,u(:,j+1),'g-o', ...
       x,u(:,j+2),'m-o', ...
       x,u(:,j+3),'b-o', ...
       x,u(:,j+4),'c-o');
  hold on;
  if md==10 
   plot(x,u(:,j+5),'g-o', ...
        x,u(:,j+6),'m-o', ...
        x,u(:,j+7),'b-o', ...
        x,u(:,j+8),'r-o', ...
        x,u(:,j+9),'g-o');

  end

  yMin = min(u(:,j     ));
  yMax = max(u(:,j+md-1));
  delta=yMax-yMin;
  if bcOption==0 
     title(sprintf('Neumann BC=D_0D_+D_-'));
  else
     title(sprintf('Neumann Problem, BC=extrapolation'));
  end
  text(.5,yMax-.01*delta,sprintf('Eigenfunctions %i-%i, n=%i',j,j+md-1,n),'Fontsize',16,'HorizontalAlignment','center' );

  grid on;
  set(gca,'XLim',[x(1)-2*h x(ng)+2*h]);

  hold off;

  name = sprintf('bceign-%i-%s-%i.eps',n,identifier,k);
  fprintf('Save %s\n',name);
  print('-depsc2',name);
  pause

  plot(x,e(:,j  ),'r-+', ...
       x,e(:,j+1),'g-+', ...
       x,e(:,j+2),'m-+', ...
       x,e(:,j+3),'b-+', ...
       x,e(:,j+4),'c-+');

  grid on;
  set(gca,'XLim',[x(1)-2*h x(ng)+2*h]);
  yMin = min(e(ia:ib,j+4));
  yMax = max(e(ia:ib,j+4));
  delta=yMax-yMin;
  yMin=yMin-delta*.2;
  yMax=yMax+delta*.2;

  set(gca,'YLim',[yMin yMax ]);

  if bcOption==0 
    title(sprintf('Neumann BC=D_0D_+D_-'));
  else
    title(sprintf('Neumann Problem, BC=extrapolation'));
  end
  text(.5,yMax-.1*delta,sprintf('Errors in eigenfunctions %i-%i, n=%i',j,j+md-1,n),'Fontsize',16,'HorizontalAlignment','center');

  name = sprintf('bceign-err-%i-%s-%i.eps',n,identifier,k);
  fprintf('Save %s\n',name);
  print('-depsc2',name);
  
  pause
end


% print -deps2 eigenfunction.eps
% print -depsc2 eigenfunction.eps

% print -depsc2 bceig.extrap80.1-5.eps
% print -depsc2 bceig.extrap80.6-10.eps

% print -depsc2 bceig.extrap160.1-5.eps
% print -depsc2 bceig.extrap160.6-10.eps
