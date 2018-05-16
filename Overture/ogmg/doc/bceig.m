%
% Determine the eigenvalues and eigenvectors of the 4th order approximation
% to the 1D eigenvalue problem 
%
%           - u_xx = lambda u   x in (0,1)
%            u=0 at x=0,1
%         
%
%     u(0)                  u(n)
%      0  1  2              n
%      +--+--+--+--+--+--+--+                
%     x=xa               xc xd              
%

clear;
clf;
set(gca,'FontSize',16);

n = 11;  % number of grid points
n = 41;  % 320; % 160; 
m = n-1; % active points are 1..m

xa=0.;
xb=1.;

h=(xb-xa)/n;
hsq=h*h;
hsq12=h*h*12.;

ng = n+3;  % include 1 ghost point on either side
x=zeros(ng,1);
ia=2;
ib=ng-1;
for i=1:ng
  x(i)=xa+(i-ia)*h;
end

%  bcOption=0 : use BC D+D-u(0) = 0 
%  bcOption=1 : use BC D+^p u(-1) = 0 
bcOption=1;  
if bcOption==0 
  identifier='eqn';
else
  identifier='extrap';
end

adim=m;
a=zeros(adim,adim);

i=1;
% a(i,i-2)=  1./hsq12;  % u(-1) = 2*u(0)-u(1)
% a(i,i-1)=-16./hsq12;  % u(0)=0
a(i,i  )= 30/hsq12;
a(i,i+1)=-16./hsq12;
a(i,i+2)=  1./hsq12;
if bcOption==0 
  a(i,i)=a(i,i)-1./hsq12;   % u(-1) = 2*u(0)-u(1)
elseif bcOption==1
  % u(-1) = 4*u(0) -6 u(1) +4*u(2) -u(3)
  a(i,i  )=a(i,i  )-6./hsq12;
  a(i,i+1)=a(i,i+1)+4./hsq12;
  a(i,i+2)=a(i,i+2)-1./hsq12;
end

i=2;
% a(i,i-2)=  1./hsq12;  % u(0)=0
a(i,i-1)=-16./hsq12;
a(i,i  )= 30./hsq12;
a(i,i+1)=-16./hsq12;
a(i,i+2)=  1./hsq12;

for i=3:m-2
  a(i,i-2)=  1./hsq12;
  a(i,i-1)=-16./hsq12;
  a(i,i  )= 30./hsq12;
  a(i,i+1)=-16./hsq12;
  a(i,i+2)=  1./hsq12;
end 

i=m;
a(i,i-2)=  1./hsq12; 
a(i,i-1)=-16./hsq12; 
a(i,i  )= 30./hsq12;
% a(i,i+1)=-16./hsq12;
% a(i,i+2)=  1./hsq12;
if bcOption==0 
  a(i,i)=a(i,i)-1./hsq12;   % u(-1) = 2*u(0)-u(1)
elseif bcOption==1
  % u(-1) = 4*u(0) -6 u(1) +4*u(2) -u(3)
  a(i,i  )=a(i,i  )-6./hsq12;
  a(i,i-1)=a(i,i-1)+4./hsq12;
  a(i,i-2)=a(i,i-2)-1./hsq12;
end

i=m-1;
a(i,i-2)=  1./hsq12;  % u(0)=0
a(i,i-1)=-16./hsq12;
a(i,i  )= 30./hsq12;
a(i,i+1)=-16./hsq12;
% a(i,i+2)=  1./hsq12;

%  a
%  pause

[v,L] = eig(a);

% Now plot the results

for i=1:adim
  lambda(i)=L(i,i);
end

% sort from smallest to largest in absolute value
[zz,iOrder] = sort(abs(lambda));

e=zeros(ng,m);
u=zeros(ng,m);
uTrue=zeros(ng,1);


% ------  Build eigenfunctions for plotting -------------
for j=1:m  % loop over eigenfunctions
  ie=iOrder(j);

  lambdaTrue= (pi*j)*(pi*j);

  relErr = abs( lambda(ie)-lambdaTrue )/lambdaTrue;

  fprintf(' j=%i ie=%i lambda(ie)=%8.2e  true=%8.2e  relErr=%8.2e err/(m^6 h^4)=%8.2e \n',...
                j,ie,lambda(ie),lambdaTrue,relErr,abs(lambda(ie)-lambdaTrue )/(j^6*h^4));

 u(ia,j)=0;
 for i=1:m  
   u(i+ia,j)=v(i,ie);
 end 
 u(ib,j)=0.;
 if bcOption==0 
   i=1;
   u(i,j)=2*u(i+1,j)-u(i+2,j);
   i=ng;
   u(i,j)=2.*u(i-1,j)-u(i-2,j);
 else
   i=1;
   u(i,j)=4*u(i+1,j)-6.*u(i+2,j)+4.*u(i+3,j)-u(i+4,j);
   i=ng;
   u(i,j)=4*u(i-1,j)-6.*u(i-2,j)+4.*u(i-3,j)-u(i-4,j);
 end

 % v(:,ie).'
 % u(:,j).'

 % arrange the eigenfunctions so they always are initially positive (easier then to compare)
 if u(ia+1,j) < 0 
   u(:,j)=u(:,j)*(-1);
 end 

 uTrue=sin(pi*j*x);

 alpha=u(:,j).'*uTrue/(u(:,j).'*u(:,j));

 u(:,j)=u(:,j)*alpha; % scale to match uTrue: min || uTrue - alpha*v || 

 e(:,j)=u(:,j)-uTrue;  % error in the eigenfunction

 u(:,j)=u(:,j)*.45+j;


% plot(xu,vu(:,j),'r-o',xv,vv(:,j),'b-x');
% title(sprintf('Eigenfunction %i, \\lambda=%8.2e, m=%i, n=%i, d=%3.2fh_1=%3.2fh_2',j,lambda(ie),m,n,d/hu,d/hv));
% pause

end


k=0;
md=10;
md=5;
for j=1:md:m
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
    title(sprintf('Dirichlet Problem BC=D_+D_-'));
  else
    title(sprintf('Dirichlet Problem BC=extrapolation',j,j+md-1,n));
  end
  text(.5,yMax-.05*delta,sprintf('Eigenfunctions %i-%i, n=%i',j,j+md-1,n),'Fontsize',16,'HorizontalAlignment','center');

  grid on;
  set(gca,'XLim',[x(1)-2*h x(ng)+2*h]);
  hold off;

  name = sprintf('bceig-%i-%s-%i.eps',n,identifier,k);
  fprintf('Save %s\n',name);
  print('-depsc2',name);

  pause

  plot(x,e(:,j  ),'r-o', ...
       x,e(:,j+1),'g-o', ...
       x,e(:,j+2),'m-o', ...
       x,e(:,j+3),'b-o', ...
       x,e(:,j+4),'c-o');

  yMin = min(e(ia:ib,j+4));
  yMax = max(e(ia:ib,j+4));
  delta=yMax-yMin;
  yMin=yMin-delta*.2;
  yMax=yMax+delta*.2;

  if bcOption==0 
    title(sprintf('Dirichlet Problem BC=D_+D_-'));
  else
    title(sprintf('Dirichlet Problem BC=extrapolation',j,j+md-1,n));
  end
  text(.5,yMax-.1*delta,sprintf('Errors in eigenfunctions %i-%i, n=%i',j,j+md-1,n),'Fontsize',16,'HorizontalAlignment','center');

  grid on;
  set(gca,'XLim',[x(1)-2*h x(ng)+2*h]);


  set(gca,'YLim',[yMin yMax ]);
  name = sprintf('bceig-err-%i-%s-%i.eps',n,identifier,k);
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
