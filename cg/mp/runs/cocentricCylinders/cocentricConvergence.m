%
% Plot errors and convergence rate for the cocentric cylinder computations
% 
clear;
clf;
set(gca,'FontSize',18);

n=3;
h=zeros(n,1); eu=zeros(n,1); eT=zeros(n,1);
h(1)=.1; h(2)=1./20.; h(3)=1./40.;

eu(1)=1.3e-3; eu(2)=3.9e-4; eu(3)=9.4e-5;

eT(1)=6.0e-3; eT(2)=1.8e-3; eT(3)=4.5e-4; 

hl=zeros(2,1); el=zeros(2,1);
hl(1)=1./15.; hl(2)=hl(1)/1.95; 
el(1)=1.2e-3; el(2)=el(1)*(hl(2)/hl(1))^2; 


loglog(h,eT,'b-x',h,eu,'r-+','MarkerSize',15);
title(sprintf('Cocentric Cylinders Heat Exchanger'));
legend('T error','v error','Location','NorthWest');
xlabel('\Delta s');

hold on;

loglog(hl,el,'k-','LineWidth',1);

xt=1./25.; yt=.85e-3; 
text(xt,yt,'Slope 2','FontSize',18);

set(gca,'XLim',[1./50.,.15]);
set(gca,'YLim',[8.0e-5,7.0e-3]);
grid on;

print('-depsc2','cocentricCylindersLogErrors.eps');
