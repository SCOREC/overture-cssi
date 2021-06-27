%
% Compute errors in the cocentric cylinder computations
% 
clear;
clf;
set(gca,'FontSize',18);

% icase=1;
icase=2;
% icase=4;

% -- read data: 
if icase == 1 
 name = 'cocentric1';
 cocentric1Fluid
elseif icase == 2 
 name = 'cocentric2';
 cocentric2fFluid
elseif icase == 4
 name = 'cocentric4';
 cocentric4aFluid
elseif icase == 8
 name = 'cocentric8';
 cocentric8Fluid
end; 

rf2 = x00;
uf2 = v0;
Tf2 = T0;

if icase == 1
 cocentric1Solid
elseif icase ==2 
 cocentric2fSolid
elseif icase == 4
 cocentric4aSolid
elseif icase == 8
 cocentric8Solid
end 
rs2 = x00;
Ts2 = T0;


% -------- compute the exact solution

ra=.5; ri=1.; rb=1.5; ya=0.; yb=1.;  

beta = -1.; % thermal expansivity times gravity
nu=.05; 
kf=.2;  % k for solid
ks=1.;  % k for fluid
Ta = 5.; % Ts2(1);
Tb = 1.; % Tf2(nbf);

% (Ti-Ta)*ks/.5 = (Tb-Ti)*kf/.5 

Ti =  (kf*log(ri/ra)*Tb+ks*log(rb/ri)*Ta)/(kf*log(ri/ra) + ks*log(rb/ri));  % T on the interface

c2 = (Tb-Ti)/log(rb/ri);
c1 = Ti-c2*log(ri);
c3= ( (c1-c2)*(ri^2-rb^2)+c2*(ri^2*log(ri)-rb^2*log(rb)) )/( log(rb/ri) );
c4= ( (c1-c2)*(rb^2*log(ri)-ri^2*log(rb))+c2*log(ri)*log(rb)*(rb^2-ri^2) )/( log(rb/ri) );

Tba=Tb-Ta; 
Tbi=Tb-Ti; 

[ma,na]=size(uf2);

n=na;
rs=zeros(n,1); Ts=zeros(n,1); rf=zeros(n,1); Tf=zeros(n,1); uf=zeros(n,1);
for i=1:n
  rr = ra + (i-1)*(ri-ra)/(n-1);
  rs(i) = rr;
  Ts(i) = Ta + log(rr/ra)*(Ti-Ta)/log(ri/ra);
  rr = ri + (i-1)*(rb-ri)/(n-1);
  rf(i) = rr;
  Tf(i) = Ti + log(rr/ri)*(Tb-Ti)/log(rb/ri);

  uf(i) = beta/(4.*nu)*( (c1-c2)*rr^2 + c2*rr^2*log(rr)+c3*log(rr)+c4);

end;

% compute errors
uErr=0.; TErr=0.;
ue=zeros(na,1); Tse=zeros(n,1); Tfe=zeros(n,1);
for i=1:na

  rr=rs2(i);
  TTrue = Ta + log(rr/ra)*(Ti-Ta)/log(ri/ra);
  TErr = max(TErr,abs(TTrue-Ts2(i)));
  Tse(i)=Ts2(i)-TTrue;


  rr = rf2(i);
  uTrue = beta/(4.*nu)*( (c1-c2)*rr^2 + c2*rr^2*log(rr)+c3*log(rr)+c4);
  ue(i)=uTrue-uf2(i);
  uErr = max(uErr,abs(ue(i)));

  TTrue = Ti + log(rr/ri)*(Tb-Ti)/log(rb/ri);
  TErr = max(TErr,abs(TTrue-Tf2(i)));
  Tfe(i)=Tf2(i)-TTrue;


end
fprintf(' *** case=%d, vErr=%8.2e, TErr=%8.2e\n',icase, uErr,TErr);

% plot(rs2,Ts2,'b-x',rs,Ts,'c-+', rf2,Tf2,'r-x',rf,Tf,'c-+','MarkerSize',10);
plot(rs2,Ts2,'bx',rs,Ts,'k-', rf2,Tf2,'r+',rf,Tf,'k-','MarkerSize',15);
title(sprintf('Cocentric Cylinders Heat Exchanger',icase));
legend('T_s','T_s (exact)','T_f','T_f (exact)');
ylabel('T');
xlabel('r');
grid on;

print('-depsc2',sprintf('%sT.eps',name));
pause;

plot(rf2,uf2,'b-x',rf,uf,'k-','MarkerSize',15);
% plot(rf2,uf2,'r+',rf,uf,'b-');
title(sprintf('Cocentric Cylinders Heat Exchanger',icase));
legend('v','v (exact)');
ylabel('v');
xlabel('r');
grid on;
print('-depsc2',sprintf('%su.eps',name));
pause;

plot(rs2,Tse,'b-x',rf2,Tfe,'r-+','MarkerSize',15);
% title(sprintf('Cocentric Cylinders Heat Exchanger icase %d, max u-err=%8.2e',icase,uErr));
title(sprintf('Cocentric Cylinders Heat Exchanger'));
legend('T_s error','T_f error');
ylabel('T (error)');
xlabel('r');
grid on;
print('-depsc2',sprintf('%sTerr.eps',name));
pause;

plot(rf2,ue,'b-x','MarkerSize',15);
% title(sprintf('Flat plates heat exchanger icase %d, max u-err=%8.2e',icase,uErr));
title(sprintf('Cocentric Cylinders Heat Exchanger'));
legend('v error','Location','North');
ylabel('v (error)');
xlabel('r');
grid on;
print('-depsc2',sprintf('%suerr.eps',name));
pause;

% set(gca,'XLim',[-1.,1.]);
% set(gca,'YLim',[-.05,.05]);
% set(gca,'XTick',-1.:.1:1.);
% set(gca,'YTick',-.05:.005:.05);
% grid on;


% print -depsc2 afmOne16ky38t0p4GlassInterfaceEnergyDensity.eps
% print -depsc2 afmOne32ky38t0p4GlassInterfaceEnergyDensity.eps

% print -depsc2 afmOne16kym25t0p5GlassInterfaceEnergyDensity.eps
% print -depsc2 afmOne32kym50t0p5GlassInterfaceEnergyDensity.eps

% if profile == 1 
%   print -depsc2 afm.profile1.eps
% else
%   print -depsc2 afm.profile2.eps
% end;
