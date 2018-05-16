%
%  compute performance and memory numbers
%

%%%%%%%%%   2nd order %%%%%%%%
n=9;
gp  =[4.2,4.2,7.3,23.,18,16.,16.,16.,125.]; % grid points M
memo=[235,255,487,206,78,308,184,117,300];
memp=[2025,2040,3280,1500,662,1,1140,607,1];
np  =[   1,   1,   1,   8, 32,8,16,32,64];

%%%%%%%%%   4th order %%%%%%%%
%n=6;
%gp  =[4.2,7.3,23.,18,16.,125.]; % grid points M
%memo=[235,751,228,93,309,790];
%memp=[3600,5500,2920,2300,2380,  1 ];
%np  =[   1,   1,   8,  32,  32, 64];


for k=1:n


  rpg = (np(k)*memo(k)*1024.*1024.)/(gp(k)*1e6*8.);
  fprintf(' Ogmg: k=%d, mem=%e MB, np=%3d, grid-points=%5.1fM, reals/grid-point = rpg=%7.1f \n',k,memo(k),np(k),gp(k),rpg);

  rpg = (np(k)*memp(k)*1024.*1024.)/(gp(k)*1e6*8.);
  fprintf(' BICG: k=%d, mem=%e MB, np=%3d, grid-points=%5.1fM, reals/grid-point = rpg=%7.1f ratio=%4.1f\n',k,memp(k),np(k),gp(k),rpg,memp(k)/memo(k));


end;
