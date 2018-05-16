
clear;
% legend('Grid: box128.order4. 04/01/14','BC: DDDDDD.','Fourth-order accurate.','2.4e+06 points. 5 levels.','CR=0.040, ECR=0.52, W[2,1]');
t=[0 1 2 3 4 5 6 7 8 ];
defect=[6.553e+04 3.016e+02 1.255e+01 4.547e-01 1.530e-02 5.930e-04 2.415e-05 1.100e-06 5.210e-08 ];


t=[0 1 2 3 4 5 6 7 8 ];
defect=[4.915e+04 1.863e+02 2.902e+00 4.059e-02 5.531e-04 9.775e-06 1.607e-07 2.428e-09 4.190e-11 ];
tn=[0 1 2 3 4 5 6 7 ];
defectn=[2.961e+01 4.565e+00 1.049e-01 1.247e-03 2.615e-05 3.313e-07 7.160e-09 1.041e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',tn,defectn,'b-x');
title('Multigrid Convergence, Box128, W[2,1]','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
% legend('BC=Dirichlet','BC=Neumann');
text(t(4),defect(4),' \leftarrow Dirichlet, CR=.040, ECR=.52','HorizontalAlignment','left','FontSize',16);
text(tn(5),defectn(5),'Neumann, CR=0.017 \rightarrow ','HorizontalAlignment','right','FontSize',16);

yMax=max(defect(1),defectn(1))*1.5;
yMin=1.e-11;
set(gca,'YLim',[yMin yMax]);

% print -depsc2 residual.box128.order4.DN.eps


%legend('Grid: box128. 03/07/17','BC: DDDDDD.','Second-order accurate.','2.4e+06 points. 5 levels.','CR=0.016, ECR=0.43, W[2,1]');
% legend('Grid: box128. 03/07/17','BC: MMMNNN.','Second-order accurate.','2.4e+06 points. 5 levels.','CR=0.017, ECR=0.44, W[2,1]');
