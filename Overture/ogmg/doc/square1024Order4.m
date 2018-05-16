clear;
% legend('Grid: square1024.order4. 04/01/14','BC: DDDD.','Fourth-order accurate.','1.1e+06 points. 5 levels.','CR=0.018, ECR=0.54, W[2,1]');
t=[0 1 2 3 4 5 6 7 8 ];
defect=[2.796e+06 8.603e+03 1.419e+02 1.755e+00 4.277e-02 7.942e-04 1.549e-05 2.831e-07 4.820e-09 ];
% legend('Grid: square1024.order4. 04/01/14','BC: NMNN.','Fourth-order accurate.','1.1e+06 points. 5 levels.','CR=0.016, ECR=0.52, W[2,1]');
tn=[0 1 2 3 4 5 6 ];
defectn=[1.974e+01 2.096e+01 2.469e-01 3.142e-03 6.763e-05 1.181e-06 2.071e-08 ];
%
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',tn,defectn,'b-x');
title('Multigrid Convergence, Square1024, W[2,1]','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
% legend('BC=Dirichlet','BC=Neumann');
text(t(4),defect(4),' \leftarrow Dirichlet, CR=.018, ECR=.54','HorizontalAlignment','left','FontSize',16);
text(tn(6),defectn(6),'Neumann, CR=.016, ECR=.52 \rightarrow ','HorizontalAlignment','right','FontSize',16);

yMax=max(defect(1),defectn(1))*1.5;
yMin=1.e-10;
set(gca,'YLim',[yMin yMax]);


% pause
% print -deps2 residual.eps


% print -depsc2 residual.square1024.order4.DN.eps
