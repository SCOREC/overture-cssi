t=[0 1 2 3 4 5 6 7 8 9 ];
defect=[2.207e+04 1.882e+01 1.538e+00 2.363e-01 2.691e-02 3.242e-03 3.865e-04 6.180e-05 5.646e-06 9.351e-07 ];
%
t2=[0 1 2 3 4 5 6 7 8 ];
defect2=[2.207e+04 5.304e+00 3.736e-01 4.424e-02 3.316e-03 3.179e-04 3.160e-05 2.859e-06 2.074e-07 ];


clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t2,defect2,'b-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
yMax=max(defect(1),defect2(1))*1.5;
yMax=1.e5;
yMin=1.e-7;
set(gca,'YLim',[yMin yMax]);
xMin=0;
xMax=9;
set(gca,'XLim',[xMin xMax]);

yPos=yMax/50.;
text((t(1)+t(10))/2.,yPos,'subAndFins2, 3.6e+06 points, 3 levels.','Fontsize',16,'HorizontalAlignment','center');

% pause
% print -deps2 residual.eps


% print -deps2 residual.subAndFins2.fmg.eps
% print -depsc2 residual.subAndFins2.fmg.eps

text(t(3),defect(3),' \leftarrow W[1,1], CR=.122, ECR=.62','HorizontalAlignment','left','FontSize',16);
text(t2(7),defect2(7),'W[2,1], CR=.087, ECR=.65 \rightarrow ','HorizontalAlignment','right','FontSize',16);

% legend('Grid: subAndFins2.hdf. 03/07/25','BC: DDDDDD+IIPPDI+IIIIDI+....','Second-order accurate.','3.6e+06 points. 3 levels.','CR=0.122, ECR=0.62, W[1,1]');
% legend('Grid: subAndFins2.hdf. 03/07/25','BC: DDDDDD+IIPPDI+IIIIDI+....','Second-order accurate.','3.6e+06 points. 3 levels.','CR=0.087, ECR=0.65, W[2,1]');
