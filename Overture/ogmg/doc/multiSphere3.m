t=[0 1 2 3 4 5 6 7 8 9 ];
defect=[3.072e+03 1.440e+01 6.224e-01 1.263e-02 6.916e-04 4.059e-05 8.990e-06 4.315e-07 8.301e-09 5.091e-10 ];
t2=[0 1 2 3 4 5 6 7 ];
defect2=[3.072e+03 6.367e-01 1.405e-02 2.659e-04 7.396e-06 2.664e-07 6.667e-09 2.825e-10 ];

% No FMG:
% t=[0 1 2 3 4 5 6 7 8 9 10 ];
% defect=[3.072e+03 4.466e+02 8.340e+00 1.567e-01 2.340e-02 3.245e-03 5.086e-05 8.781e-07 1.308e-07 1.773e-08 2.675e-10 ];
% t2=[0 1 2 3 4 5 6 7 ];
% defect2=[3.072e+03 5.490e+01 2.358e-01 5.050e-03 1.585e-04 2.952e-06 8.422e-08 1.856e-09 ];

clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t2,defect2,'b-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
yMax=max(defect(1),defect2(1))*1.5;
yMin=1.e-10;
set(gca,'YLim',[yMin yMax]);
xMin=0;
xMax=9;
set(gca,'XLim',[xMin xMax]);

yPos=yMax/50.;
text((t(1)+t(10))/2.,yPos,'multiSphere3, 5.4e+06 points, 4 levels.','Fontsize',16,'HorizontalAlignment','center');

% pause
% print -deps2 residual.eps


% print -deps2 residual.multiSphere3.fmg.eps
% print -depsc2 residual.multiSphere3.fmg.eps

% No FMG:text(t(3),defect(3),' \leftarrow W[1,1], CR=0.044, ECR=0.55','HorizontalAlignment','left','FontSize',16);
% No FMG:text(t2(7),defect2(7),'W[2,1], CR=0.018, ECR=0.63 \rightarrow ','HorizontalAlignment','right','FontSize',16);
text(t(4),defect(4),' \leftarrow W[1,1], CR=.049, ECR=.56','HorizontalAlignment','left','FontSize',16);
text(t2(6),defect2(6),'W[2,1], CR=.028, ECR=.60 \rightarrow ','HorizontalAlignment','right','FontSize',16);

% legend('Grid: multiSphere3. 03/07/23','BC: DDDDDD+IIIIDI+IIIIDI+....','Second-order accurate.','5.4e+06 points. 4 levels.','CR=0.049, ECR=0.56, W[1,1]');
% legend('Grid: multiSphere3. 03/07/23','BC: DDDDDD+IIIIDI+IIIIDI+....','Second-order accurate.','5.4e+06 points. 4 levels.','CR=0.028, ECR=0.60, W[2,1]');

% No FMG:
% legend('Grid: multiSphere3. 03/07/22','BC: DDDDDD+IIIIDI+IIIIDI+....','Second-order accurate.','5.4e+06 points. 4 levels.','CR=0.044, ECR=0.55, W[1,1]');
% W21 legend('Grid: multiSphere3. 03/07/22','BC: DDDDDD+IIIIDI+IIIIDI+....','Second-order accurate.','5.4e+06 points. 4 levels.','CR=0.018, ECR=0.63, W[2,1]');




pause

t=[0 1 2 3 4 5 6 7 ];
defect=[3.072e+03 8.832e-01 1.318e-02 1.276e-04 2.593e-06 4.275e-07 6.845e-09 4.956e-11 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: multiSphere3.','BC: DDDDDD+IIIIDI+IIIIDI+....','Second-order accurate.','5.4e+06 points. 4 levels.','CR=0.020, ECR=0.63, W[2,1]');

yMax=defect(1)*1.5;
yMin=1.e-11;
set(gca,'YLim',[yMin yMax]);

% pause
% print -deps2 residual.eps


% print -depsc2 residual.multiSphere3.fmg.eps



% 
% 
% t=[1 2 3 4 5 6 7 ];
% defect=[1.700e+00 4.545e-02 2.686e-03 5.579e-05 1.546e-06 8.339e-08 1.842e-09 ];
% clf
% set(gca,'FontSize',16);
% plot(t,defect,'r-o');
% title('Multigrid Convergence');
% ylabel('maximum residual');
% xlabel('multigrid cycle');
% set(gca,'YScale','log');
% grid on
% legend('Grid: multiSphere3.','BC: DDDDDD+IIIIDI+IIIIDI+....','Second-order accurate.','5.4e+06 points. 4 levels.','ave CR=0.033, ECR=0.61');
% % pause
% % print -deps2 residual.eps


