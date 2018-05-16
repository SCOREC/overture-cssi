t=[0 1 2 3 4 5 6 ];
defect=[6.486e+06 1.394e+01 3.539e-01 7.550e-03 1.173e-04 3.093e-06 9.493e-08 ];
t2=[0 1 2 3 4 5 6 ];
defect2=[6.486e+06 5.480e+00 6.773e-02 1.245e-03 2.142e-05 3.688e-07 6.870e-09 ];
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

yPos=yMax/50.;
text((t(1)+t(7))/2.,yPos,'shapes.bbmg5, 3.4e+06 points, 5 levels.','Fontsize',16,'HorizontalAlignment','center');

% pause
% print -deps2 residual.eps


% print -depsc2 residual.shapes.bbmg5.eps

text(t(3),defect(3),' \leftarrow F[1,1], CR=0.023, ECR=0.49','HorizontalAlignment','left','FontSize',16);
text(t2(5),defect2(5),'F[2,1], CR=0.017, ECR=0.53 \rightarrow ','HorizontalAlignment','right','FontSize',16);

% F[1,1]: legend('Grid: shapes.bbmg5. 03/07/19','BC: DDDD+PPDI+PPDI+PPDI....','Second-order accurate.','3.4e+06 points. 5 levels.','CR=0.023, ECR=0.49, F[1,1]');
% F[2,1]: legend('Grid: shapes.bbmg5. 03/07/19','BC: DDDD+PPDI+PPDI+PPDI....','Second-order accurate.','3.4e+06 points. 5 levels.','CR=0.017, ECR=0.53, F[2,1]');
