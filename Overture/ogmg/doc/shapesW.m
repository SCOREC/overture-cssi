t=[0 1 2 3 4 5 6 ];
defect=[6.486e+06 1.315e+01 3.660e-01 8.447e-03 1.368e-04 3.828e-06 1.308e-07 ];
t2=[0 1 2 3 4 5 6 ];
defect2=[6.486e+06 5.399e+00 6.471e-02 1.270e-03 2.351e-05 4.176e-07 7.102e-09 ];
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


% print -deps2 residual.shapesW.bbmg5.eps
% print -depsc2 residual.shapesW.bbmg5.eps

text(t(3),defect(3),' \leftarrow W[1,1], CR=.025, ECR=.52','HorizontalAlignment','left','FontSize',16);
text(t2(5),defect2(5),'W[2,1], CR=.017, ECR=.57 \rightarrow ','HorizontalAlignment','right','FontSize',16);

% legend('Grid: shapes.bbmg5. 04/01/13','BC: DDDD+PPDI+PPDI+PPDI....','Second-order accurate.','3.4e+06 points. 5 levels.','CR=0.025, ECR=0.52, W[1,1]');
% legend('Grid: shapes.bbmg5. 04/01/13','BC: DDDD+PPDI+PPDI+PPDI....','Second-order accurate.','3.4e+06 points. 5 levels.','CR=0.017, ECR=0.57, W[2,1]');

% F[1,1]: legend('Grid: shapes.bbmg5. 03/07/19','BC: DDDD+PPDI+PPDI+PPDI....','Second-order accurate.','3.4e+06 points. 5 levels.','CR=0.023, ECR=0.49, F[1,1]');
% F[2,1]: legend('Grid: shapes.bbmg5. 03/07/19','BC: DDDD+PPDI+PPDI+PPDI....','Second-order accurate.','3.4e+06 points. 5 levels.','CR=0.017, ECR=0.53, F[2,1]');
