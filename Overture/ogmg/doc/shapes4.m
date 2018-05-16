% legend('Grid: shapes5.order4. 04/01/14','BC: DDDD+PPDI+PPDI+PPDI....','Fourth-order accurate.','3.4e+06 points. 4 levels.','CR=0.046, ECR=0.57, W[1,1]');
t=[0 1 2 3 4 5 6 7 8 9 10 ];
defect=[8.648e+06 6.722e+04 2.086e+03 4.840e+01 1.744e+00 9.464e-02 5.115e-03 4.904e-04 2.145e-05 1.128e-06 6.632e-08 ];
%
% legend('Grid: shapes5.order4. 04/01/14','BC: DDDD+PPDI+PPDI+PPDI....','Fourth-order accurate.','3.4e+06 points. 4 levels.','CR=0.023, ECR=0.58, W[2,1]');
t2=[0 1 2 3 4 5 6 7 8 ];
defect2=[8.648e+06 1.484e+04 2.469e+02 4.273e+00 8.062e-02 2.138e-03 5.655e-05 1.501e-06 4.586e-08 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t2,defect2,'b-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
yMax=max(defect(1),defect2(1))*1.5;
yMin=1.e-8;
set(gca,'YLim',[yMin yMax]);

yPos=yMax/50.;
text((t(1)+t(11))/2.,yPos,'shapes, 3.4e+06 points, 4 levels, 4th-order.','Fontsize',16,'HorizontalAlignment','center');

% pause
% print -deps2 residual.eps


% print -depsc2 residual.shapes.order4.eps

text(t(4),defect(4),' \leftarrow W[1,1], CR=.046, ECR=.57','HorizontalAlignment','left','FontSize',16);
text(t2(7),defect2(7),'W[2,1], CR=.023, ECR=.59 \rightarrow ','HorizontalAlignment','right','FontSize',16);

