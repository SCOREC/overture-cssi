t=[0 1 2 3 4 5 6 ];
defect=[2.556e+05 5.724e-01 1.325e-02 5.837e-04 2.247e-05 8.100e-07 2.888e-08 ];
t2=[0 1 2 3 4 5 6 ];
defect2=[2.556e+05 8.843e-02 2.327e-03 5.994e-05 1.394e-06 3.626e-08 7.529e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t2,defect2,'b-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
%legend('Grid: joukowsky. 04/01/12','BC: DDDD+PPDI+PPDI.','Second-order accurate.','2.2e+06 points. 4 levels.','CR=0.035, ECR=0.54, F[1,1]');
% legend('Grid: joukowsky. 04/01/12','BC: DDDD+PPDI+PPDI.','Second-order accurate.','2.2e+06 points. 4 levels.','CR=0.024, ECR=0.59, F[2,1]');

yMax=max(defect(1),defect2(1))*1.5;
yMin=1.e-10;
set(gca,'YLim',[yMin yMax]);

yPos=yMax/50.;
text((t(1)+t(7))/2.,yPos,'Joukowsky, 2.2e+06 points, 4 levels.','Fontsize',16,'HorizontalAlignment','center');

text(t(3),defect(3),' \leftarrow W[1,1], CR=.035, ECR=.54','HorizontalAlignment','left','FontSize',16);
text(t2(5),defect2(5),'W[2,1], CR=.024, ECR=.59 \rightarrow ','HorizontalAlignment','right','FontSize',16);

% pause
% print -deps2 residual.eps


% print -deps2 residual.joukowsky.fmg.eps
% print -depsc2 residual.joukowsky.fmg.eps


