t=[0 1 2 3 4 5 6 ];
defect=[3.072e+03 1.523e+00 3.032e-02 6.249e-04 1.296e-05 2.629e-07 5.261e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: sib3.bbmg.','BC: DDDDDD+IIIIDI+IIIIDI....','Second-order accurate.','2.8e+06 points. 4 levels.','CR=0.020, ECR=0.55, W[2,1]');

yMax=defect(1)*1.5;
yMin=1.e-10;
set(gca,'YLim',[yMin yMax]);

% pause
% print -deps2 residual.eps


% print -depsc2 residual.eps


