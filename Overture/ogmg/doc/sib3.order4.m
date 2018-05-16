t=[0 1 2 3 4 5 6 7 ];
defect=[4.086e+03 5.907e-01 1.618e-02 8.415e-04 1.814e-05 1.527e-06 5.238e-08 1.978e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: sib3.order4. 04/01/13','BC: DDDDDD+IIIIDI+IIIIDI....','Fourth-order accurate.','2.8e+06 points. 3 levels.','CR=0.039, ECR=0.59, W[2,1]');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.sib3.order4.W21.eps


