t=[0 1 2 3 4 5 6 7 ];
defect=[1.589e+03 2.022e-01 1.406e-02 3.595e-04 3.173e-05 4.698e-07 2.847e-08 5.691e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: sib2.order4. 04/01/13','BC: DDDDDD+IIIIDI+IIIIDI....','Fourth-order accurate.','8.4e+05 points. 3 levels.','CR=0.038, ECR=0.60, W[2,1]');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.sib2.order4.eps


