t=[0 1 2 3 4 5 6 7 8 9 ];
defect=[2.556e+05 7.900e+01 8.602e-01 1.811e-02 1.538e-03 5.889e-05 6.933e-06 3.497e-07 4.125e-08 2.079e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: joukowsky. 04/07/14','BC: DDDD+PPDI+PPDI.','Second-order accurate.','2.2e+06 points. 4 levels.','CR=0.048, ECR=0.64, W[2,1]');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.eps


