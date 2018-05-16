t=[0 1 2 3 4 5 6 7 8 9 ];
defect=[4.000e+00 9.655e+01 4.939e+00 1.965e-01 7.536e-03 2.879e-04 1.079e-05 4.087e-07 1.537e-08 4.657e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: square1024. 04/01/27','BC: NMNN.','Second-order accurate.','1.1e+06 points. 5 levels.','CR=0.038, ECR=0.59, W[2,1]');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.eps


