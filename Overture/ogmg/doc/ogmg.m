t=[0 1 2 3 4 5 6 ];
defect=[1.974e+01 8.646e+00 8.337e-02 8.693e-04 9.725e-06 1.102e-07 2.156e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: square1024. 06/11/17','BC: NMNN.','Second-order accurate.','1.1e+06 points. 5 levels.','CR=0.012, ECR=0.50, W[2,1]');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.eps


