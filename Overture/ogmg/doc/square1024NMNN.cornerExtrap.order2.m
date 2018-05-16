t=[1 2 3 4 5 6 7 ];
defect=[1.378e+00 5.982e-02 2.613e-03 1.143e-04 4.997e-06 2.184e-07 9.959e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: square1024.','BC: NMNN.','Second-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.044, ECR=0.60');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.eps


