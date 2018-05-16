t=[1 2 3 4 5 ];
defect=[6.400e-02 5.311e-04 5.527e-06 5.278e-08 9.949e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: square1024.','BC: NMNN.','Second-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.011, ECR=0.48');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.eps


