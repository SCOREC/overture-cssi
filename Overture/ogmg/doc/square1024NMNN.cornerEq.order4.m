t=[1 2 3 4 5 ];
defect=[3.172e-01 3.594e-03 7.079e-05 1.230e-06 2.104e-08 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: square1024.order4.','BC: NMNN.','Fourth-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.015, ECR=0.51');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.eps


