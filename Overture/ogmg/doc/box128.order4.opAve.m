t=[1 2 3 4 5 6 7 ];
defect=[1.304e+01 5.400e-01 1.746e-02 6.713e-04 2.839e-05 1.232e-06 5.840e-08 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: box128.order4.','BC: DDDDDD.','Fourth-order accurate.','2.4e+06 points. 5 levels.','ave CR=0.040, ECR=0.50');
% pause
% print -deps2 residual.eps


