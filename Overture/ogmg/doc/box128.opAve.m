t=[1 2 3 4 5 6 7 ];
defect=[5.168e+00 1.094e-01 2.241e-03 5.168e-05 1.191e-06 2.801e-08 6.799e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: box128.','BC: DDDDDD.','Second-order accurate.','2.4e+06 points. 5 levels.','ave CR=0.023, ECR=0.44');
% pause
% print -deps2 residual.eps


