t=[1 2 3 4 5 6 7 8 ];
defect=[2.728e+02 8.314e+00 2.698e-01 9.930e-03 3.833e-04 1.505e-05 6.095e-07 2.507e-08 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: square1024.','BC: DDDD.','Second-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.035, ECR=0.51');
% pause
% print -deps2 residual.eps


