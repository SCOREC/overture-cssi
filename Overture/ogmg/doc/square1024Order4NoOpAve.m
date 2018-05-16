t=[1 2 3 4 5 6 7 8 9 ];
defect=[2.130e+03 1.319e+02 8.942e+00 6.353e-01 4.585e-02 3.349e-03 2.467e-04 1.830e-05 1.363e-06 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: square1024.order4.','BC: DDDD.','Fourth-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.069, ECR=0.59');
% pause
% print -deps2 residual.eps


