t=[1 2 3 4 5 6 7 ];
defect=[1.013e+02 1.597e+00 2.544e-02 3.875e-04 5.993e-06 1.471e-07 3.548e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: square1024.','BC: DDDD.','Second-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.018, ECR=0.45');
% pause
% print -deps2 residual.eps


