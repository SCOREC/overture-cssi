t=[1 2 3 4 5 6 7 ];
defect=[3.543e+01 2.052e+00 1.433e-01 9.983e-03 6.988e-04 5.084e-05 3.752e-06 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: box128.','BC: DDDDDD.','Second-order accurate.','2.4e+06 points. 5 levels.','ave CR=0.065, ECR=0.55');
% pause
% print -deps2 residual.eps


