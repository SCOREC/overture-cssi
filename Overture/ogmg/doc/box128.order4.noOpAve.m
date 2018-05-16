t=[1 2 3 4 5 6 7 ];
defect=[1.376e+02 1.303e+01 1.321e+00 1.510e-01 1.753e-02 2.029e-03 2.346e-04 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: box128.order4.','BC: DDDDDD.','Fourth-order accurate.','2.4e+06 points. 5 levels.','ave CR=0.106, ECR=0.62');
% pause
% print -deps2 residual.eps


