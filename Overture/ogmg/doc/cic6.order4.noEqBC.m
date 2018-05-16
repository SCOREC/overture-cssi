t=[1 2 3 4 5 6 7 8 9 ];
defect=[1.918e+02 1.286e+01 9.590e-01 7.099e-02 5.220e-03 3.818e-04 2.787e-05 2.058e-06 1.519e-07 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: cic6.order4.','BC: DDDD+PPDI.','Fourth-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.074, ECR=0.61');
% pause
% print -deps2 residual.eps


