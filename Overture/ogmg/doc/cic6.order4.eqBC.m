t=[1 2 3 4 5 6 7 8 9 ];
defect=[7.923e+01 4.137e+00 2.200e-01 1.187e-02 6.477e-04 3.559e-05 1.964e-06 1.083e-07 5.971e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: cic6.order4.','BC: DDDD+PPDI.','Fourth-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.054, ECR=0.58');
% pause
% print -deps2 residual.eps


