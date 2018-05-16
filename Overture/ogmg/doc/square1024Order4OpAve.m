t=[1 2 3 4 5 6 7 ];
defect=[1.225e+02 2.062e+00 3.128e-02 6.018e-04 1.183e-05 3.044e-07 7.504e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: square1024.order4.','BC: DDDD.','Fourth-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.018, ECR=0.45');
% pause
% print -deps2 residual.eps


