t=[1 2 3 4 5 6 7 8 9 ];
defect=[3.798e+00 7.930e-02 4.015e-03 3.054e-04 3.335e-05 3.589e-06 3.068e-07 3.172e-08 3.338e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: multiSphere3.','BC: DDDDDD+IIIIDI+IIIIDI+....','Second-order accurate.','5.4e+06 points. 4 levels.','ave CR=0.071, ECR=0.76');
% pause
% print -depsc2 residual.eps


