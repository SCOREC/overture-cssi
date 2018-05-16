t=[0 1 2 3 4 5 6 ];
defect=[6.486e+06 5.416e+00 6.794e-02 1.251e-03 2.139e-05 3.588e-07 7.105e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: shapes.bbmg5. 03/07/17','BC: DDDD+PPDI+PPDI+PPDI....','Second-order accurate.','3.4e+06 points. 5 levels.','CR=0.017, ECR=0.59, F[2,1]');
yMax=defect(1)*1.5;
yMin=1.e-10;
set(gca,'YLim',[yMin yMax]);

% pause
% print -deps2 residual.eps


% print -depsc2 residual.shapes.bbmg5.fmg.eps




