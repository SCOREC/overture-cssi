t=[1 2 3 4 5 6 7 8 9 ];
defect=[2.672e+00 1.138e-01 1.052e-02 3.082e-03 2.747e-04 4.460e-05 5.028e-06 8.152e-07 7.774e-08 ];
clf
set(gca,'FontSize',14);
grid off
plot(t,defect,'r-o');
title('Residuals');
ylabel('Max norm');
xlabel('iteration');
set(gca,'YScale','log');
pause
% print -deps2 residual.eps


