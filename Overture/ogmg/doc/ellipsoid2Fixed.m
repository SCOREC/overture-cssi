t=[1 2 3 4 5 6 7 8 9 10 ];
defect=[1.445e+01 4.997e+00 2.257e+00 1.109e+00 5.796e-01 3.061e-01 1.634e-01 8.801e-02 4.785e-02 2.624e-02 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: ellipsoid2.bbmg.','BC: DDDDDD+PPIIDI+IIIIDI+....','Second-order accurate.','6.5e+05 points. 3 levels.','ave CR=0.446, ECR=0.84');
% pause
% print -deps2 residual.eps


defect0=[6.269e-02 1.733e-02 7.898e-03 3.560e-03 1.698e-03 8.274e-04 4.151e-04 2.134e-04 1.123e-04 6.040e-05 ];
defect1=[1.725e+00 7.285e-01 3.250e-01 1.491e-01 7.005e-02 3.371e-02 1.661e-02 8.380e-03 4.330e-03 2.290e-03 ];
defect2=[4.030e-01 7.554e-02 2.700e-02 8.288e-03 3.045e-03 1.044e-03 3.719e-04 1.348e-04 4.595e-05 1.770e-05 ];
defect3=[3.168e-01 5.037e-02 1.511e-02 4.285e-03 1.476e-03 4.997e-04 1.686e-04 6.367e-05 1.981e-05 8.590e-06 ];
plot(t,defect0,'r-o',t,defect1,'g-x',t,defect2,'b-s',t,defect3,'c-<');
title('Residuals by component grid, fixed smoothing');
ylabel('L_2 norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
% pause
% print -deps2 residualByGrid.eps
