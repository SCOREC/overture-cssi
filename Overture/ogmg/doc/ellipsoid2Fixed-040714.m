t=[0 1 2 3 4 5 6 7 8 9 10 ];
defect=[6.635e+02 1.378e+01 2.876e+00 8.913e-01 3.717e-01 1.912e-01 1.011e-01 5.412e-02 2.928e-02 1.601e-02 8.835e-03 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: ellipsoid2.bbmg. 04/07/14','BC: DDDDDD+PPIIDI+IIIIDI+....','Second-order accurate.','6.5e+05 points. 3 levels.','CR=0.442, ECR=0.84, V[2,1]');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.eps


defect0=[2.212e-01 4.606e-03 8.669e-04 3.005e-04 1.287e-04 6.206e-05 2.995e-05 1.544e-05 7.766e-06 4.234e-06 ];
defect1=[1.412e+00 3.360e-01 1.406e-01 6.507e-02 3.069e-02 1.472e-02 7.216e-03 3.629e-03 1.879e-03 9.985e-04 ];
defect2=[1.508e+00 1.663e-01 1.928e-02 2.300e-03 2.806e-04 3.516e-05 4.975e-06 7.543e-07 1.593e-07 3.648e-08 ];
defect3=[1.582e+00 1.744e-01 2.012e-02 2.371e-03 2.848e-04 3.484e-05 4.824e-06 7.107e-07 1.527e-07 3.451e-08 ];
plot(t,defect0,'r-o',t,defect1,'g-x',t,defect2,'b-s',t,defect3,'c-<');
title('Residuals by component grid, fixed smoothing');
ylabel('L_2 norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
% pause
% print -deps2 residualByGrid.eps
% print -depsc2 residual.eps


