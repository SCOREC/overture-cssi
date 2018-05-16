t=[1 2 3 4 5 6 7 8 9 10 ];
defect=[8.062e-01 2.604e-02 9.831e-04 2.616e-05 1.555e-06 6.359e-08 2.613e-09 1.160e-10 1.086e-11 7.069e-13 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: ellipsoid2.bbmg.','BC: DDDDDD+PPIIDI+IIIIDI+....','Second-order accurate.','6.5e+05 points. 3 levels.','ave CR=0.048, ECR=0.67');
% pause
% print -deps2 residual.eps


defect0=[2.163e-02 6.831e-04 2.291e-05 7.792e-07 3.683e-08 1.578e-09 5.875e-11 2.777e-12 1.796e-13 3.007e-14 ];
defect1=[3.427e-02 9.986e-04 3.871e-05 8.231e-07 5.404e-08 1.463e-09 7.917e-11 2.227e-12 3.811e-13 6.746e-14 ];
defect2=[2.869e-02 1.106e-03 8.273e-05 9.781e-07 1.004e-07 2.690e-09 8.865e-11 2.175e-12 4.003e-13 3.886e-14 ];
defect3=[2.514e-02 1.037e-03 9.858e-05 1.323e-06 6.196e-08 1.839e-09 1.194e-10 5.496e-12 2.677e-13 3.684e-14 ];
plot(t,defect0,'r-o',t,defect1,'g-x',t,defect2,'b-s',t,defect3,'c-<');
title('Residuals by component grid, variable smoothing');
ylabel('L_2 norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
% pause
% print -deps2 residualByGrid.eps
