t=[1 2 3 4 5 6 7 8 9 10 11 ];
defect=[2.386e+01 8.674e+00 3.817e+00 1.812e+00 8.551e-01 4.039e-01 1.904e-01 8.976e-02 4.493e-02 2.347e-02 1.237e-02 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: ellipsoid2.bbmg.','BC: DDDDDD+PPIIDI+IIIIDI+IIIIDI.','Second-order accurate.','6.5e+05 points. 3 levels.','ave CR=0.429, ECR=0.83');
% pause
% print -deps2 residual.eps


defect0=[1.757e-01 2.951e-02 1.275e-02 5.380e-03 2.439e-03 1.117e-03 5.228e-04 2.483e-04 1.198e-04 5.874e-05 2.929e-05 ];
defect1=[2.943e+00 1.198e+00 5.189e-01 2.306e-01 1.044e-01 4.804e-02 2.240e-02 1.058e-02 5.059e-03 2.450e-03 1.204e-03 ];
defect2=[8.134e-01 1.219e-01 3.601e-02 1.082e-02 3.833e-03 1.303e-03 4.727e-04 1.638e-04 5.997e-05 2.074e-05 7.695e-06 ];
defect3=[7.903e-01 1.157e-01 2.787e-02 7.957e-03 2.617e-03 8.703e-04 3.095e-04 1.059e-04 3.867e-05 1.319e-05 4.933e-06 ];
plot(t,defect0,'r-o',t,defect1,'g-x',t,defect2,'b-s',t,defect3,'c-<');
title('Residuals by component grid, fixed smoothing');
ylabel('L_2 norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
% pause
% print -deps2 residualByGrid.eps
