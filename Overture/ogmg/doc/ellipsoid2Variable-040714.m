t=[0 1 2 3 4 5 6 7 8 9 10 ];
defect=[6.635e+02 1.115e+01 3.363e-01 7.058e-03 2.731e-04 8.276e-06 2.209e-07 5.454e-09 2.045e-10 1.037e-11 1.501e-12 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o');
title('Multigrid Convergence','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Grid: ellipsoid2.bbmg. 04/07/14','BC: DDDDDD+PPIIDI+IIIIDI+....','Second-order accurate.','6.5e+05 points. 3 levels.','CR=0.037, ECR=0.64, V[2,1]');
% pause
% print -deps2 residual.eps


% print -depsc2 residual.eps


defect0=[2.237e-01 5.409e-03 1.615e-04 6.113e-06 1.784e-07 5.052e-09 1.312e-10 3.998e-12 1.566e-13 4.404e-14 ];
defect1=[1.185e+00 3.579e-02 6.563e-04 1.852e-05 3.804e-07 1.030e-08 3.240e-10 1.661e-11 7.328e-13 8.427e-14 ];
defect2=[9.005e-01 7.867e-03 8.590e-05 2.880e-06 2.531e-07 9.939e-09 2.457e-10 5.578e-12 1.580e-13 8.243e-14 ];
defect3=[9.399e-01 8.887e-03 1.039e-04 4.546e-06 4.256e-07 9.374e-09 1.892e-10 5.333e-12 2.241e-13 7.318e-14 ];
plot(t,defect0,'r-o',t,defect1,'g-x',t,defect2,'b-s',t,defect3,'c-<');
title('Residuals by component grid, variable smoothing');
ylabel('L_2 norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
% pause
% print -deps2 residualByGrid.eps
% print -depsc2 residual.eps


