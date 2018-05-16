%
% Compare variable sub-smooths to fixed

t=[1 2 3 4 5 6 7 8 9 10 ];
defect=[1.445e+01 4.997e+00 2.257e+00 1.109e+00 5.796e-01 3.061e-01 1.634e-01 8.801e-02 4.785e-02 2.624e-02 ];
defect2=[8.062e-01 2.604e-02 9.831e-04 2.616e-05 1.555e-06 6.359e-08 2.613e-09 1.160e-10 1.086e-11 7.069e-13 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Ellipsoid in a Box, 6.5e+05 points');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Fixed sub-smooths','Variable sub-smooths');
text(t(3),defect(3),' \leftarrow CR=.446, ECR=.84','HorizontalAlignment','left','FontSize',16);
text(t(6),defect2(6),'CR=.048, ECR=.67 \rightarrow ','HorizontalAlignment','right','FontSize',16);

pause
print -depsc2 variableSmooth.ellipsoid2.eps
print -deps2 variableSmooth.ellipsoid2.bw.eps

% --------------------------------------------------------------

defect0=[6.269e-02 1.733e-02 7.898e-03 3.560e-03 1.698e-03 8.274e-04 4.151e-04 2.134e-04 1.123e-04 6.040e-05 ];
defect1=[1.725e+00 7.285e-01 3.250e-01 1.491e-01 7.005e-02 3.371e-02 1.661e-02 8.380e-03 4.330e-03 2.290e-03 ];
defect2=[4.030e-01 7.554e-02 2.700e-02 8.288e-03 3.045e-03 1.044e-03 3.719e-04 1.348e-04 4.595e-05 1.770e-05 ];
defect3=[3.168e-01 5.037e-02 1.511e-02 4.285e-03 1.476e-03 4.997e-04 1.686e-04 6.367e-05 1.981e-05 8.590e-06 ];
plot(t,defect0,'b-o',t,defect1,'g-x',t,defect2,'r-s',t,defect3,'m-<');
title('L_{2h} residuals by component grid, fixed smoothing');
ylabel('L_{2h} norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
pause
print -depsc2 residualByGridFixed.ellipsoid2.eps
print -deps2 residualByGridFixed.ellipsoid2.bw.eps

% --------------------------------------------------------------

defect0=[2.163e-02 6.831e-04 2.291e-05 7.792e-07 3.683e-08 1.578e-09 5.875e-11 2.777e-12 1.796e-13 3.007e-14 ];
defect1=[3.427e-02 9.986e-04 3.871e-05 8.231e-07 5.404e-08 1.463e-09 7.917e-11 2.227e-12 3.811e-13 6.746e-14 ];
defect2=[2.869e-02 1.106e-03 8.273e-05 9.781e-07 1.004e-07 2.690e-09 8.865e-11 2.175e-12 4.003e-13 3.886e-14 ];
defect3=[2.514e-02 1.037e-03 9.858e-05 1.323e-06 6.196e-08 1.839e-09 1.194e-10 5.496e-12 2.677e-13 3.684e-14 ];
plot(t,defect0,'b-o',t,defect1,'g-x',t,defect2,'r-s',t,defect3,'m-<');
title('L_{2h} residuals by component grid, variable smoothing');
ylabel('L_{2h} norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
pause
print -depsc2 residualByGridVar.ellipsoid2.eps
print -deps2 residualByGridVar.ellipsoid2.bw.eps
