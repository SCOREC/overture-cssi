t=[1 2 3 4 5 6 7 8 9 10 ];
defect=[6.516e+00 2.618e+00 1.127e+00 4.843e-01 2.091e-01 9.049e-02 3.934e-02 1.717e-02 7.533e-03 3.321e-03 ];
clf
set(gca,'FontSize',16);
grid off
plot(t,defect,'r-o');
title('Residuals');
ylabel('Max norm');
xlabel('iteration');
set(gca,'YScale','log');
pause
% print -deps2 residual.eps


defect0=[2.247e-01 5.860e-02 2.246e-02 9.033e-03 3.680e-03 1.517e-03 6.322e-04 2.668e-04 1.141e-04 4.945e-05 ];
defect1=[1.566e+00 6.200e-01 2.489e-01 1.013e-01 4.144e-02 1.713e-02 7.139e-03 3.011e-03 1.286e-03 5.562e-04 ];
defect2=[8.405e-01 1.953e-01 5.514e-02 1.812e-02 6.116e-03 2.229e-03 7.899e-04 2.965e-04 1.058e-04 3.960e-05 ];
defect3=[8.230e-01 1.867e-01 5.357e-02 1.737e-02 6.130e-03 2.206e-03 8.163e-04 3.020e-04 1.118e-04 4.126e-05 ];
plot(t,defect0,'r-o',t,defect1,'g-x',t,defect2,'b-s',t,defect3,'c-<');
title('Residuals by component grid, fixed smoothing');
ylabel('L_2 norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
pause
% print -deps2 residualByGrid.eps
print -deps2 ellipsoidFixedSmooth.eps