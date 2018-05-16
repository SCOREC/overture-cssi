t=[1 2 3 4 5 6 7 8 9 10 ];
defect=[3.250e+00 3.660e-01 3.195e-02 3.237e-03 3.651e-04 4.022e-05 6.634e-06 9.553e-07 1.045e-07 1.278e-08 ];
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


defect0=[2.690e-01 2.144e-02 1.826e-03 1.878e-04 2.061e-05 2.415e-06 3.544e-07 4.922e-08 5.808e-09 7.832e-10 ];
defect1=[3.001e-01 1.696e-02 1.646e-03 1.791e-04 2.208e-05 4.148e-06 4.257e-07 3.709e-08 5.559e-09 8.266e-10 ];
defect2=[3.882e-01 2.833e-02 2.430e-03 2.193e-04 2.830e-05 3.002e-06 4.156e-07 4.307e-08 8.889e-09 1.175e-09 ];
defect3=[3.816e-01 2.862e-02 2.445e-03 2.212e-04 2.862e-05 3.010e-06 4.103e-07 5.651e-08 7.608e-09 8.192e-10 ];
plot(t,defect0,'r-o',t,defect1,'g-x',t,defect2,'b-s',t,defect3,'c-<');
title('Residuals by component grid, variable smoothing');
ylabel('L_2 norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
pause
print -deps2 ellipsoidVarSmooth.eps