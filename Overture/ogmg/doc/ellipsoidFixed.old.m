t=[1 2 3 4 5 6 7 8 9 ];
defect=[1.261e+01 9.593e-01 3.493e-01 1.556e-01 6.934e-02 3.091e-02 1.379e-02 6.162e-03 2.757e-03 ];
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


defect0=[1.261e+01 9.593e-01 7.591e-02 6.037e-03 4.860e-04 4.211e-05 9.567e-06 4.179e-06 1.866e-06 ];
defect1=[1.793e+00 7.937e-01 3.493e-01 1.556e-01 6.934e-02 3.091e-02 1.379e-02 6.162e-03 2.757e-03 ];
defect2=[1.649e+00 2.576e-01 4.767e-02 1.149e-02 2.851e-03 7.164e-04 1.812e-04 4.589e-05 1.158e-05 ];
defect3=[1.650e+00 2.578e-01 4.772e-02 1.145e-02 2.841e-03 7.141e-04 1.806e-04 4.575e-05 1.154e-05 ];
plot(t,defect0,'r-o',t,defect1,'g-x',t,defect2,'b-s',t,defect3,'c-<');
title('Residuals by component grid, fixed smoothing');
ylabel('Max norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
pause
print -deps2 ellipsoidFixedSmooth.eps
