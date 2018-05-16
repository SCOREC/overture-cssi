t=[1 2 3 4 5 6 7 8 9 ];
defect=[5.372e-01 5.236e-02 1.201e-03 1.175e-04 8.709e-06 8.830e-07 6.335e-08 5.468e-09 3.594e-10 ];
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


defect0=[2.248e-01 1.291e-02 9.681e-04 7.124e-05 5.554e-06 4.409e-07 3.534e-08 2.819e-09 2.263e-10 ];
defect1=[5.372e-01 2.925e-02 7.145e-04 8.864e-05 7.629e-06 6.167e-07 6.335e-08 5.468e-09 3.155e-10 ];
defect2=[4.446e-01 5.236e-02 1.199e-03 1.174e-04 8.691e-06 8.806e-07 4.283e-08 4.507e-09 3.576e-10 ];
defect3=[4.445e-01 5.236e-02 1.201e-03 1.175e-04 8.709e-06 8.830e-07 4.300e-08 4.527e-09 3.594e-10 ];
plot(t,defect0,'r-o',t,defect1,'g-x',t,defect2,'b-s',t,defect3,'c-<');
title('Residuals by component grid, variable smoothing');
ylabel('Max norm');
legend('grid1','grid2','grid3','grid4');
xlabel('iteration');
set(gca,'YScale','log');
pause
print -deps2 ellipsoidVarSmooth.eps
