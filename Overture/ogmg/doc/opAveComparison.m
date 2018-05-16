t=[1 2 3 4 5 6 7  ];
defect=[2.728e+02 8.314e+00 2.698e-01 9.930e-03 3.833e-04 1.505e-05 6.095e-07  ];
defect2=[1.013e+02 1.597e+00 2.544e-02 3.875e-04 5.993e-06 1.471e-07 3.548e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Square 1024^2, second-order accuracy');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('No operator averaging','Operator averaging');
text(t(3),defect(3),' \leftarrow CR=0.035, ECR=0.51','HorizontalAlignment','left','FontSize',16);
text(t(5),defect2(5),'CR=0.018, ECR=0.45 \rightarrow ','HorizontalAlignment','right','FontSize',16);

% legend('Grid: square1024.','BC: DDDD.','Second-order accurate.','1.1e+06 points. 5 levels.','ave CR=0.018, ECR=0.45');
pause
print -depsc2 opAveComparison.eps

%  Fourth-order------------------------------------------------------

t=[1 2 3 4 5 6 7  ];
defect=[2.130e+03 1.319e+02 8.942e+00 6.353e-01 4.585e-02 3.349e-03 2.467e-04 ];
defect2=[1.225e+02 2.062e+00 3.128e-02 6.018e-04 1.183e-05 3.044e-07 7.504e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Square 1024^2, fourth-order accuracy');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('No operator averaging','Operator averaging');
text(t(3),defect(3),' \leftarrow CR=0.069, ECR=0.59','HorizontalAlignment','left','FontSize',16);
text(t(5),defect2(5),'CR=0.018, ECR=0.45 \rightarrow ','HorizontalAlignment','right','FontSize',16);

pause
print -depsc2 opAveComparison4.eps

% ------------------------------------------------------

t=[1 2 3 4 5 6 7  ];
defect=[3.543e+01 2.052e+00 1.433e-01 9.983e-03 6.988e-04 5.084e-05 3.752e-06 ];
defect2=[5.168e+00 1.094e-01 2.241e-03 5.168e-05 1.191e-06 2.801e-08 6.799e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Box 128^3, second-order accuracy');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('No operator averaging','Operator averaging');
text(t(3),defect(3),' \leftarrow CR=0.065, ECR=0.55','HorizontalAlignment','left','FontSize',16);
text(t(5),defect2(5),'CR=0.023, ECR=0.44 \rightarrow ','HorizontalAlignment','right','FontSize',16);

pause
print -depsc2 opAveComparison.box128.eps

% ------------------------------------------------------

t=[1 2 3 4 5 6 7  ];
defect=[1.376e+02 1.303e+01 1.321e+00 1.510e-01 1.753e-02 2.029e-03 2.346e-04 ];
defect2=[1.304e+01 5.400e-01 1.746e-02 6.713e-04 2.839e-05 1.232e-06 5.840e-08 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Box 128^3, fourth-order accuracy');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('No operator averaging','Operator averaging');
text(t(3),defect(3),' \leftarrow CR=0.106, ECR=0.62','HorizontalAlignment','left','FontSize',16);
text(t(5),defect2(5),'CR=0.040, ECR=0.50 \rightarrow ','HorizontalAlignment','right','FontSize',16);

pause
print -depsc2 opAveComparison4.box128.eps

