%
% Compare BC's for the 4th order method.

t=[1 2 3 4 5 6 7 8 9 ];
defect=[2.760e+03 2.059e+02 1.533e+01 1.135e+00 8.343e-02 6.101e-03 4.462e-04 3.296e-05 2.433e-06 ];
defect2=[1.279e+03 6.687e+01 3.549e+00 1.921e-01 1.049e-02 5.770e-04 3.173e-05 1.746e-06 9.592e-08 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Fourth-order accuracy, square 1024^2');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Extrapolation BC','Equation BC');
text(t(3),defect(3),' \leftarrow CR=0.074, ECR=0.60','HorizontalAlignment','left','FontSize',16);
text(t(6),defect2(6),'CR=0.054, ECR=0.56 \rightarrow ','HorizontalAlignment','right','FontSize',16);

pause
print -depsc2 fourthOrderBC.square1024.eps

% -----------------------------------------------------


t=[1 2 3 4 5 6 7 8 9 ];
defect=[1.918e+02 1.286e+01 9.590e-01 7.099e-02 5.220e-03 3.818e-04 2.787e-05 2.058e-06 1.519e-07 ];
defect2=[7.923e+01 4.137e+00 2.200e-01 1.187e-02 6.477e-04 3.559e-05 1.964e-06 1.083e-07 5.971e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Fourth-order accuracy, circle in a channel, 1.1e+06 points');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Extrapolation BC','Equation BC');
text(t(3),defect(3),' \leftarrow CR=0.074, ECR=0.61','HorizontalAlignment','left','FontSize',16);
text(t(6),defect2(6),'CR=0.054, ECR=0.58 \rightarrow ','HorizontalAlignment','right','FontSize',16);

pause
print -depsc2 fourthOrderBC.cic6.eps


% -----------------------------------------------------

t=[1 2 3 4 5 ];
defect=[2.861e-01 6.393e-03 2.336e-04 6.317e-06 2.884e-07  ];
defect2=[3.172e-01 3.594e-03 7.079e-05 1.230e-06 2.104e-08 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Fourth-order accuracy, BC=NMNN, W[2,1], square 1024^2 ');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Corner extrapolation BC','Corner Taylor-series BC');
text(t(2),defect(2),' \leftarrow CR=0.028, ECR=0.56','HorizontalAlignment','left','FontSize',16);
text(t(4),defect2(4),'CR=0.015, ECR=0.51 \rightarrow ','HorizontalAlignment','right','FontSize',16);

pause
print -depsc2 fourthOrderCornerBC.square1024.eps

% -----------------------------------------------------

t=[1 2 3 4 5 ];
defect=[1.378e+00 5.982e-02 2.613e-03 1.143e-04 4.997e-06 ];
defect2=[6.400e-02 5.311e-04 5.527e-06 5.278e-08 9.949e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Second-order accuracy, BC=NMNN, W[2,1], square 1024^2 ');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Corner extrapolation BC','Corner Taylor-series BC');
text(t(2),defect(2),' \leftarrow CR=0.044, ECR=0.60','HorizontalAlignment','left','FontSize',16);
text(t(4),defect2(4),'CR=0.011, ECR=0.48 \rightarrow ','HorizontalAlignment','right','FontSize',16);

pause
print -depsc2 secondOrderCornerBC.square1024.eps
