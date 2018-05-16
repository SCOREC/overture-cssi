t=[1 2 3 4 5 6 7  ];
defect=[3.106e+00 3.541e-02 1.474e-03 1.385e-04 2.051e-05 1.905e-06 2.940e-07 ]; 
defect2=[3.106e+00 3.541e-02 4.792e-04 8.011e-06 1.061e-07 1.754e-09 2.380e-11 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t,defect2,'b-x');
title('Circle in a channel, fourth-order accuracy');
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('No interpolation smoothing','With interpolation smoothing');
text(t(3),defect(3),' \leftarrow CR=0.064, ECR=0.68','HorizontalAlignment','left','FontSize',16);
text(t(5),defect2(5),'CR=0.014, ECR=0.53 \rightarrow ','HorizontalAlignment','right','FontSize',16);

pause
print -depsc2 interpSmoothing.cic3.order4.eps

clear;

t=[0 1 2 3 4 5 6 7 8 9 ];
% defect=[2.556e+05 6.780e+01 8.243e-01 2.209e-02 2.633e-03 1.164e-04 5.644e-06 3.231e-07 1.563e-08 1.340e-09 ];
defect=[2.556e+05 7.900e+01 8.602e-01 1.811e-02 1.538e-03 5.889e-05 6.933e-06 3.497e-07 4.125e-08 2.079e-09 ];
t2=[0 1 2 3 4 5 6 7 ];
% defect2=[2.556e+05 6.780e+01 7.198e-01 7.112e-03 2.443e-04 3.456e-06 1.065e-07 1.299e-09 ];
defect2=[2.556e+05 7.900e+01 8.602e-01 8.600e-03 1.003e-04 2.763e-06 7.299e-08 2.068e-09 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',t2,defect2,'b-x');
title('Joukowsky airfoils, 2.2e+06 points','FontSize',18);
yMax=max(defect(1),defect2(1))*2.;
yMin=min(defect(10),defect2(8))/2.;
set(gca,'YLim',[yMin yMax]);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('No interpolation smoothing','With interpolation smoothing');
text(t(6),defect(6),' \leftarrow CR=.048, ECR=.64','HorizontalAlignment','left','FontSize',16);
text(t2(6),defect2(6),'CR=.017, ECR=.56 \rightarrow ','HorizontalAlignment','right','FontSize',16);
% legend('Grid: joukowsky. 03/07/15','BC: DDDD+PPDI+PPDI.','Second-order accurate.','2.2e+06 points. 4 levels.','CR=0.009, ECR=0.59, W[2,1]');
% pause
% print -deps2 residual.eps


% legend('Grid: joukowsky. 03/07/15','BC: DDDD+PPDI+PPDI.','Second-order accurate.','2.2e+06 points. 4 levels.','CR=0.026, ECR=0.66, W[2,1]');

print -depsc2 interpSmoothing.joukowsky.order2.eps
print -deps2  interpSmoothing.joukowsky.order2.bw.eps
