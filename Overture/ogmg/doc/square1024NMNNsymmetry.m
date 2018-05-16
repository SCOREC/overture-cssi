clear;

% Neumann with BC's applied after every RB step:
t=[1 2 3 4 5 6 ];
defect=[3.320e+01 3.251e-01 3.363e-03 3.777e-05 4.283e-07 5.006e-09 ];
% Neumann with BC's applied only at the end of the RB steps
% tn=[1 2 3 4 5 6 7 8 9 ];
% defectn=[9.655e+01 4.939e+00 1.965e-01 7.536e-03 2.879e-04 1.079e-05 4.087e-07 1.537e-08 4.657e-10 ];
tn=[1 2 3 4 5 6 7 8 9 ];
defectn=[9.655e+01 4.939e+00 1.965e-01 7.536e-03 2.879e-04 1.079e-05 4.087e-07 1.537e-08 4.657e-10 ];
%
clf
set(gca,'FontSize',16);
plot(tn,defectn,'r-o', t,defect,'b-x');
title('Multigrid Convergence, Square1024, W[2,1]','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
legend('Neumann applied at end of RB smooth','Neumann applied after each RB sub-step');
text(tn(4),defectn(4),' \leftarrow CR=.038, ECR=.59','HorizontalAlignment','left','FontSize',16);
text(t(5),defect(5),'CR=.011, ECR=.49 \rightarrow ','HorizontalAlignment','right','FontSize',16);

yMax=max(defect(1),defectn(1))*1.5;
yMin=1.e-10;
set(gca,'XLim',[1 9]);
set(gca,'XLim',[1 9]);


% pause
% print -deps2 residual.eps


% print -depsc2 residual.square1024.symmetry.eps

% **************************************************************************************
% 
% t=[0 1 2 3 4 5 6 ];
% defect=[4.000e+00 3.320e+01 3.251e-01 3.363e-03 3.777e-05 4.283e-07 5.006e-09 ];
% clf
% set(gca,'FontSize',16);
% plot(t,defect,'r-o');
% title('Multigrid Convergence','FontSize',18);
% ylabel('maximum residual');
% xlabel('multigrid cycle');
% set(gca,'YScale','log');
% grid on
% legend('Grid: square1024. 04/01/27','BC: NMNN.','Second-order accurate.','1.1e+06 points. 5 levels.','CR=0.011, ECR=0.49, W[2,1]');
% % pause
% % print -deps2 residual.eps
% 
% 
% % print -depsc2 residual.eps

% **************************************************************************************
% % Neumann with BC's applied only at the end of the RB steps
% t=[0 1 2 3 4 5 6 7 8 9 ];
% defect=[4.000e+00 9.655e+01 4.939e+00 1.965e-01 7.536e-03 2.879e-04 1.079e-05 4.087e-07 1.537e-08 4.657e-10 ];
% clf
% set(gca,'FontSize',16);
% plot(t,defect,'r-o');
% title('Multigrid Convergence','FontSize',18);
% ylabel('maximum residual');
% xlabel('multigrid cycle');
% set(gca,'YScale','log');
% grid on
% legend('Grid: square1024. 04/01/27','BC: NMNN.','Second-order accurate.','1.1e+06 points. 5 levels.','CR=0.038, ECR=0.59, W[2,1]');
% pause
