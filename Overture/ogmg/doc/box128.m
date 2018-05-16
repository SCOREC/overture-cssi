
clear;

t=[0 1 2 3 4 5 6 7 8 ];
defect=[4.915e+04 1.863e+02 2.902e+00 4.059e-02 5.531e-04 9.775e-06 1.607e-07 2.428e-09 4.190e-11 ];
tn=[0 1 2 3 4 5 6 7 ];
defectn=[2.961e+01 4.565e+00 1.049e-01 1.247e-03 2.615e-05 3.313e-07 7.160e-09 1.041e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',tn,defectn,'b-x');
title('Multigrid Convergence, Box128, W[2,1]','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
% legend('BC=Dirichlet','BC=Neumann');
text(t(4),defect(4),' \leftarrow Dirichlet, CR=.016, ECR=.45','HorizontalAlignment','left','FontSize',16);
text(tn(6),defectn(6),'Neumann, CR=.017, ECR=.45 \rightarrow ','HorizontalAlignment','right','FontSize',16);

yMax=max(defect(1),defectn(1))*1.5;
yMin=1.e-11;
set(gca,'YLim',[yMin yMax]);

% print -deps2 residual.box128.DN.eps
% print -depsc2 residual.box128.DN.eps


%legend('Grid: box128. 03/07/17','BC: DDDDDD.','Second-order accurate.','2.4e+06 points. 5 levels.','CR=0.016, ECR=0.43, W[2,1]');
% legend('Grid: box128. 03/07/17','BC: MMMNNN.','Second-order accurate.','2.4e+06 points. 5 levels.','CR=0.017, ECR=0.44, W[2,1]');


% **** Dirichlet
%  $i$   & $\vert\vert\mbox{res}\vert\vert_\infty$  &  CR     &  WU    & ECR  \\   \hline 
%  $ 1$  & $ 1.9e+02$ & $0.004$ & $ 6.9$ & $0.45$ \\ 
%  $ 2$  & $ 2.9e+00$ & $0.016$ & $ 4.9$ & $0.43$ \\ 
%  $ 3$  & $ 4.1e-02$ & $0.014$ & $ 4.9$ & $0.42$ \\ 
%  $ 4$  & $ 5.5e-04$ & $0.014$ & $ 4.9$ & $0.42$ \\ 
%  $ 5$  & $ 9.8e-06$ & $0.018$ & $ 4.9$ & $0.44$ \\ 
%  $ 6$  & $ 1.6e-07$ & $0.016$ & $ 4.9$ & $0.44$ \\ 
%  $ 7$  & $ 2.4e-09$ & $0.015$ & $ 5.9$ & $0.49$ \\ 
%  $ 8$  & $ 4.2e-11$ & $0.017$ & $ 5.9$ & $0.51$ \\ 
% \hline 
% \multicolumn{5}{|c|}{Grid: box128. 04/01/28}  \\
% \multicolumn{5}{|c|}{BC: DDDDDD.}  \\
% \multicolumn{5}{|c|}{Second-order accurate.}  \\
% \multicolumn{5}{|c|}{Trigonometric solution.}  \\
% \multicolumn{5}{|c|}{W[2,1]: rb $\omega=1.12$}  \\
% \multicolumn{5}{|c|}{2.35e+06 grid-points. 5 levels.}  \\
% \multicolumn{5}{|c|}{Average CR=$0.016$, ECR=$0.45$.}  \\
% \multicolumn{5}{|c|}{time/cycle = 2.03e+00 s.}  \\


% ******************* NMNNNN
% \begin{tabular}{|c|c|c|c|c|} \hline 
%  $i$   & $\vert\vert\mbox{res}\vert\vert_\infty$  &  CR     &  WU    & ECR  \\   \hline 
%  $ 1$  & $ 2.4e+00$ & $0.082$ & $ 6.9$ & $0.70$ \\ 
%  $ 2$  & $ 5.6e-02$ & $0.023$ & $ 4.9$ & $0.47$ \\ 
%  $ 3$  & $ 6.4e-04$ & $0.011$ & $ 4.9$ & $0.41$ \\ 
%  $ 4$  & $ 1.4e-05$ & $0.022$ & $ 4.9$ & $0.46$ \\ 
%  $ 5$  & $ 1.7e-07$ & $0.012$ & $ 4.9$ & $0.41$ \\ 
%  $ 6$  & $ 4.0e-09$ & $0.023$ & $ 4.9$ & $0.47$ \\ 
%  $ 7$  & $ 6.0e-11$ & $0.015$ & $ 5.9$ & $0.49$ \\ 
% \hline 
% \multicolumn{5}{|c|}{Grid: box128. 04/01/28}  \\
% \multicolumn{5}{|c|}{BC: NMNNNN.}  \\
% \multicolumn{5}{|c|}{Second-order accurate.}  \\
% \multicolumn{5}{|c|}{Trigonometric solution.}  \\
% \multicolumn{5}{|c|}{W[2,1]: rb $\omega=1.12$}  \\
% \multicolumn{5}{|c|}{2.35e+06 grid-points. 5 levels.}  \\
% \multicolumn{5}{|c|}{Average CR=$0.017$, ECR=$0.45$.}  \\
% \multicolumn{5}{|c|}{time/cycle = 1.96e+00 s.}  \\
