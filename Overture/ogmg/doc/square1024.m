clear;

t=[0 1 2 3 4 5 6 7 8 ];
defect=[2.097e+06 3.810e+03 5.590e+01 5.197e-01 7.591e-03 8.700e-05 1.270e-06 1.624e-08 4.415e-10 ];
tn=[0 1 2 3 4 5 6];
defectn=[1.974e+01 7.317e+00 6.400e-02 5.311e-04 5.527e-06 5.278e-08 9.949e-10 ];
clf
set(gca,'FontSize',16);
plot(t,defect,'r-o',tn,defectn,'b-x');
title('Multigrid Convergence, Square1024, W[2,1]','FontSize',18);
ylabel('maximum residual');
xlabel('multigrid cycle');
set(gca,'YScale','log');
grid on
% legend('BC=Dirichlet','BC=Neumann');
text(t(4),defect(4),' \leftarrow Dirichlet, CR=.014, ECR=.51','HorizontalAlignment','left','FontSize',16);
text(tn(6),defectn(6),'Neumann, CR=.011, ECR=.49 \rightarrow ','HorizontalAlignment','right','FontSize',16);

yMax=max(defect(1),defectn(1))*1.5;
yMin=1.e-10;
set(gca,'YLim',[yMin yMax]);


% pause
% print -deps2 residual.eps


% print -deps2 residual.square1024.DN.eps
% print -depsc2 residual.square1024.DN.eps

% **********************************************************
% Dirichlet:
% \begin{tabular}{|c|c|c|c|c|} \hline 
%  $i$   & $\vert\vert\mbox{res}\vert\vert_\infty$  &  CR     &  WU    & ECR  \\   \hline 
%  $ 1$  & $ 3.8e+03$ & $0.002$ & $ 8.1$ & $0.46$ \\ 
%  $ 2$  & $ 5.6e+01$ & $0.015$ & $ 6.1$ & $0.50$ \\ 
%  $ 3$  & $ 5.2e-01$ & $0.009$ & $ 6.1$ & $0.47$ \\ 
%  $ 4$  & $ 7.6e-03$ & $0.015$ & $ 6.1$ & $0.50$ \\ 
%  $ 5$  & $ 8.7e-05$ & $0.011$ & $ 6.1$ & $0.48$ \\ 
%  $ 6$  & $ 1.3e-06$ & $0.015$ & $ 6.1$ & $0.50$ \\ 
%  $ 7$  & $ 1.6e-08$ & $0.013$ & $ 7.1$ & $0.54$ \\ 
%  $ 8$  & $ 4.4e-10$ & $0.027$ & $ 7.1$ & $0.60$ \\ 
% \hline 
% \multicolumn{5}{|c|}{Grid: square1024. 04/01/28}  \\
% \multicolumn{5}{|c|}{BC: DDDD.}  \\
% \multicolumn{5}{|c|}{Second-order accurate.}  \\
% \multicolumn{5}{|c|}{Trigonometric solution.}  \\
% \multicolumn{5}{|c|}{W[2,1]: rb $\omega=1.10$}  \\
% \multicolumn{5}{|c|}{1.06e+06 grid-points. 5 levels.}  \\
% \multicolumn{5}{|c|}{Average CR=$0.014$, ECR=$0.51$.}  \\
% \multicolumn{5}{|c|}{time/cycle = 6.07e-01 s.}  \\
% \hline 

% ****************************************************
% Neumann: NMNN

% \begin{tabular}{|c|c|c|c|c|} \hline 
%  $i$   & $\vert\vert\mbox{res}\vert\vert_\infty$  &  CR     &  WU    & ECR  \\   \hline 
%  $ 1$  & $ 7.3e+00$ & $0.371$ & $ 8.1$ & $0.88$ \\ 
%  $ 2$  & $ 6.4e-02$ & $0.009$ & $ 6.1$ & $0.46$ \\ 
%  $ 3$  & $ 5.3e-04$ & $0.008$ & $ 6.1$ & $0.46$ \\ 
%  $ 4$  & $ 5.5e-06$ & $0.010$ & $ 6.1$ & $0.47$ \\ 
%  $ 5$  & $ 5.3e-08$ & $0.010$ & $ 6.1$ & $0.47$ \\ 
%  $ 6$  & $ 9.9e-10$ & $0.019$ & $ 7.1$ & $0.57$ \\ 
% \hline 
% \multicolumn{5}{|c|}{Grid: square1024. 04/01/28}  \\
% \multicolumn{5}{|c|}{BC: NMNN.}  \\
% \multicolumn{5}{|c|}{Second-order accurate.}  \\
% \multicolumn{5}{|c|}{Trigonometric solution.}  \\
% \multicolumn{5}{|c|}{W[2,1]: rb $\omega=1.10$}  \\
% \multicolumn{5}{|c|}{1.06e+06 grid-points. 5 levels.}  \\
% \multicolumn{5}{|c|}{Average CR=$0.011$, ECR=$0.49$.}  \\
% \multicolumn{5}{|c|}{time/cycle = 6.79e-01 s.}  \\
