%
% Compute the least square fit to the errors
% 
clear;
clf;

h = [ .1 .1/2. .1/4. ];
eT = [ 5.98e-3 1.76e-3 4.47e-4 ];
ev = [ 1.33e-3 3.85e-4 9.42e-5 ];

% e = C h^sigma -> log(e) = sigma*log(h) + log(C)


degree=1;
pT = polyfit( log(h), log(eT), degree);
pv = polyfit( log(h), log(ev), degree);

fprintf(' convergence rates:   T : %5.2f,  v : %5.2f\n',pT(1),pv(1));
