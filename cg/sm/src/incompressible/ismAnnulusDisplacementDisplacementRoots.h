// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2021-10-18
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
    1.2510520647340855e+01, 1.7985120095545963e+01, 2.5104649597130992e+01, 3.0907710595982864e+01,
    // n=2, m=1,2,..,4
    1.2334522283695935e+01, 1.8035467567538374e+01, 2.5019752889578412e+01, 3.0937071838665323e+01,
    // n=3, m=1,2,..,4
    1.2176999775790313e+01, 1.8133089033841518e+01, 2.4949743731740190e+01, 3.0993818552125962e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
