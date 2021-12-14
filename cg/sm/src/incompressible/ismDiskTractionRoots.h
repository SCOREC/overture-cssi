// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2021-10-18
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
    3.0542369282271403e+00, 6.7061331941584591e+00, 9.9694678230875958e+00, 1.3170370856016123e+01,
    // n=2, m=1,2,..,4
    2.3538876347708396e+00, 4.7844442796765740e+00, 8.1766702433977513e+00, 1.1446029444951778e+01,
    // n=3, m=1,2,..,4
    3.6471263298525899e+00, 6.4781057915316590e+00, 9.6197905984024880e+00, 1.2882258435475911e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
