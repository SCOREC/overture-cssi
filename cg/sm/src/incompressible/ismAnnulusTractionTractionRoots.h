// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2021-10-18
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
    3.4346674494392022e+00, 7.7407996399021592e+00, 1.3009389736124731e+01, 1.9138855322118575e+01,
    // n=2, m=1,2,..,4
    4.7964119338586172e+00, 9.4022451914702559e+00, 1.3622383066594675e+01, 1.9446727547714004e+01,
    // n=3, m=1,2,..,4
    2.5439592523077787e+00, 6.2948319530425706e+00, 1.0674374290844092e+01, 1.4841467907331412e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
