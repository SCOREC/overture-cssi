// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2021-10-18
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
    5.1356223018406826e+00, 8.4172441403998649e+00, 1.1619841172149059e+01, 1.4795951782351261e+01,
    // n=2, m=1,2,..,4
    6.3801618959239835e+00, 9.7610231299816697e+00, 1.3015200721698434e+01, 1.6223466160318768e+01,
    // n=3, m=1,2,..,4
    7.5883424345038044e+00, 1.1064709488501185e+01, 1.4372536671617590e+01, 1.7615966049804833e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
