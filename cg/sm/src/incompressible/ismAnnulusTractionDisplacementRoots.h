// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2021-10-18
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
    5.9301137018497625e+00, 1.0041228867524906e+01, 1.6222236245157285e+01, 2.2261799618829349e+01,
    // n=2, m=1,2,..,4
    7.3209270736757316e+00, 1.0788308757534492e+01, 1.6576115170026591e+01, 2.2312959733123941e+01,
    // n=3, m=1,2,..,4
    7.6148014507882346e+00, 1.2347886231096596e+01, 1.7321954724146530e+01, 2.2621546021972616e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
