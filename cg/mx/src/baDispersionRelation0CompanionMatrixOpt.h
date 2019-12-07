// Companion matrix for BA-GDM model Nd=6, File written by BAMX/papers/bamx/matlab/baDispersionRelation.maple
// eps11, and eps22 have a 1 term GDM expansion.
// kxp2 = kx^2, kyp2=ky^2 etc.
 const int Nd=6; 
 pcmc = new std::complex<LocalReal> [Nd*Nd]; 
 if( NpBA(0,0) > 1 ){ printF("DispMatPar:ERROR:  NpBA(k1=0,k2=0) > 1. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(0,1) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=0,k2=1) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(0,2) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=0,k2=2) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(0,3) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=0,k2=3) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(0,4) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=0,k2=4) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(0,5) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=0,k2=5) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(1,0) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=1,k2=0) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(1,1) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=1,k2=1) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(1,2) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=1,k2=2) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(1,3) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=1,k2=3) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(1,4) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=1,k2=4) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(1,5) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=1,k2=5) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(2,0) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=2,k2=0) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(2,1) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=2,k2=1) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(2,2) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=2,k2=2) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(2,3) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=2,k2=3) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(2,4) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=2,k2=4) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(2,5) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=2,k2=5) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(3,0) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=3,k2=0) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(3,1) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=3,k2=1) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(3,2) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=3,k2=2) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(3,3) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=3,k2=3) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(3,4) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=3,k2=4) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(3,5) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=3,k2=5) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(4,0) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=4,k2=0) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(4,1) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=4,k2=1) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(4,2) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=4,k2=2) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(4,3) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=4,k2=3) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(4,4) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=4,k2=4) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(4,5) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=4,k2=5) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(5,0) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=5,k2=0) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(5,1) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=5,k2=1) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(5,2) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=5,k2=2) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(5,3) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=5,k2=3) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(5,4) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=5,k2=4) > 0. Not supported.\n"); OV_ABORT("error"); }
 if( NpBA(5,5) > 0 ){ printF("DispMatPar:ERROR:  NpBA(k1=5,k2=5) > 0. Not supported.\n"); OV_ABORT("error"); }
 std::complex<LocalReal> t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100;
 std::complex<LocalReal> ci1j1,ci1j2,ci1j3,ci1j4,ci1j5,ci1j6,ci2j1,ci2j2,ci2j3,ci2j4,ci2j5,ci2j6,ci3j1,ci3j2,ci3j3,ci3j4,ci3j5,ci3j6,ci4j1,ci4j2,ci4j3,ci4j4,ci4j5,ci4j6,ci5j1,ci5j2,ci5j3,ci5j4,ci5j5,ci5j6,ci6j1,ci6j2,ci6j3,ci6j4,ci6j5,ci6j6, c77;
ci1j1 = 0;
ci1j2 = 0;
ci1j3 = 0;
ci1j4 = 0;
ci1j5 = 0;
t1 = kxp2;
t2 = K0(3, 3);
t3 = t1 * t2;
t4 = K0(4, 4);
t5 = kyp2;
t6 = t4 * t5;
t7 = K0(5, 5);
t8 = kzp2;
t9 = t7 * t8;
t10 = t3 + t6 + t9;
t11 = K0(0, 0);
t12 = t1 * t11;
t13 = K0(1, 1);
t14 = t5 * t13;
t15 = K0(2, 2);
t16 = t8 * t15;
t17 = t12 + t14 + t16;
t18 = b0(0, 0, 0);
t20 = a0(0, 0, 0);
t24 = 0.1e1 / t11;
t25 = 0.1e1 / t13;
t26 = t24 * t25;
t34 = 0.1e1 / t15 / t2 / t4 / t7;
ci1j6 = -t10 * (t1 * t20 + t17 * t18) * t26 * t34;
ci2j1 = 1;
ci2j2 = 0;
ci2j3 = 0;
ci2j4 = 0;
ci2j5 = 0;
t36 = b1(0, 0, 0);
t38 = a1(0, 0, 0);
ci2j6 = -t10 * (t1 * t38 + t17 * t36) * t26 * t34;
ci3j1 = 0;
ci3j2 = 1;
ci3j3 = 0;
ci3j4 = 0;
ci3j5 = 0;
t44 = t14 + t16;
t47 = t15 * t2;
t56 = (-t11 * t44 * t7 - t47 * (t12 + t14)) * t4 - t7 * t13 * t2 * (t12 + t16);
t60 = t47 * t1;
t65 = t1 * t7 * t13 * t2;
ci3j6 = (t56 * t18 + ((-t44 * t7 - t60) * t4 - t65) * t20 - t17 * t10) * t24 * t25 * t34;
ci4j1 = 0;
ci4j2 = 0;
ci4j3 = 1;
ci4j4 = 0;
ci4j5 = 0;
ci4j6 = (t56 * t36 - ((t44 * t7 + t60) * t4 + t65) * t38) * t24 * t25 * t34;
ci5j1 = 0;
ci5j2 = 0;
ci5j3 = 0;
ci5j4 = 1;
ci5j5 = 0;
ci5j6 = (-t20 * t13 * t15 * t2 * t4 * t7 - t11 * t18 * t13 * t47 * t4 * t7 + (-t2 * (t6 + t9) * t15 - t11 * t7 * (t3 + t6)) * t13 - t4 * t11 * t15 * (t3 + t9)) * t24 * t25 * t34;
ci6j1 = 0;
ci6j2 = 0;
ci6j3 = 0;
ci6j4 = 0;
ci6j5 = 1;
ci6j6 = -t38 * t24 - t36;
 cm(0,0) = ci1j1;
 cm(0,1) = ci1j2;
 cm(0,2) = ci1j3;
 cm(0,3) = ci1j4;
 cm(0,4) = ci1j5;
 cm(0,5) = ci1j6;
 cm(1,0) = ci2j1;
 cm(1,1) = ci2j2;
 cm(1,2) = ci2j3;
 cm(1,3) = ci2j4;
 cm(1,4) = ci2j5;
 cm(1,5) = ci2j6;
 cm(2,0) = ci3j1;
 cm(2,1) = ci3j2;
 cm(2,2) = ci3j3;
 cm(2,3) = ci3j4;
 cm(2,4) = ci3j5;
 cm(2,5) = ci3j6;
 cm(3,0) = ci4j1;
 cm(3,1) = ci4j2;
 cm(3,2) = ci4j3;
 cm(3,3) = ci4j4;
 cm(3,4) = ci4j5;
 cm(3,5) = ci4j6;
 cm(4,0) = ci5j1;
 cm(4,1) = ci5j2;
 cm(4,2) = ci5j3;
 cm(4,3) = ci5j4;
 cm(4,4) = ci5j5;
 cm(4,5) = ci5j6;
 cm(5,0) = ci6j1;
 cm(5,1) = ci6j2;
 cm(5,2) = ci6j3;
 cm(5,3) = ci6j4;
 cm(5,4) = ci6j5;
 cm(5,5) = ci6j6;
