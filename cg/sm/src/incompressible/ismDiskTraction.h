// ------------------------------------------------------------------------------------------------ 
//    caseName: Exact solution for incompressible linear elasticity 
// Solid disk of radius R=1; r=0 bounded solution; r=1 Traction BC. 
// File Written by diskAnnulusSolution.mw (maple) on 2021-10-18 
// Set c3 to scale disk solution and c4 to scale the annulus solution 
// Set n and m to choose the solution, Jn(omega*r), omega=omega(n,m) 
// The eigenvalues are save in a separate include file. 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omega; 
rho=1.;
mu=1.;
omega=omegaRoot(n,m);
Real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
    ,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
    ,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
    ,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
    ,t200;
#beginMacro diskMacro(r,theta,t,u1e,u2e,pe)
t1 = omega * r;
t2 = jn(n, t1);
t3 = t2 * c3;
t4 = 0.1e1 / r;
t5 = omega * t;
t6 = cos(t5);
t7 = t6 * t4;
t8 = n * theta;
t9 = cos(t8);
t12 = sin(t5);
t13 = t12 * t4;
t14 = sin(t8);
t17 = n * n;
t19 = omega * omega;
t23 = 0.1e1 / (0.2e1 * t17 - t19 - 0.2e1 * n) * c3;
t24 = pow(r, n);
t26 = t4 * t24 * t23;
t27 = n + 0.1e1;
t28 = jn(t27, omega);
t29 = t28 * n;
t30 = t6 * omega;
t35 = t12 * omega;
t40 = jn(n, omega);
t41 = t40 * t17;
t42 = t9 * t6;
t46 = t14 * t12;
t50 = t40 * n;
t51 = t42 * t50;
t54 = t46 * t50;
u1e = 0.2e1 * t14 * t35 * t29 * t26 + 0.2e1 * t9 * t30 * t29 * t26 + t14 * t13 * t3 - 0.2e1 * t42 * t41 * t26 - 0.2e1 * t46 * t41 * t26 + t9 * t7 * t3 + 0.2e1 * t51 * t26 + 0.2e1 * t54 * t26;
t57 = t14 * t6;
t63 = jn(t27, t1);
t64 = t63 / n * omega;
t75 = t9 * t12;
u2e = -t9 * t12 * c3 * t64 + t14 * t6 * c3 * t64 - 0.2e1 * t14 * t30 * t29 * t26 + 0.2e1 * t9 * t35 * t29 * t26 + t9 * t13 * t3 - t14 * t7 * t3 + 0.2e1 * t57 * t41 * t26 - 0.2e1 * t75 * t41 * t26 - 0.2e1 * t57 * t50 * t26 + 0.2e1 * t75 * t50 * t26;
t94 = t19 * omega * t23;
t95 = t28 * t24;
t101 = t24 * t19 * t23;
t104 = t19 * t23;
t105 = t40 * t24;
pe = 0.2e1 * t42 * t105 * t104 + 0.2e1 * t46 * t105 * t104 + 0.2e1 * t42 * t95 * t94 + 0.2e1 * t46 * t95 * t94 - 0.2e1 * t51 * t101 - 0.2e1 * t54 * t101;
#endMacro

#beginMacro diskMacroVelocity(r,theta,t,v1e,v2e)
t1 = n * n;
t2 = omega * omega;
t4 = t1 - t2 / 0.2e1 - n;
t5 = omega * r;
t6 = jn(n, t5);
t8 = pow(r, n);
t9 = n + 0.1e1;
t10 = jn(t9, omega);
t12 = jn(n, omega);
t18 = t6 * t4 - n * (-t10 * omega + (n - 0.1e1) * t12) * t8;
t19 = omega * t;
t20 = cos(t19);
t21 = n * theta;
t22 = sin(t21);
t24 = sin(t19);
t25 = cos(t21);
t34 = 0.1e1 / (0.2e1 * t1 - t2 - 0.2e1 * n) * omega;
t35 = 0.1e1 / r;
v1e = 0.2e1 * t35 * t34 * c3 * (t22 * t20 - t25 * t24) * t18;
t43 = jn(t9, t5);
v2e = 0.2e1 / n * t35 * t34 * (-t43 * omega * t4 * r + n * t18) * c3 * (t25 * t20 + t22 * t24);
#endMacro

#beginMacro diskMacroAcceleration(r,theta,t,a1e,a2e)
t1 = n * n;
t2 = omega * omega;
t4 = t1 - t2 / 0.2e1 - n;
t5 = omega * r;
t6 = jn(n, t5);
t8 = pow(r, n);
t9 = n + 0.1e1;
t10 = jn(t9, omega);
t12 = jn(n, omega);
t18 = t6 * t4 - n * (-t10 * omega + (n - 0.1e1) * t12) * t8;
t19 = omega * t;
t20 = cos(t19);
t21 = n * theta;
t22 = cos(t21);
t24 = sin(t19);
t25 = sin(t21);
t34 = 0.1e1 / (0.2e1 * t1 - t2 - 0.2e1 * n) * t2;
t35 = 0.1e1 / r;
a1e = -0.2e1 * t35 * t34 * c3 * (t22 * t20 + t25 * t24) * t18;
t44 = jn(t9, t5);
a2e = 0.2e1 / n * t35 * t34 * (-t44 * omega * t4 * r + n * t18) * c3 * (t25 * t20 - t22 * t24);
#endMacro
