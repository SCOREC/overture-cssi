// File generated by Dropbox/DARPA/RPI/adePapers/adegdmi/dispersivePlaneWaveInterface3d.maple
real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150;
// -------------------------------------------------------------------------
// Need: 
//    s=sr + I*si : complex frequency
//    [ax,ay,az] = amplitude vector of incident wave, ax=axr+I*axi etc.
//    [arx,ary,arz] = amplitude vector of reflected wave, arx=arxr+I*arxi etc.
//    [atx,aty,atz] = amplitude vector of transmitted wave, atx=atxr+I*atxi etc.
//    [kx,ky,kz] = incident wave vector, kx=kxr+I*kxi etc. 
//    [krx,kry,krz] = reflected wave vector, krx=krxr+I*krxi etc. 
//    [ktx,kty,ktz] = reflected wave vector, ktx=ktxr+I*ktxi etc. 
// -------------------------------------------------------------------------
// Evaluated:                                                               
//    [Exr,Exi] ,[Eyr,Eyi], [Ezr,Ezi]        : left state                   
//    [Etxr,Etxi] ,[Etyr,Etyi], [Etzr,Etzi]  : right state                  
// -------------------------------------------------------------------------
real Exr,Eyr,Ezr, Exi,Eyi,Ezi;
real Etxr,Etyr,Etzr, Etxi,Etyi,Etzi;
t2 = exp(t * sr);
t3 = t * si;
t4 = cos(t3);
t5 = t2 * t4;
t9 = x * krxr + y * kryr + z * krzr;
t10 = cos(t9);
t15 = exp(-x * krxi - y * kryi - z * krzi);
t16 = t10 * t15;
t17 = t16 * arxr;
t19 = sin(t9);
t20 = t15 * t19;
t21 = t20 * arxi;
t27 = exp(-x * kxi - y * kyi - z * kzi);
t31 = x * kxr + y * kyr + z * kzr;
t32 = sin(t31);
t33 = t27 * t32;
t34 = t33 * axi;
t36 = cos(t31);
t37 = t27 * t36;
t38 = t37 * axr;
t40 = sin(t3);
t41 = t2 * t40;
t42 = t16 * arxi;
t44 = t20 * arxr;
t46 = t33 * axr;
t48 = t37 * axi;
Exr = t5 * t17 - t5 * t21 - t5 * t34 + t5 * t38 - t41 * t42 - t41 * t44 - t41 * t46 - t41 * t48;
t50 = t16 * aryr;
t52 = t20 * aryi;
t54 = t33 * ayi;
t56 = t37 * ayr;
t58 = t16 * aryi;
t60 = t20 * aryr;
t62 = t33 * ayr;
t64 = t37 * ayi;
Eyr = -t41 * t58 - t41 * t60 - t41 * t62 - t41 * t64 + t5 * t50 - t5 * t52 - t5 * t54 + t5 * t56;
t66 = t16 * arzr;
t68 = t20 * arzi;
t70 = t33 * azi;
t72 = t37 * azr;
t74 = t16 * arzi;
t76 = t20 * arzr;
t78 = t33 * azr;
t80 = t37 * azi;
Ezr = -t41 * t74 - t41 * t76 - t41 * t78 - t41 * t80 + t5 * t66 - t5 * t68 - t5 * t70 + t5 * t72;
Exi = t41 * t17 - t41 * t21 - t41 * t34 + t41 * t38 + t5 * t42 + t5 * t44 + t5 * t46 + t5 * t48;
Eyi = t41 * t50 - t41 * t52 - t41 * t54 + t41 * t56 + t5 * t58 + t5 * t60 + t5 * t62 + t5 * t64;
Ezi = t41 * t66 - t41 * t68 - t41 * t70 + t41 * t72 + t5 * t74 + t5 * t76 + t5 * t78 + t5 * t80;
t110 = exp(-x * ktxi - y * ktyi - z * ktzi);
t114 = x * ktxr + y * ktyr + z * ktzr;
t115 = cos(t114);
t116 = t110 * t115;
t117 = t5 * atxr;
t119 = t41 * atxi;
t121 = sin(t114);
t122 = t110 * t121;
t123 = t5 * atxi;
t125 = t41 * atxr;
Etxr = t116 * t117 - t116 * t119 - t122 * t123 - t122 * t125;
t127 = t5 * atyr;
t129 = t41 * atyi;
t131 = t5 * atyi;
t133 = t41 * atyr;
Etyr = t116 * t127 - t116 * t129 - t122 * t131 - t122 * t133;
t135 = t5 * atzr;
t137 = t41 * atzi;
t139 = t5 * atzi;
t141 = t41 * atzr;
Etzr = t116 * t135 - t116 * t137 - t122 * t139 - t122 * t141;
Etxi = t116 * t123 + t116 * t125 + t122 * t117 - t122 * t119;
Etyi = t116 * t131 + t116 * t133 + t122 * t127 - t122 * t129;
Etzi = t116 * t139 + t116 * t141 + t122 * t135 - t122 * t137;

