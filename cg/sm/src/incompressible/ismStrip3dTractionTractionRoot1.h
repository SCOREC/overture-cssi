// ----------------------------------------------------------------------------------------------------
//     caseName : Exact solution for the incompressible elasticity equations.
//  Period strip in y and z.
//    File written by ismStri3d.maple
// ----------------------------------------------------------------------------------------------------
Real rho, c, mu, omega; 
rho = 1;
c = 1;
mu = 1;
omega = 10.7649676392784327821062965544;
Real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
    ,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
    ,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
    ,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
    ,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249
    ,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299
    ,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349
    ,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399
    ,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449
    ,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499
    ,t500;
#beginMacro evalSolution(x,y,z,t, u1e,u2e,u3e,pe, v1e,v2e,v3e, a1e,a2e,a3e)
t1 = 0.888576587631673249403176198012e1 * x;
t2 = exp(t1);
t4 = 0.607681603029061635492190309752e1 * x;
t5 = cos(t4);
t7 = sin(t4);
t9 = exp(-t1);
t11 = -0.753239723406067534184872382763e-4 * t2 + 0.204598093767530818350584185509e1 * t5 - 0.211866255833191900563999356400e0 * t7 - 0.544467928195366972426789429740e0 * t9;
t12 = 0.628318530717958647692528676656e1 * z;
t13 = cos(t12);
t19 = 0.2193926127766250949726433469e-1 * t7 - 0.211866255833191900563999359800e0 * t5 + 0.779997882698741460167705417491e-5 * t2 + 0.563809658456915764794509377163e-1 * t9;
t20 = sin(t12);
t21 = t20 * t19;
t23 = 0.628318530717958647692528676656e1 * y;
t24 = cos(t23);
t26 = sin(t23);
t33 = 0.107649676392784327821062965544e2 * t;
t34 = cos(t33);
t37 = sin(t33);
t45 = t13 * t37;
t47 = t37 * t26;
u1e = t34 * (t24 * (t13 * t11 + t21) + t13 * t26 * t19 - t26 * t20 * t11) + t24 * (t37 * t20 * t11 - t13 * t37 * t19) + t45 * t26 * t11 + t47 * t21;
t51 = 0.102453787114446688167969810662e0 * t7;
t53 = 0.551541792167429340705776308962e-5 * t2 + 0.10609336515389131435859593393e-1 * t5 + t51 - 0.398673632793356429650421871652e-1 * t9;
t59 = 0.989390663484610868564140406607e0 * t7 + 0.102453787114446688167969810662e0 * t5 + 0.532620916279509775647851210991e-4 * t2 - 0.384996964165534223106925228171e0 * t9;
t60 = t20 * t59;
u2e = t34 * (t24 * (t13 * t53 + t60) + t13 * t26 * t59 - t26 * t20 * t53) + t24 * (-t13 * t37 * t59 + t37 * t20 * t53) + t45 * t26 * t53 + t47 * t60;
t83 = 0.551541792167429340705776308960e-5 * t2 + 0.106093365153891314358595930604e-1 * t5 + t51 - 0.398673632793356429650421871650e-1 * t9;
t89 = 0.989390663484610868564140406600e0 * t7 + 0.102453787114446688167969807374e0 * t5 + 0.532620916279509775647851210990e-4 * t2 - 0.384996964165534223106925228170e0 * t9;
t90 = t20 * t89;
u3e = t34 * (t24 * (t13 * t83 + t90) + t13 * t26 * t89 - t26 * t20 * t83) + t24 * (-t13 * t37 * t89 + t37 * t20 * t83) + t45 * t26 * t83 + t47 * t90;
t114 = t20 * t2;
t116 = t9 * t20;
t120 = t26 * t2;
t122 = t26 * t9;
t134 = t9 * t37;
t150 = t37 * t20;
pe = t34 * (t24 * (t13 * (-0.982344473618620248763243082755e-3 * t2 + 0.710072827687236792905290602940e1 * t9) + 0.101724137176215384777217667547e-3 * t114 - 0.735297522086965837907740907410e0 * t116) + t13 * (0.101724137176215384777217667547e-3 * t120 - 0.735297522086965837907740907410e0 * t122) + 0.982344473618620248763243082755e-3 * t20 * t120 - 0.710072827687236792905290602940e1 * t26 * t116) + t24 * (t13 * (-0.101724137176215384777217667547e-3 * t37 * t2 + 0.735297522086965837907740907410e0 * t134) - 0.982344473618620248763243082755e-3 * t37 * t114 + 0.710072827687236792905290602940e1 * t37 * t116) + t13 * (-0.982344473618620248763243082755e-3 * t37 * t120 + 0.710072827687236792905290602940e1 * t26 * t134) + 0.101724137176215384777217667547e-3 * t150 * t120 - 0.735297522086965837907740907410e0 * t150 * t122;
t159 = -0.839665196595764678524978767219e-4 * t2 + 0.228073338789939630269443895847e1 * t5 - 0.236175437683711317856477236823e0 * t7 - 0.606939272800132397554438883761e0 * t9;
t165 = -0.228073338789939630269443892187e1 * t7 + 0.220249185846552366494056465318e2 * t5 - 0.810860124708535449352960613025e-3 * t2 - 0.586117962764809884791888083821e1 * t9;
t166 = t20 * t165;
v1e = t34 * (t24 * (t13 * t159 + t166) + t13 * t26 * t165 - t26 * t20 * t159) + t24 * (-t13 * t37 * t165 + t37 * t20 * t159) + t45 * t26 * t159 + t47 * t166;
t190 = -0.573364692775175013716097731084e-3 * t2 - 0.110291170280854028050252083644e1 * t5 - 0.106507584750160537696246687098e2 * t7 + 0.414447986046241432672814397105e1 * t9;
t192 = 0.110291170280854028050252083644e1 * t7;
t196 = t192 + 0.114209164262379012481627844618e0 * t5 + 0.593732954439200783815785753958e-4 * t2 - 0.429170875565405494811044982635e0 * t9;
t197 = t20 * t196;
v2e = t34 * (t24 * (t13 * t190 + t197) + t13 * t26 * t196 - t26 * t20 * t190) + t24 * (-t13 * t37 * t196 + t37 * t20 * t190) + t45 * t26 * t190 + t47 * t197;
t221 = -0.573364692775175013716097731086e-3 * t2 - 0.110291170280854028050252080105e1 * t5 - 0.106507584750160537696246687097e2 * t7 + 0.414447986046241432672814397104e1 * t9;
t226 = t192 + 0.114209164262379012481627841038e0 * t5 + 0.593732954439200783815785753951e-4 * t2 - 0.429170875565405494811044982636e0 * t9;
t227 = t20 * t226;
v3e = t34 * (t24 * (t13 * t221 + t227) + t13 * t26 * t226 - t26 * t20 * t221) + t24 * (-t13 * t37 * t226 + t37 * t20 * t221) + t45 * t26 * t221 + t47 * t227;
t251 = 0.872888300246865845979060058975e-2 * t2 - 0.237097535821555763859928806507e3 * t5 + 0.245520211145588663285774316485e2 * t7 + 0.630954090196297983278707421990e2 * t9;
t257 = -0.254242094385757243832246831452e1 * t7 + 0.245520211145588663285774320425e2 * t5 - 0.903896866918177004550662388841e-3 * t2 - 0.653368163070060996488433067136e1 * t9;
t258 = t20 * t257;
a1e = t34 * (t24 * (t13 * t251 + t258) + t13 * t26 * t257 - t26 * t20 * t251) + t24 * (-t13 * t37 * t257 + t37 * t20 * t251) + t45 * t26 * t251 + t47 * t258;
t280 = 0.118728087897154083062681942597e2 * t7;
t282 = -0.639151604091117254913306627012e-3 * t2 - 0.122945795739354494986511003638e1 * t5 - t280 + 0.462001058718238122043635271806e1 * t9;
t288 = -0.114655070317318329149234456111e3 * t7 - 0.118728087897154083062681942597e2 * t5 - 0.617225236322957965178608717936e-2 * t2 + 0.446151915795190848609705083554e2 * t9;
t289 = t20 * t288;
a2e = t34 * (t24 * (t13 * t282 + t289) + t13 * t26 * t288 - t26 * t20 * t282) + t24 * (-t13 * t37 * t288 + t37 * t20 * t282) + t45 * t26 * t282 + t47 * t289;
t312 = -0.639151604091117254913306627010e-3 * t2 - 0.122945795739354494986510999784e1 * t5 - t280 + 0.462001058718238122043635271804e1 * t9;
t318 = -0.114655070317318329149234456110e3 * t7 - 0.118728087897154083062681938787e2 * t5 - 0.617225236322957965178608717933e-2 * t2 + 0.446151915795190848609705083552e2 * t9;
t319 = t20 * t318;
a3e = t34 * (t24 * (t13 * t312 + t319) + t13 * t26 * t318 - t26 * t20 * t312) + t24 * (-t13 * t37 * t318 + t37 * t20 * t312) + t45 * t26 * t312 + t47 * t319;
#endMacro
