! File generated by overtureFramework/cg/mx/codes/dispersion.maple
! Here is root 3 from the dispersion relation exp( i*k*x + s*t) .
      t1 = a1*ap
      t4 = b1*ck+t1*ck
      t5 = ck ** 2
      t6 = 1./t5
      t9 = t4 ** 2
      t10 = t5 ** 2
      t11 = 1./t10
      t12 = t9*t11
      t14 = a0*ap
      t15 = t14+t5+b0
      t16 = t15*t6
      t18 = a0*a1
      t19 = ap ** 2
      t24 = a1 ** 2
      t25 = t24*t19
      t26 = t5*b0
      t32 = a0 ** 2
      t33 = t32*a0
      t34 = t19*ap
      t35 = t33*t34
      t37 = t32*t19
      t40 = b1 ** 2
      t41 = t40*t5
      t46 = b0*b1
      t47 = t46*t5
      t50 = t40*t10
      t52 = t10*t5
      t58 = b0*t40
      t59 = t58*t5
      t61 = b0*t10
      t63 = b0 ** 2
      t66 = t63*t5
      t68 = t63*b0
      t70 = t24*a1
      t71 = t70*t34
      t72 = t40*b1
      t76 = t32*t24
      t77 = t19 ** 2
      t86 = a0*t24
      t91 = t24 ** 2
      t95 = t46*t10
      t98 = t40 ** 2
      t99 = t98*t10
      t102 = t40*t52
      t110 = t77*b0
      t114 = t32*a1
      t130 = t63*b1
      t131 = t130*t5
      t134 = t58*t10
      t137 = b0*t52
      t144 = 12.*t71*t72*t10-3.*t76*t77*t40*t5-0.54E2*a0*t70*t77*t47-6.*t86*t34*t40*t10+0.81E2*t91*t77*t66-0.54E2*t71*t95+0.36E2*t25*t99-3.*t25*t102+12.*t33*t24*t77*ap*b0+0.36E2*t76*t110*t5-6.*t114*t34*t72*t5-0.168E3*t86*t34*t59+0.36E2*t86*t34*b0*t10-0.66E2*t18*t19*t72*t10+0.270E3*t71*t131-0.150E3*t25*t134+12.*t25*t137+0.36E2*t1*t98*b1*t10
      t165 = t34*t63
      t169 = t18*t19
      t170 = b0*t72
      t180 = t63*t40
      t181 = t180*t5
      t184 = t63*t10
      t198 = t10 ** 2
      t201 = -0.60E2*t1*t72*t52+0.24E2*t33*a1*t110*b1+12.*t35*t41+0.36E2*t76*t77*t63+0.312E3*t114*t34*t47-3.*t37*t98*t5+0.36E2*t37*t50-0.360E3*t86*t165*t5-0.174E3*t169*t170*t5+0.552E3*t169*t95-0.60E2*t14*t99+0.36E2*t14*t102+0.321E3*t25*t181-0.396E3*t25*t184-0.192E3*t1*t170*t10+0.264E3*t1*t46*t52+12.*t98*t40*t10+0.24E2*t98*t52+12.*t40*t198
      t203 = t32 ** 2
      t223 = b0*t98
      t231 = t68*t5
      t251 = -0.48E2*t203*t77*b0+12.*t35*t58-0.192E3*t35*t26+0.72E2*t114*t165*b1+0.312E3*t37*t59-0.288E3*t37*t61+0.36E2*t86*t34*t68-0.240E3*t169*t131-0.60E2*t14*t223*t5+0.156E3*t14*t134-0.192E3*t14*t137-0.396E3*t25*t231+0.156E3*t1*t63*t72*t5+0.264E3*t1*t130*t10-0.96E2*t223*t10-0.144E3*t58*t52-0.48E2*b0*t198-0.192E3*t35*t63+0.36E2*t37*t180
      t262 = t63 ** 2
      t278 = t68*t40
      t298 = -0.552E3*t1*t68*b1*t5+0.72E2*t18*t19*t68*b1+0.24E2*t1*t262*b1+0.24E2*t63*t98*t5-0.48E2*t262*b0+0.264E3*t180*t10-0.288E3*t68*t10+0.156E3*t14*t181+0.192E3*t14*t184+0.192E3*t14*t231-0.192E3*t14*t262+0.36E2*t14*t278+12.*t25*t262+12.*t262*t40+0.192E3*t262*t5-0.144E3*t278*t5-0.192E3*t37*t66-0.288E3*t37*t68+0.192E3*t63*t52
      t301 = sqrt(t144+t201+t251+t298)
      t304 = -0.36E2*t18*t19*b1*t5-0.36E2*t1*b1*t10+0.24E2*t37*b0+12.*t301*ck+0.180E3*t1*t47+0.24E2*t14*t10-0.240E3*t14*t26-0.36E2*t14*t41+0.24E2*t14*t63+0.108E3*t25*t26+0.24E2*t37*t5+8.*t35+0.72E2*t50+8.*t52+0.72E2*t59-0.264E3*t61-0.264E3*t66+8.*t68
      t305 = t304 ** (1./3.)
      t307 = t6*t305/6.
      t321 = 2./3.*(-3.*t1*b1*t5+2.*t14*b0+2.*t14*t5+t10+0.14E2*t26+t37-3.*t41+t63)*t6/t305
      t323 = sqrt(t12/4.-2./3.*t16+t307+t321)
      t340 = sqrt(t12/2.-4./3.*t16-t307-t321-(t15*t11*t4-2.*b1/ck-t9*t4/t52/4.)/t323)
      ss = -t4*t6/4.-t323/2.+t340/2.

