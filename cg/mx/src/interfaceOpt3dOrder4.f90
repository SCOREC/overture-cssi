! This file automatically generated from interfaceOpt3dOrder4.bf90 with bpp.
! *******************************************************************************
!   Interface boundary conditions - Fourth Order and 3D
!      See also interfaceOpt.bf90
! *******************************************************************************

! These next include files will define the macros that will define the difference approximations
! The actual macro is called below
!* #Include "defineDiffNewerOrder2f.h"
!* #Include "defineDiffNewerOrder4f.h"

! These next include file will define the macros that will define the difference approximations (in op/src)
! Defines getDuDx2(u,aj,ff), getDuDxx2(u,aj,ff), getDuDx3(u,aj,ff), ...  etc. 


! ****** Dimension 2 ******
 ! getDuDx2 operation count     : additions+2*multiplications+assignments
 ! getDuDx2 optimization savings: -assignments
 ! getDuDy2 operation count     : additions+2*multiplications+assignments
 ! getDuDy2 optimization savings: -assignments
 ! getDuDxx2 operation count     : 9*multiplications+3*assignments+4*additions
 ! getDuDxx2 optimization savings: -3*assignments
 ! getDuDxy2 operation count     : 5*additions+9*multiplications+assignments
 ! getDuDxy2 optimization savings: -assignments
 ! getDuDyy2 operation count     : 9*multiplications+3*assignments+4*additions
 ! getDuDyy2 optimization savings: -3*assignments
 ! getDuDxxx2 operation count     : 25*multiplications+3*assignments+9*additions
 ! getDuDxxx2 optimization savings: 2*multiplications-3*assignments
 ! getDuDxxy2 operation count     : 33*multiplications+3*assignments+15*additions
 ! getDuDxxy2 optimization savings: 2*multiplications-3*assignments
 ! getDuDxyy2 operation count     : 32*multiplications+5*assignments+16*additions
 ! getDuDxyy2 optimization savings: additions+3*multiplications-5*assignments
 ! getDuDyyy2 operation count     : 25*multiplications+3*assignments+9*additions
 ! getDuDyyy2 optimization savings: 2*multiplications-3*assignments
 ! getDuDxxxx2 operation count     : 59*multiplications+11*assignments+24*additions
 ! getDuDxxxx2 optimization savings: 3*additions+23*multiplications-11*assignments
 ! getDuDxxxy2 operation count     : 86*multiplications+11*assignments+37*additions
 ! getDuDxxxy2 optimization savings: 2*additions+24*multiplications-11*assignments
 ! getDuDxxyy2 operation count     : 105*multiplications+16*assignments+50*additions
 ! getDuDxxyy2 optimization savings: 8*additions+37*multiplications-16*assignments
 ! getDuDxyyy2 operation count     : 89*multiplications+21*assignments+52*additions
 ! getDuDxyyy2 optimization savings: 17*additions+47*multiplications-21*assignments
 ! getDuDyyyy2 operation count     : 59*multiplications+11*assignments+24*additions
 ! getDuDyyyy2 optimization savings: 3*additions+23*multiplications-11*assignments
 ! getDuDxxxxx2 operation count     : 154*multiplications+43*assignments+82*additions
 ! getDuDxxxxx2 optimization savings: 39*additions+196*multiplications-43*assignments
 ! getDuDxxxxy2 operation count     : 207*multiplications+42*assignments+102*additions
 ! getDuDxxxxy2 optimization savings: 45*additions+211*multiplications-42*assignments
 ! getDuDxxxyy2 operation count     : 249*multiplications+55*assignments+126*additions
 ! getDuDxxxyy2 optimization savings: 49*additions+241*multiplications-55*assignments
 ! getDuDxxyyy2 operation count     : 239*multiplications+67*assignments+150*additions
 ! getDuDxxyyy2 optimization savings: 97*additions+363*multiplications-67*assignments
 ! getDuDxyyyy2 operation count     : 215*multiplications+73*assignments+150*additions
 ! getDuDxyyyy2 optimization savings: 153*additions+371*multiplications-73*assignments
 ! getDuDyyyyy2 operation count     : 154*multiplications+43*assignments+82*additions
 ! getDuDyyyyy2 optimization savings: 39*additions+196*multiplications-43*assignments
 ! getDuDxxxxxx2 operation count     : 364*multiplications+134*assignments+258*additions
 ! getDuDxxxxxx2 optimization savings: 345*additions+1337*multiplications-134*assignments
 ! getDuDxxxxxy2 operation count     : 525*multiplications+149*assignments+326*additions
 ! getDuDxxxxxy2 optimization savings: 409*additions+1502*multiplications-149*assignments
 ! getDuDxxxxyy2 operation count     : 543*multiplications+173*assignments+339*additions
 ! getDuDxxxxyy2 optimization savings: 415*additions+1552*multiplications-173*assignments
 ! getDuDxxxyyy2 operation count     : 510*multiplications+172*assignments+360*additions
 ! getDuDxxxyyy2 optimization savings: 463*additions+1755*multiplications-172*assignments
 ! getDuDxxyyyy2 operation count     : 482*multiplications+184*assignments+391*additions
 ! getDuDxxyyyy2 optimization savings: 731*additions+2241*multiplications-184*assignments
 ! getDuDxyyyyy2 operation count     : 456*multiplications+188*assignments+384*additions
 ! getDuDxyyyyy2 optimization savings: 1019*additions+2233*multiplications-188*assignments
 ! getDuDyyyyyy2 operation count     : 366*multiplications+133*assignments+258*additions
 ! getDuDyyyyyy2 optimization savings: 345*additions+1335*multiplications-133*assignments


! ****** Dimension 3 ******
 ! getDuDx3 operation count     : 2*additions+3*multiplications+assignments
 ! getDuDx3 optimization savings: -assignments
 ! getDuDy3 operation count     : 2*additions+3*multiplications+assignments
 ! getDuDy3 optimization savings: -assignments
 ! getDuDz3 operation count     : 2*additions+3*multiplications+assignments
 ! getDuDz3 optimization savings: -assignments
 ! getDuDxx3 operation count     : 18*multiplications+4*assignments+8*additions
 ! getDuDxx3 optimization savings: -4*assignments
 ! getDuDxy3 operation count     : 11*additions+18*multiplications+assignments
 ! getDuDxy3 optimization savings: -assignments
 ! getDuDyy3 operation count     : 18*multiplications+4*assignments+8*additions
 ! getDuDyy3 optimization savings: -4*assignments
 ! getDuDxz3 operation count     : 11*additions+18*multiplications+assignments
 ! getDuDxz3 optimization savings: -assignments
 ! getDuDyz3 operation count     : 11*additions+18*multiplications+assignments
 ! getDuDyz3 optimization savings: -assignments
 ! getDuDzz3 operation count     : 18*multiplications+4*assignments+8*additions
 ! getDuDzz3 optimization savings: -4*assignments
 ! getDuDxxx3 operation count     : 58*multiplications+4*assignments+21*additions
 ! getDuDxxx3 optimization savings: 6*multiplications-4*assignments
 ! getDuDxxy3 operation count     : 82*multiplications+7*assignments+38*additions
 ! getDuDxxy3 optimization savings: 9*multiplications-7*assignments
 ! getDuDxyy3 operation count     : 76*multiplications+10*assignments+41*additions
 ! getDuDxyy3 optimization savings: 6*additions+15*multiplications-10*assignments
 ! getDuDyyy3 operation count     : 58*multiplications+4*assignments+21*additions
 ! getDuDyyy3 optimization savings: 6*multiplications-4*assignments
 ! getDuDxxz3 operation count     : 82*multiplications+7*assignments+38*additions
 ! getDuDxxz3 optimization savings: 9*multiplications-7*assignments
 ! getDuDxyz3 operation count     : 50*additions+79*multiplications+4*assignments
 ! getDuDxyz3 optimization savings: 6*additions+12*multiplications-4*assignments
 ! getDuDyyz3 operation count     : 82*multiplications+7*assignments+38*additions
 ! getDuDyyz3 optimization savings: 9*multiplications-7*assignments
 ! getDuDxzz3 operation count     : 76*multiplications+10*assignments+41*additions
 ! getDuDxzz3 optimization savings: 6*additions+15*multiplications-10*assignments
 ! getDuDyzz3 operation count     : 76*multiplications+10*assignments+41*additions
 ! getDuDyzz3 optimization savings: 6*additions+15*multiplications-10*assignments
 ! getDuDzzz3 operation count     : 58*multiplications+4*assignments+21*additions
 ! getDuDzzz3 optimization savings: 6*multiplications-4*assignments
 ! getDuDxxxx3 operation count     : 161*multiplications+29*assignments+71*additions
 ! getDuDxxxx3 optimization savings: 15*additions+98*multiplications-29*assignments
 ! getDuDxxxy3 operation count     : 246*multiplications+28*assignments+110*additions
 ! getDuDxxxy3 optimization savings: 12*additions+109*multiplications-28*assignments
 ! getDuDxxyy3 operation count     : 280*multiplications+43*assignments+142*additions
 ! getDuDxxyy3 optimization savings: 46*additions+192*multiplications-43*assignments
 ! getDuDxyyy3 operation count     : 235*multiplications+49*assignments+148*additions
 ! getDuDxyyy3 optimization savings: 97*additions+225*multiplications-49*assignments
 ! getDuDyyyy3 operation count     : 161*multiplications+29*assignments+71*additions
 ! getDuDyyyy3 optimization savings: 15*additions+98*multiplications-29*assignments
 ! getDuDxxxz3 operation count     : 247*multiplications+27*assignments+110*additions
 ! getDuDxxxz3 optimization savings: 12*additions+108*multiplications-27*assignments
 ! getDuDxxyz3 operation count     : 292*multiplications+31*assignments+154*additions
 ! getDuDxxyz3 optimization savings: 46*additions+186*multiplications-31*assignments
 ! getDuDxyyz3 operation count     : 271*multiplications+166*additions+31*assignments
 ! getDuDxyyz3 optimization savings: 97*additions+207*multiplications-31*assignments
 ! getDuDyyyz3 operation count     : 247*multiplications+27*assignments+110*additions
 ! getDuDyyyz3 optimization savings: 12*additions+108*multiplications-27*assignments
 ! getDuDxxzz3 operation count     : 277*multiplications+46*assignments+142*additions
 ! getDuDxxzz3 optimization savings: 46*additions+195*multiplications-46*assignments
 ! getDuDxyzz3 operation count     : 256*multiplications+52*assignments+175*additions
 ! getDuDxyzz3 optimization savings: 115*additions+222*multiplications-52*assignments
 ! getDuDyyzz3 operation count     : 277*multiplications+46*assignments+142*additions
 ! getDuDyyzz3 optimization savings: 46*additions+195*multiplications-46*assignments
 ! getDuDxzzz3 operation count     : 235*multiplications+49*assignments+148*additions
 ! getDuDxzzz3 optimization savings: 97*additions+225*multiplications-49*assignments
 ! getDuDyzzz3 operation count     : 235*multiplications+49*assignments+148*additions
 ! getDuDyzzz3 optimization savings: 97*additions+225*multiplications-49*assignments
 ! getDuDzzzz3 operation count     : 161*multiplications+29*assignments+71*additions
 ! getDuDzzzz3 optimization savings: 15*additions+98*multiplications-29*assignments
 ! getDuDxxxxx3 operation count     : 480*multiplications+117*assignments+279*additions
 ! getDuDxxxxx3 optimization savings: 239*additions+1000*multiplications-117*assignments
 ! getDuDxxxxy3 operation count     : 644*multiplications+120*assignments+345*additions
 ! getDuDxxxxy3 optimization savings: 266*additions+1103*multiplications-120*assignments
 ! getDuDxxxyy3 operation count     : 402*additions+732*multiplications+150*assignments
 ! getDuDxxxyy3 optimization savings: 311*additions+1309*multiplications-150*assignments
 ! getDuDxxyyy3 operation count     : 685*multiplications+177*assignments+458*additions
 ! getDuDxxyyy3 optimization savings: 576*additions+1896*multiplications-177*assignments
 ! getDuDxyyyy3 operation count     : 619*multiplications+186*assignments+461*additions
 ! getDuDxyyyy3 optimization savings: 924*additions+1926*multiplications-186*assignments
 ! getDuDyyyyy3 operation count     : 479*multiplications+117*assignments+279*additions
 ! getDuDyyyyy3 optimization savings: 239*additions+1001*multiplications-117*assignments
 ! getDuDxxxxz3 operation count     : 646*multiplications+120*assignments+345*additions
 ! getDuDxxxxz3 optimization savings: 266*additions+1101*multiplications-120*assignments
 ! getDuDxxxyz3 operation count     : 833*multiplications+116*assignments+441*additions
 ! getDuDxxxyz3 optimization savings: 323*additions+1298*multiplications-116*assignments
 ! getDuDxxyyz3 operation count     : 889*multiplications+138*assignments+512*additions
 ! getDuDxxyyz3 optimization savings: 594*additions+1827*multiplications-138*assignments
 ! getDuDxyyyz3 operation count     : 754*multiplications+156*assignments+536*additions
 ! getDuDxyyyz3 optimization savings: 984*additions+1908*multiplications-156*assignments
 ! getDuDyyyyz3 operation count     : 648*multiplications+119*assignments+345*additions
 ! getDuDyyyyz3 optimization savings: 266*additions+1099*multiplications-119*assignments
 ! getDuDxxxzz3 operation count     : 730*multiplications+149*assignments+402*additions
 ! getDuDxxxzz3 optimization savings: 311*additions+1311*multiplications-149*assignments
 ! getDuDxxyzz3 operation count     : 754*multiplications+183*assignments+500*additions
 ! getDuDxxyzz3 optimization savings: 606*additions+1899*multiplications-183*assignments
 ! getDuDxyyzz3 operation count     : 524*additions+727*multiplications+174*assignments
 ! getDuDxyyzz3 optimization savings: 978*additions+1953*multiplications-174*assignments
 ! getDuDyyyzz3 operation count     : 730*multiplications+149*assignments+402*additions
 ! getDuDyyyzz3 optimization savings: 311*additions+1311*multiplications-149*assignments
 ! getDuDxxzzz3 operation count     : 679*multiplications+183*assignments+458*additions
 ! getDuDxxzzz3 optimization savings: 576*additions+1902*multiplications-183*assignments
 ! getDuDxyzzz3 operation count     : 658*multiplications+195*assignments+533*additions
 ! getDuDxyzzz3 optimization savings: 1095*additions+2004*multiplications-195*assignments
 ! getDuDyyzzz3 operation count     : 679*multiplications+183*assignments+458*additions
 ! getDuDyyzzz3 optimization savings: 576*additions+1902*multiplications-183*assignments
 ! getDuDxzzzz3 operation count     : 619*multiplications+186*assignments+461*additions
 ! getDuDxzzzz3 optimization savings: 924*additions+1926*multiplications-186*assignments
 ! getDuDyzzzz3 operation count     : 619*multiplications+186*assignments+461*additions
 ! getDuDyzzzz3 optimization savings: 924*additions+1926*multiplications-186*assignments
 ! getDuDzzzzz3 operation count     : 481*multiplications+115*assignments+279*additions
 ! getDuDzzzzz3 optimization savings: 239*additions+999*multiplications-115*assignments
 ! getDuDxxxxxx3 operation count     : 1232*multiplications+380*assignments+916*additions
 ! getDuDxxxxxx3 optimization savings: 2414*additions+8130*multiplications-380*assignments
 ! getDuDxxxxxy3 operation count     : 1742*multiplications+419*assignments+1142*additions
 ! getDuDxxxxxy3 optimization savings: 2790*additions+9123*multiplications-419*assignments
 ! getDuDxxxxyy3 operation count     : 1728*multiplications+492*assignments+1158*additions
 ! getDuDxxxxyy3 optimization savings: 2837*additions+9455*multiplications-492*assignments
 ! getDuDxxxyyy3 operation count     : 1194*additions+1600*multiplications+481*assignments
 ! getDuDxxxyyy3 optimization savings: 3128*additions+10552*multiplications-481*assignments
 ! getDuDxxyyyy3 operation count     : 1265*additions+1496*multiplications+508*assignments
 ! getDuDxxyyyy3 optimization savings: 4761*additions+13467*multiplications-508*assignments
 ! getDuDxyyyyy3 operation count     : 1421*multiplications+505*assignments+1253*additions
 ! getDuDxyyyyy3 optimization savings: 6954*additions+13503*multiplications-505*assignments
 ! getDuDyyyyyy3 operation count     : 1230*multiplications+381*assignments+916*additions
 ! getDuDyyyyyy3 optimization savings: 2414*additions+8132*multiplications-381*assignments
 ! getDuDxxxxxz3 operation count     : 1741*multiplications+417*assignments+1142*additions
 ! getDuDxxxxxz3 optimization savings: 2790*additions+9124*multiplications-417*assignments
 ! getDuDxxxxyz3 operation count     : 2098*multiplications+378*assignments+1289*additions
 ! getDuDxxxxyz3 optimization savings: 3120*additions+10006*multiplications-378*assignments
 ! getDuDxxxyyz3 operation count     : 2226*multiplications+460*assignments+1406*additions
 ! getDuDxxxyyz3 optimization savings: 3465*additions+11159*multiplications-460*assignments
 ! getDuDxxyyyz3 operation count     : 2063*multiplications+503*assignments+1531*additions
 ! getDuDxxyyyz3 optimization savings: 5227*additions+14391*multiplications-503*assignments
 ! getDuDxyyyyz3 operation count     : 1552*additions+1886*multiplications+521*assignments
 ! getDuDxyyyyz3 optimization savings: 7924*additions+14415*multiplications-521*assignments
 ! getDuDyyyyyz3 operation count     : 1745*multiplications+412*assignments+1142*additions
 ! getDuDyyyyyz3 optimization savings: 2790*additions+9120*multiplications-412*assignments
 ! getDuDxxxxzz3 operation count     : 1733*multiplications+489*assignments+1158*additions
 ! getDuDxxxxzz3 optimization savings: 2837*additions+9450*multiplications-489*assignments
 ! getDuDxxxyzz3 operation count     : 1971*multiplications+514*assignments+1338*additions
 ! getDuDxxxyzz3 optimization savings: 3398*additions+11018*multiplications-514*assignments
 ! getDuDxxyyzz3 operation count     : 2018*multiplications+554*assignments+1463*additions
 ! getDuDxxyyzz3 optimization savings: 5157*additions+14124*multiplications-554*assignments
 ! getDuDxyyyzz3 operation count     : 1502*additions+1835*multiplications+572*assignments
 ! getDuDxyyyzz3 optimization savings: 7830*additions+14196*multiplications-572*assignments
 ! getDuDyyyyzz3 operation count     : 1735*multiplications+490*assignments+1158*additions
 ! getDuDyyyyzz3 optimization savings: 2837*additions+9448*multiplications-490*assignments
 ! getDuDxxxzzz3 operation count     : 1595*multiplications+477*assignments+1194*additions
 ! getDuDxxxzzz3 optimization savings: 3128*additions+10557*multiplications-477*assignments
 ! getDuDxxyzzz3 operation count     : 1631*multiplications+526*assignments+1373*additions
 ! getDuDxxyzzz3 optimization savings: 5079*additions+13824*multiplications-526*assignments
 ! getDuDxyyzzz3 operation count     : 1415*additions+1610*multiplications+511*assignments
 ! getDuDxyyzzz3 optimization savings: 7512*additions+14169*multiplications-511*assignments
 ! getDuDyyyzzz3 operation count     : 1595*multiplications+477*assignments+1194*additions
 ! getDuDyyyzzz3 optimization savings: 3128*additions+10557*multiplications-477*assignments
 ! getDuDxxzzzz3 operation count     : 1487*multiplications+517*assignments+1265*additions
 ! getDuDxxzzzz3 optimization savings: 4761*additions+13476*multiplications-517*assignments
 ! getDuDxyzzzz3 operation count     : 1478*multiplications+538*assignments+1406*additions
 ! getDuDxyzzzz3 optimization savings: 8178*additions+14193*multiplications-538*assignments
 ! getDuDyyzzzz3 operation count     : 1487*multiplications+517*assignments+1265*additions
 ! getDuDyyzzzz3 optimization savings: 4761*additions+13476*multiplications-517*assignments
 ! getDuDxzzzzz3 operation count     : 1421*multiplications+505*assignments+1253*additions
 ! getDuDxzzzzz3 optimization savings: 6954*additions+13503*multiplications-505*assignments
 ! getDuDyzzzzz3 operation count     : 1421*multiplications+505*assignments+1253*additions
 ! getDuDyzzzzz3 optimization savings: 6954*additions+13503*multiplications-505*assignments
 ! getDuDzzzzzz3 operation count     : 1232*multiplications+382*assignments+916*additions
 ! getDuDzzzzzz3 optimization savings: 2414*additions+8130*multiplications-382*assignments

! Define 
!    defineParametricDerivativeMacros(u,dr,dx,DIM,ORDER,COMPONENTS,MAXDERIV)
!       defines -> ur2, us2, ux2, uy2, ...            (2D)
!                  ur3, us3, ut3, ux3, uy3, uz3, ...  (3D)
! This file was generated by weights.maple


! This next macro will evaluate parametric derivatives and save in temporaries
!   u is the variable name, v is the prefix for the temporaries, e.g.
!   For example, lines of the following form will be generated:
!      v = u(i1,i2,i3) 
!      vr = ur4(i1,i2,i3) 

! This next macro will evaluate parametric derivatives and save in temporaries
!   u is the variable name, v is the prefix for the temporaries, e.g.
!   For example, lines of the following form will be generated:
!      v = u(i1,i2,i3) 
!      vr = ur4(i1,i2,i3) 

! This next macro will evaluate parametric derivatives and save in temporaries
!   u is the variable name, v is the prefix for the temporaries, e.g.
!   For example, lines of the following form will be generated:
!      v = u(i1,i2,i3) 
!      vr = ur4(i1,i2,i3) 

! This next macro will evaluate x,y,z derivatives using temporaries already computed 
!   u1 is the variable name, aj the jaocbian name and v is the prefix for the temporaries
!   For example, lines of the following form will be generated:
!      getDuDx2(u1,aj,vx) 
!      getDuDxy2(u1,aj,vxy) 
!      getDuDxxx2(u1,aj,vxxx) 

! This next macro will evaluate x,y,z derivatives using temporaries already computed 
!   u1 is the variable name, aj the jaocbian name and v is the prefix for the temporaries
!   For example, lines of the following form will be generated:
!      getDuDx2(u1,aj,vx) 
!      getDuDxy2(u1,aj,vxy) 
!      getDuDxxx2(u1,aj,vxxx) 

! This next macro will evaluate x,y,z derivatives using temporaries already computed 
!   u1 is the variable name, aj the jaocbian name and v is the prefix for the temporaries
!   For example, lines of the following form will be generated:
!      getDuDx2(u1,aj,vx) 
!      getDuDxy2(u1,aj,vxy) 
!      getDuDxxx2(u1,aj,vxxx) 

! This next macro will evaluate x,y,z derivatives of the jacobian 

! u = jacobian name (rsxy), v=prefix for derivatives: vrxr, vrys, 

! defineParametricDerivativeMacros(u,dr,dx,DIM,ORDER,COMPONENTS,MAXDERIV)
! 2D, order=6, components=1
! defineParametricDerivativeMacros(u,dr,dx,DIM,ORDER,COMPONENTS,MAXDERIV)

!- defineParametricDerivativeMacros(rsxy1,dr1,dx1,3,2,2,4)
!- defineParametricDerivativeMacros(rsxy1,dr1,dx1,3,4,2,2)
!- defineParametricDerivativeMacros(u1,dr1,dx1,3,2,1,4)
!- defineParametricDerivativeMacros(u1,dr1,dx1,3,4,1,2)
!-
!- defineParametricDerivativeMacros(rsxy2,dr2,dx2,3,2,2,4)
!- defineParametricDerivativeMacros(rsxy2,dr2,dx2,3,4,2,2)
!- defineParametricDerivativeMacros(u2,dr2,dx2,3,2,1,4)
!- defineParametricDerivativeMacros(u2,dr2,dx2,3,4,1,2)

 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************

 ! We need up to 2 derivatives of u1n to order 2 -- do order 4 for TZ
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************

 ! q1 for nonlinear dispersive MLA
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************

 ! We need up to 2 derivatives of p1 to order 2-- do order 4 for TZ
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************

 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************

 ! We need up to 2 derivatives of u1n to order 2-- do order 4 for TZ
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************

 ! q2 for nonlinear dispersive MLA
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************

 ! We need up to 2 derivatives of p2 to order 2-- do order 4 for TZ
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************
 ! *************** 0 components *************
 ! *************** 1 components *************
 ! *************** 2 components *************

! ===========================================================================================
! Macro: Output some debug info for the first few time-steps 
! ===========================================================================================


! ******************************************************************************************************************
! ************* These are altered version of those from insImp.h ***************************************************
! ******************************************************************************************************************


! ==========================================================================================
!  Evaluate the Jacobian and its derivatives (parametric and spatial). 
!    rsxy   : jacobian matrix name 
!    aj     : prefix for the name of the resulting jacobian variables, 
!             e.g. ajrx, ajsy, ajrxx, ajsxy, ...
!    MAXDER : number of derivatives to evaluate.  
! ==========================================================================================

! ==========================================================================================
!  Evaluate the parametric derivatives of u.
!    u      : evaluate derivatives of this function.
!    uc     : component to evaluate
!    uu     : prefix for the name of the resulting derivatives, e.g. uur, uus, uurr, ...
!    MAXDER : number of derivatives to evaluate.  
! ==========================================================================================


! ==========================================================================================
!  Evaluate a derivative. (assumes parametric derivatives have already been evaluated)
!   DERIV   : name of the derivative. One of 
!                x,y,z,xx,xy,xz,...
!    u      : evaluate derivatives of this function.
!    uc     : component to evaluate
!    uu     : prefix for the name of the resulting derivatives (same name used with opEvalParametricDerivative) 
!    aj     : prefix for the name of the jacobian variables.
!    ud     : derivative is assigned to this variable.
! ==========================================================================================

! ******************************************************************************************************************

! loop over the boundary points



! loop over the boundary points

! loop over the boundary points with a mask. 
! Assign pts where both mask1 and mask2 are discretization pts.
! If mask1>0 and mask2<0 then we just leave the extrapolated values in u1 and u2 .

! loop over the boundary points that includes ghost points in the tangential direction

! loop over the boundary points that includes ghost points in the tangential direction.
! Assign pts where both mask1 and mask2 are discretization pts.
! If mask1>0 and mask2<0 then we just leave the extrapolated values in u1 and u2 .


! Assign pts where both mask1 and mask2 are discretization pts.
! If mask1>0 and mask2<0 then we just leave the extrapolated values in u1 and u2 .








! *********************************************************************
! ********** MACROS FOR DISPERSIVE INTERFACE CONDITIONS ***************
! *********************************************************************
!         -*- mode: F90 -*-
! *********************************************************************
! ********** MACROS FOR DISPERSIVE INTERFACE CONDITIONS ***************
!    This file is included into interface3d.bf 
! *********************************************************************


! -------------------------------------------------------------------------
! *NEW* VERSION July 16, 2019 -- use precomputed coefficients - see paper adegdmi.pdf 
!
! Macro: Evaluate DISPERSIVE forcing terms, 2nd-order accuracy 
!   This macro can be used to eval values in either domain 1 or domain 2
!
! Input:
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : This is the value of P.tt without the term involving L(E) = c^2*Delta(E)
!   beta = 1 - alphaP*Sum_k{ C_k }
! ------------------------------------------------------------------------



! -------------------------------------------------------------------------
! OLD VERSION
! Macro: Evaluate DISPERSIVE forcing terms, 2nd-order accuracy 
!   This macro can be used to eval values in either domain 1 or domain 2
!
! Input:
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : 
!   beta = 1 - alphaP*Sum_k{ C_k }
! ------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------
! Macro: Eval twilight-zone forcing for GDM equations
! Output
!  fpv(n,jv) : RHS To Pv_{n,jv} equation 
!  fev(n)    : RHS to E_{n} equation
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 2D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Eval twilight-zone forcing for GDM equations THREE-D
! Output
!  fpv(n,jv) : RHS To Pv_{n,jv} equation 
!  fev(n)    : RHS to E_{n} equation
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 3D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------
! Macro: Assign DISPERSIVE interface ghost values, DIM=2, ORDER=2, GRID=Rectangular
! 
! Here are the jump conditions (See notes in DMX_ADE)
!   [ u.x + v.y ] = 0
!   [ (1/mu)* tv,.( curl(E) ) ]
!   [ tv.( c^2*Delta(E) -alphaP*P_tt) ] = 0  --> [ tv.( beta*c^2*Delta(E) - alphaP* F) ]=0 
!   [ (1/mu)* nv.( Delta(E) ) ]=0
! 
! -------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------
! Macro:  Evaluate the RHS to the jump conditons: 2D, Order=2, Dispersive
!
! --------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------


! --------------------------------------------------------------------
! Macro: Assign  DISPERSIVE interface ghost values, DIM=2, ORDER=2, GRID=Curvilinear
! 
! Here are the jump conditions (See notes in DMX_ADE)
!   [ u.x + v.y ] = 0
!   [ (1/mu)* tv,.( curl(E) ) ]
!   [ tv.( c^2*Delta(E) -alphaP*P_tt) ] = 0  --> [ tv.( beta*c^2*Delta(E) - alphaP* F) ]=0 
!   [ (1/mu)* nv.( Delta(E) ) ]=0
! 
! -------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------
! Macro:  Evaluate the RHS to the jump conditons: 2D, Order=2, Dispersive
!
! Here are the jump conditions (See notes in DMX_ADE)
!   (1) [ div(E) ] = 0
!   (2) [ (1/mu)* nv.( Delta(E) ) ]=0
!   (3) [ (1/mu)* tv.( curl(E) ) ]=0               -->   [ (1/mu)* \nv\times( curl(E) ) ]=0
!   (4) [ tv.(c^2*Delta(E) -alphaP*P_tt) ] = 0    -->   [ \nv X ( c^2*Delta(E) -alphaP*P_tt) ] = 0 
! 
! These 6 equations can be written as 
!   [ div(E) n + (I- n n^T)( curl(E)/mu ) ] =0                                 (3 eqns)
!   [ (1/mu) n n^T Delta(E) + (I-n n^T) ( c^2*Delta(E) -alphaP*P_tt ) ] = 0    (3 eqns)
!
! --------------------------------------------------------------------------------------------



! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=2, GRID=Curvilinear
! 
!                  DISPERSIVE CASE
!
! Here are the jump conditions (See notes in DMX_ADE)
!   (1) [ div(E) ] = 0
!   (2) [ (1/mu)* nv.( Delta(E) ) ]=0
!   (3) [ (1/mu)* tv.( curl(E) ) ]=0               -->   [ (1/mu)* \nv\times( curl(E) ) ]=0
!   (4) [ tv.(c^2*Delta(E) -alphaP*P_tt) ] = 0    -->   [ \nv X ( c^2*Delta(E) -alphaP*P_tt) ] = 0 
! 
! These 6 equations can be written as 
!   [ div(E) n + (I- n n^T)( curl(E)/mu ) ] =0                                 (3 eqns)
!   [ (1/mu) n n^T Delta(E) + (I-n n^T) ( c^2*Delta(E) -alphaP*P_tt ) ] = 0    (3 eqns)
!
! An approximation to P_tt takes the form 
!   P_tt = K Delta(E) + G(E,P)
! --------------------------------------------------------------------------



! -------------------------------------------------------------------------
! Macro: Evaluate DISPERSIVE forcing terms, FOURTH-ORDER
!   This macro can be usedto eval values in either domain 1 or domain 2
!
! Input:
!   FACE : LEFT or : RIGHT
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : 
!   c2PttLEsum   : coeff of L(E) in P.tt     (second-order)
!   c4PttLEsum   : coeff of L(E) in P.tt     (fourth-order)
!   c4PttLLEsum  : coeff of L*L(E) in P.tt 
! ------------------------------------------------------------------------


! -------------------------------------------------------------------------------
! Macro: Evaluate the TZ forcings GDM FOURTH-ORDER
! -------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 2D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------
! Macro: Evaluate the GDM jump conditions in 2D, order=4
! --------------------------------------------------------------------------

! --------------------------------------------------------------------------
! Macro: Evaluate the GDM jump conditions in 2D, order=4
! Version prior to Feb 1, 2021
! --------------------------------------------------------------------------

 

! --------------------------------------------------------------------------
! Macro: Evaluate the GDM jump conditions in 2D, order=4
! Version prior to Feb 1, 2021
! --------------------------------------------------------------------------

! ========================================================================
! Macro: Getting forcing for GDM
!  Input:
!    ec : E component
!    pc : P component
!  Output:
!    fe : forcing for E ( includes factor of dt^2)
!    fpv(iv) : forcing for polarization vector iv=0,1,2,...( includes factor of dt^2)
! ========================================================================




! ===========================================================================================
! Macro:     DISPERSIVE: CURVILINEAR, 2D, ORDER 2
!    Evaluate P to second-order
! Input:
!    i1,i2,i3 : evaluate P at this point
!    u,um : u(t-dt), u(t-2*dt)
! Output:
!    evals,pv
! ===========================================================================================




! =============================================================================================
!   Evaluate the jump conditions for the GDM interface equations
! =============================================================================================





! ===========================================================================================
!  Assign P in ghost points for the fourth-order method
!
! *** THIS IS NOT USED CURRENTLY -- was created to fix a bug that was caused by a wrong alphaP
!
! ===========================================================================================

! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=2, ORDER=4, GRID=Curvilinear
!         DISPERSIVE CASE -- GDM 
! --------------------------------------------------------------------------

! ********************************************************************************
!       INTERFACE MACROS (used in interface3d.bf and interface3dOrder4.bf)
! ********************************************************************************
! -*- mode: F90; -*-
! ==================================================================================================
! Macro: compute the normal vector (an1,an2) in 2D
! GRIDTYPE (input) : curvilinear or rectangular
! ==================================================================================================

! ==================================================================================================
! Macro: compute the normal vector (an1,an2,an3) in 3D
! GRIDTYPE (input) : curvilinear or rectangular
! ==================================================================================================

! ==================================================================================================
! This macro will assign the jump conditions on the boundary
! DIM (input): number of dimensions (2 or 3)
! GRIDTYPE (input) : curvilinear or rectangular
! ==================================================================================================

! ===========================================================================================
!  Smooth P on the interface 
! ===========================================================================================

! ----------------------------------------------------------------------------------
!  Macro:
!    --- perturb a component of E by delta
! ----------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------
!  Macro:
!    --- check matrix coefficients by delta function approach ----
! Input:
!   am : matrix of coefficients is stored here 
!   evalInterfaceEquations : macro that evaluates the interface equations
! ----------------------------------------------------------------------------------



!--------------------
! get coefficients using delta function approach
!--------------------



! ----------------------------------------------------------------------------------
!  Macro:
!    --- eval matrix coefficients by delta function approach ----
! Input:
!   am : matrix of coefficients is stored here 
!   evalInterfaceEquations : macro that evaluates the interface equations
! ----------------------------------------------------------------------------------



! ----------------------------------------------------------------------------------
!  Macro:
!    --- save matrix coefficients by delta function approach ----
! Input:
!   am : matrix of coefficients is stored here 
!   evalInterfaceEquations : macro that evaluates the interface equations
! ----------------------------------------------------------------------------------


  

! *******************************************************************************
! ********** MACROS FOR NONLINEAR DISPERSIVE INTERFACE CONDITIONS ***************
! *******************************************************************************
!         -*- mode: F90 -*-
! *********************************************************************
! ********** MACROS FOR NONLINEAR INTERFACE CONDITIONS ***************
!    This file is included into interfaceOpt.bf90 
! *********************************************************************

!  +++++++ THIS FILE STARTED AS A COPY of dispersiveInterfaceMacros ++++++++
!    BE SURE TO RE-NAME MACROS AS WHEN THEY ARE CHANGED 


! -------------------------------------------------------------------------
! 
! Macro: Evaluate Nonlinear DISPERSIVE forcing terms, 2nd-order accuracy 
!   This macro can be used to eval values in either domain 1 or domain 2
!
! Input:
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : This is the value of P.tt without the term involving L(E) = c^2*Delta(E)
!   beta = 1 - alphaP*Sum_k{ C_k }
! ------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Eval twilight-zone forcing for Maxwell-Bloch equations
! Output
!  fpv(n,jv) : RHS To Pv_{n,jv} equation s
!  fev(n)    : RHS to E_{n} equation
!  fnv(n)    : RHS to N_{l} equations
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 2D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Eval twilight-zone forcing for GDM equations THREE-D
! Output
!  fpv(n,jv) : RHS To Pv_{n,jv} equation 
!  fev(n)    : RHS to E_{n} equation
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 3D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------
! Macro: Assign nonlinear DISPERSIVE interface ghost values, DIM=2, ORDER=2, GRID=Rectangular
! 
! Here are the jump conditions (See notes in DMX_ADE)
!   [ u.x + v.y ] = 0
!   [ (1/mu)* tv,.( curl(E) ) ]
!   [ tv.( c^2*Delta(E) -alphaP*P_tt) ] = 0  --> [ tv.( beta*c^2*Delta(E) - alphaP* F) ]=0 
!   [ (1/mu)* nv.( Delta(E) ) ]=0
! 
! -------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------
! Macro:  Evaluate the RHS to the jump conditions: 2D, Order=2, nonlinear Dispersive
!
! --------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------


! --------------------------------------------------------------------
! Macro: Assign NONLINEAR interface ghost values, DIM=2, ORDER=2, GRID=Curvilinear
! 
! Here are the jump conditions (See notes in DMX_ADE)
!   [ u.x + v.y ] = 0
!   [ (1/mu)* tv,.( curl(E) ) ]
!   [ tv.( c^2*Delta(E) -alphaP*P_tt) ] = 0  --> [ tv.( beta*c^2*Delta(E) - alphaP* F) ]=0 
!   [ (1/mu)* nv.( Delta(E) ) ]=0
! 
! -------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------
! Macro:  Evaluate the RHS to the jump conditions: 2D, Order=2, Dispersive
!
! Here are the jump conditions (See notes in DMX_ADE)
!   (1) [ div(E) ] = 0
!   (2) [ (1/mu)* nv.( Delta(E) ) ]=0
!   (3) [ (1/mu)* tv.( curl(E) ) ]=0               -->   [ (1/mu)* \nv\times( curl(E) ) ]=0
!   (4) [ tv.(c^2*Delta(E) -alphaP*P_tt) ] = 0    -->   [ \nv X ( c^2*Delta(E) -alphaP*P_tt) ] = 0 
! 
! These 6 equations can be written as 
!   [ div(E) n + (I- n n^T)( curl(E)/mu ) ] =0                                 (3 eqns)
!   [ (1/mu) n n^T Delta(E) + (I-n n^T) ( c^2*Delta(E) -alphaP*P_tt ) ] = 0    (3 eqns)
!
! --------------------------------------------------------------------------------------------



! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=2, GRID=Curvilinear
! 
!                  Nonlinear DISPERSIVE CASE
!
! Here are the jump conditions (See notes in DMX_ADE)
!   (1) [ div(E) ] = 0
!   (2) [ (1/mu)* nv.( Delta(E) ) ]=0
!   (3) [ (1/mu)* tv.( curl(E) ) ]=0               -->   [ (1/mu)* \nv\times( curl(E) ) ]=0
!   (4) [ tv.(c^2*Delta(E) -alphaP*P_tt) ] = 0    -->   [ \nv X ( c^2*Delta(E) -alphaP*P_tt) ] = 0 
! 
! These 6 equations can be written as 
!   [ div(E) n + (I- n n^T)( curl(E)/mu ) ] =0                                 (3 eqns)
!   [ (1/mu) n n^T Delta(E) + (I-n n^T) ( c^2*Delta(E) -alphaP*P_tt ) ] = 0    (3 eqns)
!
! An approximation to P_tt takes the form 
!   P_tt = K Delta(E) + G(E,P)
! --------------------------------------------------------------------------



! -------------------------------------------------------------------------
! Macro: Evaluate Nonlinear DISPERSIVE forcing terms, FOURTH-ORDER
!   This macro can be used to eval values in either domain 1 or domain 2
!   At this point, the first ghost lines are filled with second order accurate E (only) field
! Input:
!   FACE : LEFT or : RIGHT
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output: dispersive forcings
! ------------------------------------------------------------------------


! -------------------------------------------------------------------------------
! Macro: Evaluate the TZ forcings Nonlinear Dispersive MLA for FOURTH-ORDER codes
! -------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for nonlinear dispersive equations (4th order)
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------
! Macro: Evaluate the GDM jump conditions in 2D, order=4
! --------------------------------------------------------------------------

! ========================================================================
! Macro: Getting forcing for GDM
!  Input:
!    ec : E component
!    pc : P component
!  Output:
!    fe : forcing for E ( includes factor of dt^2)
!    fpv(jv) : forcing for polarization vector jv=0,1,2,...( includes factor of dt^2)
! ========================================================================




! ===========================================================================================
! Macro:     DISPERSIVE: CURVILINEAR, 2D, ORDER 2
!    Evaluate P to second-order
! Input:
!    i1,i2,i3 : evaluate P at this point
!    u,um : u(t-dt), u(t-2*dt)
! Output:
!    evals,pv
! ===========================================================================================




! =============================================================================================
!   Evaluate the jump conditions for the nonlinear interface equations
! =============================================================================================

! ===========================================================================================
!  Assign P in ghost points for the fourth-order method
!
! *** THIS IS NOT USED CURRENTLY -- was created to fix a bug that was caused by a wrong alphaP
!
! ===========================================================================================

! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------

! macro for dispersive forcing


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=2, ORDER=4, GRID=Curvilinear
!         NONLINEAR DISPERSIVE CASE -- MLA 
! --------------------------------------------------------------------------
! -------------------------------------------------------------------------
! Macro: Evaluate DISPERSIVE forcing terms, FOURTH-ORDER AND 3D
!   This macro can be usedto eval values in either domain 1 or domain 2
!
! Input:
!   FACE : LEFT or : RIGHT
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : 
!   c2PttLEsum   : coeff of L(E) in P.tt     (second-order)
!   c4PttLEsum   : coeff of L(E) in P.tt     (fourth-order)
!   c4PttLLEsum  : coeff of L*L(E) in P.tt 
! ------------------------------------------------------------------------

! -------------------------------------------------------------------------------
! Macro: Evaluate the TZ forcings GDM FOURTH-ORDER 3D
! -------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 3D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------
! Macro: Evaluate the MLA jump conditions in 3D, order=4
! --------------------------------------------------------------------------


! ==========================================================================================
!   Evaluate the jump conditions (including compatibility) for the MLA interface equations 
! ==========================================================================================

 ! ==========================================================================================
!         Fill in FIRST SET of 6 interface equations -- GDM/MLA -- .
! 
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers (0,1,2,3,4,5) or (6,7,8,9,10,11)
! ==========================================================================================


! ==========================================================================================
!         Fill in SECOND SET of 6 interface equations -- GDM/MLA -- .
! 
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers  (6,7,8,9,10,11)
! ==========================================================================================




! ------------------ OLD ------------------------------
! This macro will assign the jump conditions on the boundary
! DIM (input): number of dimensions (2 or 3)
! GRIDTYPE (input) : curvilinear or rectangular


! ********************************************************************************
!     Usage: setJacobianRS( aj1, r, s)
!            setJacobianRS( aj1, s, r)
! ********************************************************************************

! ***************************************************************************
! This macro will set the temp variables rx, rxx, ry, ryx, ...
! If axis=0 then
!   rx = ajrx
!   sx = ajsx
!    ...
!  else if axis=1
!    -- permute r <-> s 
!   rx = ajsx
!   sx = ajrx
!    ...
! ***************************************************************************

! ===================================================================================
!  Optimized periodic update: (only applied in serial)
!     update the periodic ghost points used by an interface on the grid face (side,axis)
! ===================================================================================

! ===================================================================================
!  Optimized periodic update:
!     update the periodic ghost points used by an interface on the grid face (side,axis)
! ===================================================================================


! ******************************************************************************
!   This next macro is called by other macros to evaluate the first and second derivatives
!   This macro assumes that opEvalJacobianDerivatives has been called
! ******************************************************************************


! *********************************************************************************
!   Evaluate derivatives for the 2nd-order 3D interface equations
! *********************************************************************************




! XXXXXXXXXXXXXXXXXXXXXXX from interface3d.bf :  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! ******************************************************************************
!   This next macro is called by evalDerivs3dOrder4
!   This macro assumes that opEvalJacobianDerivatives has been called
! ******************************************************************************


! ******************************************************************************
!   Evaluate derivatives for the 4th-order 3D interface equations
! ******************************************************************************


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! ======================================================================================
!  Evaluate the interface jump conditions when using the full centered scheme
! ======================================================================================


! =====================================================================================
! initialization step: assign first two ghost line by extrapolation
! NOTE: assign ghost points outside the ends
! =====================================================================================


! ==========================================================================================
!         Fill in 6 of the interface equations.
! 
!  The equations come in sets of 6 (composed of 2 sets of 3) that are similar in structure so we can fill in
!  6 at a time.
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers (0,1,2,3,4,5) or (6,7,8,9,10,11)
! ==========================================================================================


! =================================================================================
! Assign a given value to a set of variables
! =================================================================================


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=4, GRID=Curvilinear
!     **OLD VERSION THAT EXTRAPOLATES SECOND GHOST LINE ****
! --------------------------------------------------------------------------


 ! here are the macros from deriv.maple (file=derivMacros.h)
 ! defineMacro lapCoeff4a(is,dr,ds) ( (-2/3.*rxx*is-2/3.*ryy*is)/dr+(4/3.*rx**2+4/3.*ry**2)/dr**2 )
 ! defineMacro lapCoeff4b(is,dr,ds) ( (1/12.*rxx*is+1/12.*ryy*is)/dr+(-1/12.*rx**2-1/12.*ry**2)/dr**2 )
 ! -- these are from op/src/derivCoeff.h (generated by op/src/derivCoeff.maple)
 
 
 
 
 ! #defineMacro xLapCoeff3DOrder2b(is,dr,ds) ( -1/2.*rz**2*rx*is/dr**3-1/2.*ry**2*rx*is/dr**3-1/2.*rx**3*is/dr**3 )
 
 
 ! #defineMacro yLapCoeff3DOrder2a(is,dr,ds) ( (2*rxy*rx+ry*rxx)/dr**2+rz**2*ry*is/dr**3+(rz*sz*sy+sz*(sz*ry+rz*sy))*is/dr/ds**2+(2*ty*tx*rx+ry*tx**2)*is/dr/dt**2+(tz*(rz*ty+tz*ry)+rz*tz*ty)*is/dr/dt**2+3*ry*ty**2*is/dr/dt**2+(2*rz*ryz+rzz*ry)/dr**2-1/2.*rxxy*is/dr-1/2.*ryyy*is/dr+3*ry*ryy/dr**2+ry**3*is/dr**3+ry*rx**2*is/dr**3+(ry*sx**2+2*sy*sx*rx)*is/dr/ds**2+3*ry*sy**2*is/dr/ds**2-1/2.*ryzz*is/dr )
 
 
 ! #defineMacro yLapCoeff3DOrder2b(is,dr,ds) ( -1/2.*rz**2*ry*is/dr**3-1/2.*ry**3*is/dr**3-1/2.*ry*rx**2*is/dr**3 )
 
 
 ! #defineMacro zLapCoeff3DOrder2a(is,dr,ds) ( -1/2.*rxxz*is/dr-1/2.*ryyz*is/dr+3*rz*rzz/dr**2+rz**3*is/dr**3+(2*tz*tx*rx+rz*tx**2)*is/dr/dt**2+rz*rx**2*is/dr**3+3*rz*sz**2*is/dr/ds**2-1/2.*rzzz*is/dr+(rz*sx**2+2*sz*sx*rx)*is/dr/ds**2+3*rz*tz**2*is/dr/dt**2+(rz*sy**2+2*sz*sy*ry)*is/dr/ds**2+(2*tz*ty*ry+rz*ty**2)*is/dr/dt**2+rz*ry**2*is/dr**3+(2*rx*rxz+rz*rxx)/dr**2+(2*ry*ryz+rz*ryy)/dr**2 )
 
 
 ! #defineMacro zLapCoeff3DOrder2b(is,dr,ds) ( -1/2.*rz**3*is/dr**3-1/2.*rz*ry**2*is/dr**3-1/2.*rz*rx**2*is/dr**3 )
 


! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------


! ------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=4, GRID=Curvilinear
! ------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------

! -------------------------------------------------------------------------------
! Macro: Evaluate the TZ forcings GDM FOURTH-ORDER 3D
! -------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 3D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------
! Macro: Evaluate the GDM jump conditions in 3D, order=4
! --------------------------------------------------------------------------
 


! -------------------------------------------------------------------------
! Macro: Evaluate DISPERSIVE forcing terms, FOURTH-ORDER AND 3D
!   This macro can be usedto eval values in either domain 1 or domain 2
!
! Input:
!   FACE : LEFT or : RIGHT
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : 
!   c2PttLEsum   : coeff of L(E) in P.tt     (second-order)
!   c4PttLEsum   : coeff of L(E) in P.tt     (fourth-order)
!   c4PttLLEsum  : coeff of L*L(E) in P.tt 
! ------------------------------------------------------------------------



! ==========================================================================================
!   Evaluate the jump conditions (including compatibility) for the GDM interface equations 
! ==========================================================================================


! ==========================================================================================
!         Fill in FIRST SET of 6 interface equations -- GDM -- .
! 
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers (0,1,2,3,4,5) or (6,7,8,9,10,11)
! ==========================================================================================



! ==========================================================================================
!         Fill in SECOND SET of 6 interface equations -- GDM -- .
! 
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers  (6,7,8,9,10,11)
! ==========================================================================================



! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=4, GRID=Curvilinear
!         DISPERSIVE CASE -- GDM 
! --------------------------------------------------------------------------



!  -- this next include file holds macros for assigning interface in 3D order 2 
! eval3dJumpOrder2()
! evalInterfaceEquations23c()
! assignInterfaceGhost23c()
! initializeInterfaceVariablesMacro(LABEL)
! setIndexBoundsExtraGhost()
! resetIndexBounds()  
!         -*- mode: F90 -*-
! --------------------------------------------------------------------------
! ******** define some interface macros for 3D *********
! --------------------------------------------------------------------------


! --------------------------------------------------------------------------
!   Apply interface jump conditions, 3D order 2
! --------------------------------------------------------------------------

! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=2, GRID=Curvilinear
! 
! --------------------------------------------------------------------------


! =====================================================================================================
!  Macro:
!     Initialize various variables and loop bounds
!     Called by interface3d.bf and interface3dOrder4.bf 
! =====================================================================================================




! ==============================================================================
! Macro: Set index bounds to include extra ghost in tangential directions
!   We do this for stage I in the new 4th order scheme (for parallel) where we
!   assign the first ghost point to second-order
! ==============================================================================

! ==============================================================================
! Macro: reset index bounds 
! ==============================================================================
      


!================================================================================================
! Macro  3D ORDER=4 CURVILINEAR
!================================================================================================

        
! ====================================================================================================
!   Main macro to generate separate interface files
!   NAME: name of subroutine
!   DIM,ORDER,GRIDTYPE : number of dimensions, order of accuracy, grid-type
!   MODEL: none, gdm, mla
! ====================================================================================================
! ================================= end interface macro ====================================      






! ***************************************************************************************************
! **************************  Build different interface files ****************************************
! ***************************************************************************************************

!** Note: comment out these next lines to avoid compiling files while developing one version 


! -----------------------------------------------------------------------------------------------
! Macro to call a given interface routine.
!  Do this since the arguments lists are the same for all cases
! -----------------------------------------------------------------------------------------------
                           
! -----------------------------------------------------------------------------------------------
! Macro to call the appropriate interface routine.
! -----------------------------------------------------------------------------------------------

    
! ==============================================================================================================
!  This subroutine calls the appropriate interface routine
! ==============================================================================================================
  subroutine interfaceOpt3dOrder4( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange1, u1,u1n,u1m, wk1, mask1,rsxy1, xy1, p1,p1n,p1m, q1,q1n,q1m, boundaryCondition1, md1a,md1b,md2a,md2b,md3a,md3b,gridIndexRange2, u2,u2n,u2m, wk2, mask2,rsxy2, xy2, p2,p2n,p2m, q2,q2n,q2m, boundaryCondition2, ipar, rpar, aa2,aa4,aa8, ipvt2,ipvt4,ipvt8, ierr )

  implicit real (a-h,o-z)

  integer nd,gridType,orderOfAccuracy
  integer ipar(0:*)

  integer rectangular,curvilinear
  parameter(rectangular=0,curvilinear=1)

  ! Dispersion models
  integer noDispersion,drude,gdm
  parameter( noDispersion=0, drude=1, gdm=2 )
  integer dispersionModel1, dispersionModel2

  ! Nonlinear models
       integer noNonlinearModel,multilevelAtomic
       parameter( noNonlinearModel=0, multilevelAtomic=1 )
  integer nonlinearModel1, nonlinearModel2
  

  gridType             =ipar(18)
  orderOfAccuracy      =ipar(19)

  dispersionModel1    = ipar(44)
  dispersionModel2    = ipar(45)

  nonlinearModel1     = ipar(49)
  nonlinearModel2     = ipar(50)
  

  if( nd.eq.3 .and. orderOfAccuracy.eq.4 .and. gridType.eq.curvilinear )then

      if( dispersionModel1.eq.noDispersion .and. dispersionModel2.eq.noDispersion )then
          call interfaceMx3dOrder4c( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange1, u1,u1n,u1m, wk1, mask1,rsxy1, xy1, p1,p1n,p1m, q1,q1n,q1m, boundaryCondition1, md1a,md1b,md2a,md2b,md3a,md3b,gridIndexRange2, u2,u2n,u2m, wk2, mask2,rsxy2, xy2, p2,p2n,p2m, q2,q2n,q2m, boundaryCondition2, ipar, rpar, aa2,aa4,aa8, ipvt2,ipvt4,ipvt8, ierr )
      else if( nonlinearModel1 .ne. noNonlinearModel .or. nonlinearModel2 .ne. noNonlinearModel )then
        ! --- MLA ---
          call interfaceMxMLA3dOrder4c( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange1, u1,u1n,u1m, wk1, mask1,rsxy1, xy1, p1,p1n,p1m, q1,q1n,q1m, boundaryCondition1, md1a,md1b,md2a,md2b,md3a,md3b,gridIndexRange2, u2,u2n,u2m, wk2, mask2,rsxy2, xy2, p2,p2n,p2m, q2,q2n,q2m, boundaryCondition2, ipar, rpar, aa2,aa4,aa8, ipvt2,ipvt4,ipvt8, ierr )
      else if( dispersionModel1.ne.noDispersion .or. dispersionModel2.ne.noDispersion )then
        ! -- GDM ---
          call interfaceMxGDM3dOrder4c( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange1, u1,u1n,u1m, wk1, mask1,rsxy1, xy1, p1,p1n,p1m, q1,q1n,q1m, boundaryCondition1, md1a,md1b,md2a,md2b,md3a,md3b,gridIndexRange2, u2,u2n,u2m, wk2, mask2,rsxy2, xy2, p2,p2n,p2m, q2,q2n,q2m, boundaryCondition2, ipar, rpar, aa2,aa4,aa8, ipvt2,ipvt4,ipvt8, ierr )
      else 
        stop 9876
      end if                  

  else

    write(*,'("interfaceOpt3dOrder4: ERROR: Unexpected options")') 
    write(*,'("   nd,order,gridType=",3i3)') nd,orderOfAccuracy,gridType
    stop 5678

  end if


  return
  end                          
