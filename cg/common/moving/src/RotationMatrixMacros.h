// 2017-01-18
// This file is generated by the Maple file /BeamINS3D/formulaForGetRotationMatrix.maple


//Entries for the rotation matrix R
#define GET_R(R) \
{\
real t1 = uSlope[0];\
real t2 = 1. + t1;\
real t3 = t2 * t2;\
real t4 = uSlope[1];\
real t5 = t4 * t4;\
real t6 = uSlope[2];\
real t7 = t6 * t6;\
real t9 = sqrt(t3 + t5 + t7);\
R(0,0) = t2 / t9;\
}\
{\
real t1 = uSlope[1];\
real t2 = uSlope[0];\
real t4 = pow(1. + t2, 2.);\
real t5 = t1 * t1;\
real t6 = uSlope[2];\
real t7 = t6 * t6;\
real t9 = sqrt(t4 + t5 + t7);\
R(0,1) = -t1 / t9;\
}\
{\
real t1 = uSlope[2];\
real t2 = uSlope[0];\
real t4 = pow(1. + t2, 2.);\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = t1 * t1;\
real t9 = sqrt(t4 + t6 + t7);\
R(0,2) = -t1 / t9;\
}\
{\
real t1 = uSlope[1];\
real t2 = uSlope[0];\
real t4 = pow(1. + t2, 2.);\
real t5 = t1 * t1;\
real t6 = uSlope[2];\
real t7 = t6 * t6;\
real t9 = sqrt(t4 + t5 + t7);\
R(1,0) = t1 / t9;\
}\
{\
real t1 = uSlope[0];\
real t2 = 1. + t1;\
real t3 = t2 * t2;\
real t4 = uSlope[1];\
real t5 = t4 * t4;\
real t6 = uSlope[2];\
real t7 = t6 * t6;\
real t8 = t3 + t5 + t7;\
real t9 = sqrt(t8);\
real t11 = t2 / t9;\
R(1,1) = t11 + t7 / t8 / (1. + t11);\
}\
{\
real t1 = uSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = t1 * t1;\
real t6 = uSlope[2];\
real t7 = t6 * t6;\
real t8 = t4 + t5 + t7;\
real t11 = sqrt(t8);\
R(1,2) = -t1 / t8 * t6 / (1. + t3 / t11);\
}\
{\
real t1 = uSlope[2];\
real t2 = uSlope[0];\
real t4 = pow(1. + t2, 2.);\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = t1 * t1;\
real t9 = sqrt(t4 + t6 + t7);\
R(2,0) = t1 / t9;\
}\
{\
real t1 = uSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = t1 * t1;\
real t6 = uSlope[2];\
real t7 = t6 * t6;\
real t8 = t4 + t5 + t7;\
real t11 = sqrt(t8);\
R(2,1) = -t1 / t8 * t6 / (1. + t3 / t11);\
}\
{\
real t1 = uSlope[0];\
real t2 = 1. + t1;\
real t3 = t2 * t2;\
real t4 = uSlope[1];\
real t5 = t4 * t4;\
real t6 = uSlope[2];\
real t7 = t6 * t6;\
real t8 = t3 + t5 + t7;\
real t9 = sqrt(t8);\
real t11 = t2 / t9;\
R(2,2) = t11 + t5 / t8 / (1. + t11);\
}\


//Entries for the time derivative of the rotation matrix R
#define GET_RD(R)  \
{\
real t1 = vSlope[0];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t17 = vSlope[1];\
real t19 = vSlope[2];\
R(0,0) = t1 / t10 - t3 / t10 / t9 * (2. * t3 * t1 + 2. * t5 * t17 + 2. * t7 * t19) / 2.;\
}\
{\
real t1 = vSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t16 = vSlope[0];\
real t19 = vSlope[2];\
R(0,1) = -t1 / t10 + t5 / t10 / t9 * (2. * t5 * t1 + 2. * t3 * t16 + 2. * t7 * t19) / 2.;\
}\
{\
real t1 = vSlope[2];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t16 = vSlope[0];\
real t18 = vSlope[1];\
R(0,2) = -t1 / t10 + t7 / t10 / t9 * (2. * t7 * t1 + 2. * t3 * t16 + 2. * t5 * t18) / 2.;\
}\
{\
real t1 = vSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t16 = vSlope[0];\
real t19 = vSlope[2];\
R(1,0) = t1 / t10 - t5 / t10 / t9 * (2. * t5 * t1 + 2. * t3 * t16 + 2. * t7 * t19) / 2.;\
}\
{\
real t1 = vSlope[0];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t11 = 1. / t10;\
real t12 = t1 * t11;\
real t17 = vSlope[1];\
real t19 = vSlope[2];\
real t22 = 2. * t3 * t1 + 2. * t5 * t17 + 2. * t7 * t19;\
real t24 = t3 / t10 / t9 * t22 / 2.;\
real t25 = 1. / t9;\
real t28 = t3 * t11 + 1.;\
real t29 = 1. / t28;\
real t33 = t9 * t9;\
real t39 = t28 * t28;\
R(1,1) = t12 - t24 + 2. * t7 * t25 * t29 * t19 - t8 / t33 * t29 * t22 - t8 * t25 / t39 * (t12 - t24);\
}\
{\
real t1 = vSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = 1. / t9;\
real t12 = sqrt(t9);\
real t13 = 1. / t12;\
real t15 = t3 * t13 + 1.;\
real t16 = 1. / t15;\
real t17 = t7 * t16;\
real t19 = t9 * t9;\
real t22 = vSlope[0];\
real t25 = vSlope[2];\
real t28 = 2. * t5 * t1 + 2. * t3 * t22 + 2. * t7 * t25;\
real t31 = t5 * t10;\
real t34 = t15 * t15;\
R(1,2) = -t1 * t10 * t17 + t5 / t19 * t17 * t28 - t31 * t25 * t16 + t31 * t7 / t34 * (t22 * t13 - t3 / t12 / t9 * t28 / 2.);\
}\
{\
real t1 = vSlope[2];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t16 = vSlope[0];\
real t18 = vSlope[1];\
R(2,0) = t1 / t10 - t7 / t10 / t9 * (2. * t7 * t1 + 2. * t3 * t16 + 2. * t5 * t18) / 2.;\
}\
{\
real t1 = vSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = 1. / t9;\
real t12 = sqrt(t9);\
real t13 = 1. / t12;\
real t15 = t3 * t13 + 1.;\
real t16 = 1. / t15;\
real t17 = t7 * t16;\
real t19 = t9 * t9;\
real t22 = vSlope[0];\
real t25 = vSlope[2];\
real t28 = 2. * t5 * t1 + 2. * t3 * t22 + 2. * t7 * t25;\
real t31 = t5 * t10;\
real t34 = t15 * t15;\
R(2,1) = -t1 * t10 * t17 + t5 / t19 * t17 * t28 - t31 * t25 * t16 + t31 * t7 / t34 * (t22 * t13 - t3 / t12 / t9 * t28 / 2.);\
}\
{\
real t1 = vSlope[0];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t11 = 1. / t10;\
real t12 = t1 * t11;\
real t17 = vSlope[1];\
real t19 = vSlope[2];\
real t22 = 2. * t3 * t1 + 2. * t5 * t17 + 2. * t7 * t19;\
real t24 = t3 / t10 / t9 * t22 / 2.;\
real t25 = 1. / t9;\
real t28 = t3 * t11 + 1.;\
real t29 = 1. / t28;\
real t33 = t9 * t9;\
real t39 = t28 * t28;\
R(2,2) = t12 - t24 + 2. * t5 * t25 * t29 * t17 - t6 / t33 * t29 * t22 - t6 * t25 / t39 * (t12 - t24);\
}\


//Entries for the second time derivative of the rotation matrix R
#define GET_RDD(R)  \
{\
real t1 = aSlope[0];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t13 = vSlope[0];\
real t15 = 1. / t10 / t9;\
real t18 = vSlope[1];\
real t20 = vSlope[2];\
real t23 = 2. * t3 * t13 + 2. * t5 * t18 + 2. * t7 * t20;\
real t25 = t9 * t9;\
real t29 = t23 * t23;\
real t33 = t13 * t13;\
real t35 = t18 * t18;\
real t36 = aSlope[1];\
real t38 = t20 * t20;\
real t39 = aSlope[2];\
R(0,0) = t1 / t10 - t13 * t15 * t23 + 3. / 4. * t3 / t10 / t25 * t29 - t3 * t15 * (2. * t3 * t1 + 2. * t5 * t36 + 2. * t7 * t39 + 2. * t33 + 2. * t35 + 2. * t38) / 2.;\
}\
{\
real t1 = aSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t13 = vSlope[1];\
real t15 = 1. / t10 / t9;\
real t17 = vSlope[0];\
real t20 = vSlope[2];\
real t23 = 2. * t5 * t13 + 2. * t3 * t17 + 2. * t7 * t20;\
real t25 = t9 * t9;\
real t29 = t23 * t23;\
real t33 = t17 * t17;\
real t34 = aSlope[0];\
real t36 = t13 * t13;\
real t38 = t20 * t20;\
real t39 = aSlope[2];\
R(0,1) = -t1 / t10 + t13 * t15 * t23 - 3. / 4. * t5 / t10 / t25 * t29 + t5 * t15 * (2. * t5 * t1 + 2. * t3 * t34 + 2. * t7 * t39 + 2. * t33 + 2. * t36 + 2. * t38) / 2.;\
}\
{\
real t1 = aSlope[2];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t13 = vSlope[2];\
real t15 = 1. / t10 / t9;\
real t17 = vSlope[0];\
real t19 = vSlope[1];\
real t23 = 2. * t7 * t13 + 2. * t3 * t17 + 2. * t5 * t19;\
real t25 = t9 * t9;\
real t29 = t23 * t23;\
real t33 = t17 * t17;\
real t34 = aSlope[0];\
real t36 = t19 * t19;\
real t37 = aSlope[1];\
real t39 = t13 * t13;\
R(0,2) = -t1 / t10 + t13 * t15 * t23 - 3. / 4. * t7 / t10 / t25 * t29 + t7 * t15 * (2. * t7 * t1 + 2. * t3 * t34 + 2. * t5 * t37 + 2. * t33 + 2. * t36 + 2. * t39) / 2.;\
}\
{\
real t1 = aSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t13 = vSlope[1];\
real t15 = 1. / t10 / t9;\
real t17 = vSlope[0];\
real t20 = vSlope[2];\
real t23 = 2. * t5 * t13 + 2. * t3 * t17 + 2. * t7 * t20;\
real t25 = t9 * t9;\
real t29 = t23 * t23;\
real t33 = t17 * t17;\
real t34 = aSlope[0];\
real t36 = t13 * t13;\
real t38 = t20 * t20;\
real t39 = aSlope[2];\
R(1,0) = t1 / t10 - t13 * t15 * t23 + 3. / 4. * t5 / t10 / t25 * t29 - t5 * t15 * (2. * t5 * t1 + 2. * t3 * t34 + 2. * t7 * t39 + 2. * t33 + 2. * t36 + 2. * t38) / 2.;\
}\
{\
real t1 = aSlope[0];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t11 = 1. / t10;\
real t12 = t1 * t11;\
real t13 = vSlope[0];\
real t15 = 1. / t10 / t9;\
real t18 = vSlope[1];\
real t20 = vSlope[2];\
real t23 = 2. * t3 * t13 + 2. * t5 * t18 + 2. * t7 * t20;\
real t24 = t13 * t15 * t23;\
real t25 = t9 * t9;\
real t29 = t23 * t23;\
real t31 = 3. / 4. * t3 / t10 / t25 * t29;\
real t32 = t3 * t15;\
real t33 = t13 * t13;\
real t35 = t18 * t18;\
real t36 = aSlope[1];\
real t38 = t20 * t20;\
real t39 = aSlope[2];\
real t42 = 2. * t3 * t1 + 2. * t5 * t36 + 2. * t7 * t39 + 2. * t33 + 2. * t35 + 2. * t38;\
real t44 = t32 * t42 / 2.;\
real t45 = 1. / t9;\
real t48 = t3 * t11 + 1.;\
real t49 = 1. / t48;\
real t52 = 1. / t25;\
real t58 = t7 * t45;\
real t59 = t48 * t48;\
real t60 = 1. / t59;\
real t65 = t13 * t11 - t32 * t23 / 2.;\
real t78 = t8 * t52;\
real t85 = t8 * t45;\
real t88 = t65 * t65;\
R(1,1) = t12 - t24 + t31 - t44 + 2. * t38 * t45 * t49 - 4. * t7 * t52 * t49 * t20 * t23 - 4. * t58 * t60 * t20 * t65 + 2. * t58 * t49 * t39 + 2. * t8 / t25 / t9 * t49 * t29 + 2. * t78 * t60 * t23 * t65 - t78 * t49 * t42 + 2. * t85 / t59 / t48 * t88 - t85 * t60 * (t12 - t24 + t31 - t44);\
}\
{\
real t1 = aSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = 1. / t9;\
real t12 = sqrt(t9);\
real t13 = 1. / t12;\
real t15 = t3 * t13 + 1.;\
real t16 = 1. / t15;\
real t17 = t7 * t16;\
real t19 = vSlope[1];\
real t20 = t9 * t9;\
real t21 = 1. / t20;\
real t23 = vSlope[0];\
real t26 = vSlope[2];\
real t29 = 2. * t5 * t19 + 2. * t3 * t23 + 2. * t7 * t26;\
real t33 = t19 * t10;\
real t34 = t26 * t16;\
real t37 = t15 * t15;\
real t38 = 1. / t37;\
real t39 = t7 * t38;\
real t42 = 1. / t12 / t9;\
real t43 = t3 * t42;\
real t46 = t23 * t13 - t43 * t29 / 2.;\
real t53 = t29 * t29;\
real t57 = t5 * t21;\
real t66 = t23 * t23;\
real t67 = aSlope[0];\
real t69 = t19 * t19;\
real t71 = t26 * t26;\
real t72 = aSlope[2];\
real t75 = 2. * t5 * t1 + 2. * t3 * t67 + 2. * t7 * t72 + 2. * t66 + 2. * t69 + 2. * t71;\
real t78 = t5 * t10;\
real t88 = t46 * t46;\
R(1,2) = -t1 * t10 * t17 + 2. * t19 * t21 * t17 * t29 - 2. * t33 * t34 + 2. * t33 * t39 * t46 - 2. * t5 / t20 / t9 * t17 * t53 + 2. * t57 * t34 * t29 - 2. * t57 * t7 * t38 * t29 * t46 + t57 * t17 * t75 - t78 * t72 * t16 + 2. * t78 * t26 * t38 * t46 - 2. * t78 * t7 / t37 / t15 * t88 + t78 * t39 * (t67 * t13 - t23 * t42 * t29 + 3. / 4. * t3 / t12 / t20 * t53 - t43 * t75 / 2.);\
}\
{\
real t1 = aSlope[2];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t13 = vSlope[2];\
real t15 = 1. / t10 / t9;\
real t17 = vSlope[0];\
real t19 = vSlope[1];\
real t23 = 2. * t7 * t13 + 2. * t3 * t17 + 2. * t5 * t19;\
real t25 = t9 * t9;\
real t29 = t23 * t23;\
real t33 = t17 * t17;\
real t34 = aSlope[0];\
real t36 = t19 * t19;\
real t37 = aSlope[1];\
real t39 = t13 * t13;\
R(2,0) = t1 / t10 - t13 * t15 * t23 + 3. / 4. * t7 / t10 / t25 * t29 - t7 * t15 * (2. * t7 * t1 + 2. * t3 * t34 + 2. * t5 * t37 + 2. * t33 + 2. * t36 + 2. * t39) / 2.;\
}\
{\
real t1 = aSlope[1];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = 1. / t9;\
real t12 = sqrt(t9);\
real t13 = 1. / t12;\
real t15 = t3 * t13 + 1.;\
real t16 = 1. / t15;\
real t17 = t7 * t16;\
real t19 = vSlope[1];\
real t20 = t9 * t9;\
real t21 = 1. / t20;\
real t23 = vSlope[0];\
real t26 = vSlope[2];\
real t29 = 2. * t5 * t19 + 2. * t3 * t23 + 2. * t7 * t26;\
real t33 = t19 * t10;\
real t34 = t26 * t16;\
real t37 = t15 * t15;\
real t38 = 1. / t37;\
real t39 = t7 * t38;\
real t42 = 1. / t12 / t9;\
real t43 = t3 * t42;\
real t46 = t23 * t13 - t43 * t29 / 2.;\
real t53 = t29 * t29;\
real t57 = t5 * t21;\
real t66 = t23 * t23;\
real t67 = aSlope[0];\
real t69 = t19 * t19;\
real t71 = t26 * t26;\
real t72 = aSlope[2];\
real t75 = 2. * t5 * t1 + 2. * t3 * t67 + 2. * t7 * t72 + 2. * t66 + 2. * t69 + 2. * t71;\
real t78 = t5 * t10;\
real t88 = t46 * t46;\
R(2,1) = -t1 * t10 * t17 + 2. * t19 * t21 * t17 * t29 - 2. * t33 * t34 + 2. * t33 * t39 * t46 - 2. * t5 / t20 / t9 * t17 * t53 + 2. * t57 * t34 * t29 - 2. * t57 * t7 * t38 * t29 * t46 + t57 * t17 * t75 - t78 * t72 * t16 + 2. * t78 * t26 * t38 * t46 - 2. * t78 * t7 / t37 / t15 * t88 + t78 * t39 * (t67 * t13 - t23 * t42 * t29 + 3. / 4. * t3 / t12 / t20 * t53 - t43 * t75 / 2.);\
}\
{\
real t1 = aSlope[0];\
real t2 = uSlope[0];\
real t3 = 1. + t2;\
real t4 = t3 * t3;\
real t5 = uSlope[1];\
real t6 = t5 * t5;\
real t7 = uSlope[2];\
real t8 = t7 * t7;\
real t9 = t4 + t6 + t8;\
real t10 = sqrt(t9);\
real t11 = 1. / t10;\
real t12 = t1 * t11;\
real t13 = vSlope[0];\
real t15 = 1. / t10 / t9;\
real t18 = vSlope[1];\
real t20 = vSlope[2];\
real t23 = 2. * t3 * t13 + 2. * t5 * t18 + 2. * t7 * t20;\
real t24 = t13 * t15 * t23;\
real t25 = t9 * t9;\
real t29 = t23 * t23;\
real t31 = 3. / 4. * t3 / t10 / t25 * t29;\
real t32 = t3 * t15;\
real t33 = t13 * t13;\
real t35 = t18 * t18;\
real t36 = aSlope[1];\
real t38 = t20 * t20;\
real t39 = aSlope[2];\
real t42 = 2. * t3 * t1 + 2. * t5 * t36 + 2. * t7 * t39 + 2. * t33 + 2. * t35 + 2. * t38;\
real t44 = t32 * t42 / 2.;\
real t45 = 1. / t9;\
real t48 = t3 * t11 + 1.;\
real t49 = 1. / t48;\
real t52 = 1. / t25;\
real t58 = t5 * t45;\
real t59 = t48 * t48;\
real t60 = 1. / t59;\
real t65 = t13 * t11 - t32 * t23 / 2.;\
real t78 = t6 * t52;\
real t85 = t6 * t45;\
real t88 = t65 * t65;\
R(2,2) = t12 - t24 + t31 - t44 + 2. * t35 * t45 * t49 - 4. * t5 * t52 * t49 * t18 * t23 - 4. * t58 * t60 * t18 * t65 + 2. * t58 * t49 * t36 + 2. * t6 / t25 / t9 * t49 * t29 + 2. * t78 * t60 * t23 * t65 - t78 * t49 * t42 + 2. * t85 / t59 / t48 * t88 - t85 * t60 * (t12 - t24 + t31 - t44);\
}\
