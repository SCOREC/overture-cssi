clear;
maple('clear;');

maple('a1:=1/x^2;');
maple('a2:=-2*(1/x^2);');
maple('a3:=1/x^2;');

maple('a:=1/4;');
maple('b:=1/2;');

maple('ac1:=a*(a1+a2/2)+b*a1/2;')
maple('ac2:=a*(a2+a3+a1)+b*(a2+a1/2+a3/2);')
maple('ac3:=a*(a3+a2/2)+b*a3/2;')

% Here is the 2D laplace operator

maple('a11:=0;');
maple('a21:=1/y^2;');
maple('a31:=0;');
maple('a12:=1/x^2;');
maple('a22:=-2*(1/x^2+1/y^2);');
maple('a32:=1/x^2;');
maple('a13:=0;');
maple('a23:=1/y^2;');
maple('a33:=0;');


% Here is the 1D average 
maple('ave1:=a*(c1+c2/2)+b*c1/2;');
maple('ave2:=a*(c2+c3+c1)+b*(c2+c1/2+c3/2);');
maple('ave3:=a*(c3+c2/2)+b*c3/2;');

% *********** averaging loop ******
for i=1:10

% average in x

maple('a11c:=subs(c1=a11,c2=a21,c3=a31,ave1);');
maple('a21c:=subs(c1=a11,c2=a21,c3=a31,ave2);');
maple('a31c:=subs(c1=a11,c2=a21,c3=a31,ave3);');

maple('a12c:=subs(c1=a12,c2=a22,c3=a32,ave1);');
maple('a22c:=subs(c1=a12,c2=a22,c3=a32,ave2);');
maple('a32c:=subs(c1=a12,c2=a22,c3=a32,ave3);');

maple('a13c:=subs(c1=a13,c2=a23,c3=a33,ave1);');
maple('a23c:=subs(c1=a13,c2=a23,c3=a33,ave2);');
maple('a33c:=subs(c1=a13,c2=a23,c3=a33,ave3);');

% average in y

maple('a11:=subs(c1=a11c,c2=a12c,c3=a13c,ave1)*4;');
maple('a12:=subs(c1=a11c,c2=a12c,c3=a13c,ave2)*4;');
maple('a13:=subs(c1=a11c,c2=a12c,c3=a13c,ave3)*4;');

maple('a21:=subs(c1=a21c,c2=a22c,c3=a23c,ave1)*4;');
maple('a22:=subs(c1=a21c,c2=a22c,c3=a23c,ave2)*4;');
maple('a23:=subs(c1=a21c,c2=a22c,c3=a23c,ave3)*4;');

maple('a31:=subs(c1=a31c,c2=a32c,c3=a33c,ave1)*4;');
maple('a32:=subs(c1=a31c,c2=a32c,c3=a33c,ave2)*4;');
maple('a33:=subs(c1=a31c,c2=a32c,c3=a33c,ave3)*4;');

maple('names:=[a11,a12,a13,a21,a22,a23,a31,a32,a33]');
maple('k:=0;');
for k=1:9
  maple('k:=k+1;');
  maple('aa:=names[k];');
  maple('aa:=expand(algsubs(1/x^2=dxsqi,aa));');
  maple('aa:=expand(algsubs(1/y^2=dysqi,aa));');

  aa =maple('C(aa):');

  m = fix((k-1)/3)+1;
  n = mod(k-1,3)+1;
  
  fprintf('CC(m%i%i)=%s ',m,n,aa(11:size(aa,2)));
  if n == 3 
   fprintf('\n');
  end

end

fprintf(' ******* i=%i ********* \n',i);
pause;

end

% Here is the 3D laplace operator

maple('a111:=0;');
maple('a211:=0;');
maple('a311:=0;');
maple('a121:=0;');
maple('a221:=1/z^2;');
maple('a321:=0;');
maple('a131:=0;');
maple('a231:=0;');
maple('a331:=0;');


maple('a112:=0;');
maple('a212:=1/y^2;');
maple('a312:=0;');
maple('a122:=1/x^2;');
maple('a222:=-2*(1/x^2+1/y^2+1/z^2);');
maple('a322:=1/x^2;');
maple('a132:=0;');
maple('a232:=1/y^2;');
maple('a332:=0;');

maple('a113:=0;');
maple('a213:=0;');
maple('a313:=0;');
maple('a123:=0;');
maple('a223:=1/z^2;');
maple('a323:=0;');
maple('a133:=0;');
maple('a233:=0;');
maple('a333:=0;');

% *********** averaging loop ******
for i=1:10

% -------------- average in x ----------------------

maple('a111c:=subs(c1=a111,c2=a211,c3=a311,ave1);');
maple('a211c:=subs(c1=a111,c2=a211,c3=a311,ave2);');
maple('a311c:=subs(c1=a111,c2=a211,c3=a311,ave3);');

maple('a121c:=subs(c1=a121,c2=a221,c3=a321,ave1);');
maple('a221c:=subs(c1=a121,c2=a221,c3=a321,ave2);');
maple('a321c:=subs(c1=a121,c2=a221,c3=a321,ave3);');

maple('a131c:=subs(c1=a131,c2=a231,c3=a331,ave1);');
maple('a231c:=subs(c1=a131,c2=a231,c3=a331,ave2);');
maple('a331c:=subs(c1=a131,c2=a231,c3=a331,ave3);');


maple('a112c:=subs(c1=a112,c2=a212,c3=a312,ave1);');
maple('a212c:=subs(c1=a112,c2=a212,c3=a312,ave2);');
maple('a312c:=subs(c1=a112,c2=a212,c3=a312,ave3);');

maple('a122c:=subs(c1=a122,c2=a222,c3=a322,ave1);');
maple('a222c:=subs(c1=a122,c2=a222,c3=a322,ave2);');
maple('a322c:=subs(c1=a122,c2=a222,c3=a322,ave3);');

maple('a132c:=subs(c1=a132,c2=a232,c3=a332,ave1);');
maple('a232c:=subs(c1=a132,c2=a232,c3=a332,ave2);');
maple('a332c:=subs(c1=a132,c2=a232,c3=a332,ave3);');


maple('a113c:=subs(c1=a113,c2=a213,c3=a313,ave1);');
maple('a213c:=subs(c1=a113,c2=a213,c3=a313,ave2);');
maple('a313c:=subs(c1=a113,c2=a213,c3=a313,ave3);');

maple('a123c:=subs(c1=a123,c2=a223,c3=a323,ave1);');
maple('a223c:=subs(c1=a123,c2=a223,c3=a323,ave2);');
maple('a323c:=subs(c1=a123,c2=a223,c3=a323,ave3);');

maple('a133c:=subs(c1=a133,c2=a233,c3=a333,ave1);');
maple('a233c:=subs(c1=a133,c2=a233,c3=a333,ave2);');
maple('a333c:=subs(c1=a133,c2=a233,c3=a333,ave3);');


% -------------- average in y ------------------------

maple('a111d:=subs(c1=a111c,c2=a121c,c3=a131c,ave1);');
maple('a121d:=subs(c1=a111c,c2=a121c,c3=a131c,ave2);');
maple('a131d:=subs(c1=a111c,c2=a121c,c3=a131c,ave3);');

maple('a211d:=subs(c1=a211c,c2=a221c,c3=a231c,ave1);');
maple('a221d:=subs(c1=a211c,c2=a221c,c3=a231c,ave2);');
maple('a231d:=subs(c1=a211c,c2=a221c,c3=a231c,ave3);');

maple('a311d:=subs(c1=a311c,c2=a321c,c3=a331c,ave1);');
maple('a321d:=subs(c1=a311c,c2=a321c,c3=a331c,ave2);');
maple('a331d:=subs(c1=a311c,c2=a321c,c3=a331c,ave3);');

maple('a112d:=subs(c1=a112c,c2=a122c,c3=a132c,ave1);');
maple('a122d:=subs(c1=a112c,c2=a122c,c3=a132c,ave2);');
maple('a132d:=subs(c1=a112c,c2=a122c,c3=a132c,ave3);');

maple('a212d:=subs(c1=a212c,c2=a222c,c3=a232c,ave1);');
maple('a222d:=subs(c1=a212c,c2=a222c,c3=a232c,ave2);');
maple('a232d:=subs(c1=a212c,c2=a222c,c3=a232c,ave3);');

maple('a312d:=subs(c1=a312c,c2=a322c,c3=a332c,ave1);');
maple('a322d:=subs(c1=a312c,c2=a322c,c3=a332c,ave2);');
maple('a332d:=subs(c1=a312c,c2=a322c,c3=a332c,ave3);');

maple('a113d:=subs(c1=a113c,c2=a123c,c3=a133c,ave1);');
maple('a123d:=subs(c1=a113c,c2=a123c,c3=a133c,ave2);');
maple('a133d:=subs(c1=a113c,c2=a123c,c3=a133c,ave3);');

maple('a213d:=subs(c1=a213c,c2=a223c,c3=a233c,ave1);');
maple('a223d:=subs(c1=a213c,c2=a223c,c3=a233c,ave2);');
maple('a233d:=subs(c1=a213c,c2=a223c,c3=a233c,ave3);');

maple('a313d:=subs(c1=a313c,c2=a323c,c3=a333c,ave1);');
maple('a323d:=subs(c1=a313c,c2=a323c,c3=a333c,ave2);');
maple('a333d:=subs(c1=a313c,c2=a323c,c3=a333c,ave3);');


% --------------- average in z --------------

maple('a111:=subs(c1=a111d,c2=a112d,c3=a113d,ave1)*4;');
maple('a112:=subs(c1=a111d,c2=a112d,c3=a113d,ave2)*4;');
maple('a113:=subs(c1=a111d,c2=a112d,c3=a113d,ave3)*4;');

maple('a121:=subs(c1=a121d,c2=a122d,c3=a123d,ave1)*4;');
maple('a122:=subs(c1=a121d,c2=a122d,c3=a123d,ave2)*4;');
maple('a123:=subs(c1=a121d,c2=a122d,c3=a123d,ave3)*4;');

maple('a131:=subs(c1=a131d,c2=a132d,c3=a133d,ave1)*4;');
maple('a132:=subs(c1=a131d,c2=a132d,c3=a133d,ave2)*4;');
maple('a133:=subs(c1=a131d,c2=a132d,c3=a133d,ave3)*4;');

maple('a211:=subs(c1=a211d,c2=a212d,c3=a213d,ave1)*4;');
maple('a212:=subs(c1=a211d,c2=a212d,c3=a213d,ave2)*4;');
maple('a213:=subs(c1=a211d,c2=a212d,c3=a213d,ave3)*4;');

maple('a221:=subs(c1=a221d,c2=a222d,c3=a223d,ave1)*4;');
maple('a222:=subs(c1=a221d,c2=a222d,c3=a223d,ave2)*4;');
maple('a223:=subs(c1=a221d,c2=a222d,c3=a223d,ave3)*4;');

maple('a231:=subs(c1=a231d,c2=a232d,c3=a233d,ave1)*4;');
maple('a232:=subs(c1=a231d,c2=a232d,c3=a233d,ave2)*4;');
maple('a233:=subs(c1=a231d,c2=a232d,c3=a233d,ave3)*4;');

maple('a311:=subs(c1=a311d,c2=a312d,c3=a313d,ave1)*4;');
maple('a312:=subs(c1=a311d,c2=a312d,c3=a313d,ave2)*4;');
maple('a313:=subs(c1=a311d,c2=a312d,c3=a313d,ave3)*4;');

maple('a321:=subs(c1=a321d,c2=a322d,c3=a323d,ave1)*4;');
maple('a322:=subs(c1=a321d,c2=a322d,c3=a323d,ave2)*4;');
maple('a323:=subs(c1=a321d,c2=a322d,c3=a323d,ave3)*4;');

maple('a331:=subs(c1=a331d,c2=a332d,c3=a333d,ave1)*4;');
maple('a332:=subs(c1=a331d,c2=a332d,c3=a333d,ave2)*4;');
maple('a333:=subs(c1=a331d,c2=a332d,c3=a333d,ave3)*4;');


maple('names:=[a111,a121,a131,a211,a221,a231,a311,a321,a331,a112,a122,a132,a212,a222,a232,a312,a322,a332,a113,a123,a133,a213,a223,a233,a313,a323,a333]');
maple('k:=0;');
for k=1:27
  maple('k:=k+1;');
  maple('aa:=names[k];');
  maple('aa:=expand(algsubs(1/x^2=dxsqi,aa));');
  maple('aa:=expand(algsubs(1/y^2=dysqi,aa));');
  maple('aa:=expand(algsubs(1/z^2=dzsqi,aa));');

  aa =maple('C(aa):');

  l = mod(fix((k-1)/3),3)+1;
  m = mod((k-1),3)+1;
  n = fix((k-1)/9)+1;
  
  fprintf('CC(m%i%i%i)=    cLap*(%s); \n',l,m,n,aa(11:size(aa,2)-1));
%  if m == 3
%   fprintf('\n');
%  end

end

fprintf(' ******* i=%i *********\n ',i);
pause

end




