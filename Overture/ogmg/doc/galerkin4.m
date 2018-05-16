% averaging for fourth order approximation
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
maple('a21:=0;');
maple('a31:=-1/12/y^2;');
maple('a41:=0;');
maple('a51:=0;');

maple('a12:=0;');
maple('a22:=0;');
maple('a32:=16/12/y^2;');
maple('a42:=0;');
maple('a52:=0;');

maple('a13:=-1/12/x^2;');
maple('a23:=16/12/x^2;');
maple('a33:=-30/12*(1/x^2+1/y^2);');
maple('a43:=16/12/x^2;');
maple('a53:=-1/12/x^2;');

maple('a14:=0;');
maple('a24:=0;');
maple('a34:=16/12/y^2;');
maple('a44:=0;');
maple('a54:=0;');

maple('a15:=0;');
maple('a25:=0;');
maple('a35:=-1/12/y^2;');
maple('a45:=0;');
maple('a55:=0;');


% Here is the 1D average 
maple('ave1:=.5*a*c1;');
maple('ave2:=a*c2+b*c1+.5*a*c1+.5*(a*c3+b*c2+a*c1);');

maple('ave3:=a*c4+b*c3+a*c2+.5*(a*c3+b*c2+a*c1)+.5*(a*c5+b*c4+a*c3);');
maple('ave4:=a*c4+b*c5+.5*(a*c5+b*c4+a*c3)+.5*a*c5;');
maple('ave5:=.5*a*c5;');

% *********** averaging loop ******
% for i=1:10
for i=1:10

% average in x

maple('a11c:=subs(c1=a11,c2=a21,c3=a31,c4=a41,c5=a51,ave1);');
maple('a21c:=subs(c1=a11,c2=a21,c3=a31,c4=a41,c5=a51,ave2);');
maple('a31c:=subs(c1=a11,c2=a21,c3=a31,c4=a41,c5=a51,ave3);');
maple('a41c:=subs(c1=a11,c2=a21,c3=a31,c4=a41,c5=a51,ave4);');
maple('a51c:=subs(c1=a11,c2=a21,c3=a31,c4=a41,c5=a51,ave5);');

maple('a12c:=subs(c1=a12,c2=a22,c3=a32,c4=a42,c5=a52,ave1);');
maple('a22c:=subs(c1=a12,c2=a22,c3=a32,c4=a42,c5=a52,ave2);');
maple('a32c:=subs(c1=a12,c2=a22,c3=a32,c4=a42,c5=a52,ave3);');
maple('a42c:=subs(c1=a12,c2=a22,c3=a32,c4=a42,c5=a52,ave4);');
maple('a52c:=subs(c1=a12,c2=a22,c3=a32,c4=a42,c5=a52,ave5);');

maple('a13c:=subs(c1=a13,c2=a23,c3=a33,c4=a43,c5=a53,ave1);');
maple('a23c:=subs(c1=a13,c2=a23,c3=a33,c4=a43,c5=a53,ave2);');
maple('a33c:=subs(c1=a13,c2=a23,c3=a33,c4=a43,c5=a53,ave3);');
maple('a43c:=subs(c1=a13,c2=a23,c3=a33,c4=a43,c5=a53,ave4);');
maple('a53c:=subs(c1=a13,c2=a23,c3=a33,c4=a43,c5=a53,ave5);');

maple('a14c:=subs(c1=a14,c2=a24,c3=a34,c4=a44,c5=a54,ave1);');
maple('a24c:=subs(c1=a14,c2=a24,c3=a34,c4=a44,c5=a54,ave2);');
maple('a34c:=subs(c1=a14,c2=a24,c3=a34,c4=a44,c5=a54,ave3);');
maple('a44c:=subs(c1=a14,c2=a24,c3=a34,c4=a44,c5=a54,ave4);');
maple('a54c:=subs(c1=a14,c2=a24,c3=a34,c4=a44,c5=a54,ave5);');

maple('a15c:=subs(c1=a15,c2=a25,c3=a35,c4=a45,c5=a55,ave1);');
maple('a25c:=subs(c1=a15,c2=a25,c3=a35,c4=a45,c5=a55,ave2);');
maple('a35c:=subs(c1=a15,c2=a25,c3=a35,c4=a45,c5=a55,ave3);');
maple('a45c:=subs(c1=a15,c2=a25,c3=a35,c4=a45,c5=a55,ave4);');
maple('a55c:=subs(c1=a15,c2=a25,c3=a35,c4=a45,c5=a55,ave5);');

% average in y

maple('a11:=subs(c1=a11c,c2=a12c,c3=a13c,c4=a14c,c5=a15c,ave1)*4;');
maple('a12:=subs(c1=a11c,c2=a12c,c3=a13c,c4=a14c,c5=a15c,ave2)*4;');
maple('a13:=subs(c1=a11c,c2=a12c,c3=a13c,c4=a14c,c5=a15c,ave3)*4;');
maple('a14:=subs(c1=a11c,c2=a12c,c3=a13c,c4=a14c,c5=a15c,ave4)*4;');
maple('a15:=subs(c1=a11c,c2=a12c,c3=a13c,c4=a14c,c5=a15c,ave5)*4;');

maple('a21:=subs(c1=a21c,c2=a22c,c3=a23c,c4=a24c,c5=a25c,ave1)*4;');
maple('a22:=subs(c1=a21c,c2=a22c,c3=a23c,c4=a24c,c5=a25c,ave2)*4;');
maple('a23:=subs(c1=a21c,c2=a22c,c3=a23c,c4=a24c,c5=a25c,ave3)*4;');
maple('a24:=subs(c1=a21c,c2=a22c,c3=a23c,c4=a24c,c5=a25c,ave4)*4;');
maple('a25:=subs(c1=a21c,c2=a22c,c3=a23c,c4=a24c,c5=a25c,ave5)*4;');

maple('a31:=subs(c1=a31c,c2=a32c,c3=a33c,c4=a34c,c5=a35c,ave1)*4;');
maple('a32:=subs(c1=a31c,c2=a32c,c3=a33c,c4=a34c,c5=a35c,ave2)*4;');
maple('a33:=subs(c1=a31c,c2=a32c,c3=a33c,c4=a34c,c5=a35c,ave3)*4;');
maple('a34:=subs(c1=a31c,c2=a32c,c3=a33c,c4=a34c,c5=a35c,ave4)*4;');
maple('a35:=subs(c1=a31c,c2=a32c,c3=a33c,c4=a34c,c5=a35c,ave5)*4;');

maple('a41:=subs(c1=a41c,c2=a42c,c3=a43c,c4=a44c,c5=a45c,ave1)*4;');
maple('a42:=subs(c1=a41c,c2=a42c,c3=a43c,c4=a44c,c5=a45c,ave2)*4;');
maple('a43:=subs(c1=a41c,c2=a42c,c3=a43c,c4=a44c,c5=a45c,ave3)*4;');
maple('a44:=subs(c1=a41c,c2=a42c,c3=a43c,c4=a44c,c5=a45c,ave4)*4;');
maple('a45:=subs(c1=a41c,c2=a42c,c3=a43c,c4=a44c,c5=a45c,ave5)*4;');

maple('a51:=subs(c1=a51c,c2=a52c,c3=a53c,c4=a54c,c5=a55c,ave1)*4;');
maple('a52:=subs(c1=a51c,c2=a52c,c3=a53c,c4=a54c,c5=a55c,ave2)*4;');
maple('a53:=subs(c1=a51c,c2=a52c,c3=a53c,c4=a54c,c5=a55c,ave3)*4;');
maple('a54:=subs(c1=a51c,c2=a52c,c3=a53c,c4=a54c,c5=a55c,ave4)*4;');
maple('a55:=subs(c1=a51c,c2=a52c,c3=a53c,c4=a54c,c5=a55c,ave5)*4;');


maple('names:=[a11,a12,a13,a14,a15,a21,a22,a23,a24,a25,a31,a32,a33,a34,a35,a41,a42,a43,a44,a45,a51,a52,a53,a54,a55]');
maple('k:=0;');
for k=1:25
  maple('k:=k+1;');
  maple('aa:=names[k];');
  maple('aa:=expand(algsubs(1/x^2=dxsqi,aa));');
  maple('aa:=expand(algsubs(1/y^2=dysqi,aa));');

  aa =maple('C(aa):');

  m = fix((k-1)/5)+1;
  n = mod(k-1,5)+1;
  
  fprintf('CC(m%i%i)=  cLap*(%s);\n',m,n,aa(11:size(aa,2)-1));
%  if n == 5 
%   fprintf('\n');
%  end

end

fprintf(' ******* i=%i ********* \n',i);
pause;

end

% Here is the 3D laplace operator

maple('a[1,1,1]:=0;');
maple('a[2,1,1]:=0;');
maple('a[3,1,1]:=0;');
maple('a[4,1,1]:=0;');
maple('a[5,1,1]:=0;');

maple('a[1,2,1]:=0;');
maple('a[2,2,1]:=0;');
maple('a[3,2,1]:=0;');
maple('a[4,2,1]:=0;');
maple('a[5,2,1]:=0;');

maple('a[1,3,1]:=0;');
maple('a[2,3,1]:=0;');
maple('a[3,3,1]:=-1/12/z^2;');
maple('a[4,3,1]:=0;');
maple('a[5,3,1]:=0;');

maple('a[1,4,1]:=0;');
maple('a[2,4,1]:=0;');
maple('a[3,4,1]:=0;');
maple('a[4,4,1]:=0;');
maple('a[5,4,1]:=0;');

maple('a[1,5,1]:=0;');
maple('a[2,5,1]:=0;');
maple('a[3,5,1]:=0;');
maple('a[4,5,1]:=0;');
maple('a[5,5,1]:=0;');

maple('a[1,1,2]:=0;');
maple('a[2,1,2]:=0;');
maple('a[3,1,2]:=0;');
maple('a[4,1,2]:=0;');
maple('a[5,1,2]:=0;');

maple('a[1,2,2]:=0;');
maple('a[2,2,2]:=0;');
maple('a[3,2,2]:=0;');
maple('a[4,2,2]:=0;');
maple('a[5,2,2]:=0;');

maple('a[1,3,2]:=0;');
maple('a[2,3,2]:=0;');
maple('a[3,3,2]:=16/12/z^2;');
maple('a[4,3,2]:=0;');
maple('a[5,3,2]:=0;');

maple('a[1,4,2]:=0;');
maple('a[2,4,2]:=0;');
maple('a[3,4,2]:=0;');
maple('a[4,4,2]:=0;');
maple('a[5,4,2]:=0;');

maple('a[1,5,2]:=0;');
maple('a[2,5,2]:=0;');
maple('a[3,5,2]:=0;');
maple('a[4,5,2]:=0;');
maple('a[5,5,2]:=0;');


maple('a[1,1,3]:=0;');
maple('a[2,1,3]:=0;');
maple('a[3,1,3]:=-1/12/y^2;');
maple('a[4,1,3]:=0;');
maple('a[5,1,3]:=0;');

maple('a[1,2,3]:=0;');
maple('a[2,2,3]:=0;');
maple('a[3,2,3]:=16/12/y^2;');
maple('a[4,2,3]:=0;');
maple('a[5,2,3]:=0;');

maple('a[1,3,3]:=-1/12/x^2;');
maple('a[2,3,3]:=16/12/x^2;');
maple('a[3,3,3]:=-30/12*(1/x^2+1/y^2+1/z^2);');
maple('a[4,3,3]:=16/12/x^2;');
maple('a[5,3,3]:=-1/12/x^2;');

maple('a[1,4,3]:=0;');
maple('a[2,4,3]:=0;');
maple('a[3,4,3]:=16/12/y^2;');
maple('a[4,4,3]:=0;');
maple('a[5,4,3]:=0;');

maple('a[1,5,3]:=0;');
maple('a[2,5,3]:=0;');
maple('a[3,5,3]:=-1/12/y^2;');
maple('a[4,5,3]:=0;');
maple('a[5,5,3]:=0;');

maple('a[1,1,4]:=0;');
maple('a[2,1,4]:=0;');
maple('a[3,1,4]:=0;');
maple('a[4,1,4]:=0;');
maple('a[5,1,4]:=0;');

maple('a[1,2,4]:=0;');
maple('a[2,2,4]:=0;');
maple('a[3,2,4]:=0;');
maple('a[4,2,4]:=0;');
maple('a[5,2,4]:=0;');

maple('a[1,3,4]:=0;');
maple('a[2,3,4]:=0;');
maple('a[3,3,4]:=16/12/z^2;');
maple('a[4,3,4]:=0;');
maple('a[5,3,4]:=0;');

maple('a[1,4,4]:=0;');
maple('a[2,4,4]:=0;');
maple('a[3,4,4]:=0;');
maple('a[4,4,4]:=0;');
maple('a[5,4,4]:=0;');

maple('a[1,5,4]:=0;');
maple('a[2,5,4]:=0;');
maple('a[3,5,4]:=0;');
maple('a[4,5,4]:=0;');
maple('a[5,5,4]:=0;');

maple('a[1,1,5]:=0;');
maple('a[2,1,5]:=0;');
maple('a[3,1,5]:=0;');
maple('a[4,1,5]:=0;');
maple('a[5,1,5]:=0;');

maple('a[1,2,5]:=0;');
maple('a[2,2,5]:=0;');
maple('a[3,2,5]:=0;');
maple('a[4,2,5]:=0;');
maple('a[5,2,5]:=0;');

maple('a[1,3,5]:=0;');
maple('a[2,3,5]:=0;');
maple('a[3,3,5]:=-1/12/z^2;');
maple('a[4,3,5]:=0;');
maple('a[5,3,5]:=0;');

maple('a[1,4,5]:=0;');
maple('a[2,4,5]:=0;');
maple('a[3,4,5]:=0;');
maple('a[4,4,5]:=0;');
maple('a[5,4,5]:=0;');

maple('a[1,5,5]:=0;');
maple('a[2,5,5]:=0;');
maple('a[3,5,5]:=0;');
maple('a[4,5,5]:=0;');
maple('a[5,5,5]:=0;');

% *********** averaging loop ******
for i=1:10

% -------------- average in x ----------------------

maple('k:=0;');
for k=1:5
  maple('k:=k+1;');
  maple('n:=0;');
  for n=1:5
    maple('n:=n+1;');
    maple('ac[1,n,k]:=subs(c1=a[1,n,k],c2=a[2,n,k],c3=a[3,n,k],c4=a[4,n,k],c5=a[5,n,k],ave1);');
    maple('ac[2,n,k]:=subs(c1=a[1,n,k],c2=a[2,n,k],c3=a[3,n,k],c4=a[4,n,k],c5=a[5,n,k],ave2);');
    maple('ac[3,n,k]:=subs(c1=a[1,n,k],c2=a[2,n,k],c3=a[3,n,k],c4=a[4,n,k],c5=a[5,n,k],ave3);');
    maple('ac[4,n,k]:=subs(c1=a[1,n,k],c2=a[2,n,k],c3=a[3,n,k],c4=a[4,n,k],c5=a[5,n,k],ave4);');
    maple('ac[5,n,k]:=subs(c1=a[1,n,k],c2=a[2,n,k],c3=a[3,n,k],c4=a[4,n,k],c5=a[5,n,k],ave5);');
  end
end

% -------------- average in y ------------------------

maple('k:=0;');
for k=1:5
  maple('k:=k+1;');
  maple('n:=0;');
  for n=1:5
    maple('n:=n+1;');
    maple('ad[n,1,k]:=subs(c1=ac[n,1,k],c2=ac[n,2,k],c3=ac[n,3,k],c4=ac[n,4,k],c5=ac[n,5,k],ave1);');
    maple('ad[n,2,k]:=subs(c1=ac[n,1,k],c2=ac[n,2,k],c3=ac[n,3,k],c4=ac[n,4,k],c5=ac[n,5,k],ave2);');
    maple('ad[n,3,k]:=subs(c1=ac[n,1,k],c2=ac[n,2,k],c3=ac[n,3,k],c4=ac[n,4,k],c5=ac[n,5,k],ave3);');
    maple('ad[n,4,k]:=subs(c1=ac[n,1,k],c2=ac[n,2,k],c3=ac[n,3,k],c4=ac[n,4,k],c5=ac[n,5,k],ave4);');
    maple('ad[n,5,k]:=subs(c1=ac[n,1,k],c2=ac[n,2,k],c3=ac[n,3,k],c4=ac[n,4,k],c5=ac[n,5,k],ave5);');
  end
end

% --------------- average in z --------------

maple('k:=0;');
for k=1:5
  maple('k:=k+1;');
  maple('n:=0;');
  for n=1:5
    maple('n:=n+1;');
    maple('a[n,k,1]:=subs(c1=ad[n,k,1],c2=ad[n,k,2],c3=ad[n,k,3],c4=ad[n,k,4],c5=ad[n,k,5],ave1)*4;');
    maple('a[n,k,2]:=subs(c1=ad[n,k,1],c2=ad[n,k,2],c3=ad[n,k,3],c4=ad[n,k,4],c5=ad[n,k,5],ave2)*4;');
    maple('a[n,k,3]:=subs(c1=ad[n,k,1],c2=ad[n,k,2],c3=ad[n,k,3],c4=ad[n,k,4],c5=ad[n,k,5],ave3)*4;');
    maple('a[n,k,4]:=subs(c1=ad[n,k,1],c2=ad[n,k,2],c3=ad[n,k,3],c4=ad[n,k,4],c5=ad[n,k,5],ave4)*4;');
    maple('a[n,k,5]:=subs(c1=ad[n,k,1],c2=ad[n,k,2],c3=ad[n,k,3],c4=ad[n,k,4],c5=ad[n,k,5],ave5)*4;');
  end
end


maple('l:=0;');
for l=1:5
  maple('l:=l+1;');
  maple('k:=0;');
  for k=1:5
    maple('k:=k+1;');
    maple('n:=0;');
    for n=1:5
      maple('n:=n+1;');

      maple('aa:=expand(algsubs(1/x^2=dxsqi,a[n,k,l]));');
      maple('aa:=expand(algsubs(1/y^2=dysqi,aa));');
      maple('aa:=expand(algsubs(1/z^2=dzsqi,aa));');

      aa =maple('C(aa):');

      if n~=3 || k~=3 || l~=3 
        fprintf('CC(m%i%i%i)=   cLap*(%s); \n',n,k,l,aa(11:size(aa,2)-1));
      else
        fprintf('CC(m%i%i%i)=cI+cLap*(%s); \n',n,k,l,aa(11:size(aa,2)-1));
      end
    end
  end
end

fprintf(' ******* 3d i=%i *********\n ',i);
pause

end






