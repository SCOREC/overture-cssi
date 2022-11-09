// Macros to evaluate matrix Coefficients for CHAMp conditions.
// File written by champGenCoeff.maple

#beginMacro getChampCoeffOrder2NormalDirection0()
// ------  ORDER = 2 , NormalDirection = 0,------

const Real h2By2=h*h/2., h3By6=h2By2*h/3., h4By24=h3By6*h/4.;

Real L1IC = L1ICoeff(i1,i2,i3);
Real L1ICs = L1ICoeffs(i1,i2,i3);
Real L1rC = L1rCoeff(i1,i2,i3);
Real L1rCs = L1rCoeffs(i1,i2,i3);
Real L1sC = L1sCoeff(i1,i2,i3);
Real L1sCs = L1sCoeffs(i1,i2,i3);
Real L2IC = L2ICoeff(i1,i2,i3);
Real L2rC = L2rCoeff(i1,i2,i3);
Real L2sC = L2sCoeff(i1,i2,i3);
Real L2rrC = L2rrCoeff(i1,i2,i3);
Real L2rsC = L2rsCoeff(i1,i2,i3);
Real L2ssC = L2ssCoeff(i1,i2,i3);
Real cI  = (Sl*L1IC-bn*L2IC-bt*L1ICs)*h+Sl*h2By2*L2IC-bn*L1IC+Sl;
Real cr  = (Sl*L1rC-bn*L2rC-bt*L1rCs)*h+Sl*h2By2*L2rC-bn*L1rC;
Real cs  = ((-L1sCs-L1IC)*h-L2IC*h2By2-2)*bt+(Sl*L1sC-bn*L2sC)*h+Sl*h2By2*L2sC-bn*L1sC;
Real crr  = L2rrC*(Sl*h2By2-bn*h);
Real crs  = (-bn*L2rsC-bt*L1rC)*h+(Sl*L2rsC-bt*L2rC)*h2By2;
Real css  = (-bn*L2ssC-bt*L1sC)*h+(Sl*L2ssC-bt*L2sC)*h2By2;
#endMacro

#beginMacro getChampCoeffOrder4NormalDirection0()
// ------  ORDER = 4 , NormalDirection = 0,------

const Real h2By2=h*h/2., h3By6=h2By2*h/3., h4By24=h3By6*h/4.;

Real L1IC = L1ICoeff(i1,i2,i3);
Real L1ICs = L1ICoeffs(i1,i2,i3);
Real L1ICss = L1ICoeffss(i1,i2,i3);
Real L1ICsss = L1ICoeffsss(i1,i2,i3);
Real L1rC = L1rCoeff(i1,i2,i3);
Real L1rCs = L1rCoeffs(i1,i2,i3);
Real L1rCss = L1rCoeffss(i1,i2,i3);
Real L1rCsss = L1rCoeffsss(i1,i2,i3);
Real L1sC = L1sCoeff(i1,i2,i3);
Real L1sCs = L1sCoeffs(i1,i2,i3);
Real L1sCss = L1sCoeffss(i1,i2,i3);
Real L1sCsss = L1sCoeffsss(i1,i2,i3);
Real L2IC = L2ICoeff(i1,i2,i3);
Real L2ICs = L2ICoeffs(i1,i2,i3);
Real L2ICss = L2ICoeffss(i1,i2,i3);
Real L2rC = L2rCoeff(i1,i2,i3);
Real L2rCs = L2rCoeffs(i1,i2,i3);
Real L2rCss = L2rCoeffss(i1,i2,i3);
Real L2sC = L2sCoeff(i1,i2,i3);
Real L2sCs = L2sCoeffs(i1,i2,i3);
Real L2sCss = L2sCoeffss(i1,i2,i3);
Real L2rrC = L2rrCoeff(i1,i2,i3);
Real L2rrCs = L2rrCoeffs(i1,i2,i3);
Real L2rrCss = L2rrCoeffss(i1,i2,i3);
Real L2rsC = L2rsCoeff(i1,i2,i3);
Real L2rsCs = L2rsCoeffs(i1,i2,i3);
Real L2rsCss = L2rsCoeffss(i1,i2,i3);
Real L2ssC = L2ssCoeff(i1,i2,i3);
Real L2ssCs = L2ssCoeffs(i1,i2,i3);
Real L2ssCss = L2ssCoeffss(i1,i2,i3);
Real L3IC = L3ICoeff(i1,i2,i3);
Real L3ICs = L3ICoeffs(i1,i2,i3);
Real L3rC = L3rCoeff(i1,i2,i3);
Real L3rCs = L3rCoeffs(i1,i2,i3);
Real L3sC = L3sCoeff(i1,i2,i3);
Real L3sCs = L3sCoeffs(i1,i2,i3);
Real L3rrC = L3rrCoeff(i1,i2,i3);
Real L3rrCs = L3rrCoeffs(i1,i2,i3);
Real L3rsC = L3rsCoeff(i1,i2,i3);
Real L3rsCs = L3rsCoeffs(i1,i2,i3);
Real L3ssC = L3ssCoeff(i1,i2,i3);
Real L3ssCs = L3ssCoeffs(i1,i2,i3);
Real L3rrrC = L3rrrCoeff(i1,i2,i3);
Real L3rrrCs = L3rrrCoeffs(i1,i2,i3);
Real L3rrsC = L3rrsCoeff(i1,i2,i3);
Real L3rrsCs = L3rrsCoeffs(i1,i2,i3);
Real L3rssC = L3rssCoeff(i1,i2,i3);
Real L3rssCs = L3rssCoeffs(i1,i2,i3);
Real L3sssC = L3sssCoeff(i1,i2,i3);
Real L3sssCs = L3sssCoeffs(i1,i2,i3);
Real L4IC = L4ICoeff(i1,i2,i3);
Real L4rC = L4rCoeff(i1,i2,i3);
Real L4sC = L4sCoeff(i1,i2,i3);
Real L4rrC = L4rrCoeff(i1,i2,i3);
Real L4rsC = L4rsCoeff(i1,i2,i3);
Real L4ssC = L4ssCoeff(i1,i2,i3);
Real L4rrrC = L4rrrCoeff(i1,i2,i3);
Real L4rrsC = L4rrsCoeff(i1,i2,i3);
Real L4rssC = L4rssCoeff(i1,i2,i3);
Real L4sssC = L4sssCoeff(i1,i2,i3);
Real L4rrrrC = L4rrrrCoeff(i1,i2,i3);
Real L4rrrsC = L4rrrsCoeff(i1,i2,i3);
Real L4rrssC = L4rrssCoeff(i1,i2,i3);
Real L4rsssC = L4rsssCoeff(i1,i2,i3);
Real L4ssssC = L4ssssCoeff(i1,i2,i3);
Real cI  = Sl*(h*L1IC+h2By2*L2IC+h3By6*L3IC+h4By24*L4IC+1)+(-h*L2IC-h2By2*L3IC-h3By6*L4IC-L1IC)*bn-bt*(h*L1ICs+h2By2*L2ICs+h3By6*L3ICs);
Real cr  = Sl*(h*L1rC+h2By2*L2rC+h3By6*L3rC+h4By24*L4rC)+(-h*L2rC-h2By2*L3rC-h3By6*L4rC-L1rC)*bn-bt*(h*L1rCs+h2By2*L2rCs+h3By6*L3rCs);
Real cs  = ((-L2IC-L2sCs)*h2By2+(-L3IC-L3sCs)*h3By6+(-L1sCs-L1IC)*h-L4IC*h4By24-2)*bt+(Sl*L2sC-bn*L3sC)*h2By2+(Sl*L3sC-bn*L4sC)*h3By6+(h*L1sC+h4By24*L4sC)*Sl-bn*(h*L2sC+L1sC);
Real crr  = (Sl*L2rrC-bn*L3rrC-bt*L2rrCs)*h2By2+(Sl*L3rrC-bn*L4rrC-bt*L3rrCs)*h3By6-bn*L2rrC*h+Sl*h4By24*L4rrC;
Real crs  = ((-L2rsCs-L2rC)*h2By2+(-L3rC-L3rsCs)*h3By6-L1rC*h-L4rC*h4By24)*bt+(Sl*L2rsC-bn*L3rsC)*h2By2+(Sl*L3rsC-bn*L4rsC)*h3By6-bn*L2rsC*h+Sl*h4By24*L4rsC;
Real css  = ((-L2ssCs-L2sC)*h2By2+(-L3sC-L3ssCs)*h3By6-L1sC*h-L4sC*h4By24)*bt+(Sl*L2ssC-bn*L3ssC)*h2By2+(Sl*L3ssC-bn*L4ssC)*h3By6-bn*L2ssC*h+Sl*h4By24*L4ssC;
Real crrr  = (Sl*L3rrrC-bn*L4rrrC-bt*L3rrrCs)*h3By6+Sl*h4By24*L4rrrC-bn*h2By2*L3rrrC;
Real crrs  = ((-L3rrsCs-L3rrC)*bt+Sl*L3rrsC-L4rrsC*bn)*h3By6+(-h2By2*L2rrC-h4By24*L4rrC)*bt+Sl*h4By24*L4rrsC-bn*h2By2*L3rrsC;
Real crss  = ((-L3rssCs-L3rsC)*bt+Sl*L3rssC-L4rssC*bn)*h3By6+(-h2By2*L2rsC-h4By24*L4rsC)*bt+Sl*h4By24*L4rssC-bn*h2By2*L3rssC;
Real csss  = ((-L3sssCs-L3ssC)*bt+Sl*L3sssC-L4sssC*bn)*h3By6+(-h2By2*L2ssC-h4By24*L4ssC)*bt+Sl*h4By24*L4sssC-bn*h2By2*L3sssC;
Real crrrr  = L4rrrrC*(Sl*h4By24-bn*h3By6);
Real crrrs  = (-h3By6*L3rrrC-h4By24*L4rrrC)*bt+(Sl*h4By24-bn*h3By6)*L4rrrsC;
Real crrss  = (-h3By6*L3rrsC-h4By24*L4rrsC)*bt+(Sl*h4By24-bn*h3By6)*L4rrssC;
Real crsss  = (-h3By6*L3rssC-h4By24*L4rssC)*bt+(Sl*h4By24-bn*h3By6)*L4rsssC;
Real cssss  = (-h3By6*L3sssC-h4By24*L4sssC)*bt+(Sl*h4By24-bn*h3By6)*L4ssssC;
#endMacro

#beginMacro getChampCoeffOrder2NormalDirection1()
// ------  ORDER = 2 , NormalDirection = 1,------

const Real h2By2=h*h/2., h3By6=h2By2*h/3., h4By24=h3By6*h/4.;

Real L1IC = L1ICoeff(i1,i2,i3);
Real L1ICr = L1ICoeffr(i1,i2,i3);
Real L1rC = L1rCoeff(i1,i2,i3);
Real L1rCr = L1rCoeffr(i1,i2,i3);
Real L1sC = L1sCoeff(i1,i2,i3);
Real L1sCr = L1sCoeffr(i1,i2,i3);
Real L2IC = L2ICoeff(i1,i2,i3);
Real L2rC = L2rCoeff(i1,i2,i3);
Real L2sC = L2sCoeff(i1,i2,i3);
Real L2rrC = L2rrCoeff(i1,i2,i3);
Real L2rsC = L2rsCoeff(i1,i2,i3);
Real L2ssC = L2ssCoeff(i1,i2,i3);
Real cI  = (Sl*L1IC-bn*L2IC-bt*L1ICr)*h+Sl*h2By2*L2IC-bn*L1IC+Sl;
Real cr  = ((-L1rCr-L1IC)*h-L2IC*h2By2-2)*bt+(Sl*L1rC-bn*L2rC)*h+Sl*h2By2*L2rC-bn*L1rC;
Real cs  = (Sl*L1sC-bn*L2sC-bt*L1sCr)*h+Sl*h2By2*L2sC-bn*L1sC;
Real crr  = (-bn*L2rrC-bt*L1rC)*h+(Sl*L2rrC-bt*L2rC)*h2By2;
Real crs  = (-bn*L2rsC-bt*L1sC)*h+(Sl*L2rsC-bt*L2sC)*h2By2;
Real css  = L2ssC*(Sl*h2By2-bn*h);
#endMacro

#beginMacro getChampCoeffOrder4NormalDirection1()
// ------  ORDER = 4 , NormalDirection = 1,------

const Real h2By2=h*h/2., h3By6=h2By2*h/3., h4By24=h3By6*h/4.;

Real L1IC = L1ICoeff(i1,i2,i3);
Real L1ICr = L1ICoeffr(i1,i2,i3);
Real L1ICrr = L1ICoeffrr(i1,i2,i3);
Real L1ICrrr = L1ICoeffrrr(i1,i2,i3);
Real L1rC = L1rCoeff(i1,i2,i3);
Real L1rCr = L1rCoeffr(i1,i2,i3);
Real L1rCrr = L1rCoeffrr(i1,i2,i3);
Real L1rCrrr = L1rCoeffrrr(i1,i2,i3);
Real L1sC = L1sCoeff(i1,i2,i3);
Real L1sCr = L1sCoeffr(i1,i2,i3);
Real L1sCrr = L1sCoeffrr(i1,i2,i3);
Real L1sCrrr = L1sCoeffrrr(i1,i2,i3);
Real L2IC = L2ICoeff(i1,i2,i3);
Real L2ICr = L2ICoeffr(i1,i2,i3);
Real L2ICrr = L2ICoeffrr(i1,i2,i3);
Real L2rC = L2rCoeff(i1,i2,i3);
Real L2rCr = L2rCoeffr(i1,i2,i3);
Real L2rCrr = L2rCoeffrr(i1,i2,i3);
Real L2sC = L2sCoeff(i1,i2,i3);
Real L2sCr = L2sCoeffr(i1,i2,i3);
Real L2sCrr = L2sCoeffrr(i1,i2,i3);
Real L2rrC = L2rrCoeff(i1,i2,i3);
Real L2rrCr = L2rrCoeffr(i1,i2,i3);
Real L2rrCrr = L2rrCoeffrr(i1,i2,i3);
Real L2rsC = L2rsCoeff(i1,i2,i3);
Real L2rsCr = L2rsCoeffr(i1,i2,i3);
Real L2rsCrr = L2rsCoeffrr(i1,i2,i3);
Real L2ssC = L2ssCoeff(i1,i2,i3);
Real L2ssCr = L2ssCoeffr(i1,i2,i3);
Real L2ssCrr = L2ssCoeffrr(i1,i2,i3);
Real L3IC = L3ICoeff(i1,i2,i3);
Real L3ICr = L3ICoeffr(i1,i2,i3);
Real L3rC = L3rCoeff(i1,i2,i3);
Real L3rCr = L3rCoeffr(i1,i2,i3);
Real L3sC = L3sCoeff(i1,i2,i3);
Real L3sCr = L3sCoeffr(i1,i2,i3);
Real L3rrC = L3rrCoeff(i1,i2,i3);
Real L3rrCr = L3rrCoeffr(i1,i2,i3);
Real L3rsC = L3rsCoeff(i1,i2,i3);
Real L3rsCr = L3rsCoeffr(i1,i2,i3);
Real L3ssC = L3ssCoeff(i1,i2,i3);
Real L3ssCr = L3ssCoeffr(i1,i2,i3);
Real L3rrrC = L3rrrCoeff(i1,i2,i3);
Real L3rrrCr = L3rrrCoeffr(i1,i2,i3);
Real L3rrsC = L3rrsCoeff(i1,i2,i3);
Real L3rrsCr = L3rrsCoeffr(i1,i2,i3);
Real L3rssC = L3rssCoeff(i1,i2,i3);
Real L3rssCr = L3rssCoeffr(i1,i2,i3);
Real L3sssC = L3sssCoeff(i1,i2,i3);
Real L3sssCr = L3sssCoeffr(i1,i2,i3);
Real L4IC = L4ICoeff(i1,i2,i3);
Real L4rC = L4rCoeff(i1,i2,i3);
Real L4sC = L4sCoeff(i1,i2,i3);
Real L4rrC = L4rrCoeff(i1,i2,i3);
Real L4rsC = L4rsCoeff(i1,i2,i3);
Real L4ssC = L4ssCoeff(i1,i2,i3);
Real L4rrrC = L4rrrCoeff(i1,i2,i3);
Real L4rrsC = L4rrsCoeff(i1,i2,i3);
Real L4rssC = L4rssCoeff(i1,i2,i3);
Real L4sssC = L4sssCoeff(i1,i2,i3);
Real L4rrrrC = L4rrrrCoeff(i1,i2,i3);
Real L4rrrsC = L4rrrsCoeff(i1,i2,i3);
Real L4rrssC = L4rrssCoeff(i1,i2,i3);
Real L4rsssC = L4rsssCoeff(i1,i2,i3);
Real L4ssssC = L4ssssCoeff(i1,i2,i3);
Real cI  = Sl*(h*L1IC+h2By2*L2IC+h3By6*L3IC+h4By24*L4IC+1)+(-h*L2IC-h2By2*L3IC-h3By6*L4IC-L1IC)*bn-bt*(h*L1ICr+h2By2*L2ICr+h3By6*L3ICr);
Real cr  = ((-L2IC-L2rCr)*h2By2+(-L3IC-L3rCr)*h3By6+(-L1rCr-L1IC)*h-L4IC*h4By24-2)*bt+(Sl*L2rC-bn*L3rC)*h2By2+(Sl*L3rC-bn*L4rC)*h3By6+(h*L1rC+h4By24*L4rC)*Sl-bn*(h*L2rC+L1rC);
Real cs  = Sl*(h*L1sC+h2By2*L2sC+h3By6*L3sC+h4By24*L4sC)+(-h*L2sC-h2By2*L3sC-h3By6*L4sC-L1sC)*bn-bt*(h*L1sCr+h2By2*L2sCr+h3By6*L3sCr);
Real crr  = ((-L2rrCr-L2rC)*h2By2+(-L3rC-L3rrCr)*h3By6-L1rC*h-L4rC*h4By24)*bt+(Sl*L2rrC-bn*L3rrC)*h2By2+(Sl*L3rrC-bn*L4rrC)*h3By6-bn*L2rrC*h+Sl*h4By24*L4rrC;
Real crs  = ((-L2rsCr-L2sC)*h2By2+(-L3sC-L3rsCr)*h3By6-L1sC*h-L4sC*h4By24)*bt+(Sl*L2rsC-bn*L3rsC)*h2By2+(Sl*L3rsC-bn*L4rsC)*h3By6-bn*L2rsC*h+Sl*h4By24*L4rsC;
Real css  = (Sl*L2ssC-bn*L3ssC-bt*L2ssCr)*h2By2+(Sl*L3ssC-bn*L4ssC-bt*L3ssCr)*h3By6-bn*L2ssC*h+Sl*h4By24*L4ssC;
Real crrr  = ((-L3rrrCr-L3rrC)*bt+Sl*L3rrrC-bn*L4rrrC)*h3By6+(-h2By2*L2rrC-h4By24*L4rrC)*bt+Sl*h4By24*L4rrrC-bn*h2By2*L3rrrC;
Real crrs  = ((-L3rrsCr-L3rsC)*bt+Sl*L3rrsC-L4rrsC*bn)*h3By6+(-h2By2*L2rsC-h4By24*L4rsC)*bt+Sl*h4By24*L4rrsC-bn*h2By2*L3rrsC;
Real crss  = ((-L3rssCr-L3ssC)*bt+Sl*L3rssC-L4rssC*bn)*h3By6+(-h2By2*L2ssC-h4By24*L4ssC)*bt+Sl*h4By24*L4rssC-bn*h2By2*L3rssC;
Real csss  = (Sl*L3sssC-bn*L4sssC-bt*L3sssCr)*h3By6+Sl*h4By24*L4sssC-bn*h2By2*L3sssC;
Real crrrr  = (-h3By6*L3rrrC-h4By24*L4rrrC)*bt+L4rrrrC*(Sl*h4By24-bn*h3By6);
Real crrrs  = (-h3By6*L3rrsC-h4By24*L4rrsC)*bt+(Sl*h4By24-bn*h3By6)*L4rrrsC;
Real crrss  = (-h3By6*L3rssC-h4By24*L4rssC)*bt+(Sl*h4By24-bn*h3By6)*L4rrssC;
Real crsss  = (-h3By6*L3sssC-h4By24*L4sssC)*bt+(Sl*h4By24-bn*h3By6)*L4rsssC;
Real cssss  = (Sl*h4By24-bn*h3By6)*L4ssssC;
#endMacro
