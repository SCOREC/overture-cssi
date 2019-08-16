#defineMacro lapCoeff4a(is,dr,ds) ( (-2/3.*rxx*is-2/3.*ryy*is)/dr+(4/3.*rx**2+4/3.*ry**2)/dr**2 )

#defineMacro lapCoeff4b(is,dr,ds) ( (1/12.*rxx*is+1/12.*ryy*is)/dr+(-1/12.*rx**2-1/12.*ry**2)/dr**2 )

#defineMacro xLapCoeff4a(is,dr,ds) ( (-1/2.*rxyy*is-1/2.*rxxx*is+(sy*(ry*sx*is+sy*rx*is)+3*rx*sx**2*is+ry*sy*sx*is)/ds**2)/dr+(2*ry*rxy+3*rx*rxx+ryy*rx)/dr**2+(ry**2*rx*is+rx**3*is)/dr**3 )

#defineMacro xLapCoeff4b(is,dr,ds) ( (-1/2.*rx**3*is-1/2.*ry**2*rx*is)/dr**3 )

#defineMacro yLapCoeff4a(is,dr,ds) ( (-1/2.*ryyy*is-1/2.*rxxy*is+(3*ry*sy**2*is+ry*sx**2*is+2*sy*rx*sx*is)/ds**2)/dr+(2*rxy*rx+ry*rxx+3*ry*ryy)/dr**2+(ry**3*is+ry*rx**2*is)/dr**3 )

#defineMacro yLapCoeff4b(is,dr,ds) ( (-1/2.*ry*rx**2*is-1/2.*ry**3*is)/dr**3 )

#defineMacro lapSqCoeff4a(is,dr,ds) ( (-1/2.*ryyyy*is-rxxyy*is-1/2.*rxxxx*is+(2*sy*(2*rx*sxy*is+ry*sxx*is+2*rxy*sx*is+sy*rxx*is)+2*ry*(2*sxy*sx*is+sy*sxx*is)+2*sy*(2*rxy*sx*is+2*rx*sxy*is)+4*syy*rx*sx*is+7*ry*sy*syy*is+ryy*sy**2*is+2*ryy*sx**2*is+rxx*sx**2*is+sy*(3*ry*syy*is+3*ryy*sy*is)+sy*(2*ryy*sy*is+2*ry*syy*is)+7*rx*sxx*sx*is+sx*(2*rxx*sx*is+2*rx*sxx*is)+4*ry*sxy*sx*is+sx*(3*rx*sxx*is+3*rxx*sx*is))/ds**2)/dr+(3*rxx**2+4*rxy**2+4*rxyy*rx+4*ry*ryyy+4*ry*rxxy+3*ryy**2+4*rx*rxxx+2*ryy*rxx+(2*sy*(-2*sy*rx**2-4*ry*rx*sx)+2*ry*(-4*sy*rx*sx-2*ry*sx**2)-12*rx**2*sx**2-12*ry**2*sy**2)/ds**2)/dr**2+(6*ry**2*ryy*is+6*rxx*rx**2*is+2*ryy*rx**2*is+4*ry*rxy*rx*is+2*ry*(ry*rxx*is+2*rxy*rx*is))/dr**3+(-4*rx**4-4*ry**4-8*ry**2*rx**2)/dr**4 )

#defineMacro lapSqCoeff4b(is,dr,ds) ( (2*ry*(-rxy*rx*is-1/2.*ry*rxx*is)-ryy*rx**2*is-2*ry*rxy*rx*is-3*ry**2*ryy*is-3*rxx*rx**2*is)/dr**3+(rx**4+ry**4+2*ry**2*rx**2)/dr**4 )

#defineMacro xxxxCoeff4a(is,dr,ds) ( (-1/2.*rxxxx*is+(rxx*sx**2*is+7*rx*sxx*sx*is+sx*(2*rxx*sx*is+2*rx*sxx*is)+sx*(3*rx*sxx*is+3*rxx*sx*is))/ds**2)/dr+(4*rx*rxxx+3*rxx**2-12*sx**2*rx**2/ds**2)/dr**2+6*rxx*rx**2*is/dr**3-4*rx**4/dr**4 )

#defineMacro xxxxCoeff4b(is,dr,ds) ( -3*rxx*rx**2*is/dr**3+rx**4/dr**4 )

#defineMacro yyyyCoeff4a(is,dr,ds) ( (-1/2.*ryyyy*is+(ryy*sy**2*is+7*ry*sy*syy*is+sy*(2*ryy*sy*is+2*ry*syy*is)+sy*(3*ry*syy*is+3*ryy*sy*is))/ds**2)/dr+(4*ry*ryyy+3*ryy**2-12*sy**2*ry**2/ds**2)/dr**2+6*ryy*ry**2*is/dr**3-4*ry**4/dr**4 )

#defineMacro yyyyCoeff4b(is,dr,ds) ( -3*ryy*ry**2*is/dr**3+ry**4/dr**4 )

#defineMacro xxyyCoeff4a(is,dr,ds) ( (-1/2.*rxxyy*is+(ryy*sx**2*is+2*syy*rx*sx*is+sy*(2*rxy*sx*is+2*rx*sxy*is)+2*ry*sxy*sx*is+ry*(2*sxy*sx*is+sy*sxx*is)+sy*(2*rx*sxy*is+ry*sxx*is+2*rxy*sx*is+sy*rxx*is))/ds**2)/dr+(2*rxyy*rx+2*rxy**2+ryy*rxx+2*ry*rxxy+(ry*(-4*sy*rx*sx-2*ry*sx**2)+sy*(-2*sy*rx**2-4*ry*rx*sx))/ds**2)/dr**2+(ry*(ry*rxx*is+2*rxy*rx*is)+2*ry*rxy*rx*is+ryy*rx**2*is)/dr**3-4*ry**2*rx**2/dr**4 )

#defineMacro xxyyCoeff4b(is,dr,ds) ( (-1/2.*ryy*rx**2*is-ry*rxy*rx*is+ry*(-rxy*rx*is-1/2.*ry*rxx*is))/dr**3+ry**2*rx**2/dr**4 )

