      ! File created by Dropbox/GDM/maple/interface.maple
      d1=2.+(a1*alpha+b1)*dt

      ! --------------- Here is Et to second-order ------------

      ! Et = -1/2/dt*(2*a0*alpha*dt^2-2*b1*dt-4)/d1*E-1/2/dt*(2*b1*dt+4)/d1*Em-1/2/dt*(-b1*dt^3-2*dt^2)/d1*LE-1/2/dt*(-2*alpha*b0*dt^2-2*alpha*b1*dt)/d1*P-alpha*b1/d1*Pm-1/2/dt*(-b1*dt^3-2*dt^2)/d1*fE-dt*alpha*fP/d1
      ! Et = c2EtLE*LE + c2EtE*E + c2EtEm*Em + c2EtP*P + c2EtPm*Pm + c2EtfE*fE + c2EtfP*fP
      !   LE = c^2Delta(E) 
      c2EtLE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtE=(-a0*alpha*dt**2+b1*dt+2.)/dt/d1
      c2EtEm=(-b1*dt-2.)/dt/d1
      c2EtP=alpha*(b0*dt+b1)/d1
      c2EtPm=-1.*alpha*b1/d1
      c2EtfE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtfP=-1.*dt*alpha/d1
      ! --------------- Here is Pt to second-order ------------

      ! Pt = 1/2*(2*a0*dt^2+2*a1*dt)/dt/d1*E-a1/d1*Em+1/2*a1*dt^2/d1*LE+1/2*(2*a1*alpha*dt-2*b0*dt^2+4)/dt/d1*P+1/2*(-2*a1*alpha*dt-4)/dt/d1*Pm+1/2*a1*dt^2/d1*fE+fP*dt/d1
      ! Pt = c2PtLE*LE + c2PtE*E + c2PtEm*Em + c2PtP*P + c2PtPm*Pm + c2PtfE*fE + c2PtfP*fP
      !   LE = c^2Delta(E) 
      c2PtLE=.5*a1*dt**2/d1
      c2PtE=(1.*a0*dt+1.*a1)/d1
      c2PtEm=-1.*a1/d1
      c2PtP=(alpha*a1*dt-b0*dt**2+2.)/dt/d1
      c2PtPm=(-alpha*a1*dt-2.)/dt/d1
      c2PtfE=.5*a1*dt**2/d1
      c2PtfP=1.*dt/d1
      ! --------------- Here is Ptt to second-order ------------

      ! Ptt = a1*dt/d1*LE+((-a0*a1*alpha-a0*b1)*dt^2+d1*a0*dt+2*a1)/dt/d1*E-2*a1/dt/d1*Em+((a1*alpha*b0+b0*b1)*dt^2-d1*b0*dt-2*b1)/dt/d1*P+2*b1/dt/d1*Pm+a1*dt/d1*fE+((-a1*alpha-b1)*dt^2+dt*d1)/dt/d1*fP
      ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      c2PttLE=1.*a1*dt/d1
      c2PttE=((-1.*a1*alpha-1.*b1)*a0*dt**2+d1*a0*dt+2.*a1)/dt/d1
      c2PttEm=-2.*a1/dt/d1
      c2PttP=(b0*(a1*alpha+b1)*dt**2-1.*d1*b0*dt-2.*b1)/dt/d1
      c2PttPm=2.*b1/dt/d1
      c2PttfE=1.*a1*dt/d1
      c2PttfP=(1.*d1+(-1.*a1*alpha-1.*b1)*dt)/d1
