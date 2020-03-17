c Mie-Gruneisen EOS data for HMX products (modified, scaled)
c
c reference scales: rhoRef=.37778d0 g/cm^3
c                   pRef  = 21.3d0 GPa
c
c constants to define Gamma(rho) and Pref(rho)
      data gam0 / .43901d0 /
      data rhoc / 6.088199482d0 /
      data d2g / 4.381106011d-2 /
      data d3g / -9.196435650d-3 /
      data Ag / 6.999985625d-4 /
      data akg / -.1504017736d0 /
      data rhom / 9.264651386d0 /
      data qq / 3.93d0 /
c
c constants for the deriv of Gamma(rho)
      data d2gp / 8.762212022d-2 /   ! 2*d2g
      data d3gp / -2.758930695d-2 /  ! 3*d3g
c
c constants for the integral in the def of mu
      data e0g /  1.025036279d0 /
      data e1g / -8.722688728d-2 /
      data e2g /  3.065478550d-3 /
      data elog / 4.138032231d0 /
c
c specific heat constant
      data cvgas / 1.1120555120858d-2 /
c
c CJ gas temperature
      data tcjgas / 8.8333333333333d0 /
      data d0gas / 3.699233521d0 /
      data rhor / 4.859971412d0 /
