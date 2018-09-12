#!/usr/bin/python
#McDermott
#15 Sep 2017
#
# Calculations for compressible orifice flow
#
# Refs:
#   See my notes from 1996
#   Munson, Young, Okishi. Fundamentals of Fluid Mechanics. Wiley, 1990.

import math

HOC   = 50010. # heat of combustion [kJ/kg]

psig  = 0.0003
T_F   = 100.
C_d   = 0.85           # orifice discharge coefficient

N     = 1844             # number of holes
D_in  = 1./8.        # diameter [in]
D0_in = 8.*D_in       # upstream manifold diameter [in]

D     = D_in*2.54/100.        # fuel port diameter [m]
A     = N*math.pi*(D/2.)**2   # total flow area [m^2]

D0    = D0_in*2.54/100.
A0    = N*math.pi*(D0/2.)**2  # upstream flow area [m^2]

beta  = A/A0  # "beta ratio"

k     = 1.4    # isentropic coefficient
W     = 16.    # molecular weight
R     = 8314.5 # universal gas constant [Pa*m3/(kmol*K)]

T0    = 293. #(T_F+459.67)/1.8             # upstream absolute temperature [K]
patm  = 101325.                      # atmospheric pressure [Pa]
pcon  = 101325./14.696               # pressure units conversion factor
p0    = (psig + patm/pcon)*pcon      # upstream absolute pressure [Pa]
pb    = patm                         # downstream absolute pressure [Pa]

print('T0 [K] = '+str(T0))
print('p0 [Pa] = '+str(p0))
print('A [m2] = '+str(A))
print('beta = '+str(beta))

# determine critical pressure for choked flow

pstar = p0*(2./(k+1.))**(k/(k-1.))      # MYO (11.61)
Tstar = T0*(pstar/p0)**(k/(k-1.))       # MYO (11.58)

print('pb/p0 = '+str(pb/p0))
print('p*/p0 = '+str(pstar/p0))

if pb/p0 < pstar/p0:
    # sonic (choked)
    print('sonic')
    mdot = C_d*A*p0*math.sqrt( 2.*W/(R*T0) * (k/(k-1.)) * (1.-(2./(k+1.))) / ( ((k+1.)/2.)**(2./(k-1.)) - beta**4 ) )  # RJM notes (37)
    rho  = pstar*W/(R*Tstar)
else:
    # subsonic
    print('subsonic')
    mdot = C_d*A*p0*math.sqrt( 2.*W/(R*T0) * (k/(k-1.)) * ( 1.-(pb/p0)**((k-1.)/k) ) / ( (p0/pb)**(2./k) - beta**4 ) )  # RJM notes (39)
    rho  = pb*W/(R*T0)

print('mdot [kg/s] = '+str(mdot))
print('HOC [kJ/kg] = '+str(HOC))
print('HRR [kW] = '+str(mdot*HOC))
print('HRR [MBTU/h] = '+str(mdot*HOC*0.94783*3600/1.e6))

# determine velocity at nozzle exit

vdot = mdot/rho
vel  = vdot/A
print('vel [m/s] = '+str(vel))
