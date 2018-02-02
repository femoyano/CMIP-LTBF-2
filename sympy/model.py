# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from math import exp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
plt.ion()

"""
Created on Fri Jan  8 10:33:56 2016

@author: nano

Model
"""

# Define varialbes and functions

year = 31104000       # seconds in a year
month = 2592000       # seconds in a month
day = 86400           # seconds in a day
hour = 3600           # seconds in an hour
sec = 1               # seconds in a second!
tstep = hour


def T_resp_eq(k_ref, T, T_ref, E, R):
    # calculates the T response
    return k_ref * exp(-E/R * (1/T-1/T_ref))

# Define parameters
R = 0.008314
r_md_ref = 0.00028 / hour * tstep
r_ed_ref = 0.001 / hour * tstep
f_me = 5.6e-06 / hour * tstep
V_Dref = 1.0 / hour * tstep
D_dc0 = 8.1e-10 / sec * tstep
D_ec0 = 8.1e-11 / sec * tstep
k_ASref = 1.08e-6 / sec * tstep
k_DSref = 1.19e-10 / sec * tstep
K_Dref = 300000
f_mp = 0.5
T_ref = 293.15
E_VU = 47
E_VD = 47
E_KU = 30
E_KD = 30
E_ka = 10
E_kd = 10
E_md = 47
E_ed = 47
f_gr_ref = 0.7
f_gr_s = -0.016
pd = 2.7
psi_Rth = 15000
psi_fc = 33
d_pm = 10**-7
T = 280
clay = 0.1
silt = 0.2
sand = 0.7
ps = 0.45
z = 0.3
f_de = 0.01
M = 0.2
I_sl = 0.0848
I_ml = 0.00898

# Calculate some intermediate variables
b = 2.91 + 15.9 * clay
psi_sat = exp(6.5 - 1.3 * sand) / 1000
Rth = ps * (psi_sat / psi_Rth)**(1 / b)
fc = ps * (psi_sat / psi_fc)**(1 / b)
Md = 200 * (100 * clay)**0.6 * pd * (1 - ps)
K_D = T_resp_eq(K_Dref, T, T_ref, E_KD, R)
k_AS = T_resp_eq(k_ASref, T, T_ref, E_ka, R)
k_DS = T_resp_eq(k_DSref, T, T_ref, E_kd, R)
V_D = T_resp_eq(V_Dref, T, T_ref, E_VD, R)
r_md = T_resp_eq(r_md_ref, T, T_ref, E_md, R)
r_ed = T_resp_eq(r_ed_ref, T, T_ref, E_ed, R)
f_gr = f_gr_ref
diff_mod = (ps - Rth)**1.5 * ((M - Rth)/(ps - Rth))**2.5
D_dc = D_dc0 * diff_mod / d_pm
D_ec = D_ec0 * diff_mod / d_pm
M_fc = min(1, M / fc)


# Define fluxes

# Set flags
f_ads = 0  # adsorption desorption flag
f_mic = 0  # simulate microbial pool explicitly
f_fcs = 1  # scale Cp, Ca, ECs, M to field capacity (with max at fc)
f_sew = 1  # calculate EC and SC concentration in water


# model function for options A0M0F1S1
def f(y, t):
    Cp = y[0]
    Cd = y[1]
    Cew = y[2]
    Cem = y[3]
    # define fluxes
    F_slpc = I_sl
    F_mldc = I_ml
    F_pcdc = ((V_D * (Cew / (M * z)) * (Cp / z * M_fc)) /
              (K_D + (Cp / z * M_fc)) * z)
    F_slpc = I_sl
    F_mldc = I_ml
    F_pcdc = ((V_D * (Cew / (M * z)) * (Cp / z * M_fc)) /
              (K_D + (Cp / z * M_fc)) * z)
    F_dcrc = D_dc * (Cd - 0) * (1 - f_gr)
    F_dcpc = D_dc * (Cd - 0) * f_gr * (1 - f_de)
    F_dcem = D_dc * (Cd - 0) * f_gr * f_de
    F_emew = D_ec * (Cem - Cew)
    F_ewdc = Cew * r_ed
    F_emdc = Cem * r_ed
    # define differential equations
    dCp = F_slpc + F_dcpc - F_pcdc
    dCd = F_mldc + F_pcdc + F_ewdc + F_emdc - F_dcrc - F_dcpc - F_dcem
    dCew = F_emew - F_ewdc
    dCem = F_dcem - F_emew - F_emdc
    return [dCp, dCd, dCew, dCem]

# initial conditions
Cp0 = 1
Cd0 = 1
Cew0 = 1
Cem0 = 1
y0 = [Cp0, Cd0, Cew0, Cem0]
t = np.linspace(0, 216000, 1000)   # time grid

f(y0, t)

# solve the DEs
soln = odeint(f, y0, t)
Cp = soln[:, 0]
Cd = soln[:, 1]
Cew = soln[:, 2]
Cem = soln[:, 3]

# plot results
plt.figure()
plt.plot(t, Cp, label='Particulate C')
plt.plot(t, Cd, label='Dissolved C')
