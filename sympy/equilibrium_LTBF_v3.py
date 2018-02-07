# -*- coding: utf-8 -*-
from __future__ import division
import time
import sympy as sy
from math import exp

"""
Created on Tue Dec 29 11:56:58 2015
Last modified: Dec 2016
@author: Fernando Moyano
Script for solving the steady state equations
"""

# Model options used in this script:
diff_fun  = "hama" # Options: 'hama', 'cubic'
dec_fun   = "MM" # One of: 'MM', '2nd', '1st'
upt_fun   = "1st" # One of: 'MM', '2nd', '1st'
#flag_mmr  = 1 # microbial maintenance respiration
#flag_mic  = 1 # simulate microbial pool explicitly
#flag_fcs  = 0 # scale C_P and M to field capacity (with max at fc)
#flag_sew  = 0 # calculate C_E and C_D concentration in water
#flag_dte  = 0 # diffusivity temperature effect on/off
#flag_dce  = 0 # diffusivity carbon effect on/off

# Define varialbes
year = 31104000       # seconds in a year
month = 2592000       # seconds in a month
day = 86400           # seconds in a day
hour = 3600           # seconds in an hour
sec = 1               # seconds in a second!
tstep = hour

# Define functions

def T_resp_eq(k_ref, T, T_ref, E, R):
    return k_ref * sy.exp(-E/R * (1/T-1/T_ref))

# Define symbols
C_P, C_D, C_A, C_M = \
    sy.symbols('C_P C_D C_A C_M')


r_md, f_ug, f_ge, r_mr= (sy.symbols('r_md f_ug f_ge r_mr'))
V_D, K_D, V_U, K_U = (sy.symbols('V_D  K_D V_U K_U '))
g, M_fc, min_md = (sy.symbols('g M_fc min_md'))
M, I_sl, I_ml, mc, depth = sy.symbols('M I_sl I_ml mc depth')
Ka, k_ads, k_des, Amax = sy.symbols('Ka k_ads k_des Amax')

# Define fluxes

D_diff = g * C_D

F_slcp = I_sl
F_mlcd = I_ml

if dec_fun == "MM":
      F_cpcd = depth * V_D * C_P/depth * C_M/depth / (K_D + C_P/depth)    
if dec_fun == "2nd":
      F_cpcd = depth * V_D * C_P/depth * C_M/depth
if dec_fun == "1st":
      F_cpcd = V_D * C_P
      
if upt_fun == "MM":
      U_cd = depth * V_U * D_diff/depth * C_M/depth / (K_U + C_D/depth)
if upt_fun == "2nd":
      U_cd = depth * V_U * D_diff/depth * C_M/depth
if upt_fun == "1st":
      U_cd = V_U * D_diff
      

F_cdcm = U_cd * f_ug
F_cdcr = U_cd * (1 - f_ug)

F_cmcp = C_M * (C_M * r_md + min_md)
F_cmcr = C_M * r_mr

dC_P = F_slcp + F_cmcp - F_cpcd
dC_D = F_mlcd + F_cpcd - F_cdcr - F_cdcm
dC_M = F_cdcm - F_cmcp - F_cmcr

sol = sy.solve([dC_P, dC_D, dC_M],
                        [C_P, C_D, C_M], dict=True)

sol = sol[0]
sol_C_P = sol[C_P]
sol_C_D = sol[C_D]
sol_C_M = sol[C_M]

#%%

# Site data
clay = 0.15
sand = 0.28
silt = 0.57
ps = 0.45
I_sl_v = 0.00005
I_ml_v = 0.000005
depth_v = 0.3

# Intermediate parameter values
g_0 = 2.2 / hour * tstep
E_m = 10
E_e = 10
E_K = 89
E_V = 87
E_r = 95
K_D_ref = 50
K_U_ref = 1
k_ads_ref = 1.08e-6 / sec * tstep
k_des_ref = 1.19e-10 / sec * tstep
mc_0 = 0.00001
pd = 2700
psi_fc = 33
psi_Dth = 15000
R = 0.008314
r_md_ref = 0.0015 / hour * tstep
r_mr_ref = 0.000042
T = 288.15
T_ref = 293.15
V_D_ref = 0.35 / hour * tstep
V_U_ref = 0.09 / hour * tstep
n = 2.3 
m = 1.2
Dth = 0.06

# End parameter values
f_ug_v = 0.50
f_ge_v = 0.01 / hour * tstep
M_v = 0.2

# Calculate intermediate variables
b = 2.91 + 15.9 * clay
k_ads_v = T_resp_eq(k_ads_ref, T, T_ref, E_V, R)
k_des_v = T_resp_eq(k_des_ref, T, T_ref, E_V, R)
psi_sat = exp(6.5 - 1.3 * sand) / 1000
fc = ps * (psi_sat / psi_fc)**(1 / b)
if diff_fun == "hama":
    g_sm = (ps - Dth)**m * ((M_v - Dth)/(ps - Dth))**n
if diff_fun == "cubic":
    D_sm = M_v**3

# Calculate end variables
K_D_v = T_resp_eq(K_D_ref, T, T_ref, E_K, R)
V_D_v = T_resp_eq(V_D_ref, T, T_ref, E_V, R)
K_U_v = T_resp_eq(K_U_ref, T, T_ref, E_K, R)
V_U_v = T_resp_eq(V_U_ref, T, T_ref, E_V, R)
r_md_v = T_resp_eq(r_md_ref, T, T_ref, E_m, R)
r_mr_v  = T_resp_eq(r_mr_ref, T, T_ref, E_r, R)
g_v = g_0 * g_sm
Amax_v = 200 * (100 * clay)**0.6 * pd * (1 - ps) / 1000000 #from mg kg-1 to kg m-3
M_fc_v = 1 #sy.Min(1, M / fc)
Ka_v = k_ads/k_des
mc_v = mc_0 * pd * (1 - ps) * depth # [kgC m-2] basal microbial carbon

# Substitute variables (parameters) with values

eq_C_P = sol_C_P.subs([
    (g, g_v), (f_ug, f_ug_v), (f_ge, f_ge_v),
    (r_mr, r_mr_v),(I_ml, I_ml_v), (I_sl, I_sl_v), (K_D, K_D_v), (K_U, K_U_v), 
    (M, M_v), (mc, mc_v), (M_fc, M_fc), (r_md, r_md_v),
    (V_D, V_D_v), (V_U, V_U_v), (depth, depth_v)
    ])

eq_C_D = sol_C_D.subs([
    (g, g_v), (f_ug, f_ug_v), (f_ge, f_ge_v),
    (r_mr, r_mr_v),(I_ml, I_ml_v), (I_sl, I_sl_v), (K_D, K_D_v), (K_U, K_U_v), 
    (M, M_v), (mc, mc_v), (M_fc, M_fc), (r_md, r_md_v),
    (V_D, V_D_v), (V_U, V_U_v), (depth, depth_v)
    ])

eq_C_M = sol_C_M.subs([
    (g, g_v), (f_ug, f_ug_v), (f_ge, f_ge_v),
    (r_mr, r_mr_v),(I_ml, I_ml_v), (I_sl, I_sl_v), (K_D, K_D_v), (K_U, K_U_v), 
    (M, M_v), (mc, mc_v), (M_fc, M_fc), (r_md, r_md_v),
    (V_D, V_D_v), (V_U, V_U_v), (depth, depth_v)
    ])

eq_C = eq_C_P + eq_C_M + eq_C_D

#%% Calculate equilibrium value for adsorbed C
e_C_A = sy.Eq((k_ads/k_des), C_A / (C_D * (Amax - C_A)))   # Ka = LR / (L * R) = k_ads / k_des
sol_C_A = sy.solve(e_C_A, C_A)[0]

eq_C_A = sol_C_A.subs([(Ka, Ka_v), (C_D, eq_C_D), (Amax, Amax_v), (depth, depth_v)])
eq_C2 = eq_C + eq_C_A


#%%
file = open("python_out.txt", "a")

file.write("\n--------------\n" +
           "Time: " + time.strftime('%Y/%m/%d %H:%M:%S') + "\n\n" +
           "Options \n" + 
           "dec_fun: " + str(dec_fun) + " , upt_fun: " + str(upt_fun) +
           "\n\n" + "Solutions" + "\n\n" + 
           "C_P \n" + str(sol_C_P) + "\n\n" + 
           "C_D \n" + str(sol_C_D) + "\n\n" +
           "C_M \n" + str(sol_C_M) + "\n\n" +
           "C_A \n" + str(sol_C_A) + "\n\n" + "\n\n")
file.close()
