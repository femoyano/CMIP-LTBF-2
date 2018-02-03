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
dec_fun   = "MM" # One of: 'MM', '2nd', '1st'
#flag_mmr  = 1 # microbial maintenance respiration
#flag_mic  = 1 # simulate microbial pool explicitly

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
C_P, C_A, C_M = \
    sy.symbols('C_P C_A C_M')


r_md, f_ug, f_ge, r_mr= (sy.symbols('r_md f_ug f_ge r_mr'))
V_D, K_D = (sy.symbols('V_D  K_D'))
M_fc, min_md = (sy.symbols('M_fc min_md'))
M, I_sl, depth = sy.symbols('M I_sl depth')
Ka, k_ads, k_des, Amax = sy.symbols('Ka k_ads k_des Amax')

# Define fluxes

F_slcp = I_sl

if dec_fun == "MM":
      U_cp = depth * V_D * C_P/depth * C_M/depth / (K_D + C_P/depth)    
if dec_fun == "2nd":
      U_cp = depth * V_D * C_P/depth * C_M/depth
if dec_fun == "1st":
      U_cp = V_D * C_P
      

F_cpcm = U_cp * f_ug
F_cpcr = U_cp * (1 - f_ug)

F_cmcp = C_M * (C_M * r_md + min_md)
F_cmcr = C_M * r_mr

F_cpca  = C_P * (Amax - C_A) * k_ads
F_cacp  = C_A * k_des

dC_P = F_slcp + F_cmcp - F_cpcr - F_cpcm + F_cacp - F_cpca
dC_A = F_cpca - F_cacp
dC_M = F_cpcm - F_cmcp - F_cmcr

sol = sy.solve([dC_P, dC_A, dC_M],
                        [C_P, C_A, C_M], dict=True)

sol = sol[0]
sol_C_P = sol[C_P]
sol_C_A = sol[C_A]
sol_C_M = sol[C_M]

#%%

# Site data
clay = 0.15
sand = 0.28
silt = 0.57
ps = 0.45
I_sl_v = 0.00005
depth_v = 0.3

# Intermediate parameter values
E_m = 10
E_ad = 10
E_K = 89
E_V = 87
E_r = 95
K_D_ref = 50
k_ads_ref = 1.08e-6 / sec * tstep
k_des_ref = 1.19e-10 / sec * tstep
pd = 2700
psi_fc = 33
R = 0.008314
r_md_ref = 0.002 / hour * tstep
r_mr_ref = 0.000042
T = 288.15
T_ref = 293.15
V_D_ref = 0.35 / hour * tstep

# End parameter values
f_ug_v = 0.50
f_ge_v = 0.01 / hour * tstep
M_v = 0.2
min_md_v = 0

# Calculate end variables
K_D_v = T_resp_eq(K_D_ref, T, T_ref, E_K, R)
V_D_v = T_resp_eq(V_D_ref, T, T_ref, E_V, R)
r_md_v = T_resp_eq(r_md_ref, T, T_ref, E_m, R)
r_mr_v  = T_resp_eq(r_mr_ref, T, T_ref, E_r, R)
r_md_v = T_resp_eq(r_md_ref, T, T_ref, E_m, R)
r_mr_v  = T_resp_eq(r_mr_ref, T, T_ref, E_r, R)
k_ads_v = T_resp_eq(k_ads_ref, T, T_ref, E_ad, R)
k_des_v = T_resp_eq(k_des_ref, T, T_ref, E_ad, R)
Amax_v = 200 * (100 * clay)**0.6 * pd * (1 - ps) / 1000000 #from mg kg-1 to kg m-3

# Substitute variables (parameters) with values

eq_C_P = sol_C_P.subs([
    (f_ug, f_ug_v), (f_ge, f_ge_v), (min_md, min_md_v),
    (r_mr, r_mr_v), (I_sl, I_sl_v), (K_D, K_D_v), 
    (M, M_v), (M_fc, M_fc), (r_md, r_md_v),
    (V_D, V_D_v), (depth, depth_v),
    (Amax, Amax_v), (k_ads, k_ads_v), (k_des, k_des_v)
    ])
    
eq_C_A = sol_C_A.subs([
    (f_ug, f_ug_v), (f_ge, f_ge_v), (min_md, min_md_v),
    (r_mr, r_mr_v), (I_sl, I_sl_v), (K_D, K_D_v), 
    (M, M_v), (M_fc, M_fc), (r_md, r_md_v),
    (V_D, V_D_v), (depth, depth_v),
    (Amax, Amax_v), (k_ads, k_ads_v), (k_des, k_des_v)
    ]) 

eq_C_M = sol_C_M.subs([
    (f_ug, f_ug_v), (f_ge, f_ge_v), (min_md, min_md_v),
    (r_mr, r_mr_v), (I_sl, I_sl_v), (K_D, K_D_v), 
    (M, M_v), (M_fc, M_fc), (r_md, r_md_v),
    (V_D, V_D_v), (depth, depth_v),
    (Amax, Amax_v), (k_ads, k_ads_v), (k_des, k_des_v)
    ])

eq_C = eq_C_P + eq_C_M + eq_C_A

#%% Calculate equilibrium value for adsorbed C
e_C_A = sy.Eq((k_ads/k_des), C_A / (C_P * (Amax - C_A)))   # Ka = LR / (L * R) = k_ads / k_des
sol_C_A = sy.solve(e_C_A, C_A)[0]

eq_C_A = sol_C_A.subs([(Ka, Ka_v), (C_P, eq_C_P), (Amax, Amax_v), (depth, depth_v)])
eq_C2 = eq_C + eq_C_A


#%%
file = open("python_out.txt", "a")

file.write("\n--------------\n" +
           "Time: " + time.strftime('%Y/%m/%d %H:%M:%S') + "\n\n" +
           "Options \n" + 
           "dec_fun: " + str(dec_fun) + " , no C_D" +
           "\n\n" + "Solutions" + "\n\n" + 
           "C_P \n" + str(sol_C_P) + "\n\n" + 
           "C_M \n" + str(sol_C_M) + "\n\n" +
           "C_A \n" + str(sol_C_A) + "\n\n" + "\n\n")
file.close()
