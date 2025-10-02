# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 12:04:23 2025

@author: sarve
"""

import v2_fea2d_heat
from plotting_tools import plot_Heat2D_results
import numpy as np
import matplotlib.pyplot as plt

#########################Usage example for reading thermal problem#########################
#filename = r"X:\UIUC\Sem 4\ME 471 FEA\Homeworks\Project 3\in progress\thermal_inputs\fine_rectangle_hole.inp"
filename = r"X:\UIUC\Sem 4\ME 471 FEA\Homeworks\Project 3\in progress\thermal_inputs\coarse_rectangle.inp"
[N_NODE, N_ELEM, NNODE_ELE, ShapeOrder, Ng, N_PRE_T, 
N_FLUX_q, N_FLUX_c,COORDS, ELEM_NODE, ELEM_PROPS, ELEM_BLOAD, 
T_NODE, T_VAL,FLUX_q_ELE, FLUX_q_VAL,FLUX_c_ELE, FLUX_c_VAL] = v2_fea2d_heat.make_thermal_problem_data(filename)

'''
# Print all outputs
print('N_NODE = ', N_NODE)
print('N_ELEM = ', N_ELEM)
print('NNODE_ELE = ', NNODE_ELE)
print('ShapeOrder = ', ShapeOrder)
print('Ng = ', Ng)
print('N_PRE_T = ', N_PRE_T)
print('N_FLUX_q = ', N_FLUX_q)
print('N_FLUX_c = ', N_FLUX_c)
print('COORDS = \n', COORDS)
print('ELEM_NODE = \n', ELEM_NODE)
print('ELEM_PROPS = \n', ELEM_PROPS)
print('ELEM_BLOAD = \n', ELEM_BLOAD)
print('T_NODE = \n', T_NODE)
print('T_VAL = \n', T_VAL)
print('FLUX_q_ELE = \n', FLUX_q_ELE)
print('FLUX_q_VAL = \n', FLUX_q_VAL)
print('FLUX_c_ELE = \n', FLUX_c_ELE)
print('FLUX_c_VAL = \n', FLUX_c_VAL)
'''
##########################End of reading thermal problem##################################

#steady state

TUR, EQ_NUM, LM, K, KPP, KPF, KFP, KFF, PP, PF, UP, UF = v2_fea2d_heat.solve_steady(N_NODE, N_ELEM, NNODE_ELE, ShapeOrder, Ng, N_PRE_T, N_FLUX_q, N_FLUX_c,COORDS, ELEM_NODE, ELEM_PROPS, ELEM_BLOAD, T_NODE, T_VAL,FLUX_q_ELE, FLUX_q_VAL,FLUX_c_ELE, FLUX_c_VAL)
plot_Heat2D_results(TUR, N_ELEM, ELEM_NODE, COORDS, field='temperature', show_plot=True, filename='temperature_contour_example.png')



# transient
# theta and solution time need to be modified in fea2d_heat transient solver code. 
'''
factor = 0.1
times, TUR = v2_fea2d_heat.general_theta_solver(N_NODE, N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, 
                     T_NODE, T_VAL, FLUX_c_ELE, FLUX_c_VAL, N_FLUX_c,
                     N_FLUX_q, FLUX_q_ELE, FLUX_q_VAL, N_PRE_T,factor)

idx = np.argmin(np.abs(times - 2000))
#print(f"closest time = {times[idx100]:.3f} s at index {idx100}")

# extract the nodal temperatures at that time
T100 = TUR[:, idx]
plot_Heat2D_results(T100, N_ELEM, ELEM_NODE, COORDS, field='temperature', show_plot=True, filename='temperature_contour_example.png')

T1 = times
T51 = TUR[4, :]
T61 = TUR[5, :]
plt.figure()
plt.plot(T1, T51, '-',  label='Node 5')
plt.plot(T1, T61, '-', label='Node 6')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
#plt.ylim(125, 150) 
plt.title('Forward‐Euler, 0.1Δt_cr')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

factor = 0.9
times, TUR = v2_fea2d_heat.general_theta_solver(N_NODE, N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, 
                     T_NODE, T_VAL, FLUX_c_ELE, FLUX_c_VAL, N_FLUX_c,
                     N_FLUX_q, FLUX_q_ELE, FLUX_q_VAL, N_PRE_T,factor)
T2 = times
T52 = TUR[4, :]
T62 = TUR[5, :]


plt.figure()
plt.plot(T2, T52, '-',  label='Node 5')
plt.plot(T2, T62, '-', label='Node 6')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
#plt.ylim(125, 150) 
plt.title('Forward‐Euler, 0.9Δt_cr')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

factor = 1.1
times, TUR = v2_fea2d_heat.general_theta_solver(N_NODE, N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, 
                     T_NODE, T_VAL, FLUX_c_ELE, FLUX_c_VAL, N_FLUX_c,
                     N_FLUX_q, FLUX_q_ELE, FLUX_q_VAL, N_PRE_T,factor)
T3 = times
T53 = TUR[4, :]
T63 = TUR[5, :]

plt.figure()
plt.plot(T3, T53, '-',  label='Node 5')
plt.plot(T3, T63, '-', label='Node 6')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
#plt.ylim(100, 200) 
plt.title('Forward‐Euler, 1.1Δt_cr')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
'''
factor = 0.1
times, TUR = v2_fea2d_heat.general_theta_solver(N_NODE, N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, 
                     T_NODE, T_VAL, FLUX_c_ELE, FLUX_c_VAL, N_FLUX_c,
                     N_FLUX_q, FLUX_q_ELE, FLUX_q_VAL, N_PRE_T,factor)
'''
T2 = times
T4 = TUR[3, :]
T5 = TUR[4, :]
T455 = TUR[454, :]
plt.figure()
plt.plot(T2, T4, '-',  label='Node 5')
plt.plot(T2, T5, '-', label='Node 6')
plt.plot(T2, T455, '-', label='Node 455')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
#plt.ylim(125, 200) 
plt.title('Crank-Nicoleson time integrator, Δt = 50 s')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
'''

    #plot_Heat2D_results(T5, N_ELEM, ELEM_NODE, COORDS, field='temperature', show_plot=True, filename='temperature_contour_example.png')

##########################End of plotting temperature##################################

#To run this file, run python3 main.py from the workspace.