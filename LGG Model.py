# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 11:44:15 2021

@author: Group 47
"""
import numpy as np
from matplotlib import pyplot as plt
import timeit



# COMBUSTION CHAMBER VALUES
D_c = 100e-3 # Diameter of the Combustion Chamber
L_c = 0.20   # Length of the Combustion Chamber
P_c0 = 5e+6 # Initial pressure in the Combustion Chamber (before detonation)
gamma_c = 1.4 # Gamma for the combustion products

# PISTON VALUES
D_pis = 12e-3 
m_pis = 600e-3 # Mass of the piston
mu_pis = 0.4   # Coefficient of friction for the piston against the pump tube

# PUMP TUBE VALUES
P_t0 = 300e+3  # Initial pressure in the pump tube (ahead of the piston)
L_t0 = 0.8     # Length of the pump tube
D_t = D_pis  # Diameter of the pump tube

# RUPTURE DISK VALUES
P_rupt = 16e+6 # Pressure at which the rupture disk ruptures

# BARREL VALUES
L_b = 1   # Length of the barrel
D_b = 5e-3 # Diameter of the barrel
P_b0 = 0.1 # Initial upstream pressure in the barrel
gamma_lg = 1.4 # Gamma for the light gas

# IMPACT CHAMBER
gamma_ic = 1.4 # Gamma for the (near vaccum) gas in the impact chamber
V_ic = 1 

# PROJECTILE & SABOT VALUES
m_p = 0.01e-3
m_s = m_p
mu_s = 0.4

# SIMULATION DETAILS
delta_t = 1e-6 # Time step length

# Calculating areas from diameters
A_c = np.pi * (D_c/2)**2
A_t = np.pi * (D_t/2)**2
A_b = np.pi * (D_b/2)**2
A_pis = np.pi * (D_pis/2)**2


def DoIt(A_c = A_c, L_c = L_c, P_c0 = P_c0, gamma_c = gamma_c,\
         A_pis = A_pis, m_pis = m_pis, mu_pis = mu_pis, \
         P_t0 = P_t0, L_t0 = L_t0, A_t = A_t,\
         P_rupt = P_rupt,\
         L_b =  L_b, D_b = D_b, P_b0 = P_b0, gamma_lg = gamma_lg,\
         m_p = m_p, m_s = m_s, mu_s = mu_s,\
         gamma_ic = gamma_ic, V_ic = V_ic,\
         delta_t = delta_t):
    
    start = timeit.default_timer()  # Start Timer
    
    # DEFINING ARRAYS
    """
    any initial value that's in one of these arrays that is NOT a 
    variable will need to be replaced, I've just put in some placeholder
    values so that I can test the program.
    - Euan
    """
    P_c_array = [P_c0]
    P_t_array = [P_t0]
    P_b_array = [1e-3]
    x_pis_array = [0.01]
    x_p_array = [0.01]
    v_pis_array = [0]
    v_p_array = [0]
    a_pis_array = [0]
    a_p_array = [0]
    V_pis_array = [0]
    V_p_array = [0]
    n_array = [] # Counts the time step number
    t_array = [] # Holds the time step value
    
    # COMBUSTION CHAMBER
    P_c_array = [700e+5] # 
    
    # DISK RUPTURE BOOLEANS
    diskBroken = False
    diskJustBroken = True
    
    """
    Handles a single time step element of the piston along the pump tube.
    """
    def pistonElement():
        delta_P = P_c_array[-1] - P_t_array[-1]
        F_pressure = delta_P*A_pis
        F_fric = mu_pis*0  # change later
        F_res = F_pressure - F_fric
        a_pis = F_res/m_pis
        a_pis_array.append(a_pis)
        v_pis = v_pis_array[-1] + a_pis*delta_t
        x_pis = v_pis_array[-1]*delta_t + 0.5*a_pis*(delta_t**2) + x_pis_array[-1] 
        v_pis_array.append(v_pis)
        x_pis_array.append(x_pis)
        
        P_c = P_c_array[0]*np.power(x_pis_array[0]/x_pis_array[-1], gamma_c)
        P_c_array.append(P_c)
        P_t = P_t_array[0]*np.power((L_t0 - x_pis_array[0])/(L_t0 - x_pis_array[-1]), gamma_lg)
        P_t_array.append(P_t)
        
        if v_pis_array[-1] < 0:
            raise ValueError("Chamber pressure after detotation was insufficient to keep piston moving forwards up to disk rupture against the pressure of the lg and the piston began moving backwards. Simulation ended.")
    

    """
    Handles a single time step element of the projectile/sabot combination
    along the barrel.
    When this function is run, each time step will have two entries for the
    pump tube pressure.
    """
    def barrelElement():
        # Movement of Sabot/Projectile
        delta_P_p = P_t_array[-1] - P_b_array[-1]
        F_pressure = delta_P_p*A_b
        F_fric = mu_s*0 # change later
        F_res = F_pressure - F_fric
        a_p = F_res/(m_p + m_s)
        a_p_array.append(a_p)
        v_p = v_p_array[-1] + a_p*delta_t
        x_p = v_p_array[-1]*delta_t + 0.5*a_p*(delta_t**2) + x_p_array[-1]
        v_p_array.append(v_p)
        x_p_array.append(x_p)
        
        # Movement of Piston
        delta_P_pis = P_c_array[-1] - P_t_array[-1]
        F_pressure = delta_P_pis*A_pis
        F_fric = mu_pis*0  # change later
        F_res = F_pressure - F_fric
        a_pis = F_res/m_pis
        a_pis_array.append(a_pis)
        v_pis = v_pis_array[-1] + a_pis*delta_t
        x_pis = v_pis_array[-1]*delta_t + 0.5*a_pis*(delta_t**2) + x_pis_array[-1] 
        v_pis_array.append(v_pis)
        x_pis_array.append(x_pis)
        
        # Calculation of (changed) pressures
        P_c = P_c_array[0]*np.power(x_pis_array[0]/x_pis_array[-1], gamma_c)
        P_c_array.append(P_c)
        P_t = P_t_array[0]*np.power((A_pis*(L_t0 - x_pis_array[0]))/(A_pis*(L_t0 - x_pis_array[-1]) + A_b*x_p_array[-1]), gamma_lg)
        P_t_array.append(P_t)
        P_b = P_b_array[0]*np.power((V_ic + A_b*(L_b - x_p_array[0]))/(V_ic + A_b*(L_b - x_p_array[-1])), gamma_ic)
        P_b_array.append(P_b)
        
        
    
    i = 0
    n_array.append(i)
    t_array.append(i*delta_t)
    while x_p_array[-1] < L_b:
        if P_t_array[-1] < P_rupt and diskBroken == False:
            # Disk unbroken, just the pump tube
            pistonElement()
            #print(P_t_array[-1])
        elif diskJustBroken == True:
            # Disk broken - pump tube and sabot dynamics
            print("Barrel Reached at n = "+str(n_array[-1])+" (or "+str(t_array[-1])+"s)")
            diskBroken = True
            diskJustBroken = False
            n_disk_rupture = n_array[-1]
            barrelElement()
        else:
            # Disk broken- pump tube and sabot dynamics
            barrelElement()
        
        i += 1
        n_array.append(i)
        t_array.append(i*delta_t)
    
    """
    # Seperating out the pump tube pressure disk rupture
    P_t_barrelElement = []
    P_t_pistonElement = []
    P_t_average = P_t_array[0:n_disk_rupture]
    
    for i in range(n_disk_rupture, n_disk_rupture + len(n_array[n_disk_rupture:-1])):
        P_t_barrelElement.append(P_t_array[i])
        P_t_pistonElement.append(P_t_array[i+1])
        P_t_average.append(0.5*(P_t_array[i]+ P_t_array[i+1]))
     """   
     
    stop = timeit.default_timer()           # Stops Timer
    time_to_run = stop - start              # Calculates Code Run Time
    
    print(v_p_array[-1])
    print(L_b)
            
    return n_array, t_array, x_pis_array, x_p_array, n_disk_rupture, v_pis_array, v_p_array, time_to_run

data = DoIt()


"""
fig, ax = plt.subplots(constrained_layout=True)
ax.plot(data[0], data[2], color = "blue")
ax.tick_params(axis = 'y', labelcolor = "blue")
ax.set_xlabel("Step number")
ax.set_ylabel("Position [m]", color = "blue")
ax.axhline(y = L_t0, color = "black", ls = '--')
"""

"""
ax2 = ax.twinx()
ax2.plot(data[0],data[3], color =  "red")
ax2.tick_params(axis = 'y', labelcolor = "red")
ax2.set_ylabel("Acceleration [m/s^2]", color = "red")
"""

"""
fig, ax3 = plt.subplots(constrained_layout=True)
ax3 = ax2.twinx()
ax3.plot(data[4], data[6], color = "green")
ax3.tick_params(axis = 'y', labelcolor = "green")
ax3.set_ylabel("Velocity", color = "green")
plt.show()
"""

"""
fig, ax = plt.subplots(constrained_layout=True)
ax.plot(data[4], data[1], color = "blue")
ax.tick_params(axis = 'y', labelcolor = "blue")
ax.set_xlabel("Step number")
ax.set_ylabel("P_c (behind)", color = "blue")
ax2 = ax.twinx()
ax2.plot(data[4], data[0], color =  "red")
ax2.tick_params(axis = 'y', labelcolor = "red")
ax2.set_ylabel("P_t (infront)", color = "red")
plt.show()
"""

"""
p_pos_to_plot = []
for i in range(0, len(data[3])):
    p_pos_to_plot.append(data[3][i] + data[2][-1])

plt.plot(data[0], data[2])
plt.plot(data[0][data[4]:], p_pos_to_plot)
plt.legend(["Piston position", "Projectile Position"])
"""







""" SENSITIVITY STUDY """





res = 20
result = np.zeros(res)


D_pis_arr = np.linspace(1e-3,150e-3,res)
result_1 = np.zeros(res)
result_2 = np.zeros(res)
result_3 = np.zeros(res)
result_4 = np.zeros(res)
result_5 = np.zeros(res)
for i in range(res):
    D_t = D_pis_arr[i]
    A_t = np.pi * (D_t/2)**2
    D_pis = D_pis_arr[i]
    A_pis = np.pi * (D_pis_arr[i]/2)**2
    result_1[i] = DoIt(A_pis=A_pis, A_t=A_t, m_pis=1e-3)[6][-1]
    result_2[i] = DoIt(A_pis=A_pis, A_t=A_t, m_pis=10e-3)[6][-1]
    result_3[i] = DoIt(A_pis=A_pis, A_t=A_t, m_pis=50e-3)[6][-1]
    result_4[i] = DoIt(A_pis=A_pis, A_t=A_t, m_pis=100-3)[6][-1]
    result_5[i] = DoIt(A_pis=A_pis, A_t=A_t, m_pis=500e-3)[6][-1]

plt.figure()
plt.plot(D_pis_arr, result_1, label='m_pis = 1g')
plt.plot(D_pis_arr, result_2, label='m_pis = 10g')
plt.plot(D_pis_arr, result_3, label='m_pis = 50g')
plt.plot(D_pis_arr, result_4, label='m_pis = 100g')
plt.plot(D_pis_arr, result_5, label='m_pis = 500g')
plt.xlabel('Piston Diameter (m)')
plt.ylabel('Exit velocity (m/s)')
plt.grid()
plt.legend(bbox_to_anchor=(-0.159, 1, 1, 0), loc="lower left", ncol=3,columnspacing = 1)
plt.show()




D_pis_arr = np.linspace(1e-3,150e-3,res)
result_1 = np.zeros(res)
result_2 = np.zeros(res)
result_3 = np.zeros(res)
result_4 = np.zeros(res)
result_5 = np.zeros(res)
for i in range(res):
    D_t = D_pis_arr[i]
    A_t = np.pi * (D_t/2)**2
    D_pis = D_pis_arr[i]
    A_pis = np.pi * (D_pis_arr[i]/2)**2
    result_1[i] = DoIt(A_pis=A_pis, A_t=A_t, m_p=0.001e-3)[6][-1]
    result_2[i] = DoIt(A_pis=A_pis, A_t=A_t, m_p=0.01e-3)[6][-1]
    result_3[i] = DoIt(A_pis=A_pis, A_t=A_t, m_p=0.1e-3)[6][-1]
    result_4[i] = DoIt(A_pis=A_pis, A_t=A_t, m_p=1e-3)[6][-1]
    result_5[i] = DoIt(A_pis=A_pis, A_t=A_t, m_p=10e-3)[6][-1]

plt.figure()
plt.plot(D_pis_arr, result_1, label='m_p = 0.001g')
plt.plot(D_pis_arr, result_2, label='m_p = 0.01g')
plt.plot(D_pis_arr, result_3, label='m_p = 0.1g')
plt.plot(D_pis_arr, result_4, label='m_p = 1g')
plt.plot(D_pis_arr, result_5, label='m_p = 10g')
plt.plot(108e-3, 2481, 'x', color='blue',label='CAS ~ 19g')
plt.plot(108e-3, 2260, 'x', color='blue',label='CAS ~ 100g')
plt.plot(12.7e-3, 5700, 'x', color = 'black',label='Kent Uni H2 ~ 1g')
plt.plot(1.42*0.0254, 11000*0.3048, 'x', color = 'red', label = 'NMIMT ~ 2.46g') #https://aip.scitation.org/doi/pdf/10.1063/1.1722882
plt.plot(1.42*0.0254, 11070*0.3048, 'x', color = 'red', label = 'NMIMT ~ 3.51g')
plt.plot(1.42*0.0254, 11180*0.3048, 'x', color = 'red', label = 'NMIMT ~ 4.42g')
plt.plot(0.0508,3834, 'x',color ='green',label ='UBC ~ 5.894g') #https://open.library.ubc.ca/soa/cIRcle/collections/ubctheses/831/items/1.0098802
plt.plot(27e-3,2040, 'x', color = 'magenta',label = 'ORNL ~ 0.055g') # https://aip.scitation.org/doi/pdf/10.1063/1.1142402
plt.plot(27e-3,3000, 'x', color = 'magenta',label = 'ORNL ~ 0.29g')
plt.xlabel('Piston Diameter (m)')
plt.ylabel('Exit velocity (m/s)')
plt.grid()
plt.legend(bbox_to_anchor=(-0.159, 1, 1, 0), loc="lower left", ncol=3,columnspacing = 1)
plt.show()






m_pis_arr = np.linspace(1e-3, 600e-3, res)
for i in range(res):
    result[i] = DoIt(m_pis=m_pis_arr[i])[6][-1]
plt.figure()
plt.plot(m_pis_arr, result)
plt.xlabel('Piston Mass (kg)')
plt.ylabel('Exit velocity (m/s)')
plt.ylim(0,7000)
plt.savefig("piston_mass.png")



mu_pis_arr = np.linspace(0.1, 2, res)
for i in range(res):
    result[i] = DoIt(mu_pis=mu_pis_arr[i])[6][-1]
plt.figure()
plt.plot(mu_pis_arr, result)
plt.xlabel('Piston Friction Coefficient')
plt.ylabel('Exit Velocity (m/s)')
plt.ylim(0,7000)
plt.savefig("piston_friction_coeff.png")



P_t0_arr = np.linspace(1e+4, 1e+5, res)
for i in range(res):
    result[i] = DoIt(P_t0=P_t0_arr[i])[6][-1]
plt.figure()
plt.plot(P_t0_arr/1e+5, result)
plt.xlabel('Initial Pressure in Pump Tube (bar)')
plt.ylabel('Exit Velocity (m/s)')
plt.ylim(0,7000)
plt.savefig("initial_pump_tube_pressure.png")



L_t0_arr = np.linspace(0.3, 0.8, res)
for i in range(res):
    result[i] = DoIt(L_t0=L_t0_arr[i])[6][-1]
plt.figure()
plt.plot(L_t0_arr, result)
plt.xlabel('Pump Tube Length (m)')
plt.ylabel('Exit Velocity (m/s)')
plt.ylim(0,7000)
plt.savefig("pump_tube_length.png")



P_rupt_arr = np.linspace(1e+5, 200e+5, res)
P_rupt_results_1 = np.zeros(res)
P_rupt_results_2 = np.zeros(res)
P_rupt_results_3 = np.zeros(res)
P_rupt_results_4 = np.zeros(res)
P_rupt_results_5 = np.zeros(res)
for i in range(res):
    P_rupt_results_1[i] = DoIt(P_rupt=P_rupt_arr[i], m_p=0.001e-3)[6][-1]
    P_rupt_results_2[i] = DoIt(P_rupt=P_rupt_arr[i], m_p=0.01e-3)[6][-1]
    P_rupt_results_3[i] = DoIt(P_rupt=P_rupt_arr[i], m_p=0.1e-3)[6][-1]
    P_rupt_results_4[i] = DoIt(P_rupt=P_rupt_arr[i], m_p=1e-3)[6][-1]
    P_rupt_results_5[i] = DoIt(P_rupt=P_rupt_arr[i], m_p=10e-3)[6][-1]
plt.figure()
plt.plot(P_rupt_arr/1e+5, P_rupt_results_1, label='m_p=0.001g')
plt.plot(P_rupt_arr/1e+5, P_rupt_results_2, label='m_p=0.01g')
plt.plot(P_rupt_arr/1e+5, P_rupt_results_3, label='m_p=0.1g')
plt.plot(P_rupt_arr/1e+5, P_rupt_results_4, label='m_p=1g')
plt.plot(P_rupt_arr/1e+5, P_rupt_results_5, label='m_p=10g')
plt.plot(7000e+3/1e+5, 5700, 'x', color = 'black',label='Kent Uni H2 ~ 1g')
#plt.plot(4500e+3, 4300, 'x', color ='black',label='Kent Uni He ~ 1g')
plt.plot(45e+6/1e+5, 2481, 'x', color='blue',label='CAS ~ 19g')
plt.plot(45e+6/1e+5, 2260, 'x', color='blue',label='CAS ~ 100g')
#plt.plot(132000*6894.76/1e+5, 11000*0.3048, 'x', color = 'red', label = 'NMIMT ~ 2.46g')
#plt.plot(110000*6894.76/1e+5, 11070*0.3048, 'x', color = 'red', label = 'NMIMT ~ 3.51g')
#plt.plot(103000*6894.76/1e+5, 11180*0.3048, 'x', color = 'red', label = 'NMIMT ~ 4.42g')
plt.plot(60*6894.75729/1e+5,3834, 'x',color ='green',label ='UBC ~ 5.894g')
plt.plot(62e+5/1e+5,2040, 'x', color = 'magenta',label = 'ORNL ~ 0.055g')
plt.plot(100e+5/1e+5,3000, 'x', color = 'magenta',label = 'ORNL ~ 0.29g')
plt.xlabel('Rupture Pressure (bar)')
plt.ylabel('Exit Velocity (m/s)')
plt.grid()
plt.legend(bbox_to_anchor=(-0.159, 1, 1, 0), loc="lower left", ncol=3,columnspacing = 1)
plt.show()


P_rupt_results_2_1 = np.zeros(res)
P_rupt_results_2_2 = np.zeros(res)
P_rupt_results_2_3 = np.zeros(res)
P_rupt_results_2_4 = np.zeros(res)
P_rupt_results_2_5 = np.zeros(res)
for i in range(res):
    P_rupt_results_2_1[i] = DoIt(P_rupt=P_rupt_arr[i], m_pis=1e-3)[6][-1]
    P_rupt_results_2_2[i] = DoIt(P_rupt=P_rupt_arr[i], m_pis=10e-3)[6][-1]
    P_rupt_results_2_3[i] = DoIt(P_rupt=P_rupt_arr[i], m_pis=50e-3)[6][-1]
    P_rupt_results_2_4[i] = DoIt(P_rupt=P_rupt_arr[i], m_pis=100e-3)[6][-1]
    P_rupt_results_2_5[i] = DoIt(P_rupt=P_rupt_arr[i], m_pis=500e-3)[6][-1]
plt.figure()
plt.plot(P_rupt_arr/1e+5, P_rupt_results_2_1, label='m_pis=1g')
plt.plot(P_rupt_arr/1e+5, P_rupt_results_2_2, label='m_pis=10g')
plt.plot(P_rupt_arr/1e+5, P_rupt_results_2_3, label='m_pis=50g')
plt.plot(P_rupt_arr/1e+5, P_rupt_results_2_4, label='m_pis=100g')
plt.plot(P_rupt_arr/1e+5, P_rupt_results_2_5, label='m_pis=500g')
plt.xlabel('Rupture Pressure (bar)')
plt.ylabel('Exit Velocity (m/s)')
plt.grid()
plt.legend(bbox_to_anchor=(-0.159, 1, 1, 0), loc="lower left", ncol=3,columnspacing = 1)
plt.show()





L_b_arr = np.linspace(1, 5, res)
for i in range(res):
    result[i] = DoIt(L_b=L_b_arr[i])[6][-1]
plt.figure()
plt.plot(L_b_arr, result)
plt.xlabel('Barrel Length (m)')
plt.ylabel('Exit Velocity (m/s)')
plt.ylim(0,7000)
plt.savefig("barrel_length.png")



D_b_arr = np.linspace(1e-3, 100e-3, res)
for i in range(res):
    A_b = np.pi * (D_b_arr[i]/2)**2
    result[i] = DoIt(D_b=D_b_arr[i])[6][-1]
plt.figure()
plt.plot(D_b_arr, result)
plt.xlabel('Barrel diameter (m)')
plt.ylabel('Exit Velocity (m/s)')
plt.ylim(0,7000)
plt.savefig("barrel_diameter.png")



gamma_lg_arr = np.linspace(1.0, 1.5, res)
for i in range(res):
    result[i] = DoIt(gamma_lg=gamma_lg_arr[i])[6][-1]
plt.figure()
plt.plot(gamma_lg_arr, result)
plt.xlabel('Light Gas Specific Heat Ratio')
plt.ylabel('Exit Velocity (m/s)')
plt.ylim(0,7000)
plt.savefig("light_gas_k.png")



gamma_ic_arr = np.linspace(1.0, 1.5, res)
for i in range(res):
    result[i] = DoIt(gamma_ic=gamma_ic_arr[i])[6][-1]
plt.figure()
plt.plot(gamma_ic_arr, result)
plt.xlabel('Impact Chamber Gas Specific Heat Ratio')
plt.ylabel('Exit Velocity (m/s)')
plt.ylim(0,7000)
plt.savefig("impact_chamber_k.png")



m_p_arr = np.linspace(0.001e-3, 1e-3, res)
for i in range(res):
    result[i] = DoIt(m_p=m_p_arr[i])[6][-1]
plt.figure()
plt.plot(m_p_arr*1000, result)
plt.xlabel('Projectile Mass (g)')
plt.ylabel('Exit Velocity (m/s)')
plt.ylim(0,7000)
plt.savefig("projecitle_mass.png")



plt.figure()
plt.plot(data[3], data[6])
plt.xlabel('x position along barrel (m)')
plt.ylabel('Projectile Velocity (m/s)')



"""
CONVERGENCE PLOTS

t = np.logspace(-8, -5, 20)
v_final = np.zeros(len(t))
run_times = np.zeros(len(t))
number_steps = np.zeros(len(t))
for i in range(len(t)):
    v_final[i] = DoIt(delta_t = t[i])[6][-1]
    run_times[i] = DoIt(delta_t = t[i])[-1]

fig, ax1 = plt.subplots()
ax1.set_xlabel('Time Step (s)')
ax1.set_ylabel('Barrel Exit Velocity (m/s)')
ax1.semilogx(t, v_final)
ax1.tick_params(axis='y')

fig, ax1 = plt.subplots()
ax1.set_xlabel('Time Step (s)')
ax1.set_ylabel('Time to execute code (s)')
ax1.semilogx(t, run_times)
ax1.tick_params(axis='y')
"""

    