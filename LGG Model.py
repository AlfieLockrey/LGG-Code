# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 11:44:15 2021

@author: Group 47
"""
import numpy as np
from matplotlib import pyplot as plt

# COMBUSTION CHAMBER VALUES
D_c = 100e-3 # Diameter of the Combustion Chamber
L_c = 0.20   # Length of the Combustion Chamber
P_c0 = 5e+6  # Initial pressure in the Combustion Chamber (before detonation)
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
delta_t = 1e-5 # Time step length

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
    P_c_array = [70e+6] # 
    
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
                
    return n_array, t_array, x_pis_array, x_p_array, n_disk_rupture

data = DoIt()
"""
fig, ax = plt.subplots(constrained_layout=True)
ax.plot(data[4][0:-1], data[2], color = "blue")
ax.tick_params(axis = 'y', labelcolor = "blue")
ax.set_xlabel("Step number")
ax.set_ylabel("Position [m]", color = "blue")
ax.axhline(y = L_t0, color = "black", ls = '--')
ax2 = ax.twinx()
ax2.plot(data[4][0:-1],data[3], color =  "red")
ax2.tick_params(axis = 'y', labelcolor = "red")
ax2.set_ylabel("Acceleration [m/s^2]", color = "red")
"""
"""
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

p_pos_to_plot = []
for i in range(0, len(data[3])):
    p_pos_to_plot.append(data[3][i] + data[2][-1])

plt.plot(data[0], data[2])
plt.plot(data[0][data[4]:], p_pos_to_plot)
plt.legend(["Piston position", "Projectile Position"])