# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 11:44:15 2021

@author: Group 47
"""
import numpy as np
from matplotlib import pyplot as plt
import math

plt.rcParams['figure.figsize'] = [6.0, 4.0]
plt.rcParams['figure.dpi'] = 250


# COMBUSTION CHAMBER VALUES
D_c =30e-3                      # Diameter of the Combustion Chamber
L_c = 100e-3                    # Length of the Combustion Chamber
P0_c = 1e6                    # Initial pressure in the Combustion Chamber (before detonation)
C = 12e-3                       # Charge Mass 5 to 20 g

gamma_c = 1.2238                # γ Gamma for the combustion products

BR_exp = 0.81837                # BURN RATE Exponent
u1 = 3.102e-8                   # Burning Rate Constant
e1 = 1.27e-4                    # Propellant half web size - assumes hollow cylindrical powder
R1 = u1 / e1                    # Experimentally determined Burn Rate Coefficient
density_propel = 1.6e3          # δ Propellant Density
ForceConst_propel = 1.158e6     # λ Propellant Force Constant
CoVolume_propel = 0.8e-3        # η Propellant Co-Volum
V0_c = L_c * math.pi * (D_c/2)**2
  
# PISTON VALUES
D_pis = 0.03            # Diameter of piston
m_pis = 110e-3           # Mass of the piston
mu_static_pis = 0     # Coefficient of friction for the piston against the pump tube
mu_dynamic_pis = 0
allowPistonRearward = False

# PUMP TUBE VALUES
P0_pt = 0.7e6   # Initial pressure in the pump tube (ahead of the piston)
L0_pt = 0.5      # Length of the pump tube
D_pt = D_pis     # Diameter of the pump tube

# RUPTURE DISK VALUES
P_rupt = 65e+6  # Pressure at which the rupture disk ruptures

# BARREL VALUES
L_b = 1.9         # Length of the barrel
D_b = 12.7e-3      # Diameter of the barrel
P0_b = 0.1      # Initial upstream pressure in the barrel
gamma_lg = 1.4  # Gamma for the light gas

# IMPACT CHAMBER
gamma_ic = 1.4  # Gamma for the (near vaccum) gas in the impact chamber
V_ic = 1 

# PROJECTILE & SABOT VALUES
m_pr = 17e-3     # Projectile mass
m_sb = 0.47e-3  # Sabot Mass
mu_sb = 0.4     # Sabot Friction coeff.

# SIMULATION DETAILS
delta_t = 1e-6  # Time step length

# Calculating areas from diameters
A_c = np.pi * (D_c/2)**2
A_pt = np.pi * (D_pt/2)**2
V0_pt = A_pt * L0_pt
A_b = np.pi * (D_b/2)**2
A_pis = np.pi * (D_pis/2)**2


# | COMBUSTION |   PISTON   | RUPTURE DISK | PROJECTILE |
# |    _c      | _c/_pis/_pt|     _rup     |    _pr     |

# _c   = combustion
# _pis = piston
# _pt  = pump tube
# _pr  = projectile
# _sb  = sabot

# t = time
# T = temperature
# P = Pressure
# v = velocity
# V = Volume
# a = acceleration
# x = displacement
def DoIt(A_c = A_c, L_c = L_c, P0_c = P0_c, gamma_c = gamma_c,\
         A_pis = A_pis, mu_static_pis = mu_static_pis, mu_dynamic_pis = mu_dynamic_pis, \
         P0_pt = P0_pt, L0_pt = L0_pt, A_pt = A_pt,\
         P_rupt = P_rupt,\
         L_b =  L_b, D_b = D_b, P0_b = P0_b, gamma_lg = gamma_lg,\
         m_pr = m_pr, m_sb = m_sb, mu_sb = mu_sb,\
         gamma_ic = gamma_ic, V_ic = V_ic,\
         delta_t = delta_t):
    # DEFINING ARRAYS
    """
    any initial value that's in one of these arrays that is NOT a 
    variable will need to be replaced, I've just put in some placeholder
    values so that I can test the program.
    - Euan
    """
    global n_array, t_array, x_pis_array, x_pr_array, n_disk_rupture, P_c_array, \
            v_pis_array, v_pr_array, P_pt_array, Z_c_array, t_disk_rupture, \
            Zr_c_array, n_burnout, t_burnout
    
    P_pt_array = [P0_pt]        # PUMP TUBE Pressure Ahead of Piston
    P_b_array = [1e-3]          # BARREL PRESSURE Behind the projectile
    x_pis_array = [0]           # PISTON DISPLACEMENT from initial position
    x_pr_array = [0]            # PROJECTILE DISPLACEMENT from initial position
    v_pis_array = [0]           # PISTON VELOCITY
    v_pr_array = [0]            # PROJECTILE VELOCITY
    a_pis_array = [0]           # PISTON ACCELERATION
    a_pr_array = [0]            # PROJECTILE ACCELERATION
    Vol_c_array = [V0_c]           # PISTON/COMB. VOLUME of chamber
    Vol_pis_array = [V0_pt]         # PISTON/BARREL VOLUME behind projectile
    n_array = []                # Counts the time step number
    t_array = []                # Holds the time step value
     
    # COMBUSTION CHAMBER
    P_c_array = [P0_c] # Combustion chamber pressure array, starts with pressure for now 
   
    Zr_c_array = [0]
    Z0_br = P0_pt * (V0_c - C / density_propel) / (C * (ForceConst_propel + P0_pt * (CoVolume_propel - 1 / density_propel)))
    Z_c_array = [Z0_br]          # POWDER BURN Decimal - may need to start at z0
    
    # Event RUPTURE BOOLEANS
    global diskBroken, diskJustBroken, burnoutTF
    diskBroken = False
    diskJustBroken = True
    burnoutTF = False 

    def combustElement(S):
        """ Handles the current pressure in the combustion chamber.
            Current time drives the initial expansion and increase in pressure.
            Position of piston drives the decrease in pressure due to 
            expansion. """
        dz_dt = R1 * P_c_array[-1]**BR_exp                  # Burn Rate
        Z_cur = Z_c_array[-1] + delta_t * dz_dt             # Current Burnt decimal %
        Z_c_array.append(Z_cur)                             # Append to Powder Burn array
        Zr_c_array.append(dz_dt)                            # Append current Burn Rate
        
        if len(x_pis_array) == 1: dx = 0                    # If this is the first value then dx is 0
        else: dx = x_pis_array[-2] - x_pis_array[-1]        # Calculate the previous movement in piston to find work done
        
        S = S + dx * P_pt_array[-1]                         # Work done on piston (energy removed from gas)
        # Combustion Pressure using energy balance between internal, work done and PV
        P_c = (ForceConst_propel * C * Z_cur - (gamma_c-1)*(m_pis / 2 * v_pis_array[-1]**2 + A_c * S)) \
            / (V0_c + A_c * x_pis_array[-1] - C / density_propel - (CoVolume_propel - 1 / density_propel)* C * Z_cur)
        
        P_c_array.append(P_c)                               # Append the Combustion Pressure
        
    def burnoutElement():
        """ Handles the combustion chamber once the powder has been burnt """ 
        global n_burnout, t_burnout
        try: 
            n_burnout
        except:
            n_burnout = n_array[-1]
            t_burnout = t_array[-1]
            
        Z_c_array.append(1)
        Zr_c_array.append(0)
        P_c_array.append(P_c_array[-1])
    
    
    def pistonElement():
        """ Handles a single time step element of the piston along the
            pump tube."""    
        delta_P = P_c_array[-1] - P_pt_array[-1]                                # Identify the Pressure Differential across Piston
        F_pressure = delta_P*A_pis                                              # Identify the pressure force on piston
        
        # Check if we need to consider static or dynamic friction -------------
        if v_pis_array[-1] == 0:                                                # Static or Dynamic friction - are we stationary
            F_fric = mu_static_pis                                              # Calculate Static Friction
            F_res = (abs(F_pressure) - abs(F_fric)) * (F_pressure/abs(F_pressure))  # Calculate the resultant force with friction opposing it
            if abs(F_pressure) > F_fric:                                        # Will pressure force overcome the friction force?
                a_pis = F_res/m_pis                                             # Calculate acceleration from resultance force
            else:
                a_pis = 0                                                       # If pressure cant overcome static friction then acceleration is 0
        else:
            F_fric = mu_dynamic_pis                              
            F_res = (abs(F_pressure) - abs(F_fric)) * (F_pressure/abs(F_pressure))  # Calculate the resultant force with friction opposing it
            a_pis = F_res/m_pis                                                 # Calculate the acceleration of the piston
        
        # Check if were at the backstop ---------------------------------------
        if x_pis_array[-1] <= 0 and a_pis < 0 and allowPistonRearward == False:
            # Check if the piston is against the backstop and dont allow it to accelerate rearward
            a_pis_array.append(0)                                               # Acceleration is 0
            v_pis_array.append(0)                                               # Velocity of piston is 0
            x_pis_array.append(x_pis_array[-1])                                 # Displacement of the piston remains the same
        else:
            # Resultant force is forward
            v_pis = v_pis_array[-1] + a_pis*delta_t
            x_pis = v_pis_array[-1]*delta_t + 0.5*a_pis*(delta_t**2) + x_pis_array[-1]
            a_pis_array.append(a_pis)
            v_pis_array.append(v_pis)
            x_pis_array.append(x_pis)
        
        P_pt = P_pt_array[0]*np.power((L0_pt - x_pis_array[0])/(L0_pt - x_pis_array[-1]), gamma_lg)
        P_pt_array.append(P_pt)
        
        
    

    
    def barrelElement():
        """ Handles a single time step element of the projectile/sabot 
            combination along the barrel.
            When this function is run, each time step will have two entries 
            for the pump tube pressure. """
        # Movement of Sabot/Projectile
        delta_P_pr = P_pt_array[-1] - P_b_array[-1]
        F_pressure = delta_P_pr*A_b
        F_fric = mu_sb*0 # change later
        F_res = F_pressure - F_fric
        a_pr = F_res/(m_pr + m_sb)
        a_pr_array.append(a_pr)
        v_pr = v_pr_array[-1] + a_pr*delta_t
        x_pr = v_pr_array[-1]*delta_t + 0.5*a_pr*(delta_t**2) + x_pr_array[-1]
        v_pr_array.append(v_pr)
        x_pr_array.append(x_pr)
        P_b = P_b_array[0]*np.power((V_ic + A_b*(L_b - x_pr_array[0]))/(V_ic + A_b*(L_b - x_pr_array[-1])), gamma_ic)
        P_b_array.append(P_b)
        
    # Running the Discrete Model ----------------------------------------------
    i = 0               
    n_array.append(i)
    t_array.append(i*delta_t)
    S = 0
    while x_pr_array[-1] < L_b:
        if P_pt_array[-1] < P_rupt and diskBroken == False:
            # Disk unbroken, just the pump tube
            if Z_c_array[-1] < 1:
                combustElement(S)
            else:
                burnoutElement()
            pistonElement()
            x_pr_array.append(0)
            v_pr_array.append(0)
            a_pr_array.append(0)

        elif diskJustBroken == True:
            # Disk broken - pump tube and sabot dynamics
            diskBroken = True
            diskJustBroken = False
            n_disk_rupture = n_array[-1]
            t_disk_rupture = t_array[-1]
            
            if Z_c_array[-1] < 1: 
                combustElement(S)
            else:
                burnoutElement()
            barrelElement()
        else:
            # Disk broken- pump tube and sabot dynamics
            if Z_c_array[-1] < 1: 
                combustElement(S)
            else:
                burnoutElement()
            barrelElement()
        
        i += 1
        n_array.append(i)
        t_array.append(i*delta_t)
        
        if i*delta_t > 2:                 # If projectile hasnt left the barrel after 2 seconds then something is wrong
            n_disk_rupture = 0
            t_disk_rupture = 0
            break
            
   
    print('Projectile Exit Velocity = ', round(v_pr_array[-1],3), 'mps')
    print('Peak Piston Velocity =     ', round(max(v_pis_array),3), 'mps')
    print('Maximum P_c : ', round(max(P_c_array)/1e6,3), 'MPa')
    print('Maximum P_pt: ', round(max(P_pt_array)/1e6,3), 'MPa') 
    print('Disk Ruptured after: {}s (n = {})'.format(t_disk_rupture, n_disk_rupture))
    print('Burnout after: {}s (n = {})'.format(t_burnout, n_burnout))         
    return n_array, t_array, x_pis_array, x_pr_array, n_disk_rupture, P_c_array

data = DoIt() # Data format: n_array, t_array, x_pis_array, x_p_array, n_disk_rupture, P_c_array

# ------------------------ Displacement-Time Plots --------------------------------
fig_DT = plt.figure()               # Create Figure
fig_DT.suptitle('Displacement vs Time') # Set Figure Title
ax_DT = fig_DT.add_subplot()        # Add axes to figure
    
ax_DT.set_xlabel('Time (s)')        # Set x label
ax_DT.set_ylabel('Displacement (m)')   # Set y label

# Plot Pressures with time   
ax_DT.plot(t_array, x_pis_array,  label='Piston')                        
ax_DT.plot(t_array[0:len(x_pr_array)], x_pr_array, label='Projectile') 
ax_DT.grid()                        # Apply a grid to plot area
ax_DT.legend()                      # Enable Legends
ax_DT.axvline(t_disk_rupture, color='grey', linestyle='--')
ax_DT.axvline(t_burnout, color='red', linestyle='--')
# -----------------------------------------------------------------------------

# ------------------------ Velocity-Time Plots --------------------------------
fig_vT = plt.figure()               # Create Figure
fig_vT.suptitle('Velocity vs Time') # Set Figure Title
ax_vT = fig_vT.add_subplot()        # Add axes to figure
    
ax_vT.set_xlabel('Time (s)')        # Set x label
ax_vT.set_ylabel('Velocity (m/s)')   # Set y label
                       # Apply a grid to plot area
# Plot Velocities with time   
ax_vT.plot(t_array, v_pr_array, label='Projectile')                        
ax_vT.plot(t_array, v_pis_array, label='Piston') 
ax_vT.grid()                        # Apply a grid to plot area
ax_vT.legend()                      # Enable Legends
ax_vT.axvline(t_disk_rupture, color='grey', linestyle='--')
ax_vT.axvline(t_burnout, color='red', linestyle='--')
# -----------------------------------------------------------------------------

# ------------------------ Pressure-Time Plots --------------------------------
fig_PT = plt.figure()               # Create Figure
fig_PT.suptitle('Pressure vs Time') # Set Figure Title
ax_PT = fig_PT.add_subplot()        # Add axes to figure
ax_PT.set_yscale('log')             # Use logarithmic Y scale
    
ax_PT.set_xlabel('Time (s)')        # Set x label
ax_PT.set_ylabel('Pressure (Pa)')   # Set y label

# Plot Pressures with time   
ax_PT.plot(t_array, P_c_array,  label='Combustion')                        
ax_PT.plot(t_array, P_pt_array, label='Pump Tube') 
ax_PT.grid()                        # Apply a grid to plot area
ax_PT.legend()                      # Enable Legends
ax_PT.axhline(y=P_rupt, color='k')
ax_PT.axvline(t_disk_rupture, color='grey', linestyle='--')
ax_PT.axvline(t_burnout, color='red', linestyle='--')
# ----------------------------------------------------------------------------- 

# ---------------------------- Burn-Time Plots --------------------------------
fig_BT = plt.figure()               # Create Figure
fig_BT.suptitle('Powder Burn vs Time') # Set Figure Title
ax_BT = fig_BT.add_subplot()        # Add axes to figure
    
ax_BT.set_xlabel('Time (s)')        # Set x label
ax_BT.set_ylabel('Powder Burn (decimal %)')   # Set y label
                       # Apply a grid to plot area
# Plot Velocities with time   
ax_BT.plot(t_array[0:len(Z_c_array)], Z_c_array, label='Powder Burnt')                        
ax_BT.grid()                        # Apply a grid to plot area
ax_BT.legend()                      # Enable Legends
ax_BT.axvline(t_disk_rupture, color='grey', linestyle='--')
ax_BT.axvline(t_burnout, color='red', linestyle='--')
# ----------------------------------------------------------------------------- 

# ------------------------- BurnRate-Time Plots -------------------------------
fig_BrT = plt.figure()               # Create Figure
fig_BrT.suptitle('Powder Burn Rate vs Time') # Set Figure Title
ax_BrT = fig_BrT.add_subplot()        # Add axes to figure
    
ax_BrT.set_xlabel('Time (s)')        # Set x label
ax_BrT.set_ylabel('Powder Burn Rate(decimal % per second)')   # Set y label
                       # Apply a grid to plot area
# Plot Velocities with time   
ax_BrT.plot(t_array[0:len(Zr_c_array)], Zr_c_array, label='Powder Burn Rate')                        
ax_BrT.grid()                        # Apply a grid to plot area
ax_BrT.legend()                      # Enable Legends
ax_BrT.axvline(t_disk_rupture, color='grey', linestyle='--')
ax_BrT.axvline(t_burnout, color='red', linestyle='--')
# -----------------------------------------------------------------------------   
