# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 11:44:15 2021

To Do:
    x Add burnout element when powder is burnt and stop increasing Z_c and Zr_c
    / Add volume/pressure relations to burnout element
    o Add proper barrel-pump-tube relations to barrel element
        o Implement friction to barrel-pump-tube
        o Fix issue where projectile doesnt accelerate soon enough
    / Add some form of volume tracking to be able to plot volumes of the CVs
    o Add acceleration graphs
    o Fix where all the random global variables are and make it less janky
    o Find a reasonable value for friction
    o Add some form of rifling losses - kinetic energy losses maybe be negligible
    x Add legend for vertical lines
    o Check piston is slower than speed of sound in helium

o = to do
/ = partially done
x = done

@author: Group 47
"""
import numpy as np
from matplotlib import pyplot as plt
import math

plt.rcParams['figure.figsize'] = [6.0, 4.0]
plt.rcParams['figure.dpi'] = 250


# COMBUSTION CHAMBER VALUES
D_c = 30e-3                     # Diameter of the Combustion Chamber
L0_c = 100e-3                   # Length of the Combustion Chamber
P0_c = 1e6                      # Initial pressure in the Combustion Chamber (before detonation)
C = 12e-3                       # Charge Mass 5 to 20 g

gamma_c = 1.2238                # γ Gamma for the combustion products

BR_exp = 0.81837                # BURN RATE Exponent
u1 = 3.102e-8                   # Burning Rate Constant
e1 = 1.27e-4                    # Propellant half web size - assumes hollow cylindrical powder
R1 = u1 / e1                    # Experimentally determined Burn Rate Coefficient
density_propel = 1.6e3          # δ Propellant Density
ForceConst_propel = 1.158e6     # λ Propellant Force Constant
CoVolume_propel = 0.8e-3        # η Propellant Co-Volum
V0_c = L0_c * math.pi * (D_c / 2)**2

# PISTON VALUES
D_pis = 0.03            # Diameter of piston
m_pis = 110e-3          # Mass of the piston
mu_static_pis = 0       # Coefficient of friction for the piston against the pump tube
mu_dynamic_pis = 0
allowPistonRearward = False

# PUMP TUBE VALUES
P0_pt = 0.7e6       # Initial pressure in the pump tube (ahead of the piston)
L0_pt = 0.5         # Length of the pump tube
D_pt = D_pis        # Diameter of the pump tube

# RUPTURE DISK VALUES
P_rupt = 65e+6  # Pressure at which the rupture disk ruptures

# BARREL VALUES
L_b = 1.9           # Length of the barrel
D_b = 12.7e-3       # Diameter of the barrel
P0_b = 0.1          # Initial upstream pressure in the barrel
gamma_lg = 1.4      # Gamma for the light gas

# IMPACT CHAMBER
gamma_ic = 1.4      # Gamma for the (near vaccum) gas in the impact chamber
V_ic = 1

# PROJECTILE & SABOT VALUES
m_pr = 17e-3     # Projectile mass
m_sb = 0.47e-3  # Sabot Mass
mu_sb = 0.4     # Sabot Friction coeff.

# SIMULATION DETAILS
delta_t = 1e-7  # Time step length

# Calculating areas from diameters
A_c = np.pi * (D_c / 2)**2
A_pt = np.pi * (D_pt / 2)**2
V0_pt = A_pt * L0_pt
A_b = np.pi * (D_b / 2)**2
A_pis = np.pi * (D_pis / 2)**2

# Some global variables because im a bad programmer
t_burnout = 0
n_burnout = 0

# | COMBUSTION |   PISTON   | RUPTURE DISK | PROJECTILE |
# |    _c      | _c/_pis/_pt|     _rup     |    _pr     |

# _c   = combustion
# _pis = piston
# _pt  = pump tube
# _b   = barrel
# _bpt = barrel - pump-tube
# _pr  = projectile
# _sb  = sabot

# t = time
# T = temperature
# P = Pressure
# v = velocity
# V = Volume
# a = acceleration
# x = displacement


def DoIt(A_c=A_c, L0_c=L0_c, P0_c=P0_c, gamma_c=gamma_c, A_pis=A_pis,
         mu_static_pis=mu_static_pis, mu_dynamic_pis=mu_dynamic_pis,
         P0_pt=P0_pt, L0_pt=L0_pt, A_pt=A_pt, P_rupt=P_rupt, L_b=L_b,
         D_b=D_b, P0_b=P0_b, gamma_lg=gamma_lg, m_pr=m_pr, m_sb=m_sb,
         mu_sb=mu_sb, gamma_ic=gamma_ic, V_ic=V_ic, delta_t=delta_t):
    # DEFINING ARRAYS
    global n_array, t_array, x_pis_array, x_pr_array, n_disk_rupture, P_c_array, \
        v_pis_array, v_pr_array, P_pt_array, Z_c_array, t_disk_rupture, \
        Zr_c_array, n_burnout, t_burnout, V_c_array, V_pt_array

    P_pt_array = [P0_pt]        # PUMP TUBE Pressure Ahead of Piston
    P_b_array = [1e-3]          # BARREL PRESSURE Behind the projectile
    x_pis_array = [0]           # PISTON DISPLACEMENT from initial position
    x_pr_array = [0]            # PROJECTILE DISPLACEMENT from initial position
    v_pis_array = [0]           # PISTON VELOCITY
    v_pr_array = [0]            # PROJECTILE VELOCITY
    a_pis_array = [0]           # PISTON ACCELERATION
    a_pr_array = [0]            # PROJECTILE ACCELERATION
    V_c_array = [V0_c]          # PISTON/COMB. VOLUME of chamber
    V_pt_array = [V0_pt]        # PISTON/BARREL VOLUME behind projectile
    n_array = []                # Counts the time step number
    t_array = []                # Holds the time step value

    # COMBUSTION CHAMBER
    P_c_array = [P0_c]  # Combustion chamber pressure array, starts with pressure for now
    Zr_c_array = [0]
    Z0_br = P0_pt * (V0_c - C / density_propel) / (C * (ForceConst_propel + P0_pt * (CoVolume_propel - 1 / density_propel)))
    Z_c_array = [Z0_br]          # POWDER BURN Decimal - may need to start at z0

    # Event BOOLEANS and Values
    global diskBroken, diskJustBroken, burnoutTF, n_burnout, t_burnout
    diskBroken = False
    diskJustBroken = True
    burnoutTF = False

    def combustElement(S):
        """ Handles the current pressure in the combustion chamber.
        Current time drives the initial expansion and increase in pressure.
        Position of piston drives the decrease in pressure due to expansion.
        """
        dz_dt = R1 * P_c_array[-1]**BR_exp                  # Burn Rate
        Z_cur = Z_c_array[-1] + delta_t * dz_dt             # Current Burnt decimal %
        Z_c_array.append(Z_cur)                             # Append to Powder Burn array
        Zr_c_array.append(dz_dt)                            # Append current Burn Rate

        if len(x_pis_array) == 1:
            dx = 0                    # If this is the first value then dx is 0
        else:
            dx = x_pis_array[-2] - x_pis_array[-1]        # Calculate the previous movement in piston to find work done

        S = S + dx * P_pt_array[-1]                         # Work done on piston (energy removed from gas)
        # Combustion Pressure using energy balance between internal, work done and PV
        P_c = (ForceConst_propel * C * Z_cur - (gamma_c - 1) * (m_pis / 2 * v_pis_array[-1]**2 + A_c * S)) \
            / (V0_c + A_c * x_pis_array[-1] - C / density_propel - (CoVolume_propel - 1 / density_propel) * C * Z_cur)

        P_c_array.append(P_c)                               # Append the Combustion Pressure
        V_c = A_pis * (L0_c + x_pis_array[-1])              # Calculate the Volume of the Combustion Chamber for this time step
        V_c_array.append(V_c)

    def setBurnoutValues():
        """ One time function that sets the values for the combustion chamber
        pressures, volumes at the time and step that Z_c exceeds 1
        """
        global n_burnout, t_burnout, Vbo_c, Pbo_c
        n_burnout = n_array[-1]
        t_burnout = t_array[-1]
        Vbo_c = V_c_array[-1]
        Pbo_c = P_c_array[-1]

    def burnoutElement():
        """ Handles the combustion chamber pressures once burnout has occured
        """
        Z_c_array.append(1)
        Zr_c_array.append(0)
        V_c_old = V_c_array[-1]
        V_c_new = A_pis * (L0_c + x_pis_array[-1])              # Calculate the Volume of the Combustion Chamber for this time step
        V_c_array.append(V_c_new)

        P_c_old = P_c_array[-1]
        P_c_new = P_c_old * V_c_old / V_c_new
        P_c_array.append(P_c_new)

    def pistonElement():
        """ Handles a single time step element of the piston along the
        pump tube.
        """
        delta_P = P_c_array[-1] - P_pt_array[-1]                                # Identify the Pressure Differential across Piston
        F_pressure = delta_P * A_pis                                              # Identify the pressure force on piston

        # Check if we need to consider static or dynamic friction -------------
        if v_pis_array[-1] == 0:                                                # Static or Dynamic friction - are we stationary
            F_fric = mu_static_pis                                              # Calculate Static Friction
            F_res = (abs(F_pressure) - abs(F_fric)) * (F_pressure / abs(F_pressure))  # Calculate the resultant force with friction opposing it
            if abs(F_pressure) > F_fric:                                        # Will pressure force overcome the friction force?
                a_pis = F_res / m_pis                                             # Calculate acceleration from resultance force
            else:
                a_pis = 0                                                       # If pressure cant overcome static friction then acceleration is 0
        else:
            F_fric = mu_dynamic_pis
            F_res = (abs(F_pressure) - abs(F_fric)) * (F_pressure / abs(F_pressure))  # Calculate the resultant force with friction opposing it
            a_pis = F_res / m_pis                                                 # Calculate the acceleration of the piston

        # Check if were at the backstop ---------------------------------------
        if allowPistonRearward == False and x_pis_array[-1] <= 0 and a_pis < 0:
            # Check if the piston is against the backstop and dont allow it to accelerate rearward
            a_pis_array.append(0)                                               # Acceleration is 0
            v_pis_array.append(0)                                               # Velocity of piston is 0
            x_pis_array.append(x_pis_array[-1])                                 # Displacement of the piston remains the same
        else:
            # Resultant force is forward
            v_pis = v_pis_array[-1] + a_pis * delta_t
            x_pis = v_pis_array[-1] * delta_t + 0.5 * a_pis * (delta_t**2) + x_pis_array[-1]
            a_pis_array.append(a_pis)
            v_pis_array.append(v_pis)
            x_pis_array.append(x_pis)

        P_pt = P_pt_array[0] * np.power((L0_pt - x_pis_array[0]) / (L0_pt - x_pis_array[-1]), gamma_lg)
        P_pt_array.append(P_pt)
        V_pt = (L0_pt - x_pis_array[-1]) * A_pis
        V_pt_array.append(V_pt)

    def barrelElement():
        """ Handles the pressure and movement of the pump-tube / barrel combo
        Does not handle the pressure of the combustion chamber pre or post
        burnout
        """
        # Volume of Chamber
        V_pt = (L0_pt - x_pis_array[-1]) * A_pis
        V_b = A_b * x_pr_array[-1]
        V_bpt = V_pt + V_b
        V_pt_array.append(V_bpt)

        # Movement of Piston
        delta_P_pis = P_c_array[-1] - P_pt_array[-1]
        F_pressure = delta_P_pis * A_pis
        F_fric = 0
        F_res = F_pressure - F_fric
        a_pis = F_res / m_pis
        a_pis_array.append(a_pis)
        v_pis = v_pis_array[-1] + a_pis * delta_t
        x_pis = v_pis_array[-1] * delta_t + 0.5 * a_pis * (delta_t**2) + x_pis_array[-1]
        v_pis_array.append(v_pis)
        x_pis_array.append(x_pis)

        # Movement of Sabot/Projectile
        delta_P_pr = P_pt_array[-1] - P_b_array[-1]
        F_pressure = delta_P_pr * A_b
        F_fric = mu_sb * 0  # change later
        F_res = F_pressure - F_fric
        a_pr = F_res / (m_pr + m_sb)
        a_pr_array.append(a_pr)
        v_pr = v_pr_array[-1] + a_pr * delta_t
        x_pr = v_pr_array[-1] * delta_t + 0.5 * a_pr * (delta_t**2) + x_pr_array[-1]
        v_pr_array.append(v_pr)
        x_pr_array.append(x_pr)

        # Calculation of (changed) pressures
        P_pt = P_pt_array[0] * np.power((A_pis * (L0_pt - x_pis_array[0])) / (A_pis * (L0_pt - x_pis_array[-1]) + A_b * x_pr_array[-1]), gamma_lg)
        P_pt_array.append(P_pt)
        P_b = P_b_array[0] * np.power((V_ic + A_b * (L_b - x_pr_array[0])) / (V_ic + A_b * (L_b - x_pr_array[-1])), gamma_ic)
        P_b_array.append(P_b)

    # Running the Discrete Model ----------------------------------------------
    i = 0
    n_array.append(i)
    t_array.append(i * delta_t)
    S = 0
    while x_pr_array[-1] < L_b:
        if P_pt_array[-1] < P_rupt and diskBroken == False:
            # Disk unbroken, just the pump tube
            if Z_c_array[-1] < 1:
                combustElement(S)
            else:
                if n_burnout == 0:
                    setBurnoutValues()
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
                if n_burnout == 0:
                    setBurnoutValues()
                burnoutElement()
            barrelElement()
        else:
            # Disk broken- pump tube and sabot dynamics
            if Z_c_array[-1] < 1:
                combustElement(S)
            else:
                if n_burnout == 0:
                    setBurnoutValues()
                burnoutElement()
            barrelElement()

        i += 1
        n_array.append(i)
        t_array.append(i * delta_t)

        if i * delta_t > 2:                 # If projectile hasnt left the barrel after 2 seconds then something is wrong
            n_disk_rupture = 0
            t_disk_rupture = 0
            break

    print('Projectile Exit Velocity = ', round(v_pr_array[-1], 3), 'mps')
    print('Exit at t = {} (n = {})'.format(t_array[-1], n_array[-1]))
    print('Peak Piston Velocity =     ', round(max(v_pis_array), 3), 'mps')
    print('Maximum P_c : ', round(max(P_c_array) / 1e6, 3), 'MPa')
    print('Maximum P_pt: ', round(max(P_pt_array) / 1e6, 3), 'MPa')
    print('Disk Ruptured after: {}s (n = {})'.format(t_disk_rupture, n_disk_rupture))
    if t_burnout == 0:
        print('Burnout did not occur, projectile exited when burn ratio was {}'.format(round(Z_c_array[-1], 4)))
    else:
        print('Burnout after: {}s (n = {})'.format(t_burnout, n_burnout))
    return n_array, t_array, x_pis_array, x_pr_array, n_disk_rupture, P_c_array


data = DoIt()  # Data format: n_array, t_array, x_pis_array, x_p_array, n_disk_rupture, P_c_array

# ------------------------------- New Plot ------------------------------------
fig_Block = plt.figure()
ax_Block = fig_Block.add_subplot()
ax_Block.set_facecolor("yellow")
text_kwargs = dict(ha='center', va='center', fontsize=50, color='k')
ax_Block.text(.5, .5, 'NEW RUN', **text_kwargs)
# -----------------------------------------------------------------------------

# ---------------------------- Burn-Time Plots --------------------------------
fig_BT = plt.figure()               # Create Figure
fig_BT.suptitle('Powder Burnt vs Time')  # Set Figure Title
ax_BT = fig_BT.add_subplot()        # Add axes to figure

ax_BT.set_xlabel('Time (s)')        # Set x label
ax_BT.set_ylabel('Powder Burnt (decimal %)')   # Set y label
# Apply a grid to plot area
# Plot Velocities with time
ax_BT.plot(t_array[0:len(Z_c_array)], Z_c_array, label='Powder Burnt')
ax_BT.grid()                        # Apply a grid to plot area
ax_BT.axvline(t_disk_rupture, color='grey', linestyle='--', label='Disk Rupture')
ax_BT.axvline(t_burnout, color='red', linestyle='--', label='Burnout')
ax_BT.legend()                      # Enable Legends
# -----------------------------------------------------------------------------

# ------------------------- BurnRate-Time Plots -------------------------------
fig_BrT = plt.figure()               # Create Figure
fig_BrT.suptitle('Powder Burn Rate vs Time')  # Set Figure Title
ax_BrT = fig_BrT.add_subplot()        # Add axes to figure

ax_BrT.set_xlabel('Time (s)')        # Set x label
ax_BrT.set_ylabel('Powder Burn Rate(decimal % per second)')   # Set y label
# Apply a grid to plot area
# Plot Velocities with time
ax_BrT.plot(t_array[0:len(Zr_c_array)], Zr_c_array, label='Powder Burn Rate')
ax_BrT.grid()                        # Apply a grid to plot area
ax_BrT.axvline(t_disk_rupture, color='grey', linestyle='--', label='Disk Rupture')
ax_BrT.axvline(t_burnout, color='red', linestyle='--', label='Burnout')
ax_BrT.legend()                      # Enable Legends
# -----------------------------------------------------------------------------
"""
# ---------- Propellant Burnt vs Powder Chamber Pressure ----------------------
fig_BP = plt.figure()               # Create Figure
fig_BP.suptitle('Powder Burnt vs Powder Chamber Pressure')  # Set Figure Title
ax_BP = fig_BP.add_subplot()        # Add axes to figure

ax_BP.set_xlabel('Powder Burnt Mass Ratio (-)')        # Set x label
ax_BP.set_ylabel('Combustion Chamber Pressure (Pa)')   # Set y label
# Apply a grid to plot area
# Plot Velocities with time
ax_BP.plot(Z_c_array, P_c_array, label='C = {}g'.format(C * 1e3))
ax_BP.grid()                        # Apply a grid to plot area
ax_BP.legend()                      # Enable Legends
ax_BP.axvline(Z_c_array[n_disk_rupture], color='grey', linestyle='--', label='Disk Rupture')
# -----------------------------------------------------------------------------
"""
# ------------------------ Pressure-Time Plots --------------------------------
fig_PT = plt.figure()               # Create Figure
fig_PT.suptitle('Pressure vs Time')  # Set Figure Title
ax_PT = fig_PT.add_subplot()        # Add axes to figure
ax_PT.set_yscale('log')             # Use logarithmic Y scale

ax_PT.set_xlabel('Time (s)')        # Set x label
ax_PT.set_ylabel('Pressure (Pa)')   # Set y label

# Plot Pressures with time
ax_PT.plot(t_array, P_c_array, label='Combustion')
ax_PT.plot(t_array, P_pt_array, label='Pump Tube')
ax_PT.grid()                        # Apply a grid to plot area
ax_PT.axhline(y=P_rupt, color='k')
ax_PT.axvline(t_disk_rupture, color='grey', linestyle='--', label='Disk Rupture')
ax_PT.axvline(t_burnout, color='red', linestyle='--', label='Burnout')
ax_PT.legend()                      # Enable Legends
# -----------------------------------------------------------------------------

# ------------------------ Volume-Time Plots --------------------------------
fig_VT = plt.figure()               # Create Figure
fig_VT.suptitle('Volume vs Time')  # Set Figure Title
ax_VT = fig_VT.add_subplot()        # Add axes to figure

ax_VT.set_xlabel('Time (s)')        # Set x label
ax_VT.set_ylabel('Volume (m3)')   # Set y label

# Plot Volumes with time
ax_VT.plot(t_array, V_c_array, label='Combustion')
ax_VT.plot(t_array, V_pt_array, label='Pump Tube')
ax_VT.grid()                        # Apply a grid to plot area
ax_VT.axvline(t_disk_rupture, color='grey', linestyle='--', label='Disk Rupture')
ax_VT.axvline(t_burnout, color='red', linestyle='--', label='Burnout')
ax_VT.legend()                      # Enable Legends
# -----------------------------------------------------------------------------

# ------------------------ Displacement-Time Plots ----------------------------
fig_DT = plt.figure()               # Create Figure
fig_DT.suptitle('Displacement vs Time')     # Set Figure Title
ax_DT = fig_DT.add_subplot()        # Add axes to figure

ax_DT.set_xlabel('Time (s)')        # Set x label
ax_DT.set_ylabel('Displacement (m)')   # Set y label

# Plot Displacements with time
ax_DT.plot(t_array, x_pis_array, label='Piston')
ax_DT.plot(t_array[0:len(x_pr_array)], x_pr_array, label='Projectile')
ax_DT.grid()                        # Apply a grid to plot area
ax_DT.axvline(t_disk_rupture, color='grey', linestyle='--', label='Disk Rupture')
ax_DT.axvline(t_burnout, color='red', linestyle='--', label='Burnout')
ax_DT.legend()                      # Enable Legends
# -----------------------------------------------------------------------------
"""
# ------------------------ Velocity-Time Plots --------------------------------
fig_vT = plt.figure()               # Create Figure
fig_vT.suptitle('Velocity vs Time')  # Set Figure Title
ax_vT = fig_vT.add_subplot()        # Add axes to figure

ax_vT.set_xlabel('Time (s)')        # Set x label
ax_vT.set_ylabel('Velocity (m/s)')   # Set y label
# Apply a grid to plot area
# Plot Velocities with time
ax_vT.plot(t_array, v_pr_array, label='Projectile')
ax_vT.plot(t_array, v_pis_array, label='Piston')
ax_vT.grid()                        # Apply a grid to plot area
ax_vT.axvline(t_disk_rupture, color='grey', linestyle='--', label='Disk Rupture')
ax_vT.axvline(t_burnout, color='red', linestyle='--', label='Burnout')
ax_vT.legend()                      # Enable Legends
# -----------------------------------------------------------------------------
"""
