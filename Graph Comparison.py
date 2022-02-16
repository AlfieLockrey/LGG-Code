# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 15:40:49 2022

@author: alfie
"""

import csv
import numpy as np
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.stats import linregress
plt.rcParams['figure.figsize'] = [6, 4]
plt.rcParams['figure.dpi'] = 250


def LazyPoly(PolyCoeffs, xValues):
    power = len(PolyCoeffs) - 1
    yValues = []
    for x in xValues:
        y_sum = 0
        for pwr in range(power + 1):
            coeff = PolyCoeffs[pwr]
            y_i = coeff * x**(power - pwr)
            y_sum += y_i
        yValues.append(y_sum)
    return(yValues)


x1_data = []
y1_data = []
with open("Experimental 8g17g12.7mm.csv", "r", encoding='utf-8-sig') as file_exp:
    reader = csv.reader(file_exp)
    for i, line in enumerate(reader):
        x1 = float(line[0])
        y1 = float(line[1])
        x1_data.append(x1)
        y1_data.append(y1)
file_exp.close()

x2_data = []
y2_data = []
with open("PR17.0g-BORE12.7mm-C8.0g.csv", "r", encoding='utf-8-sig') as file_model:
    reader = csv.reader(file_model)
    for i, line in enumerate(reader):
        x2 = float(line[0])
        y2 = float(line[1])
        x2_data.append(x2)
        y2_data.append(y2)
file_model.close()

# ------------------------ Pressure-Time Plots --------------------------------
fig_PT = plt.figure()               # Create Figure
fig_PT.suptitle('Combustion Chamber Pressure vs Time: Model vs. Experimental')  # Set Figure Title
ax_PT = fig_PT.add_subplot()        # Add axes to figure
# ax_PT.set_yscale('log')             # Use logarithmic Y scale

ax_PT.set_xlabel('Time (s)')        # Set x label
ax_PT.set_ylabel('Pressure (MPa)')   # Set y label

# Plot Pressures with time
P1_c_MPa = [P / 1e6 for P in y1_data]
P2_c_MPa = [P / 1e6 for P in y2_data]
ax_PT.plot(x2_data, P2_c_MPa, label='Model')
ax_PT.plot(x1_data, P1_c_MPa, label='Experimental')

ax_PT.grid()                        # Apply a grid to plot area
ax_PT.legend()   # Enable Legends
# -----------------------------------------------------------------------------

max1_index = P1_c_MPa.index(max(P1_c_MPa))
max2_index = P2_c_MPa.index(max(P2_c_MPa))
timeshift = x1_data[max1_index] - x2_data[max2_index]

newexpXdata_long = [i - timeshift for i in x1_data]


polypower = 15
expPoly, expres, exprank, expsing, exprcond = np.polyfit(newexpXdata_long, P1_c_MPa, polypower, full=True)
modelPoly, modres, modrank, modsing, modcond = np.polyfit(x2_data, P2_c_MPa, deg=polypower, full=True)

times = np.linspace(0, min(max(x1_data), max(x2_data)), 500)

yp1 = LazyPoly(expPoly, times)
yp2 = LazyPoly(modelPoly, times)

# ------------------------ Pressure-Time Plots --------------------------------
fig_poly = plt.figure()               # Create Figure
fig_poly.suptitle('Combustion Chamber P vs T: Polynomial Comparison')  # Set Figure Title
ax_poly = fig_poly.add_subplot()        # Add axes to figure
# ax_PT.set_yscale('log')             # Use logarithmic Y scale

ax_poly.set_xlabel('Time (s)')        # Set x label
ax_poly.set_ylabel('Pressure (MPa)')   # Set y label

# Plot Pressures with time

ax_poly.plot(newexpXdata_long, P1_c_MPa, label='Experimental', color='orange')
ax_poly.plot(x2_data, P2_c_MPa, label='Model', color='blue')
ax_poly.plot(times, yp1, label='Experimental Poly', color='#ff7f0e', linestyle='--')
ax_poly.plot(times, yp2, label='Model Poly', color='#1f77b4', linestyle='--')

ax_poly.grid()                        # Apply a grid to plot area
ax_poly.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
               fancybox=True, shadow=True, ncol=4)    # Enable Legends
# -----------------------------------------------------------------------------

fig_reg = plt.figure()               # Create Figure
fig_reg.suptitle('Linear Regression of LGG Model and Experimental Data')
ax_reg = fig_reg.add_subplot()        # Add axes to figure
ax_reg.plot(yp1, yp2, label='Polynomial Data')
ax_reg.plot([0, max(yp2)], [0, max(yp2)], color='k', linestyle='--', label='x=y (R=1)')
ax_reg.grid()                        # Apply a grid to plot area
ax_reg.set_ylabel('Pressure Model (MPa)')
ax_reg.set_xlabel('Pressure Experimental (MPa)')
ax_reg.legend()

slope, intercept, r_value, p_value, std_err = linregress(yp1, yp2)
print(r_value, r_value**2)
