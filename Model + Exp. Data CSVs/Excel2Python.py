# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 12:39:37 2022

@author: alfie
"""

import pandas as pd
from matplotlib import pyplot as plt

plt.rcParams['figure.figsize'] = [6, 4]         # Sets the figure size which affects the size of text, try to keep the same where possible
plt.rcParams['figure.dpi'] = 250                # Setst the figure DPI which makes the image sharp but not a crazy large file size
# ----------------------------------------------------------------------------
locationName = 'Riad 12.7 Comparison Data.xlsx'          # Input the name and location if not right next to where you have saved the spreadsheet
sheetName = 'Sheet1'                            # Sheet name, set to Sheet1 if only one unamed sheet
columnNames = ['Charge Mass', 'E Pc', 'E Ppt', 'E Pr Vel', 'M Pc', 'M Ppt', 'M Pr Vel', 'M1 Vel', 'M2 Vel']  # List of the Column Names in the order they are in the sheet
xAxisName = 'Charge Mass (g)'              # Label for X axis
yAxisName = 'Projectile Velocity (mps)'     # Label for Y axis
# ----------------------------------------------------------------------------

sheet = pd.read_excel(locationName, sheetName)  # Reads in the full sheet

data = sheet[columnNames]                       # only takes the data from the columns that we named above


X = data[columnNames[0]].tolist()               # Takes the X and Y data from the first and second named columns
Ec = data[columnNames[1]].tolist()
Ept = data[columnNames[2]].tolist()
Ev = data[columnNames[3]].tolist()
Mc = data[columnNames[4]].tolist()
Mpt = data[columnNames[5]].tolist()
Mv = data[columnNames[6]].tolist()
M1v = data[columnNames[7]].tolist()
M2v = data[columnNames[8]].tolist()

fig = plt.figure()                              # Creates the figure object
ax = fig.add_subplot()                          # Creates the axis object
#ax2 = ax.twinx()

ax.scatter(X, Mv, marker='x', color='#1f77b4', label='Our Model', alpha=0.75)                    # Plots a scatter with markers that are xs
# ax.scatter(X, Mpt, marker='1', color='#1f77b4', label='Model Pump.', alpha=0.75)
ax.plot(X, Mv, color='#1f77b4', alpha=0.75)

ax.scatter(X, M1v, marker='x', color='darkcyan', label='Riad M1', alpha=0.75)                    # Plots a scatter with markers that are xs
# ax.scatter(X, Mpt, marker='1', color='#1f77b4', label='Model Pump.', alpha=0.75)
ax.plot(X, M1v, color='darkcyan', alpha=0.75)

ax.scatter(X, M2v, marker='x', color='mediumturquoise', label='Riad M2', alpha=0.75)                    # Plots a scatter with markers that are xs
# ax.scatter(X, Mpt, marker='1', color='#1f77b4', label='Model Pump.', alpha=0.75)
ax.plot(X, M2v, color='mediumturquoise', alpha=0.75)
# ax.plot(X, Mpt, color='#1f77b4', alpha=0.75)

ax.scatter(X, Ev, marker='x', color='orange', label='Experim. ', alpha=0.95)                    # Plots a scatter with markers that are xs
# ax.scatter(X, Ept, marker='1', color='orange', label='Exp. Pump.', alpha=0.75)                    # Plots a scatter with markers that are xs
# ax.set_ylim(50, 275)

#ax2.scatter(X, Mv, marker='o', color='green', label='Model Vel.', alpha=0.75)
#ax2.scatter(X, Ev, marker='o', color='red', label='Exp. Vel.', alpha=0.75)
#ax2.plot(X, Mv, color='green', alpha=0.75)
#ax2.set_ylim(800, 1200)


ax.set_xlabel(xAxisName)                        # Sets the x axis name to that prescribed above
ax.set_ylabel(yAxisName)                        # Sets the y axis name to that prescribed above
#ax2.set_ylabel('Projectile Velocity (mps)')
# ax.set_xscale('log')                            # Sets the x axis to logarithmic
ax.grid()                                       # Sets the axis to display a grid corresponding to the major ticks

ax.legend()
# ax2.legend()
# Saves the figure with a long descriptive name
#fig.savefig(locationName + '_' + sheetName + '(' + columnNames[0] + '_' + columnNames[1] + ').png')
fig.savefig('Riad Velocity Comparison' + '.png')
