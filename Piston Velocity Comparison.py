# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 16:04:16 2022

@author: alfie
"""

import csv
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import linregress
plt.rcParams['figure.figsize'] = [6, 4]
plt.rcParams['figure.dpi'] = 250

M1Path = 'PR17.0g-BORE12.7mm-C8.0g.csv'
E1Path = 'Model + Exp. Data CSVs/Experimental 9g17g12.7mm PISTON VELOCITY.csv'


def ReadIn(path, n=1, m=1e6):
    """ Reads in data from a 2 column csv file and returns the columns values
        as two lists """
    x_data = []
    y_data = []
    with open(path, "r", encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        for i, line in enumerate(reader):
            x = float(line[0])
            y = float(line[n]) / m
            x_data.append(x)
            y_data.append(y)
    file.close()
    return(x_data, y_data)


def Compare(modelPath, expPath, dx_e=0, dy_e=0):
    """ Compares the model and experimental data. Can shift the experimental
        in both x and y as specified"""

    mod_x, mod_y = ReadIn(modelPath, n=3, m=1)
    exp_x_unaltered, exp_y_unaltered = ReadIn(expPath, m=1)

    exp_x = [x + dx_e for x in exp_x_unaltered]
    exp_y = [y + dy_e for y in exp_y_unaltered]

    xs = [exp_x_unaltered, exp_x, mod_x, ]
    ys = [exp_y_unaltered, exp_y, mod_y]
    labels = ['Unaltered Experimental', 'Altered Experimental', 'Model']
    linestyles = ['-', '-', '-']
    colors = ['lightgrey', 'orange', 'blue']
    title = 'Piston Velocity vs Time: Model vs Experimental'
    xtitle = 'Time (s)'
    ytitle = 'Velocity (mps)'

    Plot(xs, ys, labels, linestyles, colors, title, xtitle, ytitle)

    exp_x_trimmed = []
    exp_y_trimmed = []
    for x, y in zip(exp_x, exp_y):
        if x > min(mod_x) and x < max(mod_x):
            exp_x_trimmed.append(x)
            exp_y_trimmed.append(y)

    mod_y_aligned = np.interp(exp_x_trimmed, mod_x, mod_y)
    mod_x_aligned = exp_x_trimmed

    xs = [exp_x_trimmed, mod_x_aligned]
    ys = [exp_y_trimmed, mod_y_aligned]
    labels = ['Experimental', 'Model']
    linestyles = ['-', '-']
    colors = ['orange', 'blue']
    title = 'Piston Velocity vs Time: Model vs Experimental'
    xtitle = 'Time (s)'
    ytitle = 'Piston Velocity (mps)'

    Plot(xs, ys, labels, linestyles, colors, title, xtitle, ytitle)

    exp_x_unaltered_trimmed = []
    exp_y_unaltered_trimmed = []
    for x, y in zip(exp_x_unaltered, exp_y_unaltered):
        if x > min(mod_x) and x < max(mod_x):
            exp_x_unaltered_trimmed.append(x)
            exp_y_unaltered_trimmed.append(y)

    mod_y_aligned2 = np.interp(exp_x_unaltered_trimmed, mod_x, mod_y)

    slope, intercept, r_value_unalt, p_value, std_err = linregress(exp_y_unaltered_trimmed, mod_y_aligned2)
    R2_unalt = r_value_unalt**2
    print('Unadusted Datasets have an R squared value of: ', R2_unalt)

    xs = [exp_y_unaltered_trimmed, [min(mod_y_aligned), max(mod_y_aligned)]]
    ys = [mod_y_aligned2, [min(mod_y_aligned), max(mod_y_aligned)]]
    labels = ['Unaltered Model vs Exp. (R2={})'.format(round(R2_unalt, 3)), 'x=y (R=1)']
    linestyles = ['-', '--']
    colors = ['blue', 'k']
    title = 'Piston Velocity vs Time: Model vs Experimental'
    xtitle = 'Experimental Piston Velocity (mps)'
    ytitle = 'Model Piston Velocity (mps)'

    Plot(xs, ys, labels, linestyles, colors, title, xtitle, ytitle, 2, True)

    meansquared_sum = 0
    meanvalue_sum = 0
    for exp, mod in zip(exp_y_unaltered_trimmed, mod_y_aligned2):
        squared_diff = (exp - mod)**2
        meansquared_sum += squared_diff
        meanvalue_sum += exp + mod
    meansquared = meansquared_sum / len(exp_y_unaltered_trimmed)
    meanvalue = meanvalue_sum / (2 * len(exp_y_unaltered_trimmed))
    print('Mean Squared Difference = ', meansquared)
    print('Mean value = ', meanvalue)


def Plot(xs, ys, labels, linestyles, colors, title='', xtitle='', ytitle='', n_col=4, scatter=False):
    fig, ax = plt.subplots()
    fig.suptitle(title)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

    for n in range(0, len(xs)):
        if n != len(xs) - 1 and scatter == True:
            ax.scatter(xs[n], ys[n], label=labels[n], linestyle=linestyles[n], color=colors[n], s=4)
        else:
            ax.plot(xs[n], ys[n], label=labels[n], linestyle=linestyles[n], color=colors[n])
    ax.grid()
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True,
              shadow=True, ncol=n_col)


def DeltaTPeak(modelPath, expPath):
    mod_x, mod_y = ReadIn(modelPath, n=3, m=1)
    exp_x, exp_y = ReadIn(expPath, m=1)

    mod_max = max(mod_y)
    mod_max_index = mod_y.index(mod_max)
    exp_max = max(exp_y)
    exp_max_index = exp_y.index(exp_max)
    timeshift = exp_x[exp_max_index] - mod_x[mod_max_index]
    return(timeshift)


dt = 0  # DeltaTPeak(M1Path, E1Path)
Compare(M1Path, E1Path, dx_e=-dt)


"""
plt.rcParams['figure.figsize'] = [6, 4]
plt.rcParams['figure.dpi'] = 250

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

xValues = np.linspace(0, max(x1_data), 500)

yValues = np.interp(xValues, x1_data, y1_data)

fig, ax = plt.subplots()
ax.plot(x1_data, y1_data, color='Orange', label='Input Data')
ax.plot(xValues, yValues, label='Interpolated', linestyle='--', color='red')
ax.grid(True, color='dimgray', linestyle='--', linewidth=0.5)
ax.set_axisbelow(True)
ax.set_ylabel('y')
ax.set_xlabel('x')
ax.legend()
"""
