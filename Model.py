#!/usr/bin/env python
# coding: utf-8

# Imports
import numpy as np
import matplotlib.pyplot as plt

from itertools import product

""" Function Declarations """

# Returns the Michaelis-Menten Velocity for given
# Substrate Concentration, Km, and Vmax
def mich_menten(S, k_m, v_max):
    return v_max * S / (k_m + S)

# Returns a dictionary containing kinetic data for
# all given parameter combinations
def run(params, data={}):
    for s, km, vmax in params:
        if (km, vmax) not in data:
            data[(km,vmax)] = list()
    
        data[(km, vmax)].append(mich_menten(s,km, vmax))
    
    return data

# Returns which enzyme dominated from the given
# kinetic data
# Probably, needs a minor rework as it doesn't fit its use currently.
def dominant(dataset1, dataset2):
    result = 'Indeterminant'
    
    if max(dataset1) > max(dataset2):
        result = '1'
    
    if max(dataset1) < max(dataset2):
        result = '2'
    
    return result

""" Data Simulation """
# Define parameter space
S = np.linspace(0, 5, 100)
Km = np.linspace(0.1, 1, 10)
Vmax = np.linspace(0.9, 1.1, 2)

# Get all combinations of 
# the parameter space
params = product(S, Km, Vmax)

# Generate the data
data = run(params)


""" Plot 'Em """

# Get Km and Vmax combos from the data
combos = [(km, vmax) for (km, vmax) in data]
combos = product(combos, combos)

# to keep from repeating work record km/vmax
done = list()

# to keep a count of plots produced for labeling
count = 0

# plot each combo of km/vmax
for c1, c2 in combos:
    combo = (c1, c2)
    
    # ignore enzyme with the same kinetics
    if c1 == c2:
        continue
    elif combo in done:
        # don't repeat work
        continue
    else:
        # don't repeat work
        done.append(combo)
        done.append(reversed(combo))
        
    count += 1
        
    plt.plot(S, data[c1], label=f'Enzyme 1: Km={c1[0]} Vmax={c1[1]}')
    plt.plot(S, data[c2], label=f'Enzyme 2: Km={c2[0]} Vmax={c2[1]}')
    
    plt.xlabel('Substrate Concentration')
    plt.ylabel('Reaction Velocity')
    plt.title(f'{count}: Enzyme {dominant(data[c1], data[c2])} dominates')
    plt.legend()
    
    plt.show()   
