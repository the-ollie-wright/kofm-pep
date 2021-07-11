import math
import numpy as np

import argparse
import os
from datetime import datetime

import csv

import leaders
import constructions
import functions as funs
import csvget
import copy

# Set up Parser

parser = argparse.ArgumentParser(description='Run DSMC on a Boltzmann-type Equation')

# Discretization Constants

parser.add_argument("-T", "--timesteps", type=int, default=100000, help='The number of time steps.')
parser.add_argument("-N", "--spacesteps", type=int, default=50, help='The number of space steps.')
parser.add_argument("-M", "--indivno", type=int, default=4000, help='The number of individuals in the simulations.')

# Powers for functions

parser.add_argument("-a", "--a", type=float, default=1, help='Value of a in the P function')
parser.add_argument("-b", "--b", type=float, default=2, help='Value of b in the P function')
parser.add_argument("-alph", "--alpha", type=float, default=2, help="Value of alpha in the D^2 function")
parser.add_argument("-sig", "--sigma", type=float, default=0.01, help="The variance of the eta random variable")

# Equation constants

parser.add_argument("-g", "--gamma", type=float, default=0.3, help='The diffusion constant')
parser.add_argument("-r", "--tol", type=float, default=0.1, help='Radius of Interaction')
parser.add_argument("-L", "--tauL", type=float, default=1, help='The relaxation time for the Leader species')
parser.add_argument("-i", "--taui", type=float, default=1, help='The relaxation time for the other species species')

# Region Indicator

parser.add_argument("-s", "--kkey", type=str, default='E', help='The Region Indicator')
parser.add_argument("-l", "--leadpercent", type=float, default=0.05, help='The percentage of the total individuals the leader species individuals have')


args = parser.parse_args()

# Discretization Constants

T = args.timesteps
N = args.spacesteps
M = args.indivno

# Mesh, Space and timestep caluculations and constants

start = -1.0
end = 1.0

# Powers for functions

a = args.a
b = args.b
alpha = args.alpha
sigma = args.sigma

# Equation Constants

r = args.tol
gamma = args.gamma
tauL = args.tauL
taui = args.taui

# Region Indicator and types of files

true_key = args.kkey
leaderpercent = args.leadpercent

comma = '/csvfiles'

print("%s - Running with %d time steps with %d individuals, r=%f, gamma=%f, sigma=%f" % (datetime.now(), T, M*(1+leaderpercent), r, gamma, sigma))

# Demographic Constants 
brex_no = 2  
age_no = 4  
total_no = brex_no * age_no

# Spaces and stuff

xdom, deltax = funs.setup_space_domain(start, end, N)

# Parties and boundary vectors

parties = [-0.7, -0.4, -0.1, 0.1, 0.4, 0.7]  

boundary_vec = [0.5 * (parties[i - 1] + parties[i]) for i in range(1, len(parties))]
boundary_vec.extend([-1.0, 1.0])
boundary_vec.sort()

# Characteristic Matrix 

characteristic_mat = csvget.read_characteristic(true_key, len(parties))

keylist = []

for k in range(brex_no):
    for l in range(brex_no, brex_no + age_no):
        keylist.append((k,l)) 

# Setup Leader sp[ecies

print("%s - Setting up the leader species" % (datetime.now()))

leader_vec, key_L = leaders.setup_leader_species(M, leaderpercent, true_key, parties)

# Setup follower species

boundary_vec2 = [0.5 * (parties[i - 1] + parties[i]) for i in range(1, len(parties))]
boundary_vec2.extend([-0.95, 0.95])
boundary_vec2.sort()


print("%s - Setting up the follower species" % (datetime.now()))

opin_dict, M_ = constructions.setup_follower_species(M, keylist, boundary_vec2, true_key)

opin_dict[key_L] = leader_vec
M_[key_L] = len(leader_vec)
keylist.append(key_L)
M_total = sum(M_.values())

# Creating an ordered list of probabilities to choose a list from opin_dict such that the overall probability of choosing an individual is uniform.

P_chooselist = funs.create_prob_list(M, leaderpercent, M_, tauL, taui)

# Setting up the objects used to save the last 5% of time steps

smooth_time = T * 0.05

smooth_vec = []
tobeadded = {}

for key in keylist:
    tobeadded[key] = [0 for i in range(len(xdom)-1)]

# Saving the initial conditions of the method

current_graph = {}
for key in keylist:
    current_graph[key] = [0 for i in range(len(xdom)-1)]
    for opin in opin_dict[key]:
        current_graph[key][funs.binsearch(opin, xdom)] += 1.0

for key in keylist:

    temp_folder = 'dataandsolutions/dataandsolutions' + true_key + comma + str(key)
    tf = temp_folder + '/timeis0.csv'


    try:
        os.mkdir(temp_folder)
    except OSError:
        pass

    print("%s - Creating file %s " % (datetime.now(), tf))

    with open(tf, 'w') as resultfile:
        wr = csv.writer(resultfile)
        
        wr.writerow(current_graph[key])
        wr.writerow(xdom)

# Starting the time stepping

for t in range(T):
    if t % (T/10.0) == 0 and t<T-smooth_time:
        print("%s - time is %d" % (datetime.now(), t))

    # Running interactions for this time step

    for interactions in range(int(M*(1+leaderpercent))):
        # the choose_indivs functions pops indiv0 and indiv1 thus removing them before adding back in the post interaction individuals
        indiv0, indiv1 = funs.choose_indivs(opin_dict, P_chooselist, keylist)

        new_indiv0, new_indiv1 = funs.interaction(indiv0, indiv1, a, b, alpha, sigma, r, boundary_vec, characteristic_mat, gamma)

        key0, opin0 = new_indiv0
        key1, opin1 = new_indiv1

        opin_dict[key0].append(opin0)
        opin_dict[key1].append(opin1)


    # Collecting the last 5% of time steps into a dictionary 
    if t > T - smooth_time:
        for key in keylist:
            for opin in opin_dict[key]:
                tobeadded[key][funs.binsearch(opin, xdom)] += 1.0
            
    
for key in keylist:
    # normalising the data to scale it to the density
    tobeadded[key] = [a*deltax/M_[key] for a in tobeadded[key]]

    temp_folder = 'dataandsolutions/dataandsolutions' + true_key + comma + str(key)
    tf = temp_folder + '/finaltime4.csv'


    try:
        os.mkdir(temp_folder)
    except OSError:
        pass

    print("%s - Creating file %s " % (datetime.now(), tf))

    with open(tf, 'w') as resultfile:
        wr = csv.writer(resultfile)
        
        wr.writerow(tobeadded[key])
        wr.writerow(xdom)














