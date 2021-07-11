from fenics import *
import math
import numpy as np

import argparse
import os
from datetime import datetime

import csv
import matplotlib.pyplot as plt
import matplotlib as mpl

import functions as funs
import csvget
import expressions
import leaderstuff
import constructions 



"""

 This code can be used to solve the Fokker-Planck equation derived in $PAPER$. As presented here, we have used it to simulate the 2019 UK General Election using age and self-reported Vote in
the 2016 European Referendum for the UK. In principle, using methods to derive data as described in the aforementioned paper, this code could be used to simulate other such elections with any
number of demographic groups. A small amount of code would need to be edited which has been marked in the most part.

"""

# Set up Parser

parser = argparse.ArgumentParser(description='Run FEA on a Fokker-Plank Equation')

# Discretization Constants

parser.add_argument("-T", "--totaltime", type=float, default=1, help='The total time we want to run the model for')
parser.add_argument("-k", "--timechange", type=float, default=0.01, help='The time change between timesteps')
parser.add_argument("-N", "--spacesteps", type=int, default=100, help='The number of spacesteps')

# Powers for functions

parser.add_argument("-a", "--a", type=float, default=1, help='Value of a in the P function')
parser.add_argument("-b", "--b", type=float, default=2, help='Value of b in the P function')
parser.add_argument("-alph", "--alpha", type=float, default=2, help="Value of alpha in the D^2 function")

# Equation constants

parser.add_argument("-l", "--lamda", type=float, default=1, help='The diffusion constant')
parser.add_argument("-r", "--tol", type=float, default=0.1, help='Radius of Interaction')
parser.add_argument("-L", "--tauL", type=float, default=1, help='The relaxation time for the Leader species')
parser.add_argument("-i", "--taui", type=float, default=1, help='The relaxation time for the other species species')

# Region Indicator

parser.add_argument("-s", "--kkey", type=str, default='N', help='The Region Indicator')

args = parser.parse_args()

# boundary function

def boundary(x, on_boundary):
    return on_boundary

# Discretization Constants

T = args.totaltime
dt = args.timechange
N = args.spacesteps

# Mesh, Space and timestep caluculations and constants

start = -1.0
end = 1.0
functype = 'P'
n = 4

timesteps = int(T/dt)

# Powers for functions

a = args.a
b = args.b
alpha = args.alpha

# Equation Constants

r = args.tol
l = args.lamda
tauL = args.tauL
taui = args.taui
scale = 0.08

# Region Indicator and types of files

true_key = args.kkey
pic = '/pngfiles'
comma = '/csvfiles'
path0 = 'dataandsolutions/dataandsolutions'
# path0 = 'dataandsolutions/dataandsolutionsNot/dataandsolutions'


'/data/brexitdata.csv'
# Spaces and stuff

try:
    (mesh, V, W) = funs.setup_mesh_and_space(start, end, N, n, functype, total_no)
except ValueError:
    (mesh, V, V_cg, W) = funs.setup_mesh_and_space(start, end, spacesteps, n, functype, total_no)

xdom = mesh.coordinates().transpose()[0] 
deltax = xdom[1] - xdom[0] 

# Boundary Conditions

# f_D = Constant([0.0 for i in range(total_no)])

# bc = DirichletBC(W, f_D, boundary)

# Composite Constants translated into Fenics (no minus signs included) and some basic expressions

print("%s - Running with %d time steps with lambda=%f, r=%f, tauL=%f, taui=%f" % (datetime.now(), timesteps, l, r, tauL, taui))

first_ord_const_i = Constant(dt/(2.0*taui))
first_ord_const_L = Constant(dt/(2.0*tauL))

second_ord_const_i = Constant(dt*l/(4.0*taui))
second_ord_const_L = Constant(dt*l/(4.0*tauL))

D_2 = interpolate(Expression('pow(1-x[0]*x[0],2*alph)', alph=alpha, degree=8), V)
Good = interpolate(Expression('x[0]', degree=1), V)

# Parties and boundary vectors

parties = [-0.7, -0.4, -0.1, 0.1, 0.4, 0.7]  

boundary_vec = [0.5 * (parties[i - 1] + parties[i]) for i in range(1, len(parties))]
boundary_vec.extend([-1.0, 1.0])
boundary_vec.sort()

# Demographic Constants and Characteristic Matrix 

# to run a simulation for a different set of demographics the code here needs to be updated

brex_no = 2  
age_no = 4  
total_no = brex_no * age_no

for k in range(brex_no):
    for l in range(brex_no, brex_no + age_no):
        key = (k, l)    
        key_list.append(key)

characteristic_mat = csvget.read_characteristic(true_key, len(parties), path0)



# Stuff for Leader Species

print("%s - Setting up the Leader Species" % (datetime.now()))

K_L, mass_lead, f_L = leaderstuff.setup_leader_species(V, parties, scale, r, xdom, deltax, true_key)


png_folder = path0 + str(true_key) + pic + 'L'


try:
    os.mkdir(png_folder)
except OSError:
    pass

png_filename = png_folder + '/LeaderSpecies.eps'

function_name = r'$f_L$'

plot ( f_L, color='xkcd:blue') 

plt.xlabel(r'$w$')
y = plt.ylabel(function_name)
y.set_rotation(0)
plt.yticks(rotation=45, ha='right') 
plt.savefig ( png_filename , format='eps')

print("%s - Saving picture as %s" % (datetime.now(), png_filename))
plt.close ( )

# Stuff for non-Leader Species

print("%s - Setting up the follower species and the (Bi)Linear operator(s)" % (datetime.now()))

[functions, convol_operators, masses, data] = constructions.setup_follower_species(V, W, parties, brex_no, age_no, true_key, r, deltax, xdom, boundary_vec, characteristic_mat, a, b, path0, key_list)

f, f_, f_n, v = functions 

f_tup = split(f)
f__tup = split(f_)
f_n_tup = split(f_n)
v_tup = split(v)

Kappa_jij, Kappa_iji = convol_operators
mass, mass_float, mass_sum = masses

f_n_savable = f_n.split()

# eps figures are saved for all species using MatPlotLib at time t=0

for c0 in range(int(len(key_list)/2.0)):
    c1 = c0+4
    key0 = key_list[c0]
    mass0 = mass_float[c0]
    key1 = key_list[c1]
    mass1 = mass_float[c1]

    png_folder = path0 + str(true_key) + pic + str(key0[1])

    try:
        os.mkdir(png_folder)
    except OSError:
        pass
    
    png_filename = png_folder + '/timeis0.eps'

    # Creating the correct scaling for plotting superimposed graphs on the same figure

    biggest = max(mass0, mass1)

    if biggest == mass0:
        scale_fac = mass0/mass1
        proj_f0 = project(f_n_savable[c0], V=V)
        proj_f1 = project(scale_fac*f_n_savable[c1], V=V)
    else:
        scale_fac = mass1/mass0
        proj_f0 = project(scale_fac*f_n_savable[c0], V=V)
        proj_f1 = project(f_n_savable[c1], V=V)

    proj_f0_vec = proj_f0.compute_vertex_values(mesh)
    proj_f1_vec = proj_f1.compute_vertex_values(mesh)
    largest_value = np.max(proj_f0_vec) < np.max(proj_f1_vec)

    if largest_value:
        plot ( proj_f1, color='xkcd:blue', label=str(key0) )  
        plot ( proj_f0, color='xkcd:red', linestyle='--', label=str(key1) )
    else:
        plot ( proj_f0, color='xkcd:red', linestyle='--', label=str(key0) )    
        plot ( proj_f1, color='xkcd:blue', label=str(key1) )


    function_name = r'$f_{{{time}}}$'.format(time=0)

    plt.xlabel(r'$w$')
    y = plt.ylabel(function_name)
    y.set_rotation(0)
    plt.yticks(rotation=65, ha='right')
    plt.savefig ( png_filename , format='eps' )

    print("%s - Saving picture in %s" % (datetime.now(), png_filename))
    plt.close ( )


#     demo('default', 'f_0', png_filename, xdom, up, down)


for n, key in enumerate(key_list):
    png_folder = path0 + str(true_key) + pic + str(key)
    # png_folder = 'dataandsolutions/dataandsolutionsNot/dataandsolutions' + str(true_key) + pic + str(key) 
    try:
        os.mkdir(png_folder)
    except OSError:
        pass

    png_filename = png_folder + '/timeis0.eps'
    function_name = r'$f_{{{time}}}$'.format(time=0)

    up = project(f_n_savable[n], V=V)
    # plot ( up, title = function_name )
    plot ( up)

    plt.xlabel(r'$w$')
    y = plt.ylabel(function_name)
    y.set_rotation(0)
    plt.yticks(rotation=65, ha='right')

    plt.savefig ( png_filename , format='eps')
    print("%s - Saving picture in %s" % (datetime.now(), png_filename))
    plt.close ( )



# Linear and Bilinear Operator Set up

termi_Expression = {}
termj_Expression = {}
sec_termj_Expression= {}
a = []
L = []

for i in range(len(key_list)):
    termi_Expression = sum([Kappa_jij[i][j] for j in range(len(key_list))])
    termj_Expression = sum([Kappa_iji[i][j] * f__tup[j] for j in range(len(key_list))])

    sec_termj_Expression = sum(split(f_))

    a.append(  f__tup[i] * v_tup[i] * dx - (first_ord_const_L * K_L + first_ord_const_i * termi_Expression) * f__tup[i] * dot(grad(v_tup[i]), grad(Good)) * dx \
        - first_ord_const_i * termj_Expression * dot(grad(v_tup[i]), grad(Good)) * dx \
        + (second_ord_const_L * mass_lead + second_ord_const_i * mass_sum) * dot(grad(v_tup[i]), grad(D_2 * f__tup[i])) * dx \
        + second_ord_const_i * mass[i] * dot(grad(v_tup[i]), grad(D_2 * sec_termj_Expression)) * dx)
    
    L.append( f_n_tup[i] * v_tup[i] * dx)

a = sum(a)
L = sum(L)

# Start of Time stepping

for t in range(1,timesteps+1):
    print("%s - time is %d" % (datetime.now(), t))

    # solve(a==L, f, bc)
    solve(a==L, f)

    f_n.assign(f)

    f_n_savable = f_n.split()

    Kappa_iji = []
    Kappa_jij = []

    for i in range(len(key_list)):
        datai = data[i]
        Kappa_jij.append([])
        Kappa_iji.append([])
        for j in range(len(key_list)):
            dataj = data[j]
            Kappa_jij[i].append( constructions.convol(r, deltax, xdom, boundary_vec, f_n_savable[j], datai, dataj, a, b, degree=1) )
            Kappa_iji[i].append( constructions.convol(r, deltax, xdom, boundary_vec, f_n_savable[i], dataj, datai, a, b, degree=1) )

    # Figures are saved every 10% of time steps using MatPlotLib

    if t % (timesteps/10.0) == 0:
        for n, key in enumerate(key_list):

            # eps figures are saved for all species

            png_folder = path0 + str(true_key) + pic + str(key)
            png_filename = png_folder + '/timeis' + str(dt* float(t)) + '.eps'

            function_name = r'$f_{{{time}}}$'.format(time=t)

            plot ( f_n_savable[n])
            plt.xlabel(r'$w$')
            y = plt.ylabel(function_name)
            y.set_rotation(0)
            plt.yticks(rotation=65, ha='right')
            plt.savefig ( png_filename , format='eps')
            print("%s - Saving picture as %s" % (datetime.now(), png_filename))
            plt.close ( )

        # Brexit vote superimposed on an eps figure are saved for all ages

        for c0 in range(int(len(key_list)/2.0)):
            c1 = c0+4
            key0 = key_list[c0]
            mass0 = mass_float[c0]
            key1 = key_list[c1]
            mass1 = mass_float[c1]

            png_folder = path0 + str(true_key) + pic + str(key0[1])
            try:
                os.mkdir(png_folder)
            except OSError:
                pass
            
            png_filename = png_folder + '/timeis' + str(dt*float(t)) + '.eps'

            # Creating the correct scaling for plotting superimposed graphs on the same figure

            biggest = max(mass0, mass1)
            if biggest == mass0:
                scale_fac = mass0/mass1
                proj_f0 = project(f_n_savable[c0], V=V)
                proj_f1 = project(scale_fac*f_n_savable[c1], V=V)
            else:
                scale_fac = mass1/mass0
                proj_f0 = project(scale_fac*f_n_savable[c0], V=V)
                proj_f1 = project(f_n_savable[c1], V=V)

            proj_f0_vec = proj_f0.compute_vertex_values(mesh)
            proj_f1_vec = proj_f1.compute_vertex_values(mesh)
            largest_value = np.max(proj_f0_vec) <= np.max(proj_f1_vec)

            function_name = r'$f_{{{time}}}$'.format(time=t)

            if largest_value:
                plot ( proj_f1, color='xkcd:red', linestyle='--', label=str(key1) )
                plot ( proj_f0, color='xkcd:blue', label=str(key0) )  
            else:
                plot ( proj_f0,color='xkcd:blue' , label=str(key0) )
                plot ( proj_f1, color='xkcd:red', linestyle='--', label=str(key1) )  


            plt.xlabel(r'$w$')
            y = plt.ylabel(function_name)
            y.set_rotation(0)
            plt.yticks(rotation=65, ha='right')
            plt.savefig ( png_filename , format='eps')

            print("%s - Saving picture as %s" % (datetime.now(), png_filename))
            plt.close ( )
            if t*dt >= T-DOLFIN_EPS:
                csv_folder = path0 + str(true_key) + comma + str(key0[1])
                tf = csv_folder + '/finaltime.csv'


                try:
                    os.mkdir(csv_folder)
                except OSError:
                    pass

                print("%s - Creating file %s " % (datetime.now(), tf))

                with open(tf, 'w') as resultfile:
                    wr = csv.writer(resultfile)
                    
                    wr.writerow(proj_f0_vec)
                    wr.writerow(proj_f1_vec)
                    wr.writerow(xdom)

print_tup = expressions.find_inf_masses(f_n_savable, boundary_vec, W, key_list)

for i in range(len(f_n_savable)):
    print('#--------------------------------Final masses for %s--------------------------------#\nLabour: %f\nGreen: %f\nLibdems: %f\nOther: %f\nConservative: %f\nBrexit: %f' % print_tup[i])

infl_mass = [print_tup[i][1:] for i in range(total_no)]

# Saving the full results as a csv file to be used later

with open(path0 + true_key + '/weights' + '.csv', 'w') as resultfile:
    wr = csv.writer(resultfile)
    wr.writerows(infl_mass)
