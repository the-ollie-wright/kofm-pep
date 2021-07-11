from fenics import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import csvget
import numpy as np

vec0 = []
vec1 = []

with open('dataandsolutions/dataandsolutionsE/csvfiles(1, 5)/finaltime0.csv', 'r') as data: #
    lines = data.readlines() #
    for line in lines: #
        datcsv = line.split(',') #
        for j in range(len(datcsv)):
            datcsv[j] = float(datcsv[j])
        vec0.append(datcsv)

x_0 = list(vec0[1])
y_0 = list(vec0[0])
y_0.append(0.0)

plt.bar(x_0, y_0, alpha=0.5, width=0.04)

with open('dataandsolutions/dataandsolutionsE/csvfilesfenics/finaltime.csv', 'r') as data:
    lines = data.readlines() #
    for line in lines: #
        datcsv = line.split(',') #
        for j in range(len(datcsv)):
            datcsv[j] = float(datcsv[j])
        vec1.append(datcsv)

x_1 = list(vec1[2])
y_1 = list(vec1[1])

plt.plot(x_1, y_1, color='xkcd:red', linestyle='--')

function_name = r'$f_{1}$'
plt.xlabel(r'$w$')
y_0 = plt.ylabel(function_name)
y_0.set_rotation(0)
plt.yticks(rotation=65, ha='right')

plt.savefig ( 'dataandsolutions/dataandsolutionsE/csvfiles(1, 5)/graphact.eps' , format='eps')

plt.close ( )
