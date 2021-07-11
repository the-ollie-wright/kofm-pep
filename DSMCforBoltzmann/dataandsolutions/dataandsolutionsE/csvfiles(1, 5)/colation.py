import numpy as np
from numpy import random as ran
import copy
import csv
import tempfile


matrix0 = []
matrix1 = []
matrix2 = []
matrix3 = []
matrix4 = []
matrix5 = []
matrix6 = []
matrix7 = []
matrix8 = []

with open('finaltime0.csv', 'r') as tmp:
    # reads the percentage of intended votes for each party into a nested list
    lines0 = tmp.readlines()
    for i, line in enumerate(lines0):
        datcsv0 = line.split(',')
        for n, x in enumerate(datcsv0):
            datcsv0[n] = float(x)
        matrix0.append(datcsv0)

with open('finaltime1.csv', 'r') as tmp:
    # reads the percentage of intended votes for each party into a nested list
    lines1 = tmp.readlines()
    for i, line in enumerate(lines1):
        datcsv1 = line.split(',')
        for n, x in enumerate(datcsv1):
            datcsv1[n] = float(x)
        matrix1.append(datcsv1)

with open('finaltime2.csv', 'r') as tmp:
    # reads the percentage of intended votes for each party into a nested list
    lines2 = tmp.readlines()
    for i, line in enumerate(lines2):
        datcsv2 = line.split(',')
        for n, x in enumerate(datcsv2):
            datcsv2[n] = float(x)
        matrix2.append(datcsv2)

with open('finaltime3.csv', 'r') as tmp:
    # reads the percentage of intended votes for each party into a nested list
    lines3 = tmp.readlines()
    for i, line in enumerate(lines3):
        datcsv3 = line.split(',')
        for n, x in enumerate(datcsv3):
            datcsv3[n] = float(x)
        matrix3.append(datcsv3)

with open('finaltime4.csv', 'r') as tmp:
    # reads the percentage of intended votes for each party into a nested list
    lines4 = tmp.readlines()
    for i, line in enumerate(lines4):
        datcsv4 = line.split(',')
        for n, x in enumerate(datcsv4):
            datcsv4[n] = float(x)
        matrix4.append(datcsv4)

with open('finaltime5.csv', 'r') as tmp:
    # reads the percentage of intended votes for each party into a nested list
    lines5 = tmp.readlines()
    for i, line in enumerate(lines5):
        datcsv5 = line.split(',')
        for n, x in enumerate(datcsv5):
            datcsv5[n] = float(x)
        matrix5.append(datcsv5)

with open('finaltime6.csv', 'r') as tmp:
    # reads the percentage of intended votes for each party into a nested list
    lines6 = tmp.readlines()
    for i, line in enumerate(lines6):
        datcsv6 = line.split(',')
        for n, x in enumerate(datcsv6):
            datcsv6[n] = float(x)
        matrix6.append(datcsv6)

with open('finaltime7.csv', 'r') as tmp:
    # reads the percentage of intended votes for each party into a nested list
    lines7 = tmp.readlines()
    for i, line in enumerate(lines7):
        datcsv7 = line.split(',')
        for n, x in enumerate(datcsv7):
            datcsv7[n] = float(x)
        matrix7.append(datcsv7)

with open('finaltime8.csv', 'r') as tmp:
    # reads the percentage of intended votes for each party into a nested list
    lines8 = tmp.readlines()
    for i, line in enumerate(lines8):
        datcsv8 = line.split(',')
        for n, x in enumerate(datcsv8):
            datcsv8[n] = float(x)
        matrix8.append(datcsv8)


final = [[] for i in range(len(matrix0))]

for n, demog0 in enumerate(matrix0):
    demog = []
    demog.append(demog0)
    demog.append(matrix1[n])
    demog.append(matrix2[n])
    demog.append(matrix3[n])
    demog.append(matrix4[n])
    demog.append(matrix5[n])
    demog.append(matrix6[n])
    demog.append(matrix7[n])
    demog.append(matrix8[n])
    demog = np.array(demog).transpose()
    
    for line in demog:
        final[n].append(sum(list(line))/len(list(line)))

with open('solution.csv','w') as resultfile:
    wr = csv.writer(resultfile)
    wr.writerows(final)