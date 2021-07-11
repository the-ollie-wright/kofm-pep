import numpy as np
import math

def D_2(w,k,a,b):
    return (1-(k*w)**a)**b


def setup_space_domain(start, end, N):
    deltax = (end-start)/float(N)

    xdom = [start + deltax*i for i in range(N+1)]

    return (xdom, deltax)


def Heaviside(x):
    if x<=0:
        return 0.0
    else:
        return 1.0


def binsearch(val, b):
    low = 0
    high = len(b) -1
    while low <= high:
        mid = math.ceil(((high+low)/2.0))
        if b[mid-1]<= val <=b[mid] or mid == 1 or mid == len(b) - 1:
            return mid - 1
        elif b[mid-1] > val:
            high = mid - 1
        elif b[mid] < val:
            low = mid + 1


def create_prob_list(M, leaderpercent, M_, tauL, taui):
    M_keys = M_.keys()
    
    intimediate_dict = {}

    for key in M_keys:
        if key == 'L':
            intimediate_dict[key] = float(M_[key])
        else:
            intimediate_dict[key] = float(M_[key])
    
    total_M = sum(intimediate_dict.values())

    intimediate_list = [float(value)/total_M for value in intimediate_dict.values()]

    runningtotal = 0
    P_chooselist = [0]
    for key, prob in zip(M_keys, intimediate_list):
        runningtotal += prob
        P_chooselist.append(runningtotal)

    return P_chooselist


def F(x, a, b):
    return (1 - (0.25 * (2 + x)) ** a) ** b

def P(w, v, a, b, r, boundary_vec, data0='L', data1='L'):
    
    if abs(w-v) <= r:
    
        if data0 == 'L' or data1 == 'L':
            return 1.0
        else: 
            p00 = data0[0][binsearch(v, boundary_vec)]
            p01 = data1[0][binsearch(w, boundary_vec)]
            p10 = data0[1][binsearch(v, boundary_vec)]
            p11 = data1[1][binsearch(w, boundary_vec)]

            P_char = F(p01  - p00 , a, b) * F(p11  - p10, a, b)
            return P_char
    else:
        return 0.0


def choose_indivs(opin_dict, P_chooselist, keylist):
    U0 = np.random.uniform(0,1)
    key0 = keylist[binsearch(U0, P_chooselist)]
    index0 = np.random.randint(len(opin_dict[key0]))
    opin0 = opin_dict[key0].pop(index0)
    indiv0 = [key0, opin0]

    U1 = np.random.uniform(0,1)
    key1 = keylist[binsearch(U1, P_chooselist)]
    index1 = np.random.randint(len(opin_dict[key1]))
    opin1 = opin_dict[key1].pop(index1)
    indiv1 = [key1, opin1]

    return (indiv0, indiv1)


def interaction(indiv0, indiv1, a, b, alpha, sigma, r, boundary_vec, characteristicdata, gamma):
    key0, opin0 = indiv0
    key1, opin1 = indiv1

    if key0 == 'L' and key1 == 'L':

        new_opin0 = opin0
        new_opin1 = opin1

    elif key0 == 'L':

        new_opin1 = opin1 + gamma * P(opin1, opin0, a, b, r, boundary_vec)*(opin0-opin1) + D_2(opin1,1.0, 2.0, alpha) * (2*np.random.randint(2)-1) * sigma
        new_opin0 = opin0

    elif key1 == 'L':

        new_opin0 = opin0 + gamma * P(opin0, opin1, a, b, r, boundary_vec)*(opin1-opin0) + D_2(opin0,1.0, 2.0, alpha) * (2*np.random.randint(2)-1) * sigma
        new_opin1 = opin1

    else:
        
        data0 = [characteristicdata[k] for k in key0]
        data1 = [characteristicdata[k] for k in key1]

        new_opin0 = opin0 + gamma * P(opin0, opin1, a, b, r, boundary_vec, data1, data0)*(opin1-opin0) + D_2(opin0,1.0, 2.0, alpha) * (2*np.random.randint(2)-1) * sigma
        new_opin1 = opin1 + gamma * P(opin1, opin0, a, b, r, boundary_vec, data0, data1)*(opin0-opin1) + D_2(opin1,1.0, 2.0, alpha) * (2*np.random.randint(2)-1) * sigma


    new_indiv0 = [key0, new_opin0]
    new_indiv1 = [key1, new_opin1]

    return (new_indiv0, new_indiv1)

