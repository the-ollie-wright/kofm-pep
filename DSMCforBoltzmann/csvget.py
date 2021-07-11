import numpy as np
import os

def read_data(dataset, scale=1.0):
    vec = []
    with open(dataset, 'r') as data: #
        lines = data.readlines() #
        lines.pop(0) #
        for line in lines: #
            datcsv = line.split(',') #
            datcsv = datcsv[1:]
            for j in range(len(datcsv)):
                datcsv[j] = float(datcsv[j]) / float(scale)
            vec.append(datcsv)
        vec = np.array(vec).transpose()
    return vec

def read_characteristic(key, party_no):
    # sets up the characteristic matrix. The columns represent the number in the key, i.e. (0, 2) will require column 0 and column 2 etc.
    brexdataset = 'dataandsolutions/dataandsolutions' + key + '/data/brexitdata.csv'
    agedataset = 'dataandsolutions/dataandsolutions' + key + '/data/agedata.csv'
    brex_vec = read_data(brexdataset)
    age_vec = read_data(agedataset)

    return prepare_characteristic(party_no, brex_vec, age_vec)


def prepare_characteristic(party_no, *args):
    length = [np.shape(a)[0] for a in args]
    total_length = sum(length)
    characteristicdata = np.empty([total_length, party_no])
    counter = 0
    current_length = 0
    for n, arg in enumerate(args):
        prev_length = current_length
        current_length += length[n]
        for line in arg:
            characteristicdata[counter] = line
            counter += 1

    return characteristicdata

