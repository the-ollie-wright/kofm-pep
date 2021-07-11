import numpy as np
import csvget

def setup_follower_species(M, keylist, B, true_key):
    opin_dict = {}
    M_ = {}

    for key in keylist:

        opin_vec = []

        # looking up the weighting of the initial values

        weights = csvget.read_data('dataandsolutions/dataandsolutions' + str(true_key) + '/data/initial' + str(key) + '.csv')[0]
        total_weight = sum(weights)

        # adding the number of individuals in the key's species to a dictionary for later use

        M_[key] = int(M*total_weight)

        # creating the list that use assigns an opinion to each member of the key's species

        for j in range(len(B)-1):
            M_weight = int(M*weights[j])
            for k in range(M_weight):
                indiv = np.random.uniform(B[j], B[j+1])
                opin_vec.append( indiv )

        opin_dict[key] = opin_vec


    return (opin_dict, M_)