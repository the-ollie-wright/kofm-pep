import math
import csvget
import numpy as np

def setup_leader_species(M,leaderpercent,true_key, parties):
    M_lead = leaderpercent*M

    party_weights = csvget.read_data('dataandsolutions/dataandsolutions' + str(true_key) + '/data/lastelection.csv')[0]
    
    
    M_lead_weighted = [int(M_lead*percent) for percent in party_weights]
    indiv_list = []
    for weight, party in zip(M_lead_weighted, parties):
        indiv_list.extend( [party for i in range(weight)] )
    
    keys = 'L'

    return (indiv_list, keys)