from fenics import *
import numpy as np

import functions as funs


class IndicExp(UserExpression):
    def __init__(self, a, b, total_no, **kwargs):
        self.a = a
        self.b = b
        self.total_no = total_no
        super().__init__(degree=kwargs["degree"])

    def eval(self, value, x):
        for v in range(self.total_no):
            value[v] = 1 if self.a < x[0] <= self.b else 0

    def value_shape(self):
        return (self.total_no,)


def find_inf_masses(f_n_tup, boundary_vec, vec_space, key_list):

    total_no = len(f_n_tup)
    indicator = [interpolate(IndicExp(boundary_vec[i], boundary_vec[i + 1], total_no,  degree=0), vec_space) for i in range(len(boundary_vec) - 1)]
    
    indicator = [ind.split() for ind in indicator]

    infl_mass = [[assemble(ind[0] * f_n_tup[i]*dx) for ind in indicator] for i in range(total_no)]

    print_tup = [(str(key_list[i]),) + tuple(infl_mass[i]) for i in range(total_no)]
    return print_tup
