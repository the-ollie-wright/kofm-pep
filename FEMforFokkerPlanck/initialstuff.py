from fenics import *
import expressions
import csvget
import functions as funs

class f_initialA(UserExpression):
    def __init__(self, parties, data, our_boundary, total_no, **kwargs):
        self.parties = parties
        self.data = data
        self.our_boundary = our_boundary
        self.total_no = total_no
        super().__init__(degree=kwargs["degree"])

    def eval(self, value, x):
        for i in range(self.total_no):
            value[i] = initial_funcA(x[0], self.data[i], self.parties, self.our_boundary)

    def value_shape(self):
        return (8,)


class f_initialB(UserExpression):
    def __init__(self, radius, centre, heights, total_no, **kwargs):
        self.radius = radius
        self.centre = centre
        self.total_no = total_no
        self.heights = heights
        super().__init__(degree=kwargs["degree"])

    def eval(self, value, x):
        for i in range(self.total_no):
            value[i] = initial_funcB(x[0], self.radius, self. centre, self.heights[i])


    def value_shape(self):
        return (8,)


def initial_funcA(x, data, parties, our_boundary):
    ret_vals = [data[i]*(2.0*our_boundary[i]) if abs(x-parties[i]) < our_boundary[i] else 0 for i in range(len(parties))]
    return sum(ret_vals)


def initial_funcB(x, radius, centre, height):
    return height * (1.0 / (2.0 * radius)) * 1.0 if (radius - abs(x - centre)) >= 0 else 0.0


def initialA(parties, data, our_boundary, vec_space):
    return interpolate(f_initialA(parties, data, our_boundary, total_no=8, degree=0), vec_space)



def initialB(radius, centre, data, vec_space):
    heights = [sum(dat) for dat in data]
    return interpolate(f_initialB(radius, centre , heights, total_no=8, degree=0), vec_space)
