from fenics import *
import expressions
import csvget
import functions as funs
import math

# Leader Specific Classes

class InitialStuff(UserExpression):
    def __init__(self, height, party, scale, **kwargs):
        super().__init__(degree=kwargs["degree"])
        self.height = height
        self.party = party
        self.scale = scale

    def eval(self, value, x):
        val = [mollifier(x[0],p,self.scale,h) if abs(x[0]-p) < self.scale else 0 for h, p in zip(self.height, self.party)]
        value[0] = sum(val)

    def value_shape(self):
        return ()


class K_L(UserExpression):
    def __init__(self, func, r, xdom, deltax, mass, **kwargs):
        super().__init__(degree=kwargs["degree"])
        self.func = func
        self.r = r
        self.xdom = xdom
        self.deltax = deltax
        self.mass = mass

    def eval(self, value, x):
        con_step = [(y-x[0])*self.func(y) if abs(x[0] - y) <= self.r else 0 for y in self.xdom]
        
        value[0] = self.deltax * sum(con_step)*self.mass

    def value_shape(self):
        return ()

# Leader specific Functions

def mollifier(x, place, scale, height):
    if abs(x - place) < scale:
        val = (x - place) / scale
        return height * math.exp(-val ** 2 / (1 - val ** 2))
    else:
        return 0

def read_leader_data(true_key):
    dataset = 'dataandsolutions/dataandsolutions'+ str(true_key) + '/data/lastelection.csv'
    return csvget.read_data(dataset)


# Final set up of the Leader Species stuff

def setup_leader_species(space, parties, scale, r, xdom, deltax, true_key):
    Leader_heights = read_leader_data(true_key)[0]
    print(Leader_heights)

    f_L = InitialStuff(Leader_heights, parties, scale, degree=1)
    f_L = interpolate(f_L, space)
    mass = Constant(assemble(f_L * dx))
    mass_lead = 0.05 / mass


    Kappa_L = K_L(f_L, r, xdom, deltax, mass_lead, degree=2)

    return (Kappa_L, mass_lead, f_L)
