from fenics import *
import expressions
import csvget
import functions as funs
import initialstuff

# Expressions specific to the intitial conditions and Ks

class convol(UserExpression):
    def __init__(self, r, deltax, xdom, bounds, func, datai, dataj, a, b, **kwargs):
        super().__init__(degree=kwargs["degree"])
        self.r = r
        self.deltax = deltax
        self.xdom = xdom
        self.func = func
        self.bounds = bounds
        self.datai = datai
        self.dataj = dataj
        self.a = a
        self.b = b

    def eval(self, value, x):
    
        con_step = []
        for n, y in enumerate(self.xdom):
            if abs(y-x[0]) <= self.r:
                
                # For the proposed model in D\"uring and Wright (2021)

                p00 = self.datai[0][binsearch(y, self.bounds)]
                p01 = self.dataj[0][binsearch(x[0], self.bounds)]
                p10 = self.datai[1][binsearch(y, self.bounds)]
                p11 = self.dataj[1][binsearch(x[0], self.bounds)]

                P_char = F(p01  - p00 , self.a, self.b) * F( p11  - p01 , self.a, self.b)
                con_eval = P_char * (y-x[0]) * self.func(y)
                
                # For the equivalent model in D\"uring et al. (2009)
                
                # con_eval = (y-x[0]) * self.func(y)
                                
                con_step.append(con_eval)

        value[0] = self.deltax*sum(con_step)

    def value_shape(self):
        return ()

# Functions specific to the initial conditions and Ks

def read_initial(region, demog, path0):
    region = str(region)
    demog = str(demog)
    dataset = path0 + region + '/data/initial' + demog + '.csv'
    return csvget.read_data(dataset, 100.0)

def F(x, a, b):
    return (1 - (0.25 * (2 + x)) ** a) ** b


def binsearch(val, b):
    low = 0
    high = len(b) -1
    while low <= high:
        mid = (high+low)//2
        if b[mid-1]<= val <=b[mid]:
            return mid - 1
        elif b[mid-1] > val:
            high = mid - 1
        elif b[mid] < val:
            low = mid + 1


# Final set up of the non-Leader Species stuff

def setup_follower_species(space, vec_space, parties, brex_no, age_no, true_key, r, deltax, xdom, boundary_vec, characteristicdata, a, b, path0, key_list):
    our_boundary = [0.15, 0.15, 0.1, 0.1, 0.15, 0.15]

    f = Function(vec_space)
    f_ = TrialFunction(vec_space)
    v = TestFunction(vec_space)

    init_data = []

    for key in key_list:
            init_data.append(read_initial(true_key, key, path0)[0])

    # Initial conditions as described in D\"uring and Wright (2021)

    f_n = initialstuff.initialA(parties, init_data, our_boundary, vec_space)

    # Uniform initial conditions

    # f_n = initialstuff.initialB(0.96, 0.0, init_data, vec_space)


    f_n_tup = f_n.split()

    mass_float = [assemble(f_i * dx) for f_i in f_n_tup]

    mass = [Constant(m) for m in mass_float]

    mass_sum = sum(mass)

    # Data for the non-local convolution operators (the Ks)

    data = []

    for key in key_list:
        data.append([characteristicdata[k] for k in key])

    # non-local convolution operators (the Ks) set up

    Kappa_iji = []
    Kappa_jij = []

    for i in range(len(key_list)):
        datai = data[i]
        Kappa_jij.append([])
        Kappa_iji.append([])
        for j in range(len(key_list)):
            dataj = data[j]
            Kappa_jij[i].append( convol(r, deltax, xdom, boundary_vec, f_n_tup[j], datai, dataj, a, b, degree=1) )
            Kappa_iji[i].append( convol(r, deltax, xdom, boundary_vec, f_n_tup[i], dataj, datai, a, b, degree=1) )

    functions = (f, f_, f_n, v)
    convol_operators = (Kappa_jij, Kappa_iji)
    masses = (mass, mass_float, mass_sum)

    return [functions, convol_operators, masses, key_list, data]
