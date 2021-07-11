from fenics import *
import numpy as np


def setup_mesh_and_space(start, stop, N, n, functype, total_no):
    if functype == 'DPG' or functype == 'DG':
        mesh = IntervalMesh(N, start, stop)
        V = FunctionSpace(mesh, functype, n)
        V_cg = FunctionSpace(mesh, 'CG', n)

        P1 = FiniteElement(functype, interval, n)
        element_stuff = [P1 for i in range(total_no)]
        element = MixedElement(element_stuff)

        W = FunctionSpace(mesh, element)
        return (mesh, V, W, V_cg)

    else:
        mesh = IntervalMesh(N, start, stop)
        V = FunctionSpace(mesh, functype, n)

        P1 = FiniteElement(functype, interval, n)
        element_stuff = [P1 for i in range(total_no)]
        element = MixedElement(element_stuff)

        W = FunctionSpace(mesh, element)
        return (mesh, V, W)


def Heaviside(x):
    if x<=0:
        return 0.0
    else:
        return 1.0

