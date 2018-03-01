# Christopher Iliffe Sprague
# christopher.iliffe.sprague@gmail.com

import numpy as np
from auv_dynamics import *

class AUV(object):

    def __init__(
        self,
        m   = 100,
        rho = 1000,
        cd  = 0.48,
        r   = 0.5,
        T   = 10,
        w   = 0.2,
        vf  = 1
    ):

        # fluid velocity
        t  = np.random.uniform(5*np.pi/12, 7*np.pi/12)
        p  = np.random.uniform(0, 2*np.pi)
        st = np.sin(t)
        ct = np.cos(t)
        sp = np.sin(p)
        cp = np.cos(p)
        vf = st * cp, st * sp, ct

        # parametres
        self.params = [*vf, m, rho, cd, np.pi*r**2, T, w]
        self.slb  = [-1000, -1000, 0, -10, -10, -10, -1, -1, -1, -1]
        self.sub  = [1000, 1000, 0, 10, 10, 10, 1, 1, 1, 1]
        self.ulb  = [0, -1, -1]
        self.uub  = [1, 1, 1]
        self.sdim = len(self.slb)
        self.udim = len(self.ulb)

    def lagrangian(self, u, a):
        return lagrangian(*u, a)

    def hamiltonian(self, fs, u, a):
        return hamiltonian(*fs, *u, a, *self.params)

    def pontryagin(self, fs, a, bound):
        u = pontryagin(*fs, a, *self.params).reshape(self.udim)
        if bound:
            u = [min(max(u, lb), ub) for u, lb, ub in zip(u, self.ulb, self.uub)]
        return u

    def eom_fullstate(self, fs, u):
        return eom_fullstate(*fs, *u, *self.params).reshape(self.sdim*2)

    def eom_fullstate_jac(self, fs, u):
        return eom_fullstate_jac(*fs, *u, *self.params)
