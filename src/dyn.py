# Christopher Iliffe Sprague
# christopher.iliffe.sprague@gmail.com

class System(object):

    def __init__(self, ulb, uub, sdim, params, funcs):

        self.ulb    = ulb
        self.uub    = uub
        self.udim   = len(self.ulb)
        self.sdim   = sdim
        self.params = params
        self.funcs  = funcs
        print(self.params)

    def lagrangian(self, u, a):
        return self.funcs.lagrangian(*u, a, *self.params)

    def hamiltonian(self, fs, u, a):
        return self.funcs.hamiltonian(*fs, *u, a, *self.params)

    def pontryagin(self, fs, a, bound):
        u = self.funcs.pontryagin(*fs, a, *self.params).reshape(self.udim)
        if bound:
            u = [min(max(u, lb), ub) for u, lb, ub in zip(u, self.ulb, self.uub)]
        return u

    def eom_fullstate(self, fs, u):
        return self.funcs.eom_fullstate(*fs, *u, *self.params).reshape(self.sdim*2)

    def eom_fullstate_jac(self, fs, u):
        return self.funcs.eom_fullstate_jac(*fs, *u, *self.params)
