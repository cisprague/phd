# Christopher Iliffe Sprague
# christopher.iliffe.sprague@gmail.com

from scipy.integrate import ode
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


class Leg(object):

    def __init__(self, dynamics, alpha=0, freetime=True, bound=True):

        # dynamical system
        self.dynamics = dynamics

        # homotopy parametre
        self.alpha = alpha

        # control bound
        self.bound = bound

        # infinite time horizon
        self.freetime = bool(freetime)

        # default equality constraint dimension
        if self.freetime:
            self.cdim = 11
        else:
            self.nec = 10

        # numerical integrator
        self.integrator = ode(
            lambda t, fs: self.dynamics.eom_fullstate(
                fs, self.dynamics.pontryagin(fs, self.alpha, self.bound)
            )
        )

    def recorder(self, t, fs):

        # times
        self.times = np.append(self.times, t)

        # fullstates
        self.states = np.vstack((self.states, fs))

        # actions
        self.actions = np.vstack((
            self.actions,
            self.dynamics.pontryagin(fs, self.alpha, self.bound)
        ))

    def set_times(self, t0, tf):
        self.t0 = float(t0)
        self.tf = float(tf)

    def set_states(self, s0, sf):
        self.s0 = np.array(s0, dtype=float)
        self.sf = np.array(sf, dtype=float)

    def set_costates(self, l0):
        self.l0 = np.array(l0, dtype=float)

    def set(self, t0, s0, l0, tf, sf):
        self.set_times(t0, tf)
        self.set_states(s0, sf)
        self.set_costates(l0)

    def set_params(self, alpha, bound, freetime):

        # homotopy parametre
        self.alpha = alpha

        # control bounds
        self.bound = bool(bound)

        # inifinite time horizon
        self.freetime = bool(freetime)

        # equality constraint dimension
        if self.freetime:
            self.nec = 11
        else:
            self.nec = 10

    def propagate(self, atol=1e-5, rtol=1e-5):

        # departure fullstate
        fs0 = np.hstack((self.s0, self.l0))

        # nondimensionalise
        #fs0[0:3] /= self.dynamics.upos
        #fs0[3:6] /= self.dynamics.uvel

        # reset trajectory records
        self.times = np.empty((1, 0))
        self.states = np.empty((0, 20))
        self.actions = np.empty((0, 4))

        # set integration method
        self.integrator.set_integrator("dop853", atol=atol, rtol=rtol)

        # set recorder
        self.integrator.set_solout(self.recorder)

        # set departure configuration
        #self.integrator.set_initial_value(fs0, self.t0/self.dynamics.utim)
        self.integrator.set_initial_value(fs0, self.t0)

        # numerically integrate
        #self.integrator.integrate(self.tf/self.dynamics.utim)
        self.integrator.integrate(self.tf)

    def mismatch(self, atol=1e-5, rtol=1e-5):

        # propagte trajectory
        self.propagate(atol=atol, rtol=rtol)

        # final state mismatch equality constraint
        ceq = self.states[-1, 0:10] - self.sf

        # infinite time horizon
        if self.freetime:

            # final Hamiltonian
            H = self.dynamics.hamiltonian(
                self.states[-1], self.actions[-1], self.alpha
            )

            # append equality constraint
            ceq = np.hstack((ceq, [H]))

        return ceq

    def plot_traj(self, ax=None, mark="k.-", quiver=True):

        # create new axes if not supplied
        if ax is None:
            fig = plt.figure()
            ax = fig.gca(projection="3d")

        # plot the positions
        ax.plot(*[self.states[:, dim] for dim in [0, 1, 2]], mark)

        # show thrust profile if desired
        if quiver:

            # quaternions
            qr, qx, qy, qz = [self.states[:, dim] for dim in [6, 7, 8, 9]]

            # Rene Descartes directions
            utx = 2 * (np.multiply(qx, qz) - np.multiply(qy, qr))
            uty = 2 * (np.multiply(qy, qz) - np.multiply(qx, qr))
            utz = 1 - 2 * (np.square(qx) + np.square(qy))

            # thrust magnitudes
            utx, uty, utz = [np.multiply(self.actions[:, 0], i) for i in [utx, uty, utz]]

            ax.set_aspect("equal")
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            # plot thrusts
            ax.quiver(
                *[self.states[:, dim] for dim in [0, 1, 2]],
                utx, uty, utz,
                length=0.01
            )

        return ax

    def plot_controls(self):

        # create subplots
        f, ax = plt.subplots(4, sharex=True)

        for i in range(4):
            ax[i].plot(self.actions[:, i], "k.-")

        plt.show()
