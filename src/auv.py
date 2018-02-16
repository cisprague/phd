# Christopher Iliffe Sprague
# christopher.iliffe.sprague@gmail.com

import numpy as np


class AUV(object):

    def __init__(
        self,
        m = 100,
        g = 9.82,
        rho = 1000,
        cd = 0.48,
        r = 0.5,
        T = 10,
        rot = 1,
        cv = 0
    ):

        # random fluid velocity direction
        t = np.random.uniform(5. * np.pi / 12., 7. * np.pi / 12.)
        p = np.random.uniform(0, 2 * np.pi)
        st = np.sin(t)
        ct = np.cos(t)
        sp = np.sin(p)
        cp = np.cos(p)
        v = np.array([st * cp, st * sp, ct])

        # set parametres
        self.mass = float(m)
        self.gravity = float(g)
        self.density = float(rho)
        self.cd = float(cd)
        self.thrust = float(T)
        self.spin = float(rot)
        self.current = float(cv) * v
        self.planaform = np.pi * r**2
        self.volume = 4. * self.planaform * r / 3.

        # set state and action bounds
        self.slb = [-1000, -1000, 0, -10, -10, -10, -1, -1, -1, -1]
        self.sub = [1000, 1000, 1000, 10, 10, 10, 1, 1, 1, 1]
        self.alb = [0, 0.2, -1, -1]
        self.aub = [1, 1, 1, 1]

        """
        # nondimensionalisation units
        self.umas = float(m)
        self.upos = float(1000)
        self.uvel = float(10)
        self.uacc = self.uvel**2/self.upos
        self.uare = self.upos**2
        self.uvol = self.upos**3
        self.utim = self.upos/self.uvel
        self.uden = self.umas/self.uvol
        self.ufor = self.umas*self.uacc
        self.uome = np.pi/self.utim

        # set nondimensional parametres
        self.mass      /= self.umas
        self.gravity   /= self.uacc
        self.density   /= self.uden
        self.thrust    /= self.ufor
        self.spin      /= self.uome
        self.current   /= self.uvel
        self.planaform /= self.uare
        self.volume    /= self.uvol
        """

    def eom_fullstate(self, fullstate, control):

        # extract fullstate
        x, y, z, vx, vy, vz, qr, qx, qy, qz, lx, ly, lz, lvx, lvy, lvz, lqr, lqx, lqy, lqz = fullstate

        # extract control
        ut, ub, ux, uy = control

        # free stream velocity
        vinfx = vx - self.current[0]
        vinfy = vy - self.current[1]
        vinfz = vz - self.current[2]

        # common subexpression elimination
        x0 = 1 / self.mass
        x1 = self.planaform * self.cd * self.density / 2
        x2 = self.thrust * ut
        x3 = 2 * qr
        x4 = 2 * qz
        x5 = self.spin * qy / 2
        x6 = self.spin * qz / 2
        x7 = self.spin * qr / 2
        x8 = self.spin * qx / 2
        x9 = self.planaform * self.cd * self.density * x0 / 2
        x10 = lqy * self.spin / 2
        x11 = lqz * self.spin / 2
        x12 = 2 * self.thrust * lvx * ut * x0
        x13 = 2 * self.thrust * lvy * ut * x0
        x14 = self.thrust * lvx * ut * x0
        x15 = self.thrust * lvy * ut * x0
        x16 = 4 * self.thrust * lvz * ut * x0
        x17 = lqr * self.spin / 2
        x18 = lqx * self.spin / 2

        # return state transition
        return np.array([
            vx,
            vy,
            vz,
            x0 * (-x1 * (vinfx) + x2 * (qx * x4 - qy * x3)),
            x0 * (-x1 * (vinfy) + x2 * (-qx * x3 + qy * x4)),
            x0 * (self.volume * self.density * self.gravity * ub - self.gravity *
                  self.mass - x1 * (vinfz) + x2 * (-2 * qx**2 - 2 * qy**2 + 1)),
            ux * x6 - uy * x5,
            ux * x5 + uy * x6,
            -ux * x8 + uy * x7,
            -ux * x7 - uy * x8,
            0,
            0,
            0,
            -lx + lvx * x9,
            -ly + lvy * x9,
            -lz + lvz * x9,
            qx * x13 + qy * x12 + ux * x11 - uy * x10,
            qx * x16 + ux * x10 + uy * x11 - x14 * x4 + x15 * x3,
            qy * x16 - ux * x18 + uy * x17 + x14 * x3 - x15 * x4,
            -qx * x12 - qy * x13 - ux * x17 - uy * x18
        ])

    def pontryagin(self, fullstate, alpha, bound=True):

        # extract fullstate
        x, y, z, vx, vy, vz, qr, qx, qy, qz, \
            lx, ly, lz, lvx, lvy, lvz, lqr, lqx, lqy, lqz = fullstate

        # common subexpression elimination
        x0 = 1 / (alpha - 1)
        x1 = x0 / self.mass
        x2 = self.thrust * lvz
        x3 = (alpha * self.mass) / 2
        x4 = self.thrust * lvx
        x5 = self.thrust * lvy
        x6 = x0 / 4
        x7 = 2 * alpha
        x8 = lqr * self.spin
        x9 = lqx * self.spin
        x10 = lqy * self.spin
        x11 = lqz * self.spin

        # optimal controls
        ut = x1 * (-qr * qx * x5 - qr * qy * x4 - qx**2 * x2 + qx *
                   qz * x4 - qy**2 * x2 + qy * qz * x5 + x2 / 2 + x3)
        ub = x1 * (self.volume * lvz * self.density * self.gravity / 2 + x3)
        ux = x6 * (-qr * x11 - qx * x10 + qy * x9 + qz * x8 + x7)
        uy = x6 * (qr * x10 - qx * x11 - qy * x8 + qz * x9 + x7)

        # if bounded
        if bound:
            ut = min(max(ut, self.alb[0]), self.aub[0])
            ub = min(max(ub, self.alb[1]), self.aub[1])
            ux = min(max(ux, self.alb[2]), self.aub[2])
            uy = min(max(uy, self.alb[3]), self.aub[3])

        # return controls
        return np.array([ut, ub, ux, uy])

    def hamiltonian(self, fullstate, control, alpha):

        # extract fullstate
        x, y, z, vx, vy, vz, qr, qx, qy, qz, \
            lx, ly, lz, lvx, lvy, lvz, lqr, lqx, lqy, lqz = fullstate

        # extract control
        ut, ub, ux, uy = control

        # free stream velocity
        vinfx = vx - self.current[0]
        vinfy = vy - self.current[1]
        vinfz = vz - self.current[2]

        # common subexpression elimination
        x0 = self.spin * qy / 2
        x1 = self.spin * qz / 2
        x2 = self.spin * qr / 2
        x3 = self.spin * qx / 2
        x4 = 1 / self.mass
        x5 = self.planaform * self.cd * self.density / 2
        x6 = self.thrust * ut
        x7 = 2 * qr
        x8 = 2 * qz

        return alpha * (ub + ut + ux + uy) + \
            lx * vx + ly * vy + lz * vz + \
            lqr * (ux * x1 - uy * x0) + \
            lqx * (ux * x0 + uy * x1) + \
            lqy * (-ux * x3 + uy * x2) + \
            lqz * (-ux * x2 - uy * x3) + \
            lvx * x4 * (-x5 * (vinfx) + x6 * (qx * x8 - qy * x7)) + \
            lvy * x4 * (-x5 * (vinfy) + x6 * (-qx * x7 + qy * x8)) + \
            lvz * x4 * (
                self.volume * self.density * self.gravity * ub -
                self.gravity * self.mass -
                x5 * (vinfz) +
                x6 * (-2 * qx**2 - 2 * qy**2 + 1)
            ) + (-alpha + 1) * (ub**2 + ut**2 + ux**2 + uy**2)
