# Christopher Iliffe Sprague
# christopher.iliffe.sprague@gmail.com

import numpy as np

def main():

    from auv import AUV
    from leg import Leg

    # instantiate AUV
    auv = AUV()

    # instantiate Leg
    leg = Leg(auv, alpha=0.99, bound=True)

    # departure and arrival states
    s0 = np.array([10, 10, 10, 1, 1, 1, 1, 0, 0, 0])
    sf = np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0])

    # departure and arrival times
    t0 = 0
    tf = 10000

    # initial costates
    l0 = np.ones(10)

    # set leg
    leg.set(t0, s0, l0, tf, sf)

    # propagate
    leg.propagate(atol=1e-12, rtol=1e-12)

    # plot trajectory
    leg.plot_traj()

    # plot controls
    leg.plot_controls()

    print(leg.actions)

if __name__ == "__main__":

    main()
