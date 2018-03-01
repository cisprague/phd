# Christopher Iliffe Sprague
# christopher.iliffe.sprague@gmail.com

# import resources
from auv import AUV
from leg import Leg
from problem import Problem
import pygmo_plugins_nonfree as pg7
import numpy as np
import matplotlib.pyplot as plt
import pygmo as pg

def main():

    ''' AUV testing '''

    # instantiate AUV
    auv = AUV()

    # instantiate leg
    leg = Leg(auv, alpha=0, freetime=True)

    # departure and arrival states
    s0 = np.array([10, 10, 10, 0, 0, 0, 1, 0, 0, 0], dtype=float)
    sf = np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype=float)

    # set boundary states
    leg.set_states(s0, sf)

    """
    leg.set(0, s0, np.ones(len(s0)), 1000, sf)
    leg.propagate(atol=1e-12, rtol=1e-12)
    leg.plot_traj()
    leg.plot_states()
    leg.plot_actions()
    plt.show()
    """

    # create pygmo problem
    udp = Problem(leg)
    prob = pg.problem(udp)

    # instantiate SNOPT algorithm
    algo = pg7.snopt7(True, "/usr/lib/libsnopt7_cpp.so")
    algo.set_integer_option("Major iterations limit", 4000)
    algo.set_integer_option("Iterations limit", 40000)
    algo.set_numeric_option("Major optimality tolerance", 1e-2)
    algo.set_numeric_option("Major feasibility tolerance", 1e-8)

    algo = pg.algorithm(pg.ipopt("cobyla"))
    algo.xtol_rel = 1e-7
    algo = pg.algorithm(algo)
    algo.set_verbosity(1)

    # population
    for epo in range(1):

        # instantiate popoulation
        pop = pg.population(prob, 1)

        # evolve the solution
        pop = algo.evolve(pop)

    # set problem
    udp.fitness(pop.champion_x)

    # plot trajectory
    udp.leg.plot_traj()
    udp.leg.plot_actions()
    udp.leg.plot_states()

    plt.show()



if __name__ == "__main__":

    main()
