import sys
import numpy as np
import matplotlib.pyplot as plt
import perc_objects as po
from time import time, clock
import sys
#import read_eclipse as re


def check_input_args(arguments):
    try:
        grid_size = int(arguments[0])
    except IndexError:
        po.fail('Expected 1 arguments, got %d' %\
                len(arguments))
    except ValueError:
        po.fail('expected int arguments, got %s' % \
                str(arguments))
    return 0

#def run_simulation(options = 'awesome'):
    #if options == 'awesome':
        #print 'awesome'
    #else:
        #po.fail('for shame')
    # sketch ideas of interface
    #e_cells, nx, ny, nz = re.read_eclipse()
    #perc_sim = PercSim(nx, ny, nz, e_cells = e_cells)
    #if steps == True:
        #growth_type = 'steps'
        #growth = 2
    #else:
        #growth_type = 'time'
        #growth = 2.
    #perc_sim.run_simulation(growth_type = sim_type, growth = growth)


if __name__ == '__main__':
    # profiling command
    # '$ python -m cProfile -s time __main__.py N >prof.out
    # get parameters from command line
    arguments = sys.argv[1:]
    check_input_args(arguments)
    nx = int(arguments[0])
    ny = int(arguments[0])
    nz = int(arguments[0])
    #nz = 1
    t0 = clock()
    r_max = 10000
    perc = po.Perc(nx, ny, nz, r_max = r_max)
    perc.run_simulation()
    t1 = clock()
    print "percolation time = ", t1 - t0
    if nz == 1:
        perc.plot_2d()
    else:
        perc.plot_3d()
