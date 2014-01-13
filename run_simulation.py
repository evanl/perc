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
    r_max = 1000000
    perc = po.Perc(nx, ny, nz, r_max = r_max)
    mass_inflow = 0.1418
    density = 700.
    sim_time = 11 * 365.25 
    perc.add_injection(mass_inflow, sim_time, density)
    injection = True
    #perc.make_uniform_grid()
    #perc.make_sleipner_csv()
    vol_dict, xyz_dict, poroperm_dict = perc.read_sleipner_csv()
    perc.make_sleipner_grid(vol_dict, xyz_dict, poroperm_dict)
    t0 = clock()
    perc.run_simulation(injection = injection)
    t1 = clock()
    print "percolation time = ", t1 - t0
    n, tseg = perc.get_time_index_gravseg()
    print "n, gravseg_days"
    print n, tseg
    yr_indices = perc.get_plan_year_indices()
    print yr_indices
    perc.plot_sleipner_plume()
    if nz == 1:
        perc.plot_2d()
    else:
        perc.plot_3d()
