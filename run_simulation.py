import sys
import numpy as np
import matplotlib.pyplot as plt
import perc_objects as po
#import read_eclipse as re


def check_input_args(arguments):
    try:
        grid_size = int(arguments[0])
    except IndexError:
        fail('Expected 3 arguments, got %d' %\
                len(arguments))
    except ValueError:
        fail('expected int arguments, got %s' % \
                str(arguments))
    return 0

#def run_simulation(options = 'awesome'):
    #if options == 'awesome':
        #print 'awesome'
    #else:
        #fail('for shame')
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
    # get parameters from command line
    #arguments = sys.argv[1:]
    #check_input_args(arguments)
    nx = 9
    ny = 9
    nz = 1
    perc = po.Perc(nx, ny, nz)
    perc.run_simulation()
    perc.plot_2d()
