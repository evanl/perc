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
    #nx = int(arguments[0])
    #ny = int(arguments[0])
    #nz = int(arguments[0])
    #nz = 1
    r_max = 1000000
    nx = 65
    ny = 119
    nz = 43
    mass_inflow = 0.1418
    mass_inflow = [0.0198, 0.0405, 0.0437, 0.0540, 0.0740, 0.1030, \
                  0.1390, 0.1830, 0.2370, 0.2960, 0.370, 0.421, 0.505]

    vfrac = 0.3
    density = 308.
    sim_title = 'pvar_42_030'

    sim_years = int(arguments[0])
    perc = po.Perc(nx, ny, nz, r_max = r_max, volume_fraction = vfrac)
    # TODO Run 501 years!
    print "sim_years: " + str(sim_years)
    if sim_years == 11:
        years = [1999, 2001, 2002, 2004, 2006]
    elif sim_years == 13:
        years = [1999, 2001, 2002, 2004, 2006, 2008]
    elif sim_years == 51:
        years = [2000, 2010, 2020, 2030, 2040]
    elif sim_years == 151:
        years = [2000, 2030, 2060, 2090, 2120]
    elif sim_years == 501:
        years = [2000, 2100, 2200, 2300, 2400]
    elif sim_years == 1001:
        years = [2000, 2200, 2400, 2600, 2800]
    sim_time = sim_years * 365.25 
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
    yr_indices = perc.get_plan_year_indices(years)
    print yr_indices
    print "mass balance?"
    print perc.inj.get_injected_mass() - \
            perc.inj.get_max_mass()
    perc.plot_sleipner_plume(years, sim_title = sim_title)
    #perc.plot_sleipner_cross_section(years, sec_index = 49)
    perc.plot_sleipner_thick_contact(years, gwc = False, sim_title = sim_title)
    perc.plot_sleipner_thick_contact(years, gwc = True, sim_title = sim_title)
    #perc.contour_top_boundary()
    #perc.contour_topo()
    if nz == 1:
        perc.plot_2d()
    else:
        perc.plot_3d()
