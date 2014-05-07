import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import random
import sys
import bisect
import read_eclipse as re
import eclipse_cells as ec
from time import time, clock
import csv

class Perc(object):
    def __init__(self, nx, ny, nz, r_max = 10, volume_fraction = 1.0):
        if nx >=3 and ny >=3:
            self.nx = nx
            self.ny = ny
        else:
            fail('expected nx >=3 and ny >=3, \n got \
                    nx = %d, ny= %d' %\
                    nx, ny)
        if nz == 1:
            self.nz = nz
        elif nz >=3:
            self.nz = nz
        else:
            fail('expected nz = 1 for 2d simulation or\
                    nz >=3 for 3d simulation \n \
                    got nz = %d' % nz)
        self.r_max = r_max
        self.x = {}
        self.y = {}
        self.z = {}
        self.thres_z = {}
        self.corners = {}
        self.perm = {}
        self.poro = {}
        self.volume = {}
        self.grid_values = {}
        self.fill_steps = []
        self.fill_times = []
        self.candidates = [] #keep sorted
        self.sbres = 0.2
        self.scmax = 1 - self.sbres
        self.vfrac = volume_fraction

    def add_injection(self, mass_inflow, end_time_days, \
            density):
        self.inj = self.Injection(mass_inflow, end_time_days, \
                density)

    def add_volume(self, choice):
        vol = self.vfrac * self.poro[choice] *\
                self.scmax * self.volume[choice]
        return vol

    class Injection(object):
        def __init__(self, mass_inflow, end_time_days, density):
            self.t_elapsed = 1998 * 365.25 # days
            self.q_index = 0
            self.t_end = end_time_days + self.t_elapsed
            self.rho = density
            self.massflow = mass_inflow
                                    #conversion factor
            self.mr_kg_sec = []
            self.q = []
            self.end_days = []
            for i in range(len(self.massflow)):
                self.mr_kg_sec.append(self.massflow[i] * 31.71)
                self.q.append(self.massflow[i] * 31.71 / self.rho) # -> m^3/s
                self.end_days.append((1999 + i) * 365.25)

            self.injected_mass = 0.
            self.injected_volume = 0.

            msum = 0.
            for i in range(len(self.massflow)):
                msum += self.massflow[i]
            massflow_avg = msum / float(len(self.massflow))
                
            self.max_mass = end_time_days * massflow_avg* 31.71 * 24 * 3600.
            self.max_volume = self.max_mass / self.rho

        def add_time(self, t_add):
            self.t_elapsed += t_add
            return 0
        
        def add_mass(self, vol_add):
            self.injected_volume += vol_add
            mass_add = self.rho * vol_add
            self.injected_mass += mass_add
            time_taken = vol_add / (self.q[self.q_index] * 24 * 3600) 
            # add time in days ^^^
            time_taken_1 = mass_add / (self.mr_kg_sec[self.q_index] * 24 * 3600)
            self.add_time(time_taken)
            if self.get_elapsed_time() > self.end_days[self.q_index] and \
                    self.q_index <= len(self.end_days) -1:
                self.increment_q_index()
            return 0

        def increment_q_index(self):
            self.q_index += 1
            return 0

        def get_elapsed_time(self):
            return self.t_elapsed
        
        def get_max_mass(self):
            return self.max_mass

        def get_injected_mass(self):
            return self.injected_mass

        def get_injected_volume(self):
            return self.injected_volume

        def get_density(self):
            return self.rho

        def get_mass_inflow(self):
            return self.massflow

        def get_end_time(self):
            return self.t_end

        def end_reached(self):
            if self.t_elapsed > self.t_end:
                return True
            else:
                return False

    class Corner(object):
        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z

        def get_x(self):
            return self.x

        def get_y(self):
            return self.y

        def get_z(self):
            return self.z

    def get_grid_value(self, key):
        """ returns the value for a given cell
        creates this value if it doesn't exist.
        """
        #if key not in self.grid_values:
            #self.set_grid_value(key)
        return self.grid_values[key]

    def set_grid_value(self, key, val = 'random'):
        """ sets grid value, sets value as not filled
        """
        if val == 'random':
            self.grid_values[key] = random.randint(1, self.r_max)
        else:
            self.grid_values[key] = val

    def mark_filled(self, key, time = '1'):
        """ marks grid values as filled if they are 
        within the bounds of the grid size
        """
        assert 0 <= key[0] < self.nx, \
                'i coordinate out of range(%d vs %d)' % \
                (key[0], self.nx)
        assert 0 <= key[1] < self.ny, \
                'j coordinate out of range(%d vs %d)' % \
                (key[1], self.ny)
        if self.nz == 1:
            assert key[2] == 0, 'k must equal zero'
        else:
            assert 0 <= key[2] < self.nz,\
                    'k coordinate out of range (%d vs %d)' % \
                    (key[2], self.nz)
        self.fill_steps.append(key)
        self.fill_times.append(time)

    def find_new_candidates(self):
        """ grabs neighbor cell values, inserts them into sorted list
        """
        key = self.fill_steps[-1]
        new_can = self.get_neighbor_candidates(key)
        for can in new_can:
            bisect.insort_left(self.candidates, (self.grid_values[can], can))
        return self.candidates

    def get_neighbor_candidates(self, key):
        """ checks neighbor candidates, ignores if already in list
        """
        neighbors = self.get_neighbor_keys(key)
        candidates = []
        for key in neighbors:
            if key not in self.fill_steps:
                candidates.append(key)
        return candidates

    def get_neighbor_keys(self, key):
        """ Checks six sides of neighbors for 3d case
        Checks four sides of neighbors for the 2d case
        """
        keys = []
        keys.append((key[0] - 1, key[1], key[2]))
        keys.append((key[0] + 1, key[1], key[2]))
        keys.append((key[0], key[1] - 1, key[2]))
        keys.append((key[0], key[1] + 1, key[2]))
        if self.nz != 1:
            keys.append((key[0], key[1], key[2] - 1))
            keys.append((key[0], key[1], key[2] + 1))
        return keys

    def end_criterion(self, end_type = 'boundary'):
        if end_type == 'boundary':
            if self.choice[0] in (0, self.nx-1) \
                    or self.choice[1] in (0, self.ny-1):
                print "x-y Boundary hit " 
                return True
            elif self.nz != 1 and self.choice[2] in (0, self.nz-1):
                return True
            else:
                return False
        elif end_type == 'injection':
            end_time = self.inj.get_end_time()
            elapsed = self.inj.get_elapsed_time()
            if elapsed > end_time:
                print "end criterion"
                print "time elapsed: " + str(elapsed)
                print " end time:    " + str(end_time)
                return True
            elif self.end_criterion(end_type = 'boundary'):
                return True
            else:
                return False

    def run_simulation(self, injection = False):
        """ fills grid. If no initial value is specified, picks
        i, j, k == nx/2, ny/2, nz/2
        """
        if injection == True:
            end_type = 'injection'
        else:
            end_type = 'boundary'
        print "PERCOLATING........"
        step_count = 0
        while True:
            step_count +=1
            self.candidates = self.find_new_candidates()
            assert self.candidates, 'no fillable cells found'
            self.choice = self.percolate()
            time = step_count
            if injection == True:
                volume_filled = self.add_volume(self.choice)
                self.inj.add_mass(volume_filled)
                time = self.inj.get_elapsed_time()

            self.mark_filled(self.choice, time = time)
            if self.end_criterion(end_type = end_type):
                print "Number of Cells filled: " + \
                        str(len(self.fill_steps))
                print "mass in system        : " + \
                        str(self.inj.get_injected_mass())
                print "maximum mass          : " + \
                        str(self.inj.get_max_mass())
                break
        return 0

    def percolate(self):
        choice = self.candidates[0][1]
        #print choice, '{:.3e}'.format(self.grid_values[choice]),\
                #" runner up -> ", self.candidates[1][1], \
                #'{:.3e}'.format(self.grid_values[self.candidates[1][1]]),\
                #" end", '{:.3e}'.format(self.grid_values[self.candidates[-1][1]])
                
        self.candidates.remove(self.candidates[0])
        return choice

    def make_uniform_grid(self):
        print "making uniform grid"
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    key = (i, j, k)
                    self.set_grid_value(key)
                    self.x[key] = i
                    self.y[key] = j
                    self.z[key] = k
        if len(self.fill_steps) == 0:
            init_key = (self.nx/2, self.ny/2, self.nz/2)
            self.mark_filled(init_key)
        print "grid with: (nx, ny, nz) = ", \
                (self.nx, self.ny, self.nz), " made!"
        return 0

    def make_sleipner_grid(self, vol_dict, xyz_dict, poroperm_dict):
        """ sets :
            self.x
            self.y
            self.z
            self.poro
            self.perm
        """
        t0 = clock()
        print "making Sleipner grid"
        self.nx = 65
        self.ny = 119
        self.nz = 43
        base_elev = xyz_dict[(32, 77, 34)][2]
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    key = (i, j, k)
                    vol = vol_dict[key]
                    x = xyz_dict[key][0]
                    y = xyz_dict[key][1]
                    z = xyz_dict[key][2]
                    poro = poroperm_dict[key][0]
                    perm = poroperm_dict[key][1]
                    self.x[key] = x
                    self.y[key] = y
                    self.z[key] = z
                    self.thres_z[key] = base_elev - z
                    if j <=49:
                        boost = 0.0
                        self.z[key] += boost
                        self.thres_z[key] += boost
                    self.volume[key] = vol
                    self.poro[key] = poro
                    self.perm[key] = perm
                    #if perm > 0.1:
                        #self.perm[key] = 2000.
                    val = self.perc_threshold(key) + 1. * pow(10,5.)
                    #if i == 32 and j == 77:
                        #print '{:d}, {:.3e}, {:.3e}'.format(k, perm, val)
                    self.set_grid_value(key, val = val)
        if len(self.fill_steps) == 0:
            init_key = (32, 77, 34)
            self.mark_filled(init_key, time = 1998. * 365.25 )
        print "grid with: (nx, ny, nz) = ", \
                (self.nx, self.ny, self.nz), " made in "
        print clock() - t0, " seconds"
        return 0

    def contour_topo(self):
        fig = plt.figure(figsize = (9., 12.))
        ax = fig.add_subplot(111)
        x = []
        y = []
        elev = []
        for i in range(65):
            b2 = []
            b3 = []
            blank = []
            #if i >= 35 and i < 50:
            for j in range(119):
                        #if j >= 45 and j < 75:
                b2.append(self.x[(i, j, 2)])
                b3.append(self.y[(i, j, 2)])
                blank.append(self.z[(i, j, 2)])
            elev.append(blank)
            x.append(b2)
            y.append(b3)
        xp = np.asarray(x)
        yp = np.asarray(y)
        elp = np.asarray(elev)
        N = 10
        c = ax.contourf(xp, yp, elp, N)
        cb = plt.colorbar(c, format='%.2f')
        cb.set_ticks(np.linspace(np.amin(elp), np.amax(elp), N))
        cb.set_label('elev [m]')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        #plt.savefig('topo.png')

    def perc_threshold(self, key):
        # TODO
        # ioannidis et al 1996. 
        c = 0.186
        c = pow(10.,8.)
        sigma = 0.045
        #c = 1.
        #sigma = 1.
        pcd = c * sigma * \
                pow(self.perm[key] / self.poro[key], -1/2.)
        rho_b = 1019.
        g = 9.81
        delta_rho = rho_b - self.inj.get_density() 
        pgrav = delta_rho * g * (self.thres_z[key])
        if key[0] == 32 and key[1] == 77:
            a = 1
            print "k, pcd, pgrav", "perm"
            print '{:d}, {:3e}, {:3e}, {:3e}'.format(key[2], pcd, \
                    pgrav, pcd + pgrav)
        return pcd + pgrav

    def get_time_index_gravseg(self):
        time_days = 0.
        n = 0
        time_indices = []
        for i in range(1, len(self.fill_steps)):
            key = self.fill_steps[i]
            key0 = self.fill_steps[i-1]
            if key[2] == 2 and key0[2] != 2:
                time_indices.append(i)
        n = time_indices[0]
        time_days = self.fill_times[n] 
        return n, time_days

    def get_plan_year_indices(self, years):
        yr_indices = []
        for year in years:
            yr_days = (year) * 365.25 
            for n in range(0, len(self.fill_times)):
                yr_ind = 0
                if n > 0 and \
                    self.fill_times[n] > yr_days and \
                    self.fill_times[n-1] < yr_days:
                    yr_ind = n
                    yr_indices.append(yr_ind)
        return yr_indices

    def plot_sleipner_thick_contact(self, years, gwc = False, sim_title = ''):
        if gwc == True:
            tc_str = 'contact'
        else:
            tc_str = 'thickness'
        yr_indices = self.get_plan_year_indices(years)
        size = 14
        font = {'size' : size}
        matplotlib.rc('font', **font)
        fig = plt.figure(figsize=(10.0, 2.5), dpi = 960)
        middle = len(years) * 10
        pos = 100 + middle
        for n in range(len(yr_indices)):
            pos +=1
            ax = fig.add_subplot(pos)
            xf = []
            yf = []
            kf = []
            for i in range(self.nx):
                tempx = []
                tempy = []
                tempk = []
                for j in range(self.ny):
                    x = self.x[(i, j, 0)]
                    y = self.y[(i, j, 0)]
                    tn = yr_indices[n]
                    thick, contact = self.get_thick_contact(i, j, tn)
                    tempx.append(x)
                    tempy.append(y)
                    if gwc == True:
                        tempk.append(contact)
                    else:
                        tempk.append(thick)
                xf.append(tempx)
                yf.append(tempy)
                kf.append(tempk)
            xp = np.asarray(xf)
            yp = np.asarray(yf)
            kp = np.asarray(kf)
            N = 10
            contour_label = False
            ax_label = False
            c = ax.contourf(xp, yp, kp, N)
            if n == len(years) - 1:
                fig.subplots_adjust(right=0.84)
                cb_axes = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                cb = fig.colorbar(c, cax = cb_axes, format = '%.2f')
                cb.set_ticks(np.linspace(np.amin(kp), np.amax(kp), N))
                cb.set_label(tc_str + ': [m]')
            if n != 0:
                ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_title(str(years[n]))
            ax.axis([0, 3000, 0, 6000])
        plt.savefig(sim_title + '_' +  tc_str + '.pdf', fmt = 'pdf')
        plt.clf()
        return 0

    def plot_sleipner_plume(self, years, sim_title = 'sleipner_perc'):
        yr_indices = self.get_plan_year_indices(years)
        size = 14
        font = {'size' : size}
        matplotlib.rc('font', **font)
        fig = plt.figure(figsize=(16.0, 5), dpi=960)
        middle = len(years) * 10
        pos = 100 + middle
        for i in range(len(yr_indices)):
            pos +=1
            ax = fig.add_subplot(pos)
            xf = []
            yf = []
            kf = []
            for n in range(yr_indices[i]):
                key = self.fill_steps[n]
                #if key[0] >= 35 and key[0] < 50:
                    #if key[1] >= 45 and key[1] < 75: 
                xf.append(self.x[key])
                yf.append(self.y[key])
                kf.append(key[2])
                if 50 == key[1]:
                    key1 = (key[0], key[1]-1, key[2])
            xp = np.asarray(xf)
            yp = np.asarray(yf)
            sc = ax.scatter(xp, yp, s=20, c=kf)
            ax.set_title(str(years[i]))
            ax.axis([0, 3000, 0, 6000])
            ax.xaxis.set_ticks(np.arange(0, 3000, 1500))
            if i != 0:
                ax.set_yticklabels([])
            #elif i == 5:
                #cb_axes = self.fig.add_axes([0.85, 0.15, 0.05, 0.7])
                #fig.colorbar(sc, cax = cb_axes)
        plt.savefig(sim_title + '_plume.pdf', fmt = 'pdf')
        plt.clf()
        return 0

    def plot_sleipner_cross_section(self, years, sec_index = 32):
        yr_indices = self.get_plan_year_indices(years)
        size = 14
        font = {'size' : size}
        matplotlib.rc('font', **font)
        fig = plt.figure(figsize=(16.0, 5))
        pos = 150
        top = []
        bot = []
        ybound = []
        for key in self.x.keys():
            if key[0] == sec_index:
                t, b = self.get_boundary_zs(key[0], key[1])
                top.append(t)
                bot.append(b)
                ybound.append(self.y[key])
        for i in range(len(yr_indices)):
            pos +=1
            ax = fig.add_subplot(pos)
            yf = []
            zf = []
            for n in range(yr_indices[i]):
                key = self.fill_steps[n]
                if key[0] == sec_index:
                    yf.append(self.y[key])
                    zf.append(self.z[key])
            yp = np.asarray(yf)
            zp = np.asarray(zf)
            tp = np.asarray(top)
            bp = np.asarray(bot)
            yb = np.asarray(ybound)
            tl = ax.scatter(yb, tp, s=5, c='r')
            bl = ax.scatter(yb, bp, s=5, c='g')
            sc = ax.scatter(yp, zp, s=10)
            ax.set_title(str(years[i]))
            ax.axis([0, 6000, -815, -800])
            ax.xaxis.set_ticks(np.arange(0, 6000, 1500))
            if i != 0:
                ax.set_yticklabels([])
        plt.savefig('sleipner_cross_section.png')
        plt.clf()
        return 0
    
    def contour_top_boundary(self):
        top = []
        x = []
        y = []
        top = []
        for i in range(self.nx):
            xinter = []
            yinter = []
            tinter = []
            for j in range(self.ny):
                key = (i, j, 2)
                xinter.append(self.x[key])
                yinter.append(self.y[key])
                tinter.append(self.z[key])
            x.append(xinter)
            y.append(yinter)
            top.append(tinter)
        xp = np.asarray(x)
        yp = np.asarray(y)
        tp = np.asarray(top)
        fig = plt.figure(figsize=(8.5,11))
        ax = fig.add_subplot(111)
        N = 50
        cs_val = ax.contour(xp, yp, tp, N)
        cb_val = plt.colorbar(cs_val, shrink = 0.8,\
                extend='both')
        cb_val.set_label('Top Boundary [z]')
        fig.savefig('top_boundary.png', bbox_inches='tight', format='png')
        return 0


    def make_scatter_plan_t0_tn(self, t0, tn):
        n = tn-t0
        x = np.zeros(n)
        y = np.zeros(n)
        for n in range(t0, tn):
            key = self.fill_steps[n]
            x[n] = self.x[key]
            y[n] = self.y[key]
        return x, y
       
    def plot_2d(self, uniform_grid = True):
        print "PLOTTING..........."
        f = plt.figure()
        ax = f.add_subplot(111)

        # make base grid of cells
        if uniform_grid == True:
            pts = []
            xs = []
            ys = []
            for i in [0, self.nx-1]:
                for j in [0, self.ny-1]:
                    key = (i, j, 0)
                    xs.append(self.x[key])
                    ys.append(self.y[key])
            xp = np.asarray(xs)
            yp = np.asarray(ys)
            ax.scatter(xp, yp, s=30, c='w', marker='s')

        # go through steps and figure out times
        xf = []
        yf = []
        tf = []
        tmin = self.fill_times[0]
        tmax = self.fill_times[-1]
        for i in range(0, len(self.fill_steps)):
            key = self.fill_steps[i]
            xf.append(self.x[key])
            yf.append(self.y[key])
            tf.append(self.fill_times[i])

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        xfp = np.asarray(xf)
        yfp = np.asarray(yf)
        cm = plt.get_cmap('bone_r')
        sc = ax.scatter(xfp, yfp, c = tf, vmin=tmin, vmax=tmax, s = 300, cmap=cm)
        plt.colorbar(sc)
        plt.savefig('sleipner_2d.png')
        #plt.show()

    def get_boundary_zs(self, i, j):
        for k in range(1, self.nz):
            key0 = (i, j, k-1)
            key1 = (i, j, k)
            if self.perm[key0] < 1. and self.perm[key1] > 1.:
                ztop = self.z[key1]
            elif self.perm[key0] > 1. and self.perm[key1] < 1.:
                zbot = self.z[key0]
        return ztop, zbot

    def get_thick_contact(self, i, j, time_index):
        column = []
        for key in self.fill_steps[:time_index]:
            if key[0] == i and key[1] == j:
                column.append(self.z[key])
        column.sort()
        if len(column) == 0:
            thick = 0.
            contact = -812.
        else:
            thick = column[-1] - column[0] + 0.52
            contact = column[0]
            if contact < -812.:
                contact = -812.
        return thick, contact

    def plot_3d(self, uniform_grid = True):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        if uniform_grid == True:
            pts = []
            xs = []
            ys = []
            zs = []
            for i in [0, self.nx-1]:
                for j in [0, self.ny-1]:
                    for k in [0, self.nz-1]:
                        key = (i, j, k)
                        xs.append(self.x[key])
                        ys.append(self.y[key])
                        zs.append(self.z[key])
            xp = np.asarray(xs)
            yp = np.asarray(ys)
            zp = np.asarray(zs)
            ax.scatter(xp, yp, zp, s=30, c='w', marker='s')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        xf = [] 
        yf = []
        zf = []
        tf = []
        tmin = self.fill_times[0]
        tmax = self.fill_times[-1]
        for i in range(0, len(self.fill_steps)):
            key = self.fill_steps[i]
            xf.append(self.x[key])
            yf.append(self.y[key])
            zf.append(self.z[key])
            tf.append(self.fill_times[i])
        xfp = np.asarray(xf)
        yfp = np.asarray(yf)
        zfp = np.asarray(zf)
        cm = plt.get_cmap('bone_r')
        sc = ax.scatter(xfp, yfp, zfp, \
                c = tf, vmin=tmin, vmax=tmax, s = 300, cmap=cm)
        plt.colorbar(sc)
        #plt.show()
        return 0

    def make_sleipner_csv(self):
        e_cells, nx, ny, nz = re.read_eclipse()
        f = open('sl_data.csv', 'w')
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    key = (i, j, k)
                    ind = self.e_cell_index(i, j, k)
                    oc = e_cells[ind].getCorners()
                    corners = []
                    for c in oc:
                        x, y = c.getXY()
                        # FLIPPING ALL ZS IN THIS. 
                        z = - c.getZ()
                        nc = self.Corner(x, y, z)
                        corners.append(nc)
                    self.corners[key] = corners
                    x = self.get_x_centroid(corners)
                    y = self.get_y_centroid(corners)
                    z = self.get_z_centroid(corners)
                    poro = e_cells[ind].getPorosity()
                    perm = e_cells[ind].getXPermeability()
                    volume = self.get_volume(x, y, z, corners)
                    vol_s = str(volume)
                    x_s = str(x)
                    y_s = str(y)
                    z_s = str(z)
                    poro_s = str(poro)
                    perm_s = str(perm)
                    f.write(', '.join([str(i), str(j), str(k), \
                            vol_s, x_s, y_s, z_s, poro_s, perm_s]))
                    f.write('\n')
        f.close()
        return 0

    def read_sleipner_csv(self):
        with open('sl_data.csv', 'rb') as csvfile:
            vol_dict = {}
            xyz_dict = {}
            poroperm_dict = {}
            rd = csv.reader(csvfile, delimiter = ',')
            for row in rd:
                key = (int(row[0]), int(row[1]), int(row[2]))
                vol_dict[key] = float(row[3])
                xyz_dict[key] = (float(row[4]), float(row[5]), float(row[6]))
                poroperm_dict[key] = (float(row[7]), float(row[8]))
        csvfile.close()
        return vol_dict, xyz_dict, poroperm_dict

    def e_cell_index(self, i, j, k):
        nx = 65
        ny = 119
        return i + nx * j + nx * ny * k

    def get_x_centroid(self, corners):
        count = 0.
        sum_c = 0.
        for c in corners:
            count += 1.
            sum_c += c.get_x()
        return sum_c / count

    def get_y_centroid(self, corners):
        count = 0.
        sum_c = 0.
        for c in corners:
            count += 1.
            sum_c += c.get_y()
        return sum_c / count

    def get_z_centroid(self, corners):
        count = 0.
        sum_c = 0.
        for c in corners:
            count += 1.
            sum_c += c.get_z()
        return sum_c / count
  
    def get_dx(self, eleme, direc):
        """ returns the length of a grid cell in a particular direction.
        dir is either 1, 2 or 3 for x, y and z directions.
        i, j and k are the indices
        """
        if direc == 1 :
            corners = self.corners[eleme]
            dx = corners[0].get_x() - corners[1].get_x()
            return dx
        elif direc == 2 :
            corners = self.corners[eleme]
            dy = corners[0].get_y() - corners[2].get_y()
            return dy
        elif direc == 3 :
            z1 = abs(e_cells[self.e_cell_index(i,j,k)].getTopZ() - \
                    e_cells[self.e_cell_index(i,j,k)].getBottomZ())
            return z1
        else:
            raise Exception("Invalid direction, \n" + \
                    " Please specify 1, 2 or 3.\n")

    def get_volume(self, x, y, z, corners):
        """ uses the equation for volume of an orientable polyhedron
            V = 1/3 \sum_i x_i \dot n^hat_i A_i
        """
        face_map = ['west', 'south', 'east', 'north', 'bot', 'top']

        v_sum = 0.0
        for face in face_map:
            a = self.get_area(corners, face)
            centroid = self.get_face_center(x, y, z, corners, face)
            cent = np.asarray(centroid)
            vec = self.get_normal_vector(x, y, z, corners, face)
            v_sum += np.dot(cent, vec) * a

        vol = 1./3. * v_sum
        return vol

    def get_area(self, corners, face):
        """ returns the area of a cell face, east, west, etc
        """ 
        if face == 'west':
            x1 = corners[2].get_y()
            x2 = corners[0].get_y()
            y1 = corners[2].get_z()
            y2 = corners[0].get_z()
            y3 = corners[6].get_z()
            y4 = corners[4].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'south':
            x1 = corners[2].get_x()
            x2 = corners[3].get_x()
            y1 = corners[2].get_z()
            y2 = corners[3].get_z()
            y3 = corners[6].get_z()
            y4 = corners[7].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'east':
            x1 = corners[3].get_y()
            x2 = corners[1].get_y()
            y1 = corners[3].get_z()
            y2 = corners[1].get_z()
            y3 = corners[7].get_z()
            y4 = corners[5].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'north':
            x1 = corners[0].get_x()
            x2 = corners[1].get_x()
            y1 = corners[0].get_z()
            y2 = corners[1].get_z()
            y3 = corners[4].get_z()
            y4 = corners[5].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'bot':
            nc = [corners[6], corners[7], corners[4], corners[5]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0],2.) + pow(c[1],2.) + 1)
            x1 = corners[2].get_x()
            x2 = corners[3].get_x()
            y1 = corners[2].get_y()
            y2 = corners[0].get_y()
            area = mag * ((x2 * y2 - x1 * y2) - (x2 * y1 - x1 * y1))
        elif face == 'top':
            nc = [corners[2], corners[3], corners[0], corners[1]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0],2.) + pow(c[1],2.) + 1)
            x1 = corners[6].get_x()
            x2 = corners[7].get_x()
            y1 = corners[6].get_y()
            y2 = corners[4].get_y()
            area = mag * ((x2 * y2 - x1 * y2) - (x2 * y1 - x1 * y1))
        else:
            raise Exception("Invalid Face, please specify" +  \
                    "one of the six faces in face_map \n\n")

        return area

    def get_face_center(self, xc, yc, zc, corners, face):
        """ center vector location relative to polyhedron center
        """
        if face == 'west':
            nc = [corners[0], corners[2], corners[4], corners[6]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'south':
            nc = [corners[2], corners[3], corners[6], corners[7]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
            a = 2
        elif face == 'east':
            nc = [corners[3], corners[1], corners[7], corners[5]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'north':
            nc = [corners[0], corners[1], corners[4], corners[5]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'bot':
            nc = [corners[6], corners[7], corners[4], corners[5]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'top':
            nc = [corners[2], corners[3], corners[0], corners[1]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        else:
            raise Exception("Invalid Face, please specify" +  \
                    "one of the six faces in face_map \n\n")

        vec = [xf - xc, yf - yc, zf - zc]
        return vec

    def get_normal_vector(self, x, y, z, corners, face):
        """ gets normal vector of face
        """
        if face == 'west':
            vec = [-1., 0., 0.]
        elif face == 'south':
            vec = [0., -1., 0.]
        elif face == 'east':
            vec = [1., 0., 0.]
        elif face == 'north':
            vec = [0., 1., 0.]
        elif face == 'bot':
            nc = [corners[6], corners[7], corners[4], corners[5]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0], 2.) + pow(c[1],2.) + 1)
            vec = [c[0]/mag, c[1]/mag, -1./mag]
        elif face == 'top':
            nc = [corners[2], corners[3], corners[0], corners[1]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0], 2.) + pow(c[1],2.) + 1)
            vec = [-c[0]/mag, -c[1]/mag, 1./mag]
        else:
            raise Exception("Invalid Face, please specify" +  \
                    "one of the six faces in face_map \n\n")
        return vec

    def fit_plane(self, corners):
        """ takes four corner points and fits a plane least squares to them
            returns in form z = c[0] x + c[1] y + c[2]
        """
        x = []
        y = []
        z = []
        for c in corners:
            x.append(c.get_x())
            y.append(c.get_y())
            z.append(c.get_z())
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        A = np.column_stack((x, y, np.ones(x.size)))
        c, resid, rank, sigma = np.linalg.lstsq(A, z)
        return c, resid, rank, sigma

    def get_area_side(self, x1, x2, y1, y2, y3, y4):
        h = x2 - x1
        b1 = y4 - y2
        b2 = y3 - y1
        return 0.5 * h * (b1 + b2)

def fail(msg):
    '''print error and quit'''
    print >> sys.stderr, msg
    sys.exit(1)
