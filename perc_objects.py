import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import sys
import bisect

class Perc(object):
    def __init__(self, nx, ny, nz, r_max = 10):
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
        self.grid_values = {}
        self.filled = {}
        self.fill_steps = []
        self.candidates = [] #keep sorted

    def get_grid_value(self, key):
        """ returns the value for a given cell
        creates this value if it doesn't exist.
        """
        #if key not in self.grid_values:
            #self.set_grid_value(key)
        return self.grid_values[key]

    def set_grid_value(self, key):
        """ sets grid value, sets value as not filled
        """
        self.grid_values[key] = random.randint(1, self.r_max)

    def mark_filled(self, key):
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
        self.filled[key] = True

    def find_new_candidates(self):
        """ grabs neighbor cell values, inserts them into sorted list
        """
        k = self.fill_steps[-1]
        new_can = self.get_neighbor_candidates(k)
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
        if self.choice[0] in (0, self.nx-1) \
                or self.choice[1] in (0, self.ny-1):
            return True
        elif self.nz != 1 and self.choice[2] in (0, self.nz-1):
            return True
        else:
            return False

    def run_simulation(self):
        """ fills grid. If no initial value is specified, picks
        i, j, k == nx/2, ny/2, nz/2
        """
        # populating grid values 
        self.make_uniform_grid()
        print "PERCOLATING........"
        step_count = 0
        while True:
            step_count +=1
            self.candidates = self.find_new_candidates()
            assert self.candidates, 'no fillable cells found'
            self.choice = self.percolate()
            self.mark_filled(self.choice)
            self.mark_filled(self.choice)
            self.fill_steps.append(self.choice)
            if self.end_criterion():
                break
        return 0

    def percolate(self):
        choice = self.candidates[0][1]
        self.candidates.remove(self.candidates[0])
        return choice

    def make_uniform_grid(self):
        print "making uniform grid"
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    self.set_grid_value((i, j, k))
        if len(self.filled) == 0:
            init_key = (self.nx/2, self.ny/2, self.nz/2)
            self.mark_filled(init_key)
            self.fill_steps.append(init_key)
        print "grid with: (nx, ny, nz) = ", \
                (self.nx, self.ny, self.nz), " made!"
        return 0

    def plot_2d(self, uniform_grid = True):
        print "PLOTTING..........."
        f = plt.figure()
        ax = f.add_subplot(111)

        # make base grid of cells
        if uniform_grid == True:
            pts = []
            xs = []
            ys = []
            for x in [0, self.nx-1]:
                for y in [0, self.ny-1]:
                    xs.append(x)
                    ys.append(y)
            xp = np.asarray(xs)
            yp = np.asarray(ys)
            ax.scatter(xp, yp, s=30, c='w', marker='s')

        # go through steps and figure out times
        xf = []
        yf = []
        tf = []
        tmin = 0
        tmax = len(self.fill_steps)
        for i in range(tmin, tmax):
            xf.append(self.fill_steps[i][0])
            yf.append(self.fill_steps[i][1])
            tf.append(i)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        xfp = np.asarray(xf)
        yfp = np.asarray(yf)
        cm = plt.get_cmap('bone_r')
        sc = ax.scatter(xfp, yfp, c = tf, vmin=tmin, vmax=tmax, s = 300, cmap=cm)
        plt.colorbar(sc)
        plt.show()
    def plot_3d(self, uniform_grid = True):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')

        if uniform_grid == True:
            pts = []
            xs = []
            ys = []
            zs = []
            for x in [0, self.nx-1]:
                for y in [0, self.ny-1]:
                    for z in [0, self.nz-1]:
                        xs.append(x)
                        ys.append(y)
                        zs.append(z)
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
        tmin = 0
        tmax = len(self.fill_steps)
        for i in range(tmin, tmax):
            xf.append(self.fill_steps[i][0])
            yf.append(self.fill_steps[i][1])
            zf.append(self.fill_steps[i][2])
            tf.append(i)
        xfp = np.asarray(xf)
        yfp = np.asarray(yf)
        zfp = np.asarray(zf)
        cm = plt.get_cmap('bone_r')
        sc = ax.scatter(xfp, yfp, zfp, \
                c = tf, vmin=tmin, vmax=tmax, s = 300, cmap=cm)
        plt.colorbar(sc)
        plt.show()

def fail(msg):
    '''print error and quit'''
    print >> sys.stderr, msg
    sys.exit(1)
