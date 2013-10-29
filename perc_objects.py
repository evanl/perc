import numpy as np
import matplotlib.pyplot as plt
import random

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
        elif self.nz >=3:
            self.nz = nz
        else:
            fail('expected nz = 1 for 2d simulation or\
                    nz >=3 for 3d simulation \n \
                    got nz = %d' % nz)
        self.r_max = r_max
        self.grid = {}
        self.filled = {}
        self.fill_steps =[]

    def get_grid_value(self, key):
        """ returns the value for a given cell
        creates this value if it doesn't exist.
        """
        if key not in self.grid:
            self.set_grid_value(key)
        return self.grid[key]

    def set_grid_value(self, key):
        """ sets grid value, sets value as not filled
        """
        self.grid[key] = random.randint(1, self.r_max)

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
            assert 0 <= key[2] < nz,\
                    'k coordinate out of range (%d vs %d)' % \
                    (key[2], self.nz)
        self.filled[key] = True

    def find_candidates(self):
        candidates = {}
        for k in self.filled.keys():
            new_can = self.get_neighbor_candidates(k)
            candidates.update(new_can)
        return candidates

    def get_neighbor_candidates(self, key):
        """ checks neighbor candidates based on minimum value
        """
        neighbors = self.get_neighbor_keys(key)
        candidates = {}
        min_val = self.r_max
        for key in neighbors:
            self.get_grid_value(key)
            if key in self.filled.keys():
                continue
            if self.grid[key] == min_val:
                candidates[key] = self.grid[key]
            elif self.grid[key] < min_val:
                candidates = {}
                min_val = self.grid[key]
                candidates[key] = self.grid[key]
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


    def run_simulation(self):
        """ fills grid. If no initial value is specified, picks
        i, j, k == nx/2, ny/2, nz/2
        """
        if len(self.filled) == 0:
            init_key = (self.nx/2, self.ny/2, self.nz/2)
            self.set_grid_value(init_key)
            self.mark_filled(init_key)
            self.fill_steps.append(init_key)

        step_count = 0
        while True:
            step_count +=1
            candidates = self.find_candidates()
            assert candidates, 'no fillable cells found'
            self.choice = random.choice(candidates.keys())
            self.mark_filled(self.choice)
            self.fill_steps.append(self.choice)
            if self.end_criterion():
                break
        return 0

    def end_criterion(self, end_type = 'boundary'):
        if self.choice[0] in (0, self.nx-1) \
                or self.choice[1] in (0, self.ny-1):
            return True
        elif self.nz != 1 and self.choice[3] in (0, self.nz-1):
            return True
        else:
            return False

    def plot_2d(self):
        f = plt.figure()
        ax = f.add_subplot(111)

        # make base grid of cells
        pts = []
        xs = []
        ys = []
        for x in range(0,self.nx):
            for y in range(0, self.ny):
                xs.append(x)
                ys.append(y)
        xp = np.asarray(xs)
        yp = np.asarray(ys)
        cells = np.asarray
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

        xfp = np.asarray(xf)
        yfp = np.asarray(yf)
        cm = plt.get_cmap('bone_r')
        sc = ax.scatter(xfp, yfp, c = tf, vmin=tmin, vmax=tmax, s = 300, cmap=cm)
        plt.colorbar(sc)
        plt.show()

def fail(msg):
    '''print error and quit'''
    print >> sys.stderr, msg
    sys.exit(1)
