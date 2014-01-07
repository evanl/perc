import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import sys
import bisect
import read_eclipse as re
import eclipse_cells as ec
from time import time, clock
import csv

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
        self.x = {}
        self.y = {}
        self.z = {}
        self.corners = {}
        self.perm = {}
        self.poro = {}
        self.volume = {}
        self.grid_values = {}
        self.fill_steps = []
        self.candidates = [] #keep sorted

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
        self.fill_steps.append(key)

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
        print "PERCOLATING........"
        step_count = 0
        while True:
            step_count +=1
            self.candidates = self.find_new_candidates()
            assert self.candidates, 'no fillable cells found'
            self.choice = self.percolate()
            self.mark_filled(self.choice)
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
        t0 = clock()
        print "making Sleipner grid"
        self.nx = 65
        self.ny = 119
        self.nz = 43
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
                    self.volume[key] = vol
                    self.poro[key] = poro
                    self.perm[key] = perm
                    #TODO
                    # USE REAL THRESHOLD INSTEAD OF RANDOM
                    self.set_grid_value(key)
        if len(self.fill_steps) == 0:
            init_key = (self.nx/2, self.ny/2, self.nz/2)
            self.mark_filled(init_key)
        print "grid with: (nx, ny, nz) = ", \
                (self.nx, self.ny, self.nz), " made in "
        print clock() - t0, " seconds"
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
        tmin = 0
        tmax = len(self.fill_steps)
        for i in range(tmin, tmax):
            key = self.fill_steps[i]
            xf.append(self.x[key])
            yf.append(self.y[key])
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
        tmin = 0
        tmax = len(self.fill_steps)
        for i in range(tmin, tmax):
            key = self.fill_steps[i]
            xf.append(self.x[key])
            yf.append(self.y[key])
            zf.append(self.z[key])
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
