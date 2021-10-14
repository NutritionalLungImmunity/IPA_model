import numpy as np


class Diffuse():

    """DIFFUSION IN CLOSED BOUNDARY CONDITIONS"""

    def __init__(self, f, pde_factor, delta_t):
        self.f = f
        self.pde_factor = pde_factor
        self.delta_t = delta_t

    def solver(self, space, molecule, index):
        xbin = len(space)
        ybin = len(space[0])
        zbin = len(space[0][0])

        import timeit

        D = [[[0 for _ in range(zbin)] for _ in range(ybin)] for _ in range(xbin)]
        next = [[[0 for _ in range(zbin)] for _ in range(ybin)] for _ in range(xbin)]

        C_x = [0 for _ in range(xbin)]
        C_y = [0 for _ in range(ybin)]
        C_z = [0 for _ in range(zbin)]

        D_x = [0 for _ in range(xbin)]
        D_y = [0 for _ in range(ybin)]
        D_z = [0 for _ in range(zbin)]

        A  = -(3 - 2 * self.f)
        B  = self.pde_factor/self.delta_t + 2 * (3 - 2 * self.f)
        Bc = self.pde_factor/self.delta_t + (3 - 2 * self.f)
        C  = -(3 - 2 * self.f)
        E  = self.pde_factor/self.delta_t - 4 * self.f

        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    #print(space[x][y][z].molecules[molecule].values[0](index))
                    D[x][y][z] = E * space[x][y][z].molecules[molecule].values[0]
                    if y == 0 and z == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (space[x][y + 1][z].molecules[molecule].values[index] +
                                                            space[x][y][z + 1].molecules[molecule].values[index] +
                                                            2 * space[x][y][z].molecules[molecule].values[index])
                    elif y == 0 and z == zbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (space[x][y + 1][z].molecules[molecule].values[index] +
                                                            space[x][y][z - 1].molecules[molecule].values[index] +
                                                            2 * space[x][y][z].molecules[molecule].values[index])
                    elif y == ybin - 1 and z == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (space[x][y - 1][z].molecules[molecule].values[index] +
                                                            space[x][y][z + 1].molecules[molecule].values[index] +
                                                            2 * space[x][y][z].molecules[molecule].values[index])
                    elif y == ybin - 1 and z == zbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (space[x][y - 1][z].molecules[molecule].values[index] +
                                                            space[x][y][z - 1].molecules[molecule].values[index] +
                                                            2 * space[x][y][z].molecules[molecule].values[index])
                    elif y == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (space[x][y + 1][z].molecules[molecule].values[index] +
                                                            space[x][y][z - 1].molecules[molecule].values[index] +
                                                            space[x][y][z + 1].molecules[molecule].values[index] +
                                                            space[x][y][z].molecules[molecule].values[index])
                    elif y == ybin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (space[x][y - 1][z].molecules[molecule].values[index] +
                                                            space[x][y][z - 1].molecules[molecule].values[index] +
                                                            space[x][y][z + 1].molecules[molecule].values[index] +
                                                            space[x][y][z].molecules[molecule].values[index])
                    elif z == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (space[x][y - 1][z].molecules[molecule].values[index] +
                                                            space[x][y + 1][z].molecules[molecule].values[index] +
                                                            space[x][y][z + 1].molecules[molecule].values[index] +
                                                            space[x][y][z].molecules[molecule].values[index])
                    elif z == zbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (space[x][y + 1][z].molecules[molecule].values[index] +
                                                            space[x][y - 1][z].molecules[molecule].values[index] +
                                                            space[x][y][z - 1].molecules[molecule].values[index] +
                                                            space[x][y][z].molecules[molecule].values[index])
                    else:
                        D[x][y][z] = D[x][y][z] + self.f * (space[x][y - 1][z].molecules[molecule].values[index] +
                                                            space[x][y + 1][z].molecules[molecule].values[index] +
                                                            space[x][y][z - 1].molecules[molecule].values[index] +
                                                            space[x][y][z + 1].molecules[molecule].values[index])

        for y in range(ybin):
            for z in range(zbin):
                for x in range(xbin):
                    if x == 0:
                        C_x[x] = C / Bc
                    #elif x == (xbin - 1):
                    #    C_x[x] = C / (Bc - A * C_x[x - 1])
                    else:
                        C_x[x] = C / (B - A * C_x[x - 1])
                for x in range(xbin):
                    if x == 0:
                        D_x[x] = D[x][y][z] / Bc
                    elif x == xbin - 1:
                        D_x[x] = (D[x][y][z] - A * D_x[x - 1]) / (Bc - A * C_x[x - 1])
                    else:
                        D_x[x] = (D[x][y][z] - A * D_x[x - 1]) / (B - A * C_x[x - 1])
                for x in reversed(range(xbin)):
                    if x == xbin - 1:
                        next[x][y][z] = D_x[x]
                    else:
                        next[x][y][z] = D_x[x] - C_x[x] * next[x + 1][y][z]

        #tic = timeit.default_timer()
        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    D[x][y][z] = E * next[x][y][z];
                    if x == 0 and z == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x+1][y][z] + next[x][y][z+1] + 2 * next[x][y][z])
                    elif x == 0 and z == zbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x+1][y][z] + next[x][y][z-1] + 2 * next[x][y][z])
                    elif x == xbin - 1 and z == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x-1][y][z] + next[x][y][z+1] + 2 * next[x][y][z])
                    elif x == xbin - 1 and z == zbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x-1][y][z] + next[x][y][z-1] + 2 * next[x][y][z])
                    elif x == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x+1][y][z] + next[x][y][z-1] + next[x][y][z+1] + next[x][y][z])
                    elif x == xbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x-1][y][z] + next[x][y][z-1] + next[x][y][z+1] + next[x][y][z])
                    elif z == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x-1][y][z] + next[x+1][y][z] + next[x][y][z+1] + next[x][y][z])
                    elif z == zbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x-1][y][z] + next[x+1][y][z] + next[x][y][z-1] + next[x][y][z])
                    else:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x-1][y][z] + next[x+1][y][z] + next[x][y][z-1] + next[x][y][z+1]);
        #print("# ", (timeit.default_timer() - tic))

        for x in range(xbin):
            for z in range(zbin):
                for y in range(ybin):
                    if y == 0:
                        C_y[y] = C / Bc
                    #elif y == ybin:
                    #    C_y[y] = C / (Bc - A * C_y[y - 1])
                    else:
                        C_y[y] = C / (B - A * C_y[y - 1])
                for y in range(ybin):
                    if y == 0:
                        D_y[y] = D[x][y][z] / Bc
                    elif y == ybin - 1:
                        D_y[y] = (D[x][y][z] - A * D_y[y - 1]) / (Bc - A * C_y[y - 1])
                    else:
                        D_y[y] = (D[x][y][z] - A * D_y[y - 1]) / (B - A * C_y[y - 1])

                for y in reversed(range(ybin)):
                    if y == ybin - 1:
                        next[x][y][z] = D_y[y]
                    else:
                        next[x][y][z] = D_y[y] - C_y[y] * next[x][y + 1][z]

        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    D[x][y][z] = E * next[x][y][z]
                    if y == 0 and x == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x][y+1][z] + next[x+1][y][z] + 2 * next[x][y][z])
                    elif y == 0 and x == xbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x][y+1][z] + next[x-1][y][z] + 2 * next[x][y][z])
                    elif y == ybin - 1 and x == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x][y-1][z] + next[x+1][y][z] + 2 * next[x][y][z])
                    elif y == ybin - 1 and x == xbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x][y-1][z] + next[x-1][y][z] + 2 * next[x][y][z])
                    elif y == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x][y+1][z] + next[x-1][y][z] + next[x+1][y][z] + next[x][y][z])
                    elif y == ybin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x][y-1][z] + next[x-1][y][z] + next[x+1][y][z] + next[x][y][z])
                    elif x == 0:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x][y-1][z] + next[x][y+1][z] + next[x+1][y][z] + next[x][y][z])
                    elif x == xbin - 1:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x][y-1][z] + next[x][y+1][z] + next[x-1][y][z] + next[x][y][z])
                    else:
                        D[x][y][z] = D[x][y][z] + self.f * (next[x][y-1][z] + next[x][y+1][z] + next[x-1][y][z] + next[x+1][y][z])

        for y in range(ybin):
            for x in range(xbin):
                for z in range(zbin):
                    if z == 0:
                        C_z[z] = C / Bc
                    #elif z == zbin:
                    #    C_z[z] = C / (Bc - A * C_z[z - 1])
                    else:
                        C_z[z] = C / (B - A * C_z[z - 1])
                for z in range(zbin):
                    if z == 0:
                        D_z[z] = D[x][y][z] / Bc
                    elif z == zbin - 1:
                        D_z[z] = (D[x][y][z] - A * D_z[z - 1]) / (Bc - A * C_z[z - 1])
                    else:
                        D_z[z] = (D[x][y][z] - A * D_z[z - 1]) / (B - A * C_z[z - 1])

                for z in reversed(range(zbin)):
                    if z == zbin - 1:
                        next[x][y][z] = D_z[z]
                    else:
                        next[x][y][z] = D_z[z] - C_z[z] * next[x][y][z + 1]

        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    space[x][y][z].set_molecule_qtty(molecule, index, next[x][y][z])
                    #print(space[x][y][z].molecules[molecule].values[0](index))

class ADI():

    def __init__(self, lbd, iter):
        self.iter = iter
        self.lbd = lbd
        self.X = None
        self.Y = None
        self.Z = None

    def create_open_matrixes(self, xbin, ybin, zbin):
        self.X = np.zeros([xbin, xbin])
        self.Y = np.zeros([ybin, ybin])
        self.Z = np.zeros([zbin, zbin])

        for i in range(xbin):
            for j in range(xbin):
                if i == j:
                    self.X[i,j] = - 2
                elif i == j + 1 or i == j - 1:
                    self.X[i, j] = 1
        self.X = np.linalg.inv((-self.lbd*self.X + np.identity(xbin)))

        for i in range(ybin):
            for j in range(ybin):
                if i == j:
                    self.Y[i, j] = - 2
                elif i == j + 1 or i == j - 1:
                    self.Y[i, j] = 1
        self.Y = np.linalg.inv((-self.lbd*self.Y + np.identity(ybin)))

        for i in range(zbin):
            for j in range(zbin):
                if i == j:
                    self.Z[i, j] = - 2
                elif i == j + 1 or i == j - 1:
                    self.Z[i, j] = 1
        self.Z = np.linalg.inv((-self.lbd*self.Z + np.identity(zbin)))

    def create_closed_matrixes(self, xbin, ybin, zbin):
        self.X = np.zeros([xbin, xbin])
        self.Y = np.zeros([ybin, ybin])
        self.Z = np.zeros([zbin, zbin])

        for i in range(xbin):
            for j in range(xbin):
                if i == j:
                    self.X[i,j] = - 2
                elif i == j + 1 or i == j - 1:
                    self.X[i, j] = 1
        self.X[0, 0] = -1
        self.X[0, 1] = 1
        self.X[xbin-1, xbin-1] = -1
        self.X[xbin-1, xbin - 2] = 1
        self.X = np.linalg.inv((-self.lbd*self.X + np.identity(xbin)))

        for i in range(ybin):
            for j in range(ybin):
                if i == j:
                    self.Y[i, j] = - 2
                elif i == j + 1 or i == j - 1:
                    self.Y[i, j] = 1
        self.Y[0, 0] = -1
        self.Y[0, 1] = 1
        self.Y[ybin - 1, ybin - 1] = -1
        self.Y[ybin - 1, ybin - 2] = 1
        self.Y = np.linalg.inv((-self.lbd*self.Y + np.identity(ybin)))

        for i in range(zbin):
            for j in range(zbin):
                if i == j:
                    self.Z[i, j] = - 2
                elif i == j + 1 or i == j - 1:
                    self.Z[i, j] = 1
        self.Z[0, 0] = -1
        self.Z[0, 1] = 1
        self.Z[zbin - 1, zbin - 1] = -1
        self.Z[zbin - 1, zbin - 2] = 1
        self.Z = np.linalg.inv((-self.lbd*self.Z + np.identity(zbin)))


    def create_periodic_matrixes(self, xbin, ybin, zbin):
        self.X = np.zeros([xbin, xbin])
        self.Y = np.zeros([ybin, ybin])
        self.Z = np.zeros([zbin, zbin])

        for i in range(xbin):
            for j in range(xbin):
                if i == j:
                    self.X[i,j] = - 2
                elif i == j + 1 or i == j - 1:
                    self.X[i, j] = 1
        self.X[0, xbin-1] = 1
        self.X[xbin-1, 0] = 1
        self.X = np.linalg.inv((-self.lbd*self.X + np.identity(xbin)))

        for i in range(ybin):
            for j in range(ybin):
                if i == j:
                    self.Y[i, j] = - 2
                elif i == j + 1 or i == j - 1:
                    self.Y[i, j] = 1
        self.Y[0, ybin - 1] = 1
        self.Y[ybin - 1, 0] = 1
        self.Y = np.linalg.inv((-self.lbd*self.Y + np.identity(ybin)))

        for i in range(zbin):
            for j in range(zbin):
                if i == j:
                    self.Z[i, j] = - 2
                elif i == j + 1 or i == j - 1:
                    self.Z[i, j] = 1
        self.Z[0, zbin - 1] = 1
        self.Z[zbin - 1, 0] = 1
        self.Z = np.linalg.inv((-self.lbd*self.Z + np.identity(zbin)))
        """
        file = open("/Users/henriquedeassis/mat.csv", "w")
        strg = ""
        for i in range(zbin):
            for j in range(zbin):
                strg += str(self.Z[i,j]) + ","
            strg += "\n"
        file.write(strg)
        file.close()
        """

    def solver(self, grid, molecule, index):
        xbin = len(grid)
        ybin = len(grid[0])
        zbin = len(grid[0][0])

        M = np.zeros([xbin, ybin, zbin])
        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    M[x,y,z] = grid[x][y][z].molecules[molecule].values[index]

        for _ in range(self.iter):
            for y in range(ybin):
                for z in range(zbin):
                    M[:,y,z] = self.X.dot(M[:,y,z])
            for x in range(xbin):
                for z in range(zbin):
                    M[x, :, z] = self.Y.dot(M[x, :, z])
            for y in range(ybin):
                for x in range(xbin):
                    M[x,y,:] = self.Z.dot(M[x,y,:])

        for x in range(xbin):
            for y in range(ybin):
                for z in range(zbin):
                    grid[x][y][z].molecules[molecule].values[index] = M[x,y,z] #if M[x,y,z] > 0 else 0